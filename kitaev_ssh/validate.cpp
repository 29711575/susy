#include "free_projector.h"
#include "pure_projector_fast.h" // regression for the shared Pfaffian phase fix
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>

using namespace kitaev_ssh;

namespace {
struct Metric { double maximum = 0., tolerance = 0.; int comparisons = 0; };
std::map<std::string,Metric> metrics;
void check(const std::string &name, double error, double tolerance) {
    auto &metric = metrics[name];
    metric.maximum = std::max(metric.maximum,error); metric.tolerance = tolerance; ++metric.comparisons;
    if (!std::isfinite(error) || error>tolerance) {
        std::ostringstream s; s<<std::setprecision(17)<<name<<": error="<<error<<" tolerance="<<tolerance;
        throw std::runtime_error(s.str());
    }
}
void require(bool condition, const std::string &name) { check(name,condition ? 0. : 1.,0.); }
double maxError(const MatType &a, const MatType &b) {
    if (a.rows()!=b.rows() || a.cols()!=b.cols()) return std::numeric_limits<double>::infinity();
    return (a-b).cwiseAbs().maxCoeff();
}
template<class F> void rejects(F f) {
    bool rejected = false;
    try { f(); } catch (const std::exception &) { rejected = true; }
    require(rejected,"invalid_input_rejection");
}

// Independent Fock oracle. Build operators in occupation-bit states, not by
// calling stage1_dense.h or converting the tested generator back into H.
std::vector<MatType> fockMajoranas(int modes) {
    const int dimension = 1<<modes;
    std::vector<MatType> c(2*modes,MatType::Zero(dimension,dimension));
    for (int j = 0; j<modes; ++j) for (int n = 0; n<dimension; ++n) {
        const int lower = n&((1<<j)-1);
        const double sign = __builtin_popcount(unsigned(lower))%2 ? -1. : 1.;
        c[j](n^(1<<j),n) = sign;
        c[modes+j](n^(1<<j),n) = DataType(0,(n&(1<<j)) ? sign : -sign);
    }
    return c;
}
struct Edge { int a,b; double t; };
std::vector<Edge> referenceEdges(int lx, int ly, Couplings j) {
    const int nc = lx*ly;
    std::vector<Edge> edges;
    for (int x = 0; x<lx; ++x) for (int y = 0; y<ly; ++y) {
        const int a = x*ly+y;
        edges.push_back({a,nc+((x+lx-1)%lx)*ly+y,j.Jx});
        edges.push_back({a,nc+x*ly+(y+ly-1)%ly,j.Jy});
        edges.push_back({a,nc+a,j.Jz});
    }
    return edges;
}
MatType fockHamiltonian(const std::vector<Edge> &edges, const std::vector<MatType> &c) {
    MatType h = MatType::Zero(c[0].rows(),c[0].cols());
    for (const auto &edge : edges) h += DataType(0,edge.t)*c[edge.a]*c[edge.b];
    return h;
}
MatType trialGenerator(int dimension, const std::vector<Edge> &edges) {
    MatType k = MatType::Zero(dimension,dimension);
    for (const auto &edge : edges) {
        k(edge.a,edge.b) += DataType(0,-2*edge.t);
        k(edge.b,edge.a) += DataType(0,2*edge.t);
    }
    return k;
}
MatType fockGreen(const cVecType &psi, const std::vector<MatType> &c) {
    MatType g = MatType::Zero(c.size(),c.size());
    for (int i = 0; i<int(c.size()); ++i) for (int j = 0; j<int(c.size()); ++j)
        if (i!=j) g(i,j) = -(psi.adjoint()*c[i]*c[j]*psi)(0);
    return g;
}
double fockParity(const cVecType &psi) {
    double p = 0.;
    for (int n = 0; n<psi.size(); ++n) p += (__builtin_popcount(unsigned(n))%2 ? -1 : 1)*std::norm(psi(n));
    return p;
}
MatType referenceExponential(const MatType &h, double tau) {
    Eigen::SelfAdjointEigenSolver<MatType> es(h);
    require(es.info()==Eigen::Success,"reference_eigensolver");
    return es.eigenvectors()*(-tau*es.eigenvalues().array()).exp().matrix().asDiagonal()*es.eigenvectors().adjoint();
}

void singleBond() {
    const double t = .73, dt = .19;
    // cA=sigma_x, cB=-sigma_y; f=(cA-i*cB)/2.
    MatType ca(2,2), cb(2,2);
    ca<<0,1,1,0; cb<<0,DataType(0,1),DataType(0,-1),0;
    const MatType h = DataType(0,t)*ca*cb;
    MatType expectedH = MatType::Zero(2,2); expectedH(0,0)=t; expectedH(1,1)=-t;
    check("single_bond_levels",maxError(h,expectedH),1e-14);
    MatType k = MatType::Zero(2,2); k(0,1)=DataType(0,-2*t); k(1,0)=-k(0,1);
    const auto trial = GaussianTrialState::fromMajoranaHamiltonian(k);
    cVecType occupied(2); occupied<<0,1;
    check("single_bond_ground_green",maxError(trial.G_T,fockGreen(occupied,{ca,cb})),2e-14);
    require(trial.fermionParity()==-1,"single_bond_matter_parity");
    const MatType b = expm(k,-dt), u = referenceExponential(h,dt), inverse = referenceExponential(h,-dt);
    const std::array<MatType,2> c{ca,cb};
    for (int i = 0; i<2; ++i) {
        MatType transformed = b(i,0)*ca+b(i,1)*cb;
        check("single_bond_operator_evolution",maxError(u*c[i]*inverse,transformed),2e-14);
    }
    const MatType expectedB = std::cosh(2*t*dt)*MatType::Identity(2,2)-std::sinh(2*t*dt)*k/(2*t);
    check("single_bond_factor_two",maxError(b,expectedB),2e-14);
    const auto anchor = pfaffianForSignOfProductWithStatus(MatType::Zero(2,2),trial.G_T);
    require(anchor.ok(),"single_bond_pfaffian_status");
    check("single_bond_raw_pfaffian_anchor",std::abs(anchor.value+DataType(1)),0.);
    for (int factors : {0,1,2,3}) {
        const std::vector<PureProjectorSlice> slices(factors,PureProjectorSlice(b,1.,"one_bond"));
        const auto weight = PureProjectorWeightEvaluator().evaluate(trial,slices);
        const auto stable = pureProjectorStableReferenceWeight(trial,slices,PureFastOptions{});
        require(weight.ok() && stable.ok(),"single_bond_pfaffian_status");
        const double exact = (occupied.adjoint()*referenceExponential(h,factors*dt)*occupied)(0).real();
        check("single_bond_weight",std::abs(weight.weight-exact),2e-14);
        check("single_bond_stable_weight",std::abs(stable.weight-exact),2e-14);
        check("single_bond_log_weight",std::abs(weight.log_abs_weight-factors*dt*t),2e-14);
    }
}

void geometry(int lx, int ly) {
    HoneycombLattice lattice(lx,ly);
    require(lattice.sites()==2*lx*ly && lattice.bondCount()==3*lx*ly,"geometry_counts");
    std::vector<std::array<int,3>> degree(lattice.sites(),{0,0,0});
    std::vector<int> occurrences(lattice.bondCount(),0);
    std::set<std::pair<int,int>> endpoints;
    int wrapsX = 0, wrapsY = 0;
    for (int cell = 0; cell<lattice.cells(); ++cell) {
        const auto r = lattice.cellCoordinate(cell);
        require(lattice.cellIndex(r.x,r.y)==cell,"geometry_indices");
    }
    for (int site = 0; site<lattice.sites(); ++site) {
        const auto r = lattice.siteCoordinate(site);
        require(lattice.siteIndex(r.x,r.y,r.sublattice)==site,"geometry_indices");
    }
    for (int i = 0; i<2; ++i) for (int j = 0; j<2; ++j) {
        const auto a = lattice.primitive(i), b = lattice.reciprocal(j);
        check("reciprocal_vectors",std::abs(a.x*b.x+a.y*b.y-(i==j ? 2*pi : 0)),3e-15);
    }
    for (const auto &b : lattice.bonds()) {
        require(lattice.bondIndex(b.cell.x,b.cell.y,b.type)==b.index,"geometry_indices");
        const int mu = static_cast<int>(b.type);
        ++degree[b.a][mu]; ++degree[b.b][mu]; endpoints.insert({b.a,b.b});
        const auto a = lattice.position(b.a), target = lattice.position(b.b);
        const auto a1 = lattice.primitive(0), a2 = lattice.primitive(1);
        check("geometry_coordinates",std::abs(target.x+b.boundary_wrap.x*lx*a1.x+b.boundary_wrap.y*ly*a2.x-a.x-b.displacement.x),1e-14);
        check("geometry_coordinates",std::abs(target.y+b.boundary_wrap.x*lx*a1.y+b.boundary_wrap.y*ly*a2.y-a.y-b.displacement.y),1e-14);
        check("bond_length",std::abs(std::hypot(b.displacement.x,b.displacement.y)-1.),1e-15);
        wrapsX += b.boundary_wrap.x!=0; wrapsY += b.boundary_wrap.y!=0;
    }
    require(int(endpoints.size())==lattice.bondCount(),"geometry_unique_bonds");
    require(wrapsX==ly && wrapsY==lx,"geometry_boundary_counts");
    for (const auto &d : degree) for (int n : d) require(n==1,"geometry_degree");
    for (int mu = 0; mu<3; ++mu) require(int(lattice.bonds(static_cast<BondType>(mu)).size())==lattice.cells(),"geometry_type_counts");
    require(int(lattice.plaquettes().size())==lattice.cells(),"geometry_plaquette_count");
    for (const auto &p : lattice.plaquettes()) {
        require(p.steps.size()==6 && lattice.pathEnd(p)==p.start,"geometry_plaquette_closure");
        std::set<int> vertices;
        double dx=0,dy=0;
        for (const auto &s : p.steps) {
            const auto &b = lattice.bond(s.bond);
            ++occurrences[s.bond]; vertices.insert(s.direction==1 ? b.a : b.b);
            dx += s.direction*b.displacement.x; dy += s.direction*b.displacement.y;
        }
        require(vertices.size()==6,"geometry_hexagon_vertices");
        check("geometry_plaquette_displacement",std::hypot(dx,dy),1e-14);
    }
    for (int n : occurrences) require(n==2,"geometry_plaquette_incidence");
    for (int wx : {-1,1}) for (int wy : {-1,1}) {
        const auto u = lattice.zeroFluxGauge(wx,wy);
        for (int p = 0; p<lattice.cells(); ++p) require(lattice.plaquetteFlux(p,u)==1,"zero_flux_sectors");
        for (int axis = 0; axis<2; ++axis) {
            const auto path = lattice.windingPath(axis);
            require(lattice.pathEnd(path)==path.start,"winding_closure");
            require(lattice.wilsonLoop(axis,u)==(axis==0 ? wx : wy),"wilson_sectors");
            double dx=0,dy=0;
            for (const auto &s : path.steps) { const auto &b=lattice.bond(s.bond); dx+=s.direction*b.displacement.x; dy+=s.direction*b.displacement.y; }
            const auto a = lattice.primitive(axis); const int length = axis==0 ? lx : ly;
            check("winding_displacement",std::hypot(dx+length*a.x,dy+length*a.y),2e-14);
        }
    }
}

void exactDiagonalization() {
    const auto c = fockMajoranas(4);
    for (int i = 0; i<8; ++i) for (int j = 0; j<8; ++j)
        check("fock_majorana_algebra",maxError(c[i]*c[j]+c[j]*c[i],(i==j ? 2. : 0.)*MatType::Identity(16,16)),0.);
    const Couplings couplings;
    const KitaevSSHModel model(HoneycombLattice(2,2),couplings);
    const auto edges = referenceEdges(2,2,couplings);
    const MatType h = fockHamiltonian(edges,c);
    Eigen::SelfAdjointEigenSolver<MatType> target(h);
    require(target.info()==Eigen::Success,"reference_eigensolver");
    // Inhomogeneous same-sublattice bilinears give a noncommuting Gaussian
    // trial with the same parity and nonzero target-ground-state overlap.
    auto otherEdges = edges; otherEdges.push_back({0,1,.17}); otherEdges.push_back({4,6,.11});
    for (bool noncommuting : {false,true}) {
        const auto &initialEdges = noncommuting ? otherEdges : edges;
        const MatType ht = fockHamiltonian(initialEdges,c);
        Eigen::SelfAdjointEigenSolver<MatType> initial(ht);
        require(initial.info()==Eigen::Success,"reference_eigensolver");
        const cVecType psi0 = initial.eigenvectors().col(0);
        const auto trial = noncommuting ? GaussianTrialState::fromMajoranaHamiltonian(trialGenerator(8,initialEdges)) : freeGroundState(model);
        check("ed_trial_green",maxError(trial.G_T,fockGreen(psi0,c)),1e-12);
        check("ed_trial_parity",std::abs(trial.fermionParity()-fockParity(psi0)),1e-12);
        const double overlap = std::norm((target.eigenvectors().col(0).adjoint()*psi0)(0));
        if (noncommuting) {
            require((h*ht-ht*h).norm()>1e-3,"ed_noncommuting_trial");
            require(overlap>1e-3 && overlap<1-1e-6,"ed_trial_overlap");
            std::cout<<"DETAIL ed_noncommuting_ground_overlap="<<overlap<<" commutator_norm="<<(h*ht-ht*h).norm()<<'\n';
        }
        double previousEnergy = std::numeric_limits<double>::infinity();
        for (double tauPerBranch : {.05,.2,.8,2.}) {
            const cVecType evolved = referenceExponential(h,tauPerBranch)*psi0;
            const double norm = evolved.squaredNorm();
            const cVecType psi = evolved/std::sqrt(norm);
            const MatType expectedGreen = fockGreen(psi,c);
            const double expectedEnergy = (psi.adjoint()*h*psi)(0).real();
            for (int block : {1,3,8}) {
                const auto projected = projectStatic(model,trial,tauPerBranch,.025,block);
                check("ed_projector_green",maxError(projected.green,expectedGreen),2e-10);
                check("ed_projector_energy",std::abs(projected.energy-expectedEnergy),2e-10);
                check("ed_projector_log_norm",std::abs(projected.log_norm_squared-std::log(norm)),2e-10);
                check("ed_projector_norm_relative",std::abs(std::expm1(projected.log_norm_squared-std::log(norm))),2e-10);
            }
            require(expectedEnergy<=previousEnergy+1e-12,"ed_projection_energy_descent"); previousEnergy=expectedEnergy;
            // Existing sequential Pfaffian lift, independent of the adapter's
            // positive U^dagger U norm. Keep the raw-product oracle short.
            if (tauPerBranch<=.2) {
                PureProjectorWeightOptions options; options.mode=PureProjectorWeightMode::GenericComplex;
                const auto weight=PureProjectorWeightEvaluator(options).evaluate(trial,staticContour(model,tauPerBranch,.025));
                require(weight.ok(),"ed_pfaffian_status");
                check("ed_pfaffian_weight_relative",std::abs(weight.weight/norm-DataType(1)),2e-10);
                check("ed_pfaffian_phase",std::abs(weight.complex_phase-DataType(1)),2e-12);
            }
        }
    }
}

// Cell-periodic Bloch basis exp(i*(k1*x+k2*y)), no basis displacement phase.
DataType blochF(double k1, double k2, Couplings j = {}) {
    return j.Jz+j.Jx*std::exp(DataType(0,-k1))+j.Jy*std::exp(DataType(0,-k2));
}
struct BlochReference { std::vector<double> spectrum; MatType green; double energy=0.; };
BlochReference blochReference(int lx, int ly, Couplings j) {
    const int nc=lx*ly;
    BlochReference result; result.green=MatType::Zero(2*nc,2*nc);
    for (int nx=0; nx<lx; ++nx) for (int ny=0; ny<ly; ++ny) {
        const double k1=2*pi*nx/lx,k2=2*pi*ny/ly;
        const auto f=blochF(k1,k2,j); const double magnitude=std::abs(f);
        require(magnitude>1e-10,"bloch_finite_grid_nonzero");
        result.spectrum.push_back(-2*magnitude); result.spectrum.push_back(2*magnitude); result.energy-=magnitude;
        for (int x=0; x<lx; ++x) for (int y=0; y<ly; ++y)
            for (int xx=0; xx<lx; ++xx) for (int yy=0; yy<ly; ++yy) {
                const DataType phase=std::exp(DataType(0,k1*(x-xx)+k2*(y-yy)))/double(nc);
                result.green(x*ly+y,nc+xx*ly+yy) += phase*DataType(0,-1)*f/magnitude;
                result.green(nc+x*ly+y,xx*ly+yy) += phase*DataType(0,1)*std::conj(f)/magnitude;
            }
    }
    std::sort(result.spectrum.begin(),result.spectrum.end());
    return result;
}

void freeLattice(int length) {
    const HoneycombLattice lattice(length,length);
    const KitaevSSHModel model(lattice);
    const auto k=model.frameworkGenerator();
    const auto a=model.realAntisymmetricA();
    check("generator_hermitian",maxError(k,k.adjoint()),0.);
    check("generator_skew",maxError(k,-k.transpose()),0.);
    check("real_A_skew",(a+a.transpose()).cwiseAbs().maxCoeff(),0.);
    const auto reference=blochReference(length,length,model.couplings());
    Eigen::SelfAdjointEigenSolver<MatType> es(k);
    require(es.info()==Eigen::Success,"model_eigensolver");
    for (int n=0; n<model.majoranas(); ++n) {
        check("bloch_spectrum",std::abs(es.eigenvalues()(n)-reference.spectrum[n]),2e-12);
        check("spectrum_pairing",std::abs(es.eigenvalues()(n)+es.eigenvalues()(model.majoranas()-1-n)),2e-12);
    }
    const auto trial=freeGroundState(model);
    require(trial.Phi.rows()==2*length*length && trial.Phi.cols()==length*length,"trial_dimensions");
    check("bloch_ground_green",maxError(trial.G_T,reference.green),2e-12);
    std::cout<<"DETAIL L="<<length<<" sites="<<lattice.sites()<<" bonds="<<lattice.bondCount()
             <<" matter_parity="<<trial.fermionParity()<<" E0="<<reference.energy
             <<" min_positive_energy="<<es.eigenvalues()(lattice.cells())<<'\n';
    for (double tauPerBranch : {.1,.8,2.}) {
        const auto projected=projectStatic(model,trial,tauPerBranch,.05,4);
        check("ground_projector_green",maxError(projected.green,reference.green),2e-10);
        check("ground_projector_energy",std::abs(projected.energy-reference.energy),2e-10);
        check("ground_projector_log_norm",std::abs(projected.log_norm_squared+2*tauPerBranch*reference.energy),2e-10);
        check("projector_green_residual",projected.green_residual,1e-10);
        check("projector_solve_residual",projected.solve_residual,1e-10);
    }
    std::mt19937_64 rng(20260915+length); std::uniform_real_distribution<double> displacement(-3.,3.);
    std::vector<double> x(lattice.bondCount()); for (auto &v:x) v=displacement(rng);
    const KitaevSSHModel decoupled(lattice,{},x);
    check("g0_X_matrix_independence",maxError(k,decoupled.frameworkGenerator()),0.);
    const auto original=projectStatic(model,trial,.4,.05,4), changed=projectStatic(decoupled,trial,.4,.05,4);
    check("g0_X_observable_independence",maxError(original.green,changed.green),0.);
    check("g0_X_observable_independence",std::abs(original.energy-changed.energy),0.);
    check("g0_X_observable_independence",std::abs(original.log_norm_squared-changed.log_norm_squared),0.);
    for (int sample=0; sample<3; ++sample) {
        std::vector<int> d(model.majoranas()); MatType D=MatType::Zero(model.majoranas(),model.majoranas());
        int determinant=1;
        for (int i=0; i<model.majoranas(); ++i) { d[i]=(rng()&1) ? 1 : -1; D(i,i)=d[i]; determinant*=d[i]; }
        const auto transformed=model.gaugeTransformed(d);
        check("gauge_generator",maxError(transformed.frameworkGenerator(),D*k*D),0.);
        const auto trialD=GaussianTrialState::fromPhi(D*trial.Phi);
        require(trialD.fermionParity()==determinant*trial.fermionParity(),"gauge_matter_parity");
        const auto projected=projectStatic(transformed,trialD,.4,.05,4);
        check("gauge_green",maxError(projected.green,D*original.green*D),2e-10);
        check("gauge_energy",std::abs(projected.energy-original.energy),2e-10);
        check("gauge_log_norm",std::abs(projected.log_norm_squared-original.log_norm_squared),2e-10);
        std::vector<BondPath> paths;
        for (const auto &bond:lattice.bonds()) { paths.push_back({bond.a,{{bond.index,1}}}); paths.push_back({bond.b,{{bond.index,-1}}}); }
        auto winding=lattice.windingPath(0); winding.steps.push_back({lattice.bondIndex(0,0,BondType::Z),1}); paths.push_back(winding);
        for (const auto &p:lattice.plaquettes()) { auto open=p; open.steps.resize(3); paths.push_back(open); }
        for (const auto &path:paths)
            check("gauge_string",std::abs(gaugeStringGreen(model,original.green,path)-gaugeStringGreen(transformed,projected.green,path)),2e-10);
        for (int axis=0; axis<2; ++axis) require(lattice.wilsonLoop(axis,transformed.gauge())==1,"gauge_wilson_invariance");
    }
}

void staticInputsAndZeroModes() {
    HoneycombLattice lattice(2,2);
    std::vector<double> x(lattice.bondCount()); for (int n=0; n<int(x.size()); ++n) x[n]=.1*(n-3);
    Couplings j; j.g=.37;
    auto u=lattice.zeroFluxGauge(-1,1);
    const KitaevSSHModel model(lattice,j,x,u);
    auto expected=referenceEdges(2,2,j);
    for (int n=0; n<int(expected.size()); ++n) expected[n].t=(expected[n].t+j.g*x[n])*u[n];
    check("static_X_generator",maxError(model.frameworkGenerator(),trialGenerator(8,expected)),0.);
    const auto c=fockMajoranas(4); const auto h=fockHamiltonian(expected,c);
    Eigen::SelfAdjointEigenSolver<MatType> es(h);
    check("static_X_trial_green",maxError(freeGroundState(model).G_T,fockGreen(es.eigenvectors().col(0),c)),2e-12);
    rejects([&]{HoneycombLattice invalid(1,4);});
    rejects([&]{lattice.siteCoordinate(-1);}); rejects([&]{lattice.cellIndex(2,0);});
    rejects([&]{lattice.gaugeProduct({0,{{0,-1}}},u);});
    rejects([&]{lattice.gaugeProduct({0,{{0,0}}},u);});
    rejects([&]{KitaevSSHModel invalid(lattice,j,{1.,2.});});
    auto badX=x; badX[0]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{KitaevSSHModel invalid(lattice,j,badX);});
    auto badU=u; badU[0]=0; rejects([&]{KitaevSSHModel invalid(lattice,j,x,badU);});
    badU=u; badU[0]*=-1; rejects([&]{KitaevSSHModel invalid(lattice,j,x,badU);});
    auto badJ=j; badJ.g=std::numeric_limits<double>::infinity(); rejects([&]{KitaevSSHModel invalid(lattice,badJ);});
    rejects([&]{model.propagator(0.);});
    rejects([&]{staticContour(model,.31,.1);});
    MatType realSym=MatType::Identity(8,8); rejects([&]{KitaevSSHModel::validateGenerator(realSym,8);});
    // Isotropic L=6 PBC includes the exact Dirac nodes. Failure is deliberate.
    Couplings isotropic; isotropic.Jz=1.;
    const KitaevSSHModel zero(HoneycombLattice(6,6),isotropic);
    bool diagnosed=false;
    try { freeGroundState(zero); } catch (const GaussianTrialStateError &e) {
        diagnosed=e.code()==GaussianTrialStateStatus::zero_mode_ambiguous;
    }
    require(diagnosed,"exact_zero_mode_diagnostic");
}

void diracNodes() {
    const double q=std::acos(-.3), kx=2*q/std::sqrt(3.);
    auto energy=[](double x,double y) { return 2*std::abs(blochF(std::sqrt(3.)*x/2+1.5*y,-std::sqrt(3.)*x/2+1.5*y)); };
    const double vx=2*std::sqrt(3.)*std::sin(q),vy=1.8;
    for (double sign : {-1.,1.}) {
        check("dirac_node_energy",energy(sign*kx,0.),2e-15);
        for (double h : {1e-3,5e-4,1e-4}) {
            const double sx=(energy(sign*kx+h,0.)+energy(sign*kx-h,0.))/(2*h);
            const double sy=(energy(sign*kx,h)+energy(sign*kx,-h))/(2*h);
            // Symmetric secants have O(h^2) truncation. Analytic bounds at
            // h=1e-3 are below 2e-7 relative; smaller h checks linear scaling.
            check("dirac_x_slope_relative",std::abs(sx/vx-1),3e-7);
            check("dirac_y_slope_relative",std::abs(sy/vy-1),3e-7);
        }
    }
    std::cout<<"DETAIL dirac_k1="<<q<<" dirac_k2="<<-q<<" cartesian_kx="<<kx
             <<" cartesian_ky=0 vx="<<vx<<" vy="<<vy<<'\n';
}
} // namespace

int main() {
    std::cout<<std::setprecision(17);
    try {
        singleBond(); geometry(2,2); geometry(3,4); geometry(4,4); geometry(6,6);
        exactDiagonalization(); staticInputsAndZeroModes(); freeLattice(4); freeLattice(6); diracNodes();
        for (const auto &item:metrics) std::cout<<"PASS "<<item.first<<" max_error="<<item.second.maximum
            <<" tolerance="<<item.second.tolerance<<" comparisons="<<item.second.comparisons<<'\n';
        std::cout<<"PASS P0-P1 fixed-gauge matter-sector baseline; no phonon sampling or spin projection\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr<<"FAIL kitaev_ssh_validate: "<<e.what()<<'\n'; return 1;
    }
}
