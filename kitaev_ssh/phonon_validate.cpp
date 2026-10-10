#include "phonon_projector.h"
#include "pure_projector_protocol.h"
#include <iostream>
#include <map>
#include <sstream>

using namespace kitaev_ssh;
namespace {
struct Metric { double error=0.,tolerance=0.; };
std::map<std::string,Metric> metrics;
void check(const std::string &name,double error,double tolerance) {
    auto &m=metrics[name]; m.error=std::max(m.error,error); m.tolerance=tolerance;
    if (!std::isfinite(error) || error>tolerance) {
        std::ostringstream s; s<<std::setprecision(17)<<name<<" error="<<error<<" tolerance="<<tolerance;
        throw std::runtime_error(s.str());
    }
}
void require(bool ok,const std::string &name) { check(name,ok?0.:1.,0.); }
template<class F> void rejects(F f) { bool failed=false; try { f(); } catch(const std::exception&) { failed=true; } require(failed,"input_transaction_rejection"); }
double error(const MatType &a,const MatType &b) { return (a-b).cwiseAbs().maxCoeff(); }

// Recover Q from implemented quadratic action by exact polarization, not by
// copying the production link coefficients into the reference matrix.
Eigen::MatrixXd precision(const PhononPath &zero) {
    const int n=zero.nodes(); Eigen::MatrixXd q=Eigen::MatrixXd::Zero(n,n);
    for (int i=0;i<n;++i) q(i,i)=2*zero.withCoordinate(i,0,1.).action();
    for (int i=0;i<n;++i) for (int j=i+1;j<n;++j) {
        q(i,j)=q(j,i)=zero.withCoordinate(i,0,1.).withCoordinate(j,0,1.).action()-(q(i,i)+q(j,j))/2;
    }
    return q;
}
Eigen::MatrixXd covariance(const PhononPath &p) {
    Eigen::MatrixXd c(p.nodes(),p.nodes());
    for (int i=0;i<p.nodes();++i) for (int j=0;j<p.nodes();++j)
        c(i,j)=std::exp(-p.parameters().omega*p.dt()*std::abs(i-j))/(2*p.parameters().mass*p.parameters().omega);
    return c;
}
double gaussianLogIntegral(const Eigen::MatrixXd &q) {
    Eigen::LLT<Eigen::MatrixXd> llt(q); require(llt.info()==Eigen::Success,"positive_phonon_precision");
    const Eigen::MatrixXd l=llt.matrixL();
    return q.rows()/2.*std::log(2*pi)-l.diagonal().array().log().sum();
}
void freePhonons() {
    for (int m:{12,24,48}) {
        const PhononPath p(1,m,2.4/m,{1.7,.8});
        const auto q=precision(p), c=q.llt().solve(Eigen::MatrixXd::Identity(p.nodes(),p.nodes())).eval();
        check("free_covariance",(c-covariance(p)).cwiseAbs().maxCoeff(),2e-12);
        check("free_log_norm",std::abs(p.logNormalization()+gaussianLogIntegral(q)+.8*2.4/2),5e-12);
        const double alpha=1.7*.8,z=.8*p.dt();
        for (int j=0;j<p.nodes();++j) {
            check("free_precision_diagonal",std::abs(q(j,j)-alpha*((j==0||j==m)?1+1/std::tanh(z):2/std::tanh(z))),2e-12);
            if(j<m) check("free_precision_neighbor",std::abs(q(j,j+1)+alpha/std::sinh(z)),2e-12);
        }
    }
    std::mt19937_64 rng(20260916);
    PhononPath blank(3,12,.2,{1.7,.8});
    const auto path=blank.freeSample(rng);
    for (int j=0;j<path.nodes();++j) for(int b=0;b<path.bonds();++b) {
        const double value=path.coordinate(j,b)+.31;
        check("local_phonon_log_ratio",std::abs(path.localLogRatio(j,b,value)-(path.withCoordinate(j,b,value).logWeight()-path.logWeight())),2e-12);
    }
    std::stringstream snapshot; path.save(snapshot); const auto restored=PhononPath::load(snapshot);
    check("path_snapshot_weight",std::abs(restored.logWeight()-path.logWeight()),0.);
    for (int j=0;j<path.nodes();++j) for(int b=0;b<path.bonds();++b)
        check("path_snapshot_coordinates",std::abs(restored.coordinate(j,b)-path.coordinate(j,b)),0.);
    // AR sampler checks have statistical, not roundoff, tolerances: six SE.
    const int samples=100000; double mean=0.,var0=0.,varMid=0.,corr=0.,cross=0.;
    for(int s=0;s<samples;++s) {
        const auto x=blank.freeSample(rng); const double a=x.coordinate(0,0),b=x.coordinate(6,0);
        mean+=b; var0+=a*a; varMid+=b*b; corr+=a*b; cross+=b*x.coordinate(6,1);
    }
    const double v=1/(2*1.7*.8), c=v*std::exp(-.8*1.2);
    check("free_sampler_mean_SE",std::abs(mean/samples)/std::sqrt(v/samples),6.);
    check("free_sampler_variance_SE",std::max(std::abs(var0/samples-v),std::abs(varMid/samples-v))/(v*std::sqrt(2./samples)),6.);
    check("free_sampler_two_time_SE",std::abs(corr/samples-c)/std::sqrt((v*v+c*c)/samples),6.);
    check("free_sampler_cross_bond_SE",std::abs(cross/samples)/(v/std::sqrt(double(samples))),6.);
    require(std::isfinite(PhononPath(1,2,1000.,{1.,1.}).logWeight()),"large_kernel_log_finite");
    require(std::isfinite(PhononPath(1,2,1e-12,{1.,1.}).logWeight()),"small_kernel_log_finite");
    rejects([]{PhononParameters p(0.,1.);}); rejects([]{PhononPath p(1,3,.1,{1.,1.});});
    rejects([&]{path.withCoordinate(0,0,std::numeric_limits<double>::infinity());});
    rejects([]{std::stringstream in("KITAEV_SSH_PHONON_PATH 1\n1 2 .1 1 1\n0");PhononPath::load(in);});
}

void compare(const CoupledPhononPath::Snapshot &a,const CoupledPhononPath::Snapshot &b) {
    for(int f=0;f<int(a.factors.size());++f) check("transaction_factors",error(MatType(a.factors[f].matrix),MatType(b.factors[f].matrix)),0.);
    for(int j=0;j<a.path.nodes();++j) for(int n=0;n<a.path.bonds();++n)
        check("transaction_coordinates",std::abs(a.path.coordinate(j,n)-b.path.coordinate(j,n)),0.);
    check("transaction_log_weight",std::abs(a.weight.logAbsTotal()-b.weight.logAbsTotal()),2e-12);
    check("transaction_phase",std::abs(a.weight.matter_phase-b.weight.matter_phase),2e-12);
    check("transaction_green",error(a.weight.center_green,b.weight.center_green),2e-12);
}
// Independent 2x2 honeycomb Fock oracle for a non-mirrored, noncommuting path.
std::vector<MatType> fockMajoranas() {
    std::vector<MatType> c(8,MatType::Zero(16,16));
    for(int a=0;a<4;++a) for(int n=0;n<16;++n) {
        const double sign=__builtin_popcount(unsigned(n&((1<<a)-1)))%2?-1.:1.;
        c[a](n^(1<<a),n)=sign;
        c[4+a](n^(1<<a),n)=DataType(0,(n&(1<<a))?sign:-sign);
    }
    return c;
}
MatType fockH(const std::vector<MatType> &c,const std::vector<double> &x,double g) {
    MatType h=MatType::Zero(16,16);
    for(int a=0;a<2;++a) for(int b=0;b<2;++b) {
        int i=2*a+b;
        h+=DataType(0,1+g*x[3*i])*c[i]*c[4+2*(1-a)+b];
        h+=DataType(0,1+g*x[3*i+1])*c[i]*c[4+2*a+1-b];
        h+=DataType(0,.5+g*x[3*i+2])*c[i]*c[4+i];
    }
    return h;
}
MatType fockExp(const MatType &h,double t) {
    Eigen::SelfAdjointEigenSolver<MatType> es(h); require(es.info()==Eigen::Success,"fock_eigensolver");
    return es.eigenvectors()*(-t*es.eigenvalues().array()).exp().matrix().asDiagonal()*es.eigenvectors().adjoint();
}
void coupledMapping() {
    std::mt19937_64 rng(1729);
    const KitaevSSHModel model(HoneycombLattice(2,2),p2Couplings(.4));
    const auto trial=freeGroundState(model); const auto generator=honeycombGenerator(model);
    const auto path=PhononPath(12,6,.07,{1.7,.8}).freeSample(rng);
    CoupledPhononPath state(path,trial,generator);
    const auto before=state.state();
    for(int f=0;f<2*path.intervals();++f) {
        const int j=factorNode(f,path.intervals());
        check("mapping_half_steps",error(MatType(before.factors[f].matrix),nodeHalfStep(path,j,8,generator)),0.);
    }
    require(factorNode(path.intervals()-1,path.intervals())==path.center() && factorNode(path.intervals(),path.intervals())==path.center(),"center_two_half_steps");
    require(path.node(1)!=path.node(path.intervals()-1),"independent_bra_ket_nodes");
    for(int j=0;j<path.nodes();++j) {
        auto proposal=state.propose(j,2,-3.); // negative effective hopping allowed
        compare(proposal.candidate(),state.fullReference(path.withCoordinate(j,2,-3.)));
        check("transaction_local_phonon",std::abs(proposal.localPhononLogRatio()-(proposal.candidate().weight.log_phonon-before.weight.log_phonon)),2e-12);
        state.reject(proposal); compare(state.state(),before);
        rejects([&]{state.accept(proposal);});
    }
    auto accepted=state.propose(path.center(),0,path.coordinate(path.center(),0)+.2);
    auto stale=state.propose(0,1,path.coordinate(0,1)-.1);
    const auto expected=accepted.candidate(); state.accept(accepted); compare(state.state(),expected);
    rejects([&]{state.accept(stale);}); rejects([&]{state.propose(0,0,std::numeric_limits<double>::quiet_NaN());});
    compare(state.state(),expected); compare(state.state(),state.fullReference(state.state().path));
    CoupledPhononPath other(path,trial,generator);
    auto foreign=other.propose(0,0,.2); rejects([&]{state.accept(foreign);});
    // Inject a failure after private factor construction: live caches must survive.
    MajoranaGenerator failing=[generator](const std::vector<double> &x) {
        if(x[0]==123.) throw std::runtime_error("injected generator failure"); return generator(x);
    };
    CoupledPhononPath guarded(path,trial,failing); const auto guardBefore=guarded.state();
    rejects([&]{guarded.propose(1,0,123.);}); compare(guarded.state(),guardBefore);

    const auto c=fockMajoranas();
    Eigen::SelfAdjointEigenSolver<MatType> initial(fockH(c,std::vector<double>(12,0.),.4));
    const cVecType psi=initial.eigenvectors().col(0);
    cVecType right=psi; Eigen::RowVectorXcd left=psi.adjoint();
    for(int l=1;l<=path.center();++l) right=fockExp(fockH(c,path.node(l),.4),path.dt()/2)*fockExp(fockH(c,path.node(l-1),.4),path.dt()/2)*right;
    for(int l=path.intervals();l>path.center();--l) left=left*fockExp(fockH(c,path.node(l),.4),path.dt()/2)*fockExp(fockH(c,path.node(l-1),.4),path.dt()/2);
    const DataType amplitude=(left*right)(0);
    check("fock_path_log_matter",std::abs(before.weight.log_abs_matter-std::log(std::abs(amplitude))),2e-12);
    check("fock_path_phase",std::abs(before.weight.matter_phase-amplitude/std::abs(amplitude)),2e-12);
    MatType green=MatType::Zero(8,8);
    for(int i=0;i<8;++i) for(int j=0;j<8;++j) if(i!=j) green(i,j)=-(left*c[i]*c[j]*right)(0)/amplitude;
    check("fock_center_green",error(before.weight.center_green,green),2e-12);
    const auto wrong=pureProjectorEndpointRebuild(trial,before.factors,path.intervals()-1,4,PureProjectorOptions{});
    require(wrong.ok() && error(wrong.green.adjoint(),green)>1e-5,"center_insertion_sensitive");
    // Known scalar lifts test phase plumbing, not a claim that this particular
    // honeycomb path has a sign problem. Neither sign nor complex phase erased.
    for (DataType phase:{DataType(-1),std::exp(DataType(0,.37))}) {
        auto factors=before.factors; factors[0].eta=phase;
        const auto weight=evaluatePath(path,trial,factors,4);
        check("retained_scalar_phase",std::abs(weight.matter_phase-before.weight.matter_phase*phase),2e-12);
        check("retained_scalar_magnitude",std::abs(weight.log_abs_matter-before.weight.log_abs_matter),2e-12);
    }
}

std::vector<PureProjectorSlice> midpointFactors(const KitaevSSHModel &model,
                                                const PhononPath &path,
                                                const PureImaginaryTimeProtocol &protocol) {
    require(protocol.slices()==path.intervals(),"midpoint_slice_count");
    std::vector<PureProjectorSlice> factors;
    factors.reserve(2*path.intervals());
    for (int link=0; link<path.intervals(); ++link) {
        const double g=protocol.midpointValue(link);
        MatType leftGenerator=model.frameworkGeneratorWith(path.node(link),g);
        MatType rightGenerator=model.frameworkGeneratorWith(path.node(link+1),g);
        const MatType left=expm(leftGenerator,-path.dt()/2);
        const MatType right=expm(rightGenerator,-path.dt()/2);
        factors.push_back(PureProjectorSlice::fromValidatedDense(
            left,1.,"driven_midpoint:left:"+std::to_string(link)));
        factors.push_back(PureProjectorSlice::fromValidatedDense(
            right,1.,"driven_midpoint:right:"+std::to_string(link)));
    }
    return factors;
}

// Canonical midpoint protocol reference. Each link owns one midpoint g and
// both endpoint half steps; no node-level merge, cache reuse, or recovery
// reconstruction is used here. The exact Fock product is the oracle.
void midpointDrivenReference() {
    const auto c=fockMajoranas();
    const KitaevSSHModel model(HoneycombLattice(2,2),p2Couplings(0.));
    const auto trial=freeGroundState(model);
    const std::vector<std::pair<double,double>> protocols={
        {0.,.4},{.1,.2},{.25,.1},{.05,.6},
        {0.,.8},{.15,.3},{.3,.05},{.02,.9},{.4,0.}};
    for (int sample=0; sample<int(protocols.size()); ++sample) {
        const PhononPath path(12,4,.1,{1.7,.8});
        PhononPath configured=path;
        for (int j=0;j<configured.nodes();++j) for (int b=0;b<configured.bonds();++b) {
            const double value=.013*(sample+1)*(j+1)*(b+2)-.021*((j+sample)%3);
            configured=configured.withCoordinate(j,b,value);
        }
        const auto protocol=(sample==int(protocols.size())-1)
            ?PureImaginaryTimeProtocol::constant(.4,path.totalTau(),path.dt())
            :PureImaginaryTimeProtocol::linearRamp(
                protocols[sample].first,protocols[sample].second,path.totalTau(),path.dt());
        const auto factors=midpointFactors(model,configured,protocol);
        double largestAdjacentCommutator=0.;
        for (int link=1; link<configured.intervals(); ++link) {
            const auto previous=model.frameworkGeneratorWith(configured.node(link),
                                                              protocol.midpointValue(link-1));
            const auto next=model.frameworkGeneratorWith(configured.node(link),
                                                         protocol.midpointValue(link));
            if (sample!=int(protocols.size())-1)
                require((previous-next).norm()>1e-8,"midpoint_node_context_distinct");
            const MatType hPrevious=fockH(c,configured.node(link),protocol.midpointValue(link-1));
            const MatType hNext=fockH(c,configured.node(link),protocol.midpointValue(link));
            largestAdjacentCommutator=std::max(largestAdjacentCommutator,
                (hPrevious*hNext-hNext*hPrevious).norm());
        }
        if (sample!=int(protocols.size())-1)
            require(largestAdjacentCommutator>1e-8,"midpoint_noncommuting_slices");
        const auto measured=evaluatePath(configured,trial,factors,4);
        if (sample==int(protocols.size())-1) {
            const KitaevSSHModel fixedModel(HoneycombLattice(2,2),p2Couplings(.4));
            const auto fixedFactors=pathFactors(configured,8,honeycombGenerator(fixedModel));
            for (int f=0;f<int(factors.size());++f)
                check("midpoint_constant_factor_regression",
                      error(MatType(factors[f].matrix),MatType(fixedFactors[f].matrix)),0.);
            const auto fixedWeight=evaluatePath(configured,trial,fixedFactors,4);
            check("midpoint_constant_log_regression",
                  std::abs(measured.log_abs_matter-fixedWeight.log_abs_matter),2e-12);
            check("midpoint_constant_phase_regression",
                  std::abs(measured.matter_phase-fixedWeight.matter_phase),2e-12);
            check("midpoint_constant_green_regression",
                  error(measured.center_green,fixedWeight.center_green),2e-12);
        }
        Eigen::SelfAdjointEigenSolver<MatType> initial(
            fockH(c,std::vector<double>(12,0.),0.));
        require(initial.info()==Eigen::Success,"midpoint_fock_initialization");
        const cVecType psi=initial.eigenvectors().col(0);
        cVecType right=psi; Eigen::RowVectorXcd left=psi.adjoint();
        for (int link=1; link<=configured.center(); ++link) {
            const double g=protocol.midpointValue(link-1);
            right=fockExp(fockH(c,configured.node(link),g),configured.dt()/2)*
                  fockExp(fockH(c,configured.node(link-1),g),configured.dt()/2)*right;
        }
        for (int link=configured.intervals(); link>configured.center(); --link) {
            const double g=protocol.midpointValue(link-1);
            left=left*fockExp(fockH(c,configured.node(link),g),configured.dt()/2)*
                 fockExp(fockH(c,configured.node(link-1),g),configured.dt()/2);
        }
        const DataType amplitude=(left*right)(0);
        const DataType phase=amplitude/std::abs(amplitude);
        check("midpoint_fock_log_matter",std::abs(measured.log_abs_matter-std::log(std::abs(amplitude))),2e-12);
        check("midpoint_fock_phase",std::abs(measured.matter_phase-phase),2e-12);
        MatType green=MatType::Zero(8,8);
        for (int i=0;i<8;++i) for (int j=0;j<8;++j) if(i!=j)
            green(i,j)=-(left*c[i]*c[j]*right)(0)/amplitude;
        check("midpoint_fock_center_green",error(measured.center_green,green),2e-12);
    }
}

void coupledSingleBond() {
    const double mass=1.7,omega=.8,g=.6,J=.5,tau=2.4;
    MajoranaGenerator generator=[&](const std::vector<double> &x) {
        MatType k=MatType::Zero(2,2); k(0,1)=DataType(0,-2*(J+g*x.at(0))); k(1,0)=-k(0,1); return k;
    };
    for(int eta:{-1,1}) {
        const auto trial=GaussianTrialState::fromMajoranaHamiltonian(-double(eta)*generator({0.}));
        double previousZ=0.,previousMean=0.;
        for(int m:{12,24,48}) {
            const PhononPath path(1,m,tau/m,{mass,omega});
            const auto q=precision(path); const Eigen::MatrixXd c=q.llt().solve(Eigen::MatrixXd::Identity(path.nodes(),path.nodes()));
            auto eval=[&](const PhononPath &p) { return evaluatePath(p,trial,pathFactors(p,2,generator),4); };
            const auto base=eval(path); Eigen::VectorXd source(path.nodes()),w(path.nodes());
            for(int j=0;j<path.nodes();++j) {
                // Read the linear source from the actual signed Majorana overlap.
                const auto shifted=eval(path.withCoordinate(j,0,1.));
                source(j)=shifted.log_abs_matter-base.log_abs_matter;
                w(j)=path.dt()*((j==0||j==m)?.5:1.);
                check("single_bond_source",std::abs(source(j)+eta*g*w(j)),2e-12);
                check("single_bond_phase",std::abs(shifted.matter_phase-DataType(1)),2e-12);
            }
            const double logZ=path.logNormalization()+gaussianLogIntegral(q)+base.log_abs_matter+.5*source.dot(c*source);
            const double mean=(c*source)(path.center());
            const auto analyticC=covariance(path);
            const double discrete=-eta*J*tau-omega*tau/2+g*g/2*w.dot(analyticC*w);
            const double discreteMean=-eta*g*(analyticC*w)(path.center());
            check("coupled_discrete_log_norm",std::abs(logZ-discrete),5e-12);
            check("coupled_discrete_center_mean",std::abs(mean-discreteMean),5e-12);
            check("coupled_center_variance",std::abs(c(path.center(),path.center())-1/(2*mass*omega)),2e-12);
            const double exact=-eta*J*tau-omega*tau/2+g*g/(2*mass*omega*omega*omega)*(omega*tau-1+std::exp(-omega*tau));
            const double exactMean=-eta*g/(mass*omega*omega)*(-std::expm1(-omega*tau/2));
            const double ez=std::abs(logZ-exact),em=std::abs(mean-exactMean);
            if(m>12) {
                check("coupled_log_norm_order_two",std::abs(previousZ/ez-4),.05);
                check("coupled_center_mean_order_two",std::abs(previousMean/em-4),.01);
            }
            std::cout<<"CONVERGENCE eta="<<eta<<" h="<<path.dt()<<" logZ_error="<<ez<<" center_mean_error="<<em<<'\n';
            previousZ=ez; previousMean=em;
        }
    }
}
void complexDirection() {
    const auto c=fockMajoranas();
    const KitaevSSHModel model(HoneycombLattice(2,2),p2Couplings(.4));
    const auto trial=freeGroundState(model);
    const auto base=honeycombGenerator(model);
    MajoranaGenerator generator=[base](const std::vector<double> &x) {
        MatType k=base(x);
        k(0,1)+=DataType(0,-1.4*x[0]); k(1,0)-=DataType(0,-1.4*x[0]);
        k(0,2)+=DataType(0,-.6*x[1]); k(2,0)-=DataType(0,-.6*x[1]);
        return k;
    };
    auto h=[&](const std::vector<double> &x) -> MatType {
        return fockH(c,x,.4)+DataType(0,.7*x[0])*c[0]*c[1]+DataType(0,.3*x[1])*c[0]*c[2];
    };
    std::mt19937_64 rng(1997);
    const auto path=PhononPath(12,6,.11,{1.7,.8},std::vector<double>(12,.1378676470588235)).freeSample(rng);
    const auto measured=evaluatePath(path,trial,pathFactors(path,8,generator),4);
    Eigen::SelfAdjointEigenSolver<MatType> es(h(std::vector<double>(12,0.)));
    const cVecType psi=es.eigenvectors().col(0); cVecType right=psi; Eigen::RowVectorXcd left=psi.adjoint();
    for(int l=1;l<=path.center();++l) right=fockExp(h(path.node(l)),path.dt()/2)*fockExp(h(path.node(l-1)),path.dt()/2)*right;
    for(int l=path.intervals();l>path.center();--l) left=left*fockExp(h(path.node(l)),path.dt()/2)*fockExp(h(path.node(l-1)),path.dt()/2);
    const DataType amp=(left*right)(0),phase=amp/std::abs(amp);
    std::cout<<"REFERENCE complex_fock_phase="<<phase<<" physical_path_phase="<<measured.matter_phase<<'\n';
    require(std::abs(phase.imag())>1e-5,"complex_phase_sensitive");
    check("complex_fock_log_matter",std::abs(measured.log_abs_matter-std::log(std::abs(amp))),2e-12);
    check("complex_fock_phase",std::abs(measured.matter_phase-phase),2e-12);
    MatType green=MatType::Zero(8,8);
    for(int i=0;i<8;++i) for(int j=0;j<8;++j) if(i!=j) green(i,j)=-(left*c[i]*c[j]*right)(0)/amp;
    check("complex_fock_center_green",error(measured.center_green,green),2e-12);
    for(const auto options:{PhononWeightOptions{false,false},PhononWeightOptions{true,false},
                            PhononWeightOptions{false,true},PhononWeightOptions{true,true}}) {
        CoupledPhononPath state(path,trial,generator,4,options);
        const auto &weight=state.state().weight;
        check("optimized_complex_fock_log",std::abs(weight.log_abs_matter-std::log(std::abs(amp))),2e-12);
        check("optimized_complex_fock_phase",std::abs(weight.matter_phase-phase),2e-12);
        check("optimized_complex_fock_green",error(weight.center_green,green),2e-12);
    }
}
}
int main() {
    std::cout<<std::setprecision(17);
    try {
        freePhonons(); coupledMapping(); midpointDrivenReference(); complexDirection(); coupledSingleBond();
        for(const auto &m:metrics) std::cout<<"PASS "<<m.first<<" max_error="<<m.second.error<<" tolerance="<<m.second.tolerance<<'\n';
        return 0;
    } catch(const std::exception &e) { std::cerr<<"FAIL P2: "<<e.what()<<'\n'; return 1; }
}
