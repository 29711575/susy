#include "fermion_observables.h"
#include "vbs_observables.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace kitaev_ssh;
namespace {
using Complex=std::complex<double>;

void require(bool condition,const std::string &name) {
    if(!condition) throw std::runtime_error("sheared_pbc_validate: "+name);
}
void check(const std::string &name,double value,double tolerance) {
    if(!std::isfinite(value) || value>tolerance)
        throw std::runtime_error("sheared_pbc_validate: "+name+" error="+std::to_string(value));
}
double positionNorm(Position p) { return std::hypot(p.x,p.y); }
Position add(Position a,Position b) { return {a.x+b.x,a.y+b.y}; }
Position scale(double x,Position a) { return {x*a.x,x*a.y}; }
Position subtract(Position a,Position b) { return {a.x-b.x,a.y-b.y}; }
double distance(Position a,Position b) { return positionNorm(subtract(a,b)); }
bool sameCell(Cell a,Cell b) { return a.x==b.x && a.y==b.y; }

void geometryChecks(const HoneycombLattice &lattice) {
    const int lx=lattice.Lx(),ly=lattice.Ly(),s=lattice.shear();
    require(lattice.torusMatrix()[0][0]==lx && lattice.torusMatrix()[0][1]==s &&
            lattice.torusMatrix()[1][0]==0 && lattice.torusMatrix()[1][1]==ly,"torus_matrix");
    for(std::int64_t x=-2*lx-3;x<=2*lx+3;++x) for(std::int64_t y=-2*ly-3;y<=2*ly+3;++y) {
        const auto reduced=lattice.reduceCell(x,y);
        require(reduced.representative.x>=0 && reduced.representative.x<lx &&
                reduced.representative.y>=0 && reduced.representative.y<ly,"reduction_representative");
        const std::int64_t reconstructed=static_cast<std::int64_t>(reduced.representative.x)+
            static_cast<std::int64_t>(reduced.winding.x)*lx+
            static_cast<std::int64_t>(reduced.winding.y)*s;
        const std::int64_t reconstructedY=static_cast<std::int64_t>(reduced.representative.y)+
            static_cast<std::int64_t>(reduced.winding.y)*ly;
        require(reconstructed==x && reconstructedY==y,"reduction_identity");
        require(lattice.cellIndexWrapped(x+lx,y)==lattice.cellIndexWrapped(x,y),"x_period");
        require(lattice.cellIndexWrapped(x,y+ly)==lattice.cellIndexWrapped(x-s,y),"y_sheared_period");
    }
    for(int x=0;x<lx;++x) for(int y=0;y<ly;++y) {
        const Cell origin{x,y};
        require(sameCell(lattice.translatedCell(origin,lx,0),origin),"translated_x_period");
        require(sameCell(lattice.translatedCell(origin,0,ly),
                          lattice.translatedCell(origin,-s,0)),"translated_y_period");
        for(const auto d1 : std::vector<Cell>{{-2,1},{1,-2},{lx+1,ly-1}})
            for(const auto d2 : std::vector<Cell>{{2,-1},{-lx,ly+2},{0,0}}) {
                const auto lhs=lattice.translatedCell(lattice.translatedCell(origin,d1.x,d1.y),d2.x,d2.y);
                const auto rhs=lattice.translatedCell(origin,
                    static_cast<std::int64_t>(d1.x)+d2.x,static_cast<std::int64_t>(d1.y)+d2.y);
                require(sameCell(lhs,rhs),"translation_commutativity");
            }
    }

    std::vector<std::array<int,3>> degree(lattice.sites(),{0,0,0});
    std::vector<int> incidence(lattice.bondCount(),0);
    std::set<std::pair<int,int>> edges;
    const auto a1=lattice.primitive(0),a2=lattice.primitive(1);
    for(const auto &bond:lattice.bonds()) {
        const int mu=static_cast<int>(bond.type);
        ++degree[bond.a][mu]; ++degree[bond.b][mu];
        edges.insert({std::min(bond.a,bond.b),std::max(bond.a,bond.b)});
        const auto pa=lattice.position(bond.a),pb=lattice.position(bond.b);
        const auto t1=lattice.periodVector(0),t2=lattice.periodVector(1);
        const auto unwrapped=add(pb,add(scale(bond.boundary_wrap.x,t1),scale(bond.boundary_wrap.y,t2)));
        check("bond_displacement",distance(subtract(unwrapped,pa),bond.displacement),2e-14);
        check("bond_length",std::abs(positionNorm(bond.displacement)-1.),2e-14);
        require(lattice.bondIndex(bond.cell.x,bond.cell.y,bond.type)==bond.index,"bond_index");
    }
    require(static_cast<int>(edges.size())==lattice.bondCount(),"unique_bonds");
    for(const auto &d:degree) for(int value:d) require(value==1,"bond_degree");
    for(int mu=0;mu<3;++mu) require(static_cast<int>(lattice.bonds(static_cast<BondType>(mu)).size())==lattice.cells(),"bond_type_count");

    const auto &vertical=lattice.bond(lattice.bondIndex(0,0,BondType::Y));
    const auto endpoint=lattice.siteCoordinate(vertical.b);
    require(endpoint.sublattice==Sublattice::B && endpoint.x==s && endpoint.y==ly-1,"sheared_corner_y_bond");

    require(static_cast<int>(lattice.plaquettes().size())==lattice.cells(),"plaquette_count");
    for(const auto &plaquette:lattice.plaquettes()) {
        require(plaquette.steps.size()==6 && lattice.pathEnd(plaquette)==plaquette.start,"plaquette_closure");
        std::set<int> vertices;
        Position displacement{0.,0.};
        for(const auto &step:plaquette.steps) {
            const auto &bond=lattice.bond(step.bond);
            ++incidence[step.bond];
            vertices.insert(step.direction==1?bond.a:bond.b);
            displacement=add(displacement,scale(step.direction,bond.displacement));
        }
        require(vertices.size()==6,"plaquette_vertices");
        check("plaquette_displacement",positionNorm(displacement),2e-14);
    }
    for(int value:incidence) require(value==2,"plaquette_incidence");

    for(int wx:{-1,1}) for(int wy:{-1,1}) {
        const auto gauge=lattice.zeroFluxGauge(wx,wy);
        for(int p=0;p<lattice.cells();++p) require(lattice.plaquetteFlux(p,gauge)==1,"zero_flux_all_sectors");
        for(int axis=0;axis<2;++axis) {
            const auto path=lattice.windingPath(axis);
            require(lattice.pathEnd(path)==path.start,"winding_closure");
            require(lattice.wilsonLoop(axis,gauge)==(axis==0?wx:wy),"wilson_sector");
            Position displacement{0.,0.};
            for(const auto &step:path.steps)
                displacement=add(displacement,scale(step.direction,lattice.bond(step.bond).displacement));
            check("winding_displacement",positionNorm(add(displacement,lattice.periodVector(axis))),2e-14);
        }
    }
}

void momentumChecks(const HoneycombLattice &lattice) {
    const auto grid=allowedVBSMomenta(lattice);
    require(static_cast<int>(grid.size())==lattice.cells(),"momentum_count");
    for(std::size_t p=0;p<grid.size();++p) {
        check("momentum_T1_character",std::abs(lattice.Lx()*grid[p].u-std::round(lattice.Lx()*grid[p].u)),2e-13);
        check("momentum_T2_character",std::abs(lattice.shear()*grid[p].u+lattice.Ly()*grid[p].v-
                                                   std::round(lattice.shear()*grid[p].u+lattice.Ly()*grid[p].v)),2e-13);
        for(std::size_t q=0;q<p;++q) {
            const double du=fractionalDelta(grid[p].u,grid[q].u),dv=fractionalDelta(grid[p].v,grid[q].v);
            require(std::abs(du)>1e-10 || std::abs(dv)>1e-10,"momentum_unique");
        }
    }
    for(const auto &p:grid) for(const auto &q:grid) {
        Complex sum=0.;
        for(int r=0;r<lattice.cells();++r) {
            const auto c=lattice.cellCoordinate(r);
            sum+=std::polar(1.,2.*pi*((p.u-q.u)*c.x+(p.v-q.v)*c.y));
        }
        check("character_orthogonality",std::abs(sum-(p.index==q.index?double(lattice.cells()):Complex(0.))),2e-11);
    }
}

std::vector<double> deterministicField(const HoneycombLattice &lattice) {
    std::vector<double> x(3*lattice.cells());
    for(int r=0;r<lattice.cells();++r) {
        const auto c=lattice.cellCoordinate(r);
        for(int mu=0;mu<3;++mu)
            x[3*r+mu]=std::sin(.173*(r+1)*(mu+1))+.031*(c.x+2*c.y+mu);
    }
    return x;
}

void fourierChecks(const HoneycombLattice &lattice) {
    const auto x=deterministicField(lattice);
    const VBSObservables all(lattice);
    const auto staged=all.fourier(x),direct=all.directFourier(x);
    require(staged.size()==direct.size(),"fourier_dimension");
    double error=0.;
    for(std::size_t q=0;q<staged.size();++q) for(int mu=0;mu<3;++mu) error=std::max(error,std::abs(staged[q][mu]-direct[q][mu]));
    check("staged_fft_direct",error,3e-10);
    const auto inverse=all.inverseStructureFactor(staged),realspace=all.directRealSpaceCorrelation(x);
    double inverseError=0.;
    for(int d=0;d<lattice.cells();++d) for(int mu=0;mu<3;++mu) for(int nu=0;nu<3;++nu)
        inverseError=std::max(inverseError,std::abs(inverse[d](mu,nu)-realspace[d](mu,nu)));
    check("inverse_quotient_correlation",inverseError,3e-10);

    const auto target=makeQPoint("target",5./12.,7./12.,lattice);
    const auto offgrid=makeQPoint("offgrid",.137,.271,lattice);
    const VBSObservables selected(lattice,VBSMode::Selected,{target,offgrid});
    const auto values=selected.selectedFourier(x);
    const auto phasesTarget=qPointPhases(lattice,target),phasesOff=qPointPhases(lattice,offgrid);
    for(std::size_t j=0;j<values.size();++j) {
        const auto &phases=j==0?phasesTarget:phasesOff;
        for(int mu=0;mu<3;++mu) {
            Complex expected=0.;
            for(int r=0;r<lattice.cells();++r) expected+=phases[r]*x[3*r+mu];
            check("selected_direct_sum",std::abs(values[j][mu]-expected),3e-10);
        }
    }
    for(const auto &p:allowedVBSMomenta(lattice)) {
        const auto q=makeQPoint("grid",p.u,p.v,lattice);
        if(p.index==0) check("zero_window",std::abs(q.window-Complex(1.)),2e-12);
        else check("nonzero_grid_window",std::abs(q.window),2e-11);
    }
    const auto modGrid=allowedVBSMomenta(lattice);
    const auto mod=modGrid.at(std::min(5,lattice.cells()-1));
    const auto modMinus=minusQPoint(makeQPoint("mod",mod.u,mod.v,lattice),lattice);
    std::vector<double> sinusoid(3*lattice.cells(),0.);
    for(int r=0;r<lattice.cells();++r) {
        const auto c=lattice.cellCoordinate(r);
        sinusoid[3*r]=std::cos(2.*pi*(mod.u*c.x+mod.v*c.y));
    }
    const auto modFourier=all.directFourier(sinusoid);
    check("known_modulation_peak",std::abs(std::abs(modFourier[mod.index][0])-lattice.cells()/2.),3e-10);
    check("known_modulation_conjugate_peak",std::abs(std::abs(modFourier[modMinus.n1*lattice.Ly()+modMinus.n2][0])-lattice.cells()/2.),3e-10);
    for(const auto &p:modGrid) if(p.index!=mod.index && p.index!=modMinus.n1*lattice.Ly()+modMinus.n2)
        check("known_modulation_background",std::abs(modFourier[p.index][0]),3e-10);
}

Complex directWindow(const HoneycombLattice &lattice,double u,double v) {
    Complex sum=0.;
    for(int r=0;r<lattice.cells();++r) {
        const auto c=lattice.cellCoordinate(r);
        sum+=std::polar(1.,-2.*pi*(u*c.x+v*c.y));
    }
    return sum/double(lattice.cells());
}

void targetAndMinusChecks() {
    const HoneycombLattice rectangular(12,3,0),sheared(12,3,3),general(5,3,1);
    const double bareU=0.4195693767448338,bareV=0.5804306232551661;
    const auto bare=makeQPoint("bare",bareU,bareV,sheared);
    require(!bare.on_grid,"bare_target_off_grid");
    require(std::abs(sheared.Lx()*bare.u-std::round(sheared.Lx()*bare.u))>1e-6,
            "bare_target_failed_T1_mismatch");
    require(std::abs(sheared.shear()*bare.u+sheared.Ly()*bare.v-
                     std::round(sheared.shear()*bare.u+sheared.Ly()*bare.v))<2e-12,
            "bare_target_T2_condition");
    check("bare_target_direct_window",std::abs(bare.window-directWindow(sheared,bareU,bareV)),2e-13);
    const auto qRect=makeQPoint("target",5./12.,7./12.,rectangular);
    const auto qShear=makeQPoint("target",5./12.,7./12.,sheared);
    require(!qRect.on_grid && qShear.on_grid,"target_grid_status");
    require(qShear.n1==5 && qShear.n2==0 && qShear.n1*sheared.Ly()+qShear.n2==15,"target_grid_label");
    check("rational_rectangular_direct_window",std::abs(qRect.window-directWindow(rectangular,5./12.,7./12.)),2e-13);
    check("rational_sheared_direct_window",std::abs(qShear.window-directWindow(sheared,5./12.,7./12.)),2e-13);
    check("rational_window_zero",std::abs(qRect.window),2e-11);
    check("rational_window_geometry_independent",std::abs(qRect.window-qShear.window),2e-13);
    const auto targetMinus=minusQPoint(qShear,sheared);
    require(targetMinus.on_grid && targetMinus.n1==7 && targetMinus.n2==0,"target_minus_lookup");
    const auto generic=allowedVBSMomenta(general).at(1*general.Ly());
    const auto genericPoint=makeQPoint("generic",generic.u,generic.v,general);
    const auto genericMinus=minusQPoint(genericPoint,general);
    require(genericMinus.n1==general.Lx()-1 && genericMinus.n2==(general.shear()%general.Ly()),"sheared_minus_not_independent");
}

Complex blochFactor(const HoneycombLattice &,const VBSGridPoint &q,Couplings j) {
    return j.Jz+j.Jx*std::polar(1.,-2.*pi*q.u)+j.Jy*std::polar(1.,-2.*pi*q.v);
}

void spectrumAndA3Checks(const HoneycombLattice &lattice) {
    Couplings j; j.Jx=1.;j.Jy=1.;j.Jz=.5;j.g=.23;
    const KitaevSSHModel model(lattice,j);
    const auto k=model.frameworkGenerator();
    const auto a=model.realAntisymmetricA();
    check("A_real_skew",(a+a.transpose()).norm(),2e-13);
    check("K_hermitian",(k-k.adjoint()).norm(),2e-13);
    check("K_transpose_skew",(k+k.transpose()).norm(),2e-13);
    std::vector<double> x(lattice.bondCount());
    for(int b=0;b<lattice.bondCount();++b) x[b]=.07*std::sin(.31*(b+1));
    const KitaevSSHModel rebuilt(lattice,j,x);
    check("A3_rebuild",(rebuilt.frameworkGenerator()-model.frameworkGeneratorWith(x)).norm(),2e-12);

    Couplings freeJ=j; freeJ.g=0.;
    const KitaevSSHModel freeModel(lattice,freeJ);
    Eigen::SelfAdjointEigenSolver<MatType> solver(freeModel.frameworkGenerator());
    require(solver.info()==Eigen::Success,"free_spectrum_solver");
    std::vector<double> expected;
    for(const auto &q:allowedVBSMomenta(lattice)) {
        const double magnitude=std::abs(blochFactor(lattice,q,freeJ));
        require(magnitude>1e-10,"free_spectrum_no_zero_mode");
        expected.push_back(-2.*magnitude); expected.push_back(2.*magnitude);
    }
    std::sort(expected.begin(),expected.end());
    require(expected.size()==static_cast<std::size_t>(solver.eigenvalues().size()),"free_spectrum_dimension");
    for(int n=0;n<solver.eigenvalues().size();++n) check("free_spectrum_allowed_grid",std::abs(solver.eigenvalues()(n)-expected[n]),3e-11);

    for(int wx:{-1,1}) for(int wy:{-1,1}) {
        const auto gauge=lattice.zeroFluxGauge(wx,wy);
        const KitaevSSHModel sector(lattice,freeJ,{},gauge);
        require(sector.lattice().wilsonLoop(0,gauge)==wx && sector.lattice().wilsonLoop(1,gauge)==wy,"matter_wilson_sector");
        check("sector_generator_hermitian",(sector.frameworkGenerator()-sector.frameworkGenerator().adjoint()).norm(),2e-13);
    }
}

void fermionChecks(const HoneycombLattice &lattice) {
    Couplings j; j.Jx=1.;j.Jy=1.;j.Jz=.5;j.g=.2;
    const KitaevSSHModel model(lattice,j);
    const FermionObservables observables(model);
    for(const auto &channel:observables.channels()) for(int x=0;x<lattice.Lx();++x) for(int y=0;y<lattice.Ly();++y) {
        const auto path=observables.pathFor(x,y,channel);
        require(lattice.pathEnd(path)==lattice.siteIndex(lattice.translatedCell({x,y},channel.dx,channel.dy).x,
                                                           lattice.translatedCell({x,y},channel.dx,channel.dy).y,
                                                           static_cast<Sublattice>(channel.beta)),"fermion_endpoint");
    }
    MatType green=MatType::Zero(lattice.sites(),lattice.sites());
    for(int i=0;i<green.rows();++i) for(int j2=0;j2<green.cols();++j2)
        green(i,j2)=Complex(.003*(i+1)*(j2+2),.002*(i-j2));
    std::vector<int> siteGauge(lattice.sites());
    for(int i=0;i<lattice.sites();++i) siteGauge[i]=(i%3==0)?-1:1;
    const auto transformed=model.gaugeTransformed(siteGauge);
    const FermionObservables transformedObservables(transformed);
    MatType transformedGreen=green;
    for(int i=0;i<green.rows();++i) for(int j2=0;j2<green.cols();++j2)
        transformedGreen(i,j2)=double(siteGauge[i]*siteGauge[j2])*green(i,j2);
    const auto before=observables.measure(green),after=transformedObservables.measure(transformedGreen);
    require(before.size()==after.size(),"fermion_measure_dimension");
    for(std::size_t i=0;i<before.size();++i) if(observables.channels()[i].string)
        check("fermion_gauge_string_invariance",std::abs(before[i]-after[i]),2e-12);
}

} // namespace

int main() {
    try {
        targetAndMinusChecks();
        for(const auto spec : std::vector<std::array<int,3>>{{5,3,1},{5,3,4},{12,3,3}}) {
            const HoneycombLattice lattice(spec[0],spec[1],spec[2]);
            geometryChecks(lattice);
            momentumChecks(lattice);
            fourierChecks(lattice);
            spectrumAndA3Checks(lattice);
            fermionChecks(lattice);
        }
        std::cout<<"sheared_pbc_validate PASS (geometry, gauge, spectrum, FFT, correlations and fermion strings).\n";
        return 0;
    } catch(const std::exception &e) {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
