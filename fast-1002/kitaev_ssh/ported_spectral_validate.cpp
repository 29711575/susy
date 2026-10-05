#include "phonon_mc.h"
#include <array>
#include <iomanip>
#include <iostream>

using namespace kitaev_ssh;

namespace {
double maxExponential=0.,maxCayley=0.,maxLogWeight=0.,maxPhase=0.,maxGreen=0.,maxRatio=0.;
int matrixCases=0,replayProposals=0;

void require(bool condition,const char *message) {
    if(!condition) throw std::runtime_error(message);
}

void matrixCheck(MatType k) {
    BipartiteSpectrum s;
    require(s.compute(k),"bipartite decomposition failed");
    const auto oracle=selfAdjointSpectral(k);
    require(oracle.ok,"reference eigendecomposition failed");
    for(double a:{.05,.1}) {
        const cVecType e=(oracle.D.array()*(-a)).exp().cast<DataType>();
        const cVecType t=(oracle.D.array()*(a/2.)).tanh().cast<DataType>();
        const MatType b=spectralRecombination(oracle,e),c=spectralRecombination(oracle,t);
        MatType reducedB,reducedC;s.exponential(a,reducedB);s.cayley(a,reducedC);
        maxExponential=std::max(maxExponential,(reducedB-b).norm()/std::max(1.,b.norm()));
        maxCayley=std::max(maxCayley,(reducedC-c).norm()/std::max(1.,c.norm()));
    }
    ++matrixCases;
}

void compareMeasured(CoupledPhononPath &a,CoupledPhononPath &b) {
    a.ensureMeasured(); b.ensureMeasured();
    require(a.configurationHash()==b.configurationHash(),"configuration identities differ");
    const auto &wa=a.state().weight,&wb=b.state().weight;
    maxLogWeight=std::max(maxLogWeight,std::abs(wa.log_abs_matter-wb.log_abs_matter));
    maxPhase=std::max(maxPhase,std::abs(wa.matter_phase-wb.matter_phase));
    maxGreen=std::max(maxGreen,(wa.center_green-wb.center_green).cwiseAbs().maxCoeff());
}

void replay(int L,int intervals,double lambda) {
    const KitaevSSHModel model(HoneycombLattice(L,L),p2Couplings(std::sqrt(lambda)));
    std::mt19937_64 rng(93000+L+intervals);
    const auto path=PhononPath(model.lattice().bondCount(),intervals,.1,{100.,.1}).freeSample(rng);
    const auto trial=freeGroundState(model);
    const auto fast=honeycombGenerator(model);
    const MajoranaGenerator reference(
        [model](const std::vector<double> &x){return model.frameworkGeneratorWith(x);},
        [model](const std::vector<double> &x,MatType &out){model.frameworkGeneratorWith(x,out);},false);
    const PhononWeightOptions options{true,true,true};
    CoupledPhononPath a(path,trial,fast,4,options),b(path,trial,reference,4,options);
    PhononSampler sa(a,rng,PhononProposal::BondRefresh),sb(b,rng,PhononProposal::BondRefresh);
    compareMeasured(a,b);
    for(int n=0;n<4;++n) {
        const int bond=(17*n+3)%path.bonds();
        const auto drawn=PhononPath(1,intervals,.1,{100.,.1}).freeSample(rng);
        std::vector<double> values(path.nodes());
        for(int j=0;j<path.nodes();++j) values[j]=drawn.coordinate(j,0);
        const double u=(n%2==0)?.1:.9;
        const auto ma=sa.attemptBond(bond,values,u),mb=sb.attemptBond(bond,values,u);
        require(ma.accepted==mb.accepted,"accept/reject decisions differ");
        maxRatio=std::max(maxRatio,std::abs(ma.log_ratio-mb.log_ratio));
        compareMeasured(a,b); ++replayProposals;
    }
    // A phase-only consumer must keep the same signed weight and expose an
    // invalid Green flag, rather than manufacturing a Green or unit phase.
    PhononWeightOptions phaseOnly=options; phaseOnly.compute_measurement_green=false;
    CoupledPhononPath noGreen(a.state().path,trial,fast,4,phaseOnly);
    noGreen.ensureMeasured();
    require(noGreen.state().weight.phase_valid&&!noGreen.state().weight.green_valid,
            "phase-only validity flags");
    require(noGreen.state().weight.center_green.size()==0,"phase-only Green populated");
    maxPhase=std::max(maxPhase,std::abs(noGreen.state().weight.matter_phase-a.state().weight.matter_phase));
    const auto canonical=noGreen.uncachedReference(noGreen.state().path,true,false);
    require(canonical.weight.phase_valid&&!canonical.weight.green_valid&&canonical.weight.center_green.size()==0,"canonical phase-only validity");
    maxPhase=std::max(maxPhase,std::abs(canonical.weight.matter_phase-a.state().weight.matter_phase));
}
}

int main(int argc,char **) {
    try {
        std::mt19937_64 rng(93001); std::normal_distribution<double> normal;
        for(int L:{2,4,6,9,12}) for(int shear:{0,1}) {
            const HoneycombLattice lattice(L,L,shear);
            std::vector<double> x(lattice.bondCount());
            for(double &v:x) v=.4*normal(rng);
            const KitaevSSHModel model(lattice,p2Couplings(std::sqrt(2.)),x);
            matrixCheck(model.frameworkGenerator());
            std::vector<int> gauge(model.majoranas());
            for(int i=0;i<model.majoranas();++i) gauge[i]=(i%3==0)?-1:1;
            matrixCheck(model.gaugeTransformed(gauge).frameworkGenerator());
        }
        // Zero, degenerate and rank-deficient hopping blocks.
        matrixCheck(MatType::Zero(8,8));
        MatType degenerate=MatType::Zero(8,8);
        for(int i=0;i<3;++i){degenerate(i,4+i)=DataType(0.,-2.);degenerate(4+i,i)=DataType(0.,2.);}
        matrixCheck(degenerate);
        BipartiteSpectrum rejected;
        MatType nonbipartite=degenerate;
        nonbipartite(0,1)=DataType(0.,-1e-16);nonbipartite(1,0)=DataType(0.,1e-16);
        require(!rejected.compute(nonbipartite),"small non-bipartite term was discarded");
        replay(4,12,1.3);
        if(argc>1) {replay(9,48,1.6); replay(12,100,2.);}
        require(maxExponential<2e-12&&maxCayley<2e-12,"matrix-function mismatch");
        require(maxLogWeight<2e-9&&maxRatio<2e-9,"weight or acceptance-ratio mismatch");
        require(maxPhase<2e-10&&maxGreen<2e-9,"signed phase or center-Green mismatch");
        std::cout<<std::setprecision(17)<<"{\"matrix_cases\":"<<matrixCases
            <<",\"replay_proposals\":"<<replayProposals<<",\"max_exponential_relative_error\":"<<maxExponential
            <<",\"max_cayley_relative_error\":"<<maxCayley<<",\"max_log_weight_error\":"<<maxLogWeight
            <<",\"max_log_ratio_error\":"<<maxRatio<<",\"max_phase_error\":"<<maxPhase
            <<",\"max_green_entry_error\":"<<maxGreen<<"}\nBIPARTITE_SPECTRAL_VALIDATE_OK\n";
    } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}
