#include "phonon_mc.h"
#include <iomanip>
#include <iostream>
using namespace kitaev_ssh;
namespace {
void require(bool b,const char *s){if(!b)throw std::runtime_error(s);}
double maxLog=0.,maxRatio=0.,maxPhase=0.,maxGreen=0.;int cases=0,proposals=0;
PureProjectorSlice fromReal(const Eigen::MatrixXd &b) {
    const int n=b.rows()/2;MatType value(b.rows(),b.cols());
    for(int j=0;j<n;++j)for(int i=0;i<n;++i) {
        value(i,j)=b(i,j);value(n+i,n+j)=b(n+i,n+j);
        value(i,n+j)=DataType(0.,-b(i,n+j));value(n+i,j)=DataType(0.,b(n+i,j));
    }
    return PureProjectorSlice::fromValidatedDense(value,1.,"real-validation");
}
void compare(const GaussianTrialState &trial,const std::vector<PureProjectorSlice> &slices,
             bool expectReal) {
    PureFastOptions a,b;a.real_magnitude=true;b.accelerated_linear_algebra=false;
    PureStableWeightWorkspace wa,wb;
    const auto fast=pureProjectorMagnitudeWeight(trial,slices,a,&wa);
    const auto oracle=pureProjectorMagnitudeWeight(trial,slices,b,&wb);
    require(fast.status==oracle.status,"real/reference status");
    require(fast.first_failing_slice==oracle.first_failing_slice,"real/reference failing slice");
    require(fast.overlap_rank==oracle.overlap_rank,"real/reference rank");
    if(fast.ok())maxLog=std::max(maxLog,std::abs(fast.log_abs_weight-oracle.log_abs_weight));
    if(expectReal)require(wa.real_magnitude.calls==1&&wa.real_magnitude.fallbacks==0,"real dispatch missing");
    else require(wa.real_magnitude.fallbacks==1,"complex fallback missing");
    ++cases;
}
void replay(int L,int shear) {
    KitaevSSHModel model(HoneycombLattice(L,L,shear),p2Couplings(.8));
    std::mt19937_64 rng(10024000+L+shear);
    auto path=PhononPath(model.lattice().bondCount(),8,.1,{100.,.1}).freeSample(rng);
    const auto trial=freeGroundState(model);require(bool(trial.real_magnitude_trial),"real trial certificate");
    const auto generator=honeycombGenerator(model);
    auto prepared=prepareOwnedPath(path,trial.Phi.rows(),generator,false);
    compare(trial,prepared.weightOrder,L>=4);
    PureStableWeightWorkspace persistent;PureFastOptions opts;opts.real_magnitude=true;
    pureProjectorMagnitudeWeight(trial,prepared.weightOrder,opts,&persistent);
    pureProjectorMagnitudeWeight(trial,prepared.weightOrder,opts,&persistent);
    if(L>=4)require(persistent.real_magnitude.cache_hits>=prepared.weightOrder.size(),"immutable real cache reuse");
    auto generic=GaussianTrialState::fromPhi(trial.Phi);compare(generic,prepared.weightOrder,false);
    auto mutated=trial;mutated.Phi.col(0)*=DataType(0.,1.);compare(mutated,prepared.weightOrder,false);
    auto changed=prepared.weightOrder;MatType extra=changed[0].matrix;
    extra(0,1)+=DataType(0.,1e-30);extra(1,0)-=DataType(0.,1e-30);
    changed[0]=PureProjectorSlice::fromValidatedDense(extra,1.,"tiny-extra-coupling");compare(trial,changed,false);
    if(L==4) {
        compare(trial,{fromReal(Eigen::MatrixXd::Zero(2*L*L,2*L*L))},false);
        compare(trial,{fromReal(1e-20*Eigen::MatrixXd::Identity(2*L*L,2*L*L))},false);
        const auto &q=trial.real_magnitude_trial->basis;
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(q);
        Eigen::MatrixXd full=qr.householderQ()*Eigen::MatrixXd::Identity(q.rows(),q.rows());
        compare(trial,{fromReal((full.rightCols(q.cols())*q.transpose()).eval())},false);
    }
    CoupledPhononPath fast(path,trial,generator,4,{true,true,true});
    MajoranaGenerator genericGenerator([model](const std::vector<double>&x){return model.frameworkGeneratorWith(x);});
    CoupledPhononPath reference(path,trial,genericGenerator,4,{true,true,true});
    PhononSampler sa(fast,rng,PhononProposal::BondRefresh),sb(reference,rng,PhononProposal::BondRefresh);
    for(int p=0;p<4;++p) {
        auto drawn=PhononPath(1,path.intervals(),path.dt(),path.parameters()).freeSample(rng);
        std::vector<double> values(path.nodes());for(int j=0;j<path.nodes();++j)values[j]=drawn.coordinate(j,0);
        const auto ma=sa.attemptBond((p*13+1)%path.bonds(),values,p%2?.91:.13);
        const auto mb=sb.attemptBond((p*13+1)%path.bonds(),values,p%2?.91:.13);
        require(ma.accepted==mb.accepted,"real/reference Metropolis decision");
        maxRatio=std::max(maxRatio,std::abs(ma.log_ratio-mb.log_ratio));
        fast.ensureMeasured();reference.ensureMeasured();
        const auto &a=fast.state().weight,&b=reference.state().weight;
        require(fast.configurationHash()==reference.configurationHash(),"real/reference path");
        maxLog=std::max(maxLog,std::abs(a.log_abs_matter-b.log_abs_matter));
        maxPhase=std::max(maxPhase,std::abs(a.matter_phase-b.matter_phase));
        maxGreen=std::max(maxGreen,(a.center_green-b.center_green).cwiseAbs().maxCoeff());
        ++proposals;
    }
}
}
int main(){try{
    for(int L:{2,4,6,9,12})for(int shear:{0,1})replay(L,shear);
    require(maxLog<2e-9&&maxRatio<2e-9&&maxPhase<2e-10&&maxGreen<2e-9,"real propagation numerical mismatch");
    std::cout<<std::setprecision(17)<<"{\"cases\":"<<cases<<",\"proposals\":"<<proposals
        <<",\"max_log_error\":"<<maxLog<<",\"max_ratio_error\":"<<maxRatio
        <<",\"max_phase_error\":"<<maxPhase<<",\"max_green_error\":"<<maxGreen
        <<"}\nREAL_MAGNITUDE_VALIDATE_OK\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
