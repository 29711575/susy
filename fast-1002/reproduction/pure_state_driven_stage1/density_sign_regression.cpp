#include "stage1_common.h"
#include "pure_projector_failure_snapshot.h"
using namespace driven_stage1;

static PureFastRatioResult replay(const GaussianTrialState &trial,
        const PureFastConfiguration &cfg,const PureFastProposal &proposal,int block) {
    PureFastOptions options;options.stabilized_local_updates=true;
    // No channel-specific positivity contract: the signed algebra must work.
    require(!options.positive_hs1_ratio_contract,"unexpected positive contract");
    PureProjectorOptions stackOptions;stackOptions.accelerated_linear_algebra=true;
    stackOptions.green_residual_tolerance=options.local_green_rebuild_tolerance;
    PureProjectorStackManager before(trial,cfg.slices,block,stackOptions,
        proposal.index,true,false,true,true,true);
    PureProjectorStackManager after(trial,cfg.slices,block,stackOptions,
        proposal.index+1,true,false,true,true,true);
    require(before.ok()&&after.ok(),"snapshot cut reconstruction failed");
    const auto result=pureFastLocalRatioFromGreen(trial,cfg.slices,proposal.index,
        proposal.new_factor,proposal.new_eta,before.phiRight(),after.greenView(),options);
    require(result.ok(),"snapshot ratio untrusted");
    return result;
}

int main(int argc,char**argv){try{
    require(argc==2,"density_sign_regression FIXTURE_DIRECTORY");
    double maximumError=0;
    for(auto fixture:std::vector<std::pair<int,double>>{
            {901,0.7699066766710073},{902,0.7453937320445516}}){
        const auto saved=pure_failure_snapshot::read(std::string(argv[1])+"/seed"+
            std::to_string(fixture.first)+".bin");
        require(saved.observed_ratio.real()<0,"fixture no longer contains original sign failure");
        ModelParameters p;p.L=4;p.hs=1;p.dt=.05;
        const auto trial=makeInitialState(p,0,0,-1);
        require((trial.Phi-saved.phi).norm()<1e-12,"snapshot trial mismatch");
        const auto forward=replay(trial,saved.configuration,saved.proposal,saved.block);
        // Independent Fock propagation of the exact frozen factors gives the
        // signed numbers above; both old and candidate weights are positive.
        require(std::abs(forward.ratio-fixture.second)<1e-10,"signed Fock ratio mismatch");
        maximumError=std::max(maximumError,std::abs(forward.ratio-fixture.second));
        auto cfg=saved.configuration;const int index=saved.proposal.index;
        auto reverse=saved.proposal;
        reverse.new_factor=cfg.slices[index].matrix;reverse.new_hs=cfg.hs_fields[index];
        reverse.new_eta=cfg.slices[index].eta;
        cfg.slices[index].matrix=saved.proposal.new_factor;
        cfg.slices[index].eta=saved.proposal.new_eta;cfg.hs_fields[index]=saved.proposal.new_hs;
        const auto backward=replay(trial,cfg,reverse,saved.block);
        require(std::abs(forward.ratio*backward.ratio-DataType(1))<1e-10,"signed reciprocity mismatch");
        std::cout<<std::setprecision(17)<<"seed="<<fixture.first<<" Fock="<<fixture.second
            <<" forward="<<forward.ratio<<" reverse="<<backward.ratio<<'\n';
    }
    // Symmetrizing arbitrary bad data is not a recovery: reject when it
    // changes the determinant that fixes the overlap magnitude.
    MatType c=MatType::Zero(2,2),g=MatType::Zero(2,2);
    c(0,1)=.3;c(1,0)=-.3;g(0,1)=.2;g(1,0)=-.2;
    require(pureSkewProductSign(c,g).ok(),"valid skew contraction rejected");
    g(0,0)=2.;g(1,1)=3.;
    require(!pureSkewProductSign(c,g).ok(),"determinant-changing symmetric corruption accepted");
    std::cout<<"PASS frozen proposals=2 reciprocal_flips=2 max_Fock_error="<<maximumError<<'\n';
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
