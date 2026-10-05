#include "stage1_dense.h"
#include "pure_projector_failure_snapshot.h"
#include <iostream>
using namespace driven_stage1;

DataType fockWeight(const ModelParameters &p,const PureImaginaryTimeProtocol &protocol,
                    const cVecType &trial,const PureFastConfiguration &cfg) {
    const auto gamma=gammaMatrices(p.L);
    cVecType state=trial;
    for(int k=int(cfg.slices.size())-1;k>=0;--k)
        state=cfg.slices[k].eta*denseFactorAtFlattenedLocation(p,protocol,gamma,
            cfg.locations[k],cfg.hs_fields[k])*state;
    return trial.dot(state);
}

PureFastRatioResult localRatio(const GaussianTrialState &trial,const PureFastConfiguration &cfg,
                               const PureFastProposal &proposal) {
    PureProjectorStackManager stack(trial,cfg.slices,4,PureProjectorOptions(),proposal.index);
    require(stack.ok(),"local stack failed");
    const MatType right=stack.phiRight();
    require(stack.moveToCut(proposal.index+1),"local cut failed");
    const auto green=stack.green();
    return pureFastLocalRatioAtCut(trial,cfg.slices,proposal.index,proposal.new_factor,
                                  proposal.new_eta,right,stack.phiLeft(),PureFastOptions(),&green);
}

void close(DataType a,DataType b,const char *message){
    require(std::abs(a-b)<1e-9*std::max(1.,std::abs(b)),message);
}
ModelParameters model(int L,int hs){ModelParameters p;p.L=L;p.hs=hs;p.dt=.05;p.analytic_hs_factors=true;return p;}

void negative(){
    auto p=model(6,0);const auto protocol=PureImaginaryTimeProtocol::constant(6,.2,p.dt);
    const auto trial=makeInitialState(p,0,0,-1);
    const auto dense=denseInitialState(p,gammaMatrices(p.L),trial,0,-1);
    for(auto fixture:std::vector<std::pair<int,int>>{{3,4},{1,50}}){
        std::mt19937_64 rng(fixture.first);auto cfg=makeDrivenContour(p,protocol,rng).configuration;
        const auto initialHash=pureFastConfigurationHash(cfg);
        PureFastOptions options;options.stabilized_local_updates=true;
        PureProjectorFastWalker walker(trial,cfg,4,PureFastRunMode::FastStrict,options);
        const auto originalGreen=walker.measurementGreen().green;
        const double originalLog=walker.currentWeight().log_abs_weight;const int originalSign=walker.z2Sign();
        DataType forward;
        for(int direction=0;direction<2;++direction){
            const auto proposal=flipProposal(p,protocol,cfg,fixture.second,0.);
            auto candidate=cfg;candidate.slices[proposal.index].matrix=proposal.new_factor;
            candidate.hs_fields[proposal.index]=proposal.new_hs;
            const auto direct=fockWeight(p,protocol,dense.vector,candidate)/fockWeight(p,protocol,dense.vector,cfg);
            const auto local=localRatio(trial,cfg,proposal);
            const auto center=pureProjectorMpEndpointSameProposal(trial,cfg.slices,proposal.index,proposal.new_factor,1.);
            require(local.ok()&&center.ok(),"HS0 ratio reference untrusted: "+center.message);
            require(direct.real()<0&&local.ratio.real()<0&&center.ratio_z2==-1,"negative HS0 sign lost");
            close(local.ratio,direct,"HS0 local/Fock mismatch");close(center.ratio,direct,"HS0 center/Fock mismatch");
            if(direction==0)forward=center.ratio;else close(forward*center.ratio,1.,"flip reciprocity failed");
            const auto accepted=walker.propose(proposal);require(accepted.accepted&&!accepted.terminated,"HS0 forced flip failed");
            std::cout<<"HS0 seed="<<fixture.first<<" index="<<fixture.second<<" direction="<<direction
                <<" Fock="<<direct<<" local="<<local.ratio<<" center="<<center.ratio<<'\n';
            cfg=std::move(candidate);
        }
        require(walker.configurationHash()==initialHash&&walker.z2Sign()==originalSign,"flip-back hash/sign not restored");
        require(std::abs(walker.currentWeight().log_abs_weight-originalLog)<1e-9,"flip-back weight not restored");
        require(relativeError(walker.measurementGreen().green,originalGreen)<1e-9,"flip-back Green not restored");
    }
}

void positive(){
    int configs=0,flips=0,centers=0;double maxError=0;
    for(int L:{4,6})for(int seed=1;seed<=4;++seed){
        auto p=model(L,1);const auto protocol=PureImaginaryTimeProtocol::constant(6,.2,p.dt);
        const auto trial=makeInitialState(p,0,0,-1);require(positiveHs1RatioContract(p,trial),"positive proof not enabled");
        const auto dense=denseInitialState(p,gammaMatrices(L),trial,0,-1);std::mt19937_64 rng(seed);
        const auto cfg=makeDrivenContour(p,protocol,rng).configuration;
        const auto weight=fockWeight(p,protocol,dense.vector,cfg);
        const auto ref=pureProjectorStableReferenceWeight(trial,cfg.slices,PureFastOptions());
        require(weight.real()>0&&ref.ok()&&ref.z2_sign==1,"HS1 fixed independent weight not positive");++configs;
        const auto indices=proposalIndices(cfg);
        for(int draw=0;draw<4;++draw){
            const int half=int(indices.size()/2);
            const int slot=draw<2?draw*half+std::uniform_int_distribution<int>(0,half-1)(rng):
                std::uniform_int_distribution<int>(0,int(indices.size())-1)(rng);
            auto proposal=flipProposal(p,protocol,cfg,indices[slot],.5);auto candidate=cfg;
            candidate.slices[proposal.index].matrix=proposal.new_factor;candidate.hs_fields[proposal.index]=proposal.new_hs;
            const auto candidateWeight=fockWeight(p,protocol,dense.vector,candidate);
            require(candidateWeight.real()>0,"HS1 candidate weight not positive");
            const auto direct=candidateWeight/weight;const auto local=localRatio(trial,cfg,proposal);
            require(local.ok()&&local.ratio.real()>0,"HS1 local ratio not positive");
            require(std::abs(std::log(direct.real()))<=proposal.positive_ratio_log_bound+1e-10,"Fock HS1 ratio outside bound");
            require(std::abs(std::log(local.ratio.real()))<=proposal.positive_ratio_log_bound+1e-10,"local HS1 ratio outside bound");
            close(local.ratio,direct,"HS1 local/Fock mismatch");maxError=std::max(maxError,std::abs(local.ratio-direct));++flips;
            if(seed==1&&draw==0){
                const auto center=pureProjectorMpEndpointSameProposal(trial,cfg.slices,proposal.index,proposal.new_factor,1.);
                require(center.ok()&&center.ratio_z2==1,"HS1 center untrusted");
                close(center.ratio,direct,"HS1 center/Fock mismatch");++centers;
            }
        }
    }
    std::cout<<"HS1 configurations="<<configs<<" positive_candidates="<<flips<<" bounded_ratios="<<flips
        <<" centers="<<centers<<" max_ratio_error="<<maxError<<'\n';
}

void eta(){
    auto p=model(4,1);const auto trial=makeInitialState(p,0,0,-1);const MatType factor=localHsFactor(p,6,0,0,1);
    const auto unit=pureProjectorStableReferenceWeight(trial,{{factor,1.,"eta"}},PureFastOptions());
    for(double eta:{2.,.5,-2.}){
        const auto weight=pureProjectorStableReferenceWeight(trial,{{factor,eta,"eta"}},PureFastOptions());
        require(unit.ok()&&weight.ok(),"eta reference failed");
        const DataType ratio=std::exp(weight.log_abs_weight-unit.log_abs_weight)*weight.complex_phase/unit.complex_phase;
        close(ratio,eta,"eta magnitude counted more than once");std::cout<<"eta="<<eta<<" ratio="<<ratio<<'\n';
    }
}

void unknown(){
    auto p=model(4,1);const auto trial=makeInitialState(p,0,0,-1);
    PureFastConfiguration cfg;cfg.slices.emplace_back(MatType::Identity(8,8),DataType(0,1),"unknown");
    cfg.hs_fields.push_back(0);cfg.locations.push_back({});
    const auto weight=pureProjectorStableReferenceWeight(trial,cfg.slices,PureFastOptions());
    require(weight.ok()&&weight.z2_sign==0,"unknown fixture does not produce unresolved reference");
    bool rejected=false;try{PureProjectorFastWalker w(trial,cfg,4,PureFastRunMode::FastStrict);}
    catch(const std::runtime_error &){rejected=true;}require(rejected,"unknown RealZ2 silently initialized");
    PureFastOptions complex;complex.weight_mode=PureProjectorWeightMode::GenericComplex;
    PureProjectorFastWalker generic(trial,cfg,4,PureFastRunMode::FastStrict,complex);
    require(generic.z2Sign()==0,"complex-mode unknown sign converted to +1");
    std::mt19937_64 rng(3);const auto protocol=PureImaginaryTimeProtocol::constant(6,.2,p.dt);
    auto independent=makeDrivenContour(p,protocol,rng).configuration;
    rejected=false;try{PureProjectorFastWalker bad(trial,independent,4,PureFastRunMode::FastStrict,
        PureFastOptions(),PureFastInitializationPolicy::MirroredTheoremZ2Plus);}
    catch(const std::runtime_error &){rejected=true;}require(rejected,"unproved mirror accepted");
    auto mirrored=makeMirroredDrivenContour(p,protocol,rng).configuration;
    PureProjectorFastWalker good(trial,mirrored,4,PureFastRunMode::FastStrict,PureFastOptions(),
        PureFastInitializationPolicy::MirroredTheoremZ2Plus);require(good.z2Sign()==1,"proved mirror not +1");
    // The endpoint is identity, but the boundary contains an unresolved Cayley
    // chart. Endpoint determinants must not silently manufacture its spin lift.
    MatType turn=MatType::Identity(8,8);turn(0,0)=-1;turn(1,1)=-1;
    const MatType identity=MatType::Identity(8,8);
    std::vector<PureProjectorSlice> singularChart={{turn,1.,"turn"},{turn,1.,"undo"},
        {identity,1.,"identity"},{identity,1.,"identity"}};
    const auto unresolved=pureProjectorMpEndpointSameProposal(trial,singularChart,0,turn,2.);
    require(!unresolved.ok()&&!unresolved.ratio_converged&&unresolved.ratio_z2==0,"unresolved spin lift certified");
    std::cout<<"unknown RealZ2 rejected; GenericComplex sign=0; unproved mirror rejected; proved mirror sign=+1; unresolved spin lift rejected\n";
}

void anomalies(){
    auto p=model(4,1);const auto trial=makeInitialState(p,0,0,-1);
    const auto protocol=PureImaginaryTimeProtocol::constant(6,.1,p.dt);std::mt19937_64 rng(3);
    const auto cfg=makeMirroredDrivenContour(p,protocol,rng).configuration;
    for(bool underflow:{true,false}){
        PureFastOptions options;options.stabilized_local_updates=true;options.positive_hs1_ratio_contract=!underflow;
        PureProjectorFastWalker walker(trial,cfg,4,PureFastRunMode::FastStrict,options,
            PureFastInitializationPolicy::MirroredTheoremZ2Plus);
        auto proposal=flipProposal(p,protocol,cfg,proposalIndices(cfg).front(),.5);
        if(underflow)proposal.new_eta=1e-200; // Nonzero scalar whose double square underflows.
        else proposal.positive_ratio_log_bound=0; // Inject an inconsistent bound to exercise recovery.
        const auto result=walker.propose(proposal);
        require(result.terminated&&!result.accepted&&!result.ratio.ok(),"ratio anomaly became normal rejection");
        require(walker.diagnostics().local_double_retries==1,"ratio anomaly must retry exactly once");
        require(walker.configurationHash()==pureFastConfigurationHash(cfg),"failure changed configuration");
        const std::string path=std::string("hs1_review_")+(underflow?"zero_":"bound_")+
            std::to_string(std::chrono::system_clock::now().time_since_epoch().count())+".bin";
        pure_failure_snapshot::write(path,trial,walker,proposal,result,4,rng);
        const auto saved=pure_failure_snapshot::read(path);
        require(saved.hash==walker.configurationHash()&&saved.proposal.index==proposal.index&&
            saved.proposal.uniform==proposal.uniform&&saved.proposal.new_eta==proposal.new_eta&&
            saved.proposal.positive_ratio_log_bound==proposal.positive_ratio_log_bound,"failure proposal roundtrip failed");
        require((saved.phi-trial.Phi).norm()==0&&(saved.proposal.new_factor-proposal.new_factor).norm()==0,"failure matrices changed");
        std::ostringstream rngText;rngText<<rng;require(saved.rng==rngText.str(),"failure RNG changed");
        std::cout<<(underflow?"zero":"bound")<<" retried once and failed closed; snapshot="<<path<<'\n';
    }
}

int main(int argc,char **argv){try{
    std::cout<<std::setprecision(17);const std::string which=argc>1?argv[1]:"all";
    if(which=="all"||which=="negative")negative();
    if(which=="all"||which=="positive")positive();
    if(which=="all"||which=="eta")eta();
    if(which=="all"||which=="unknown")unknown();
    if(which=="all"||which=="anomalies")anomalies();
    std::cout<<"PASS "<<which<<'\n';
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
