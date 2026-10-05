#include "phonon_mc.h"
#include <iostream>
using namespace kitaev_ssh;
int main() {
    try {
        double maxLog=0,maxPhase=0,maxGreen=0,maxRatio=0,maxAcceptance=0;
        int accepts=0,rejects=0;
        const KitaevSSHModel model(HoneycombLattice(4,4),p2Couplings(.3));
        std::mt19937_64 rng(7901);
        auto path=PhononPath(48,64,.1,{1.7,.8},std::vector<double>(48,.1378676470588235)).freeSample(rng);
        auto trial=freeGroundState(model);auto generator=honeycombGenerator(model);
        auto compare=[&](const PathWeight &a,const PathWeight &b) {
            maxLog=std::max(maxLog,std::abs(a.log_abs_matter-b.log_abs_matter));
            maxPhase=std::max(maxPhase,std::abs(a.matter_phase-b.matter_phase));
            maxGreen=std::max(maxGreen,(a.center_green-b.center_green).cwiseAbs().maxCoeff());
        };
        auto require=[&](bool condition,const char *message) {
            if(!condition)throw std::runtime_error(message);
        };
        // Eager bond proposals certify their signed weight but intentionally
        // defer the center Green. Rebuild that measurement before comparing it;
        // never reduce an empty matrix or bypass the original Green tolerance.
        auto measuredWeight=[&](const CoupledPhononPath::Snapshot &snapshot,
                                const GaussianTrialState &stateTrial) {
            auto weight=snapshot.weight;
            require(weight.magnitude_valid&&weight.phase_valid,"uncertified signed weight");
            if(!weight.green_valid) {
                const auto center=pureProjectorEndpointRebuild(stateTrial,snapshot.factors,
                    snapshot.path.intervals(),4,PureProjectorOptions{});
                require(center.ok(),"center Green certification failed");
                require(std::abs(weight.log_abs_matter-center.log_abs_weight)<=
                    1e-9*std::max(1.,std::abs(weight.log_abs_matter)),"center magnitude mismatch");
                weight.center_green=center.green.adjoint().eval();
                weight.green_valid=true;
            }
            require(weight.center_green.rows()==stateTrial.Phi.rows()&&
                    weight.center_green.cols()==stateTrial.Phi.rows(),"center Green dimensions");
            return weight;
        };
        for(auto options:{PhononWeightOptions{true,false},PhononWeightOptions{false,true},
                          PhononWeightOptions{true,true},PhononWeightOptions{true,true,false,true}}) {
            CoupledPhononPath ref(path,trial,generator,4),candidate(path,trial,generator,4,options);
            compare(measuredWeight(ref.state(),trial),measuredWeight(candidate.state(),trial));
            for(int n=0;n<24;++n) {
                const int bond=n%48,node=n%3==0?32:(n%3==1?0:64);
                auto before=ref.state();
                std::vector<double> values(65);
                const auto drawn=PhononPath(1,64,.1,{1.7,.8},{.1378676470588235}).freeSample(rng);
                for(int j=0;j<65;++j)values[j]=drawn.coordinate(j,0);
                auto a=n%2?ref.proposeBondPath(bond,values):ref.propose(node,bond,values[node]);
                auto b=n%2?candidate.proposeBondPath(bond,values):candidate.propose(node,bond,values[node]);
                compare(measuredWeight(a.candidate(),trial),measuredWeight(b.candidate(),trial));
                const double ra=a.candidate().weight.log_abs_matter-ref.state().weight.log_abs_matter;
                const double rb=b.candidate().weight.log_abs_matter-candidate.state().weight.log_abs_matter;
                maxRatio=std::max(maxRatio,std::abs(ra-rb));
                maxAcceptance=std::max(maxAcceptance,std::abs(std::exp(std::min(0.,ra))-std::exp(std::min(0.,rb))));
                const double u=(double(rng()>>12)+.5)*0x1p-52;
                if((std::log(u)<std::min(0.,ra))!=(std::log(u)<std::min(0.,rb)))throw std::runtime_error("decision mismatch");
                // Exercise both transaction operations independent of a random
                // path's high acceptance. MH decision agreement checked above.
                if(n%3) {ref.accept(a);candidate.accept(b);++accepts;}
                else {ref.reject(a);candidate.reject(b);++rejects;
                    if((ref.state().weight.center_green-before.weight.center_green).norm()!=0)throw std::runtime_error("rollback");}
                compare(measuredWeight(ref.state(),trial),measuredWeight(candidate.state(),trial));
            }
            auto before=candidate.state();
            try { candidate.propose(32,0,std::numeric_limits<double>::infinity());throw std::logic_error("missing failure"); }
            catch(const std::invalid_argument&) {}
            compare(measuredWeight(before,trial),measuredWeight(candidate.state(),trial));
        }
        auto prefixReplay=[&](int lx,int ly,int shear) {
            const KitaevSSHModel smallModel(HoneycombLattice(lx,ly,shear),p2Couplings(.3));
            std::mt19937_64 localRng(std::uint64_t(12000+shear));
            const int intervals=8;
            const auto smallPath=PhononPath(smallModel.lattice().bondCount(),intervals,.2,
                {1.7,.8},std::vector<double>(smallModel.lattice().bondCount(),.1378676470588235)).freeSample(localRng);
            const auto smallTrial=freeGroundState(smallModel);
            const auto smallGenerator=honeycombGenerator(smallModel);
            const PhononWeightOptions uncachedOptions{true,true,false,false};
            const PhononWeightOptions cachedOptions{true,true,false,true,true};
            CoupledPhononPath uncached(smallPath,smallTrial,smallGenerator,4,uncachedOptions);
            CoupledPhononPath cached(smallPath,smallTrial,smallGenerator,4,cachedOptions);
            CoupledPhononPath prefixOnly(smallPath,smallTrial,smallGenerator,4,{true,true,false,true});
            auto compareState=[&](const CoupledPhononPath::Snapshot &a,
                                  const CoupledPhononPath::Snapshot &b) {
                compare(measuredWeight(a,smallTrial),measuredWeight(b,smallTrial));
                require(a.configuration_hash==b.configuration_hash,"prefix configuration mismatch");
                require(a.weight.magnitude_valid==b.weight.magnitude_valid&&
                        a.weight.phase_valid==b.weight.phase_valid&&
                        a.weight.green_valid==b.weight.green_valid,"prefix validity mismatch");
            };
            require(cached.state().signed_prefix_cache&&
                    cached.state().signed_prefix_cache->checkpoints.size()==
                        cached.state().weight_order.size()+1,"prefix cache was not built");
            compareState(uncached.state(),cached.state());
            const std::vector<int> nodes{intervals,intervals/2,0,intervals/2,intervals,0};
            for(std::size_t step=0;step<nodes.size();++step) {
                const int node=nodes[step],bond=int(step)%smallModel.lattice().bondCount();
                const double value=smallPath.coordinate(node,bond)+.017*double(step+1);
                const auto beforeCache=cached.state().signed_prefix_cache;
                const auto beforeCayley=cached.state().lu_cayley_cache;
                auto a=uncached.propose(node,bond,value);
                auto b=cached.propose(node,bond,value);
                auto c=prefixOnly.propose(node,bond,value);
                compareState(b.candidate(),c.candidate());
                require(b.candidate().signed_prefix_cache&&
                        b.candidate().signed_prefix_cache->checkpoints.size()==
                            b.candidate().weight_order.size()+1,"candidate prefix cache missing");
                compareState(a.candidate(),b.candidate());
                const auto candidateCache=b.candidate().signed_prefix_cache;
                const auto candidateCayley=b.candidate().lu_cayley_cache;
                require(bool(candidateCayley),"missing candidate Cayley cache");
                if(step%3==0) {
                    uncached.accept(a);cached.accept(b);prefixOnly.accept(c);
                    require(cached.state().lu_cayley_cache==candidateCayley,"Cayley accept ownership");
                    require(cached.state().signed_prefix_cache==candidateCache,"accepted cache not committed");
                } else {
                    uncached.reject(a);cached.reject(b);prefixOnly.reject(c);
                    require(cached.state().lu_cayley_cache==beforeCayley,"Cayley reject ownership");
                    require(cached.state().signed_prefix_cache==beforeCache,"rejected cache was committed");
                }
                compareState(uncached.state(),cached.state());
            }
            const auto before=cached.state().lu_cayley_cache;
            try {cached.propose(0,0,std::numeric_limits<double>::infinity());require(false,"missing failure");}
            catch(const std::invalid_argument&) {}
            require(cached.state().lu_cayley_cache==before,"failed proposal changed cache");
            // An invalid checkpoint must rebuild the correct hashes, and that
            // rebuilt cache must actually resume on the following evaluation.
            const auto &slices=cached.state().weight_order;
            PureFastOptions opts;
            auto sameResult=[&](const PureProjectorWeightResult &a,const PureProjectorWeightResult &b) {
                require(a.ok()&&b.ok()&&a.status==b.status,"core status mismatch");
                require(std::abs(a.log_abs_weight-b.log_abs_weight)<2e-10,"core log mismatch");
                require(std::abs(a.complex_phase-b.complex_phase)<2e-10,"core phase mismatch");
                require((a.green-b.green).norm()<2e-11,"core Green mismatch");
            };
            PureProjectorSignedPrefixCache good;
            const auto reference=pureProjectorStableReferenceWeight(
                smallTrial,slices,opts,nullptr,nullptr,nullptr,nullptr,&good);
            const int cut=3;
            for(int corruption=0;corruption<2;++corruption) {
                auto bad=good;
                auto invalid=std::make_shared<PureProjectorSignedPrefixCheckpoint>(*good.checkpoints[cut]);
                if(corruption==0)invalid->prefix_hash^=1;
                else invalid->q(0,0)=DataType(std::numeric_limits<double>::quiet_NaN(),0);
                bad.checkpoints[cut]=invalid;
                PureProjectorSignedPrefixCache rebuilt;
                const auto resumes=signedPrefixResumeCount,steps=signedPropagationSliceCount;
                sameResult(reference,pureProjectorStableReferenceWeight(
                    smallTrial,slices,opts,nullptr,nullptr,nullptr,&bad,&rebuilt,cut));
                require(signedPrefixResumeCount==resumes,"invalid checkpoint resumed");
                require(signedPropagationSliceCount-steps==slices.size(),"fallback did not fully propagate");
                require(rebuilt.checkpoints.size()==good.checkpoints.size(),"rebuild size");
                for(std::size_t i=0;i<good.checkpoints.size();++i)
                    require(rebuilt.checkpoints[i]->prefix_hash==good.checkpoints[i]->prefix_hash,"rebuilt hash mismatch");
                sameResult(reference,pureProjectorStableReferenceWeight(
                    smallTrial,slices,opts,nullptr,nullptr,nullptr,&rebuilt,nullptr,cut));
                require(signedPrefixResumeCount==resumes+1,"rebuilt checkpoint did not resume");
            }
            PureProjectorLUCayleyCache records;
            auto solve=[&](const std::vector<PureProjectorSlice> &f,const PureFastOptions &o) {
                return pureProjectorStableReferenceWeight(
                    smallTrial,f,o,nullptr,nullptr,nullptr,nullptr,nullptr,-1,&records);
            };
            const auto cold=cayleyLUSolveCount;
            sameResult(reference,solve(slices,opts));
            require(cayleyLUSolveCount>cold,"cold cache skipped LU");
            const auto warm=cayleyLUSolveCount,hits=cayleyCacheHitCount;
            sameResult(reference,solve(slices,opts));
            require(cayleyLUSolveCount==warm&&cayleyCacheHitCount>hits,"warm cache did not hit");
            auto copied=slices;
            copied[cut].matrix=MatType(slices[cut].matrix); // equal values, distinct storage
            sameResult(reference,solve(copied,opts));
            require(cayleyLUSolveCount==warm+1,"unverified storage did not fall back");
            auto changedOptions=opts;changedOptions.residual_tolerance*=.5;
            const auto invalidated=cayleyLUSolveCount;
            sameResult(reference,solve(copied,changedOptions));
            require(cayleyLUSolveCount>invalidated,"option change did not invalidate");
            records[cut].reset();
            const auto missing=cayleyLUSolveCount;
            sameResult(reference,solve(copied,changedOptions));
            require(cayleyLUSolveCount==missing+1,"missing record did not use LU");
            auto singular=copied;
            singular[cut].matrix=-MatType::Identity(smallTrial.Phi.rows(),smallTrial.Phi.rows());
            const auto failed=solve(singular,opts);
            const auto failedRef=pureProjectorStableReferenceWeight(smallTrial,singular,opts);
            require(!failed.ok()&&failed.status==failedRef.status&&failed.first_failing_slice==failedRef.first_failing_slice,
                    "unavailable Cayley changed failure status");
            CoupledPhononPath deferred(smallPath,smallTrial,smallGenerator,4,{true,true,true,true,true});
            const auto deferredSolves=cayleyLUSolveCount,deferredHits=cayleyCacheHitCount;
            auto deferredProposal=deferred.propose(0,0,smallPath.coordinate(0,0)+.02);
            require(!deferredProposal.candidate().lu_cayley_cache&&!deferredProposal.candidate().weight.phase_valid&&
                    !deferredProposal.candidate().weight.green_valid,"deferred built phase/cache");
            require(cayleyLUSolveCount==deferredSolves&&cayleyCacheHitCount==deferredHits,"deferred Cayley work");
        };
        prefixReplay(3,2,0);
        prefixReplay(3,2,1);
        std::cout<<std::setprecision(17)<<"{\"max_log_abs_error\":"<<maxLog<<",\"max_phase_error\":"<<maxPhase
                 <<",\"max_green_error\":"<<maxGreen<<",\"max_log_ratio_error\":"<<maxRatio
                 <<",\"max_acceptance_error\":"<<maxAcceptance<<",\"accept_transactions\":"<<accepts
                 <<",\"reject_transactions\":"<<rejects<<"}\n";
        // Accumulated double-precision error across 128 conditioned factors;
        // core trust thresholds remain unchanged and apply to every call.
        if(maxLog>2e-10 || maxPhase>2e-10 || maxRatio>4e-10 || maxAcceptance>4e-10 || maxGreen>2e-12)
            throw std::runtime_error("optimization replay tolerance exceeded");
    } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}
