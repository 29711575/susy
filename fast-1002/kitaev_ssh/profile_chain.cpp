#define KITAEV_SSH_PROFILE
#include "phonon_mc.h"
#include <iostream>
#include <string>
using namespace kitaev_ssh;
int main(int argc,char **argv) {
    try {
        const std::string proposalName=argc>1?argv[1]:"bond-refresh";
        const bool deferred=argc>2?std::stoi(argv[2])!=0:false;
        const int intervals=argc>3?std::stoi(argv[3]):64;
        const int sweeps=argc>4?std::stoi(argv[4]):30;
        const bool prefixCache=argc>5?std::stoi(argv[5])!=0:false;
        const bool cayleyCache=argc>6?std::stoi(argv[6])!=0:false;
        if(proposalName!="local"&&proposalName!="bond-refresh") throw std::invalid_argument("proposal");
        if(intervals<=0||sweeps<=0) throw std::invalid_argument("intervals/sweeps");
        const KitaevSSHModel model(HoneycombLattice(4,4),p2Couplings(.3));
        std::mt19937_64 rng(7701);
        const auto path=PhononPath(48,intervals,.1,{1.7,.8},std::vector<double>(48,.1378676470588235)).freeSample(rng);
        CoupledPhononPath state(path,freeGroundState(model),honeycombGenerator(model),4,{true,true,deferred,prefixCache,cayleyCache});
        PhononSampler sampler(state,rng,proposalName=="local"?PhononProposal::LocalConditional:PhononProposal::BondRefresh);
        PathProfile profile;
        CoreProfile coreProfile;
        activePathProfile=&profile;
        activeCoreProfile=&coreProfile;
        const auto start=std::chrono::steady_clock::now();
        for(int sweep=0;sweep<sweeps;++sweep) sampler.sweep();
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        activePathProfile=nullptr;
        activeCoreProfile=nullptr;
        std::size_t prefixBytes=0,prefixCheckpoints=0;
        if(state.state().signed_prefix_cache) {
            prefixCheckpoints=state.state().signed_prefix_cache->checkpoints.size();
            for(const auto &checkpoint:state.state().signed_prefix_cache->checkpoints)
                if(checkpoint)prefixBytes+=sizeof(*checkpoint)+
                    std::size_t(checkpoint->q.size())*sizeof(DataType);
        }
        std::cout<<std::setprecision(17)<<"{\"seed\":7701,\"proposal\":\""<<proposalName<<"\",\"deferred\":"<<std::boolalpha<<deferred
                 <<",\"intervals\":"<<intervals<<",\"sweeps\":"<<sweeps<<",\"tau\":"<<intervals*.1<<",\"dt\":0.1,\"elapsed_seconds\":"<<elapsed
                 <<",\"acceptance\":"<<sampler.counts().acceptance()<<",\"proposals\":"<<sampler.counts().proposed<<",\"timings\":{";
        bool first=true;
        for(const auto &entry:profile.seconds) {
            if(!first)std::cout<<',';first=false;
            std::cout<<'"'<<entry.first<<"\":{\"seconds\":"<<entry.second<<",\"calls\":"<<profile.calls[entry.first]<<'}';
        }
        std::cout<<"},\"core_timings\":{";
        first=true;
        for(const auto &entry:coreProfile.seconds) {
            if(!first)std::cout<<',';first=false;
            std::cout<<'"'<<entry.first<<"\":{\"seconds\":"<<entry.second<<",\"calls\":"<<coreProfile.calls[entry.first]<<'}';
        }
        std::cout<<"},\"final_log_abs_matter\":"<<state.state().weight.log_abs_matter
                 <<",\"final_phase\":["<<state.state().weight.matter_phase.real()<<','<<state.state().weight.matter_phase.imag()<<']'
                 <<",\"cayley_cache\":"<<cayleyCache
                 <<",\"cayley_lu_calls\":"<<cayleyLUSolveCount
                 <<",\"cayley_cache_hits\":"<<cayleyCacheHitCount
                 <<",\"cayley_cache_fallbacks\":"<<cayleyCacheFallbackCount
                 <<",\"signed_factor_steps\":"<<signedPropagationSliceCount
                 <<",\"cayley_bytes_estimate\":"<<(state.state().lu_cayley_cache?pureProjectorCayleyBytes(*state.state().lu_cayley_cache):0)
                 <<",\"prefix_cache\":"<<prefixCache<<",\"prefix_resumes\":"<<signedPrefixResumeCount
                 <<",\"prefix_skipped_slices\":"<<signedPrefixSkippedSliceCount
                 <<",\"prefix_checkpoints\":"<<prefixCheckpoints<<",\"prefix_bytes_estimate\":"<<prefixBytes<<"}\n";
    } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}
