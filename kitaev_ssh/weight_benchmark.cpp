#ifdef BENCH_PROFILE
#define KITAEV_SSH_PROFILE
#endif
#include "phonon_mc.h"
#include "vbs_observables.h"
#include <iostream>
using namespace kitaev_ssh;
int main(int argc,char **argv) {
    try {
        if(argc!=2)throw std::invalid_argument("mode 0=reference 1=svd 2=merge 3=combined");
        const int mode=std::stoi(argv[1]); if(mode<0||mode>3)throw std::invalid_argument("mode");
        const KitaevSSHModel model(HoneycombLattice(4,4),p2Couplings(.3));
        std::mt19937_64 rng(7701);
        const auto initial=PhononPath(48,64,.1,{1.7,.8},std::vector<double>(48,.1378676470588235)).freeSample(rng);
        CoupledPhononPath state(initial,freeGroundState(model),honeycombGenerator(model),4,{bool(mode&1),bool(mode&2)});
        PhononSampler sampler(state,rng,PhononProposal::BondRefresh);
        VBSObservables vbs(model.lattice());auto names=basicObservableNames();
        const auto extra=vbs.names();names.insert(names.end(),extra.begin(),extra.end());
        PhaseStatistics statistics(names,5);PathProfile profile;
#ifdef BENCH_PROFILE
        activePathProfile=&profile;
#endif
        const auto begin=std::chrono::steady_clock::now();
        sample(sampler,{0,30,1,5},statistics,[&](const auto &s) {
            auto values=basicObservables(model,s);const auto more=vbs.measure(s.path.node(s.path.center()));
            values.insert(values.end(),more.begin(),more.end());return values;
        },[](const auto&){});
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        activePathProfile=nullptr;
        std::cout<<std::setprecision(17)<<"{\"mode\":"<<mode<<",\"sweeps\":30,\"seed\":7701,\"elapsed_seconds\":"<<elapsed
                 <<",\"seconds_per_sweep\":"<<elapsed/30<<",\"accepted\":"<<sampler.counts().accepted
                 <<",\"final_log_abs\":"<<state.state().weight.log_abs_matter
                 <<",\"final_phase\":["<<state.state().weight.matter_phase.real()<<','<<state.state().weight.matter_phase.imag()
                 <<"],\"timings\":{";
        bool first=true;for(const auto &v:profile.seconds){if(!first)std::cout<<',';first=false;std::cout<<'"'<<v.first<<"\":"<<v.second;}
        std::cout<<"}}\n";
    }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
