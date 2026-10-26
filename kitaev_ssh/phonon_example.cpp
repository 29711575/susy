#include "phonon_projector.h"
#include "physical_matter_trial.h"
#include <fstream>
#include <iostream>
#include <map>
#include <set>

using namespace kitaev_ssh;
int main(int argc,char **argv) {
    try {
        std::map<std::string,std::string> args;
        const std::set<std::string> allowed{"Lx","Ly","intervals","dt","mass","omega-ph","g","seed","save","load","trial-dx","trial-dy","trial-dz"};
        for(int i=1;i<argc;++i) {
            std::string key=argv[i];
            if(key=="--help") {
                std::cout<<"P2 free-path initialization and forced accept/reject demo, NOT coupled MC.\n"
                    <<"--Lx 2 --Ly 2 --intervals 8 --dt 0.05 --mass 4 --omega-ph 0.5 --g 0.3 --seed 42\n"
                    <<"--trial-dx 0 --trial-dy 0 --trial-dz 0 (uniform endpoint shifts).\n"
                    <<"--save FILE / --load FILE (coordinates, parameters and trial shifts only). J=(1,1,0.5).\n";
                return 0;
            }
            if(key.rfind("--",0)!=0 || !allowed.count(key.substr(2)) || args.count(key.substr(2)) || i+1>=argc)
                throw std::invalid_argument("unknown, repeated or incomplete argument: "+key);
            args[key.substr(2)]=argv[++i];
        }
        auto real=[&](const std::string &key,double fallback) {
            if(!args.count(key)) return fallback;
            std::size_t end; const double v=std::stod(args.at(key),&end);
            if(end!=args.at(key).size() || !std::isfinite(v)) throw std::invalid_argument("invalid --"+key);
            return v;
        };
        auto integer=[&](const std::string &key,int fallback) {
            double v=real(key,fallback);
            if(v!=std::floor(v) || v<0 || v>std::numeric_limits<int>::max()) throw std::invalid_argument("invalid --"+key);
            return int(v);
        };
        if(args.count("load")) for(const std::string key:{"intervals","dt","mass","omega-ph","trial-dx","trial-dy","trial-dz"})
            if(args.count(key)) throw std::invalid_argument("--load conflicts with --"+key);
        const double g=real("g",.3);
        const KitaevSSHModel model(HoneycombLattice(integer("Lx",2),integer("Ly",2)),p2Couplings(g));
        const double d[3]={real("trial-dx",0),real("trial-dy",0),real("trial-dz",0)};
        std::vector<double> shifts(model.lattice().bondCount());
        for(const auto &b:model.lattice().bonds()) shifts[b.index]=d[static_cast<int>(b.type)];
        std::mt19937_64 rng(integer("seed",42));
        auto initial=[&]() {
            if(args.count("load")) { std::ifstream in(args.at("load")); return PhononPath::load(in); }
            return PhononPath(model.lattice().bondCount(),integer("intervals",8),real("dt",.05),
                              {real("mass",4.),real("omega-ph",.5)},shifts).freeSample(rng);
        };
        CoupledPhononPath state(initial(),selectFreeMatterTrial(model).trial,honeycombGenerator(model));
        auto rejected=state.propose(state.state().path.center(),0,state.state().path.coordinate(state.state().path.center(),0)+.1);
        state.reject(rejected);
        auto accepted=state.propose(0,1,state.state().path.coordinate(0,1)-.15);
        const double logRatio=accepted.logAbsRatio(); const auto phaseRatio=accepted.phaseRatio();
        state.accept(accepted);
        if(args.count("save")) { std::ofstream out(args.at("save")); state.state().path.save(out); }
        const auto &s=state.state(); const auto &p=s.path.parameters();
        std::cout<<std::setprecision(17)<<"{\"scope\":\"fixed_gauge_matter_sector\",\"initialization\":\""
            <<(args.count("load")?"loaded_phonon_path":"exact_free_phonon_path")<<"\",\"coupled_mc\":false,\"forced_transaction_demo\":true"
            <<",\"Jx\":1,\"Jy\":1,\"Jz\":0.5,\"g\":"<<g<<",\"mass\":"<<p.mass<<",\"omega_ph\":"<<p.omega
            <<",\"k\":"<<p.spring()<<",\"lambda\":"<<p.lambda(g)<<",\"intervals\":"<<s.path.intervals()<<",\"nodes\":"<<s.path.nodes()
            <<",\"tau\":"<<s.path.totalTau()<<",\"center_factor_cut\":"<<s.path.intervals()
            <<",\"log_phonon\":"<<s.weight.log_phonon<<",\"log_abs_matter\":"<<s.weight.log_abs_matter
            <<",\"matter_phase\":["<<s.weight.matter_phase.real()<<','<<s.weight.matter_phase.imag()<<']'
            <<",\"proposal_log_abs_ratio\":"<<logRatio<<",\"proposal_phase_ratio\":["<<phaseRatio.real()<<','<<phaseRatio.imag()<<']'
            <<",\"phonon_trial\":\"equal_shifted_vacuum_endpoints\",\"trial_shift_by_bond\":[";
        for(int b=0;b<s.path.bonds();++b){if(b)std::cout<<',';std::cout<<s.path.trialShift(b);}
        std::cout<<"]}\n";
        return 0;
    } catch(const std::exception &e) { std::cerr<<"kitaev_ssh_phonon_example: "<<e.what()<<'\n'; return 1; }
}
