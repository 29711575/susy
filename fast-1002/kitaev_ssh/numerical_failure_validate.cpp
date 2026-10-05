// Injection is compiled into this validator only, never the production binary.
#define KITAEV_SSH_TEST_FAILURE_INJECTION
#define main coupled_mc_test_main
#include "coupled_mc.cpp"
#undef main
#include <filesystem>

namespace {
void demand(bool ok,const std::string &what) {
    if(!ok) throw std::runtime_error("numerical failure validator: "+what);
}
void equalSnapshots(const CoupledPhononPath::Snapshot &a,const CoupledPhononPath::Snapshot &b) {
    demand(a.configuration_hash==b.configuration_hash,"RNG stream / accepted coordinates changed");
    demand(std::abs(a.weight.log_abs_matter-b.weight.log_abs_matter)<2e-10,"weight changed");
    demand(a.weight.phase_valid==b.weight.phase_valid && a.weight.green_valid==b.weight.green_valid,"validity changed");
    if(a.weight.phase_valid)demand(std::abs(a.weight.matter_phase-b.weight.matter_phase)<2e-10,"phase changed");
    if(a.weight.green_valid)demand((a.weight.center_green-b.weight.center_green).norm()<2e-10,"Green changed");
}
void transactions(PhononWeightOptions options,PhononProposal mode) {
    KitaevSSHModel model(HoneycombLattice(3,2,1),p2Couplings(.3));
    std::mt19937_64 rng(19311);
    auto path=PhononPath(model.lattice().bondCount(),4,.2,{1.7,.8}).freeSample(rng);
    auto trial=freeGroundState(model); auto generator=honeycombGenerator(model);
    CoupledPhononPath normal(path,trial,generator,4,options), injected(path,trial,generator,4,options);
    std::vector<NumericalEvent> events;
    PhononSampler a(normal,rng,mode), b(injected,rng,mode,[&](const auto &e,const auto&){events.push_back(e);});
    a.setStage("warmup");b.setStage("warmup");
    int first=0,retries=0;
    injected.failureInjection=[&](const char *op,bool reference){
        if(std::string(op)=="measurement")return;
        if(reference)++retries;
        else if(first++==0)throw std::runtime_error("overlap_untrusted injected first evaluation");
    };
    // Entire sweeps, including shuffle, normal-distribution cache and all u's,
    // must remain identical after recovery. A later sweep detects extra draws.
    for(int i=0;i<2;++i){a.sweep();b.sweep();equalSnapshots(normal.state(),injected.state());}
    demand(retries==1&&events.size()==1&&events[0].recovered,"one recovery event");
    demand(a.counts().accepted==b.counts().accepted&&a.counts().proposed==b.counts().proposed,"acceptance counts");
    demand(b.counts().first_failures==1&&b.counts().reference_recovered==1&&b.counts().reference_failed==0&&b.counts().untrusted_rejected==0,"recovery counters");
    injected.failureInjection={};a.ensureMeasured();b.ensureMeasured();equalSnapshots(normal.state(),injected.state());
    const auto obsA=basicObservables(model,normal.state()),obsB=basicObservables(model,injected.state());
    for(std::size_t i=0;i<obsA.size();++i)demand(std::abs(obsA[i]-obsB[i])<2e-10,"measurements");

    // A recovered candidate must still undergo MH rejection when its ratio
    // demands it, rather than accepting merely because reconstruction worked.
    bool rejected=false;
    for(double delta:{-1.,1.,-3.,3.}) {
        const double value=injected.state().path.coordinate(2,0)+delta;
        auto p=injected.propose(2,0,value,true);
        const double ratio=p.candidate().weight.log_abs_matter-injected.state().weight.log_abs_matter;
        injected.reject(p);
        if(ratio>=-1e-6)continue;
        const auto saved=injected.state();const auto c=b.counts();
        injected.failureInjection=[](const char*,bool reference){if(!reference)throw std::runtime_error("overlap_untrusted rejected recovery");};
        const auto move=b.attempt(2,0,value,std::nextafter(1.,0.));
        demand(!move.accepted&&std::abs(move.log_ratio-ratio)<2e-10,"recovered MH rejection");
        demand(b.counts().metropolis_rejected==c.metropolis_rejected+1,"ordinary rejection count");
        equalSnapshots(saved,injected.state());
        demand(saved.signed_prefix_cache==injected.state().signed_prefix_cache,"recovered rejection cache");
        injected.failureInjection={};rejected=true;break;
    }
    demand(rejected,"negative-ratio test candidate unavailable");

    const auto before=injected.state();const auto counts=b.counts();
    injected.failureInjection=[](const char*,bool){throw std::runtime_error("green_residual_exceeded injected persistent failure");};
    bool failed=false;
    try {
        if(mode==PhononProposal::LocalConditional)b.attempt(2,0,before.path.coordinate(2,0)+.1,.25);
        else {std::vector<double> values(before.path.nodes(),.1);b.attemptBond(0,values,.25);}
    }catch(const std::runtime_error&){failed=true;}
    demand(failed,"persistent failure must terminate");equalSnapshots(before,injected.state());
    demand(before.signed_prefix_cache==injected.state().signed_prefix_cache,"failure committed cache");
    demand(b.counts().proposed==counts.proposed+1&&b.counts().accepted==counts.accepted&&
        b.counts().metropolis_rejected==counts.metropolis_rejected&&b.counts().reference_failed==1&&
        b.counts().first_failures==counts.first_failures+1&&events.size()==3&&!events.back().recovered,"fatal counters");
    if(options.defer_phase_measurement) {
        injected.failureInjection={};b.attempt(2,0,before.path.coordinate(2,0),.25);
        const auto deferred=injected.state();const auto measuredCounts=b.counts();
        injected.failureInjection=[](const char *op,bool){if(std::string(op)=="measurement")throw std::runtime_error("solve_residual_exceeded measurement");};
        PhaseStatistics stats({"x"},1);bool measured=false;failed=false;
        try {sample(b,{0,1,1,1},stats,[&](const auto&){measured=true;return std::vector<Complex>{1.};},[](const auto&){});}
        catch(const std::runtime_error&){failed=true;}
        demand(failed&&!measured&&stats.blocks().empty(),"failed measurement was included/skipped");
        demand(b.counts().reference_failed==measuredCounts.reference_failed+1,"measurement failure counted");
        const auto beforeMeasure=injected.state();
        injected.failureInjection=[](const char *op,bool reference){if(std::string(op)=="measurement"&&!reference)throw std::runtime_error("overlap_untrusted measurement recovery");};
        b.ensureMeasured();demand(injected.state().weight.phase_valid&&injected.state().weight.green_valid,"measurement recovery");
        demand(beforeMeasure.configuration_hash==injected.configurationHash(),"measurement moved path");
    }
}
void cliCase(const std::filesystem::path &dir,const std::string &name,const std::string &operation,bool persistent,bool deferred=false) {
    int fired=0;
    testNumericalFailure=[&](const char *op,bool reference){
        if(operation!=op)return;
        // Initial search must be allowed a second fresh draw after both paths
        // fail for the first draw. Production operations must instead abort.
        if(operation=="initialization") {if((persistent&&fired<2)||(!persistent&&fired<1)){++fired;throw std::runtime_error("overlap_untrusted injected initial");}return;}
        if((!reference&&fired==0)||persistent){++fired;throw std::runtime_error("solve_residual_exceeded injected CLI");}
    };
    std::vector<std::string> args={"numerical_failure_validate","--Lx","3","--Ly","2","--shear","1",
        "--intervals","2","--warmup","0","--sweeps","2","--measure-every","1","--block-samples","1",
        "--seed","19311","--defer-phase",deferred?"1":"0","--output",(dir/name).string()};
    std::vector<char*> argv;for(auto &arg:args)argv.push_back(arg.data());
    const int code=coupled_mc_test_main(int(argv.size()),argv.data());testNumericalFailure={};
    demand(code==(persistent&&operation!="initialization"?1:0),"CLI exit code "+name);
    demand(fired>0,"CLI injection not reached");
}
}
int main(int argc,char **argv) {
    try {
        transactions({false,false,false,false},PhononProposal::LocalConditional);
        transactions({true,true,false,true},PhononProposal::LocalConditional);
        transactions({true,true,true,false},PhononProposal::LocalConditional);
        transactions({true,true,false,false},PhononProposal::BondRefresh);
        if(argc==2) {
            const std::filesystem::path dir=argv[1];std::filesystem::create_directories(dir);
            cliCase(dir,"recover","local",false);
            cliCase(dir,"abort","local",true);
            cliCase(dir,"measure_recover","measurement",false,true);
            cliCase(dir,"measure_abort","measurement",true,true);
            cliCase(dir,"initial_recover","initialization",false);
            cliCase(dir,"initial_redraw","initialization",true);
        }
        std::cout<<"PASS numerical failure injection: recovery, abort, counters, RNG, transactions, deferred measurement\n";
        return 0;
    }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
