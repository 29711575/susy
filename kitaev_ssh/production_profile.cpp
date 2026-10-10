// Run the actual production CLI, including warmup, all selected measurements
// and output. The normal MC executable has no profiling clocks or registries.
#define KITAEV_SSH_PROFILE
#define main kitaevProductionMain
#include "coupled_mc.cpp"
#undef main

int main(int argc,char **argv) {
    std::string prefix="kitaev_ssh_mc";
    for(int i=1;i+1<argc;++i)
        if(std::string(argv[i])=="--output") prefix=argv[i+1];
    const std::string file=prefix+"_profile.json";
    if(outputPathExists(file)) {
        std::cerr<<"profile output exists; choose another --output prefix\n";
        return 1;
    }
    PathProfile path;
    CoreProfile core;
    activePathProfile=&path;
    activeCoreProfile=&core;
    const int rc=kitaevProductionMain(argc,argv);
    activePathProfile=nullptr;
    activeCoreProfile=nullptr;
    if(rc!=0) return rc;
    std::ofstream out(file);
    out<<std::setprecision(17)<<"{\"schema\":\"production-profile-v1\","
        <<"\"scope\":\"initialization, warmup and measurement; actual production CLI\","
        <<"\"timers_are_inclusive\":true,\"path\":{";
    auto write=[&](const auto &profile) {
        bool comma=false;
        for(const auto &entry:profile.seconds) {
            if(comma) out<<',';
            comma=true;
            jsonString(out,entry.first);
            out<<":{\"seconds\":"<<entry.second<<",\"calls\":"<<profile.calls.at(entry.first)<<'}';
        }
    };
    write(path);
    out<<"},\"core\":{";
    write(core);
    out<<"}}\n";
    out.close();
    if(!out) {
        std::cerr<<"production profile output write failed\n";
        return 1;
    }
    return 0;
}
