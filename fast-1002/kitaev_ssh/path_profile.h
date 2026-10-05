#ifndef KITAEV_SSH_PATH_PROFILE_H
#define KITAEV_SSH_PATH_PROFILE_H
#include <chrono>
#include <map>
#include <string>
namespace kitaev_ssh {
// Enabled only in the dedicated profiling executable. Times are inclusive
// within each named operation; core-internal allocations remain in that core.
struct PathProfile { std::map<std::string,double> seconds; std::map<std::string,unsigned long long> calls; };
inline thread_local PathProfile *activePathProfile=nullptr;
template<class F> auto profilePath(const char *name,F operation) -> decltype(operation()) {
#ifdef KITAEV_SSH_PROFILE
    if(activePathProfile) {
        const auto begin=std::chrono::steady_clock::now();
        auto result=operation();
        activePathProfile->seconds[name]+=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        ++activePathProfile->calls[name];
        return result;
    }
#endif
    return operation();
}
}
#endif
