#ifndef QMC_CORE_PROFILE_H
#define QMC_CORE_PROFILE_H
// Fine-grained, compile-time-gated timers for the numerical-core hot path.
//
// This header deliberately lives under inc/ (global namespace) so that the
// core projector headers can be instrumented without depending on the
// kitaev_ssh application layer. It mirrors kitaev_ssh/path_profile.h: named
// inclusive timers plus call counts in a thread-local registry.
//
// When KITAEV_SSH_PROFILE is not defined, profileCore() reduces to a direct
// call with no clock access, no map lookup and no name evaluation at runtime;
// production binaries therefore have zero instrumentation overhead.
#include <chrono>
#include <map>
#include <string>
#include <type_traits>

struct CoreProfile {
    std::map<std::string,double> seconds;
    std::map<std::string,unsigned long long> calls;
};
inline thread_local CoreProfile *activeCoreProfile=nullptr;

template<class F>
auto profileCore(const char *name, F operation) -> decltype(operation()) {
#ifdef KITAEV_SSH_PROFILE
    if (activeCoreProfile) {
        const auto begin=std::chrono::steady_clock::now();
        if constexpr (std::is_same_v<decltype(operation()),void>) {
            operation();
            activeCoreProfile->seconds[name]+=
                std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
            ++activeCoreProfile->calls[name];
            return;
        } else {
            auto result=operation();
            activeCoreProfile->seconds[name]+=
                std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
            ++activeCoreProfile->calls[name];
            return result;
        }
    }
#endif
    (void)name;
    return operation();
}

// Optional byte-watermark probe, only used by the profiling executable.
#ifdef KITAEV_SSH_PROFILE
#include <cstddef>
#include <fstream>
#include <string>
#include <unistd.h>
namespace qmc_core_profile_detail {
inline std::size_t residentBytes() {
    std::ifstream status("/proc/self/status");
    std::string key; long value=0; std::string unit;
    while(status>>key) {
        if(key=="VmRSS:"){ status>>value>>unit; return std::size_t(value)*1024; }
        std::string rest; std::getline(status,rest);
    }
    return 0;
}
}
#endif

#endif
