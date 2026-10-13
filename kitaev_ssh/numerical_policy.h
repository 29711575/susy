#ifndef KITAEV_SSH_NUMERICAL_POLICY_H
#define KITAEV_SSH_NUMERICAL_POLICY_H

namespace kitaev_ssh {
// Numerical implementation only: sampling, sources, trial states and strict
// certification/failure gates are deliberately outside this policy.
// Keep the independent uncached recovery implementation in phonon_projector.h.
struct NumericalPolicy {
    const char *name;
    const char *weight_engine;
    int defer_phase, prefix_cache, cayley_cache, center_green_fast;
};

inline constexpr NumericalPolicy productionNumerics{
    "production-v1", "combined", 1, 1, 1, 1};
// Diagnostic configuration, not a claim to disable every shared-library
// optimization. Full independent recovery still uses uncachedReference.
inline constexpr NumericalPolicy referenceNumerics{
    "reference-v1", "reference", 0, 0, 0, 0};
} // namespace kitaev_ssh

#endif
