#ifndef PURE_PROJECTOR_DRIVEN_H
#define PURE_PROJECTOR_DRIVEN_H

#include "pure_projector_fast.h"

#include <algorithm>
#include <complex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct PureDrivenContourResult {
    PureFastConfiguration configuration;
    int ket_factor_count = 0;
    int bra_factor_count = 0;
    bool bra_ket_hs_independent = true;
    std::string contour_order =
        "bra_strict_reverse_adjoint_then_ket_canonical_action";
};

inline PureDrivenContourResult pureProjectorIndependentBraKetContour(
        PureFastConfiguration ketCanonicalAction,
        PureFastConfiguration braCanonicalAction) {
    if (ketCanonicalAction.slices.empty() || braCanonicalAction.slices.empty() ||
        ketCanonicalAction.slices.size() != ketCanonicalAction.hs_fields.size() ||
        ketCanonicalAction.slices.size() != ketCanonicalAction.locations.size() ||
        braCanonicalAction.slices.size() != braCanonicalAction.hs_fields.size() ||
        braCanonicalAction.slices.size() != braCanonicalAction.locations.size() ||
        ketCanonicalAction.slices.size() != braCanonicalAction.slices.size())
        throw std::invalid_argument("invalid independent bra/ket branch factors");

    PureDrivenContourResult result;
    result.ket_factor_count = int(ketCanonicalAction.slices.size());
    result.bra_factor_count = int(braCanonicalAction.slices.size());
    result.configuration = std::move(ketCanonicalAction);
    for (int index = int(braCanonicalAction.slices.size()) - 1;
         index >= 0; --index) {
        const PureProjectorSlice &source = braCanonicalAction.slices[index];
        PureSliceLocation location = braCanonicalAction.locations[index];
        location.branch = PureBranch::Bra;
        result.configuration.slices.emplace_back(
            source.matrix.adjoint().eval(), std::conj(source.eta),
            "bra_adjoint:" + source.label);
        result.configuration.hs_fields.push_back(
            braCanonicalAction.hs_fields[index]);
        result.configuration.locations.push_back(location);
    }
    // The Majorana-to-Fock representation reverses the complete flattened
    // list.  Swap only the two completed blocks so its dense action is
    // C_bra^dagger C_ket, while preserving both branch-internal orders.
    std::rotate(result.configuration.slices.begin(),
                result.configuration.slices.begin() + result.ket_factor_count,
                result.configuration.slices.end());
    std::rotate(result.configuration.hs_fields.begin(),
                result.configuration.hs_fields.begin() + result.ket_factor_count,
                result.configuration.hs_fields.end());
    std::rotate(result.configuration.locations.begin(),
                result.configuration.locations.begin() + result.ket_factor_count,
                result.configuration.locations.end());
    return result;
}

#endif
