#ifndef PURE_STATE_DRIVEN_STAGE1_COMMON_H
#define PURE_STATE_DRIVEN_STAGE1_COMMON_H

#include "kitaevChain.h"
#include "pure_projector_driven.h"
#include "pure_projector_observables.h"
#include "pure_projector_protocol.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace driven_stage1 {

constexpr double kPi = 3.141592653589793238462643383279502884;

inline void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

struct ModelParameters {
    int L = 4;
    int boundary = 0;
    int hs = 0;
    double dt = 0.1;
    double t = 1.0;
    double delta = 1.0;
    double mu = 0.0;
    bool analytic_hs_factors = false;
};

inline MatType kineticGenerator(const ModelParameters &p, double t,
                                double delta, double mu) {
    auto component = [&](double d, double m) {
        SpinlessTvChainUtils config(
            p.L, 1.0, 0.0, 2, p.boundary, d, m, p.hs);
        MatType generator = MatType::Zero(2 * p.L, 2 * p.L);
        config.KineticGenerator(generator);
        return generator;
    };
    const MatType base = component(0.0, 0.0);
    return t * base + delta * (component(1.0, 0.0) - base) +
           mu * (component(0.0, 1.0) - base);
}

inline MatType exponential(MatType generator, double scale = 1.0) {
    return expm(generator, scale);
}

inline MatType localHsGenerator(const ModelParameters &p, double interaction,
                                int layer, int auxiliary, int sigma) {
    SpinlessTvChainUtils config(
        p.L, p.dt, interaction, 2, p.boundary, p.delta, p.mu, p.hs);
    MatType generator = MatType::Zero(2 * p.L, 2 * p.L);
    const double lambda = std::acosh(std::exp(0.5 * interaction * p.dt));
    int a, b, c, d;
    config.aux2MajoranaIdx(auxiliary, 0, layer, a, b);
    config.aux2MajoranaIdx(auxiliary, 1, layer, c, d);
    const DataType value(0.0, lambda * sigma);
    if (p.hs == 0) {
        generator(a, b) = value;
        generator(b, a) = -value;
        generator(c, d) = value;
        generator(d, c) = -value;
    } else {
        generator(a, c) = value;
        generator(c, a) = -value;
        generator(b, d) = -value;
        generator(d, b) = value;
    }
    return generator;
}

inline MatType localHsFactor(const ModelParameters &p, double interaction,
                             int layer, int auxiliary, int sigma) {
    const MatType generator=localHsGenerator(p,interaction,layer,auxiliary,sigma);
    if(!p.analytic_hs_factors)return exponential(generator);
    // Each HS generator consists of two disjoint 2 x 2 blocks with H^2=lambda^2 I.
    // Keep the same HS decomposition and signs, evaluating its exponential exactly.
    const double lambda=std::acosh(std::exp(0.5*interaction*p.dt));
    MatType factor=MatType::Identity(2*p.L,2*p.L);
    if(lambda==0)return factor;
    const double scale=std::sinh(lambda)/lambda,diagonal=std::cosh(lambda);
    for(int c=0;c<generator.cols();++c)for(int r=0;r<generator.rows();++r)
        if(generator(r,c)!=DataType(0,0)){
            factor(r,c)=scale*generator(r,c);factor(r,r)=diagonal;
        }
    return factor;
}

inline GaussianTrialState makeInitialState(
        const ModelParameters &p, double trialMu = 0.0,
        double edgeSplitting = 0.0, int expectedParity = 0) {
    MatType h = kineticGenerator(p, 1.0, 1.0, trialMu);
    if (edgeSplitting != 0.0) {
        SpinlessTvChainUtils coordinates(
            p.L, p.dt, 0.0, 2, p.boundary, 1.0, trialMu, p.hs);
        const int left = coordinates.majoranaCoord2Idx(0, 1);
        const int right = coordinates.majoranaCoord2Idx(p.L - 1, 1);
        const DataType value(0.0, edgeSplitting);
        h(left, right) += value;
        h(right, left) -= value;
    }
    GaussianTrialState state = GaussianTrialState::fromMajoranaHamiltonian(h);
    if (expectedParity != 0 && state.fermionParity() != expectedParity)
        throw std::invalid_argument(
            "Phi_0 parity does not match the explicit parity policy");
    return state;
}

// Exact occupation product state, in gamma=(c+c^dagger, i(c-c^dagger))
// coordinates. No finite trial interaction or auxiliary-field ordering is used.
inline GaussianTrialState makeCdwInitialState(
        const ModelParameters &p, int occupiedSublattice = 0,
        int expectedParity = 0) {
    require(p.L >= 2 && p.L % 2 == 0, "CDW initial state requires even L");
    require(occupiedSublattice == 0 || occupiedSublattice == 1,
            "invalid CDW occupied sublattice");
    MatType phi = MatType::Zero(2 * p.L, p.L);
    const double a = 1.0 / std::sqrt(2.0);
    for (int j = 0; j < p.L; ++j) {
        phi(j,j) = a;
        phi(p.L+j,j) = DataType(0, j % 2 == occupiedSublattice ? -a : a);
    }
    auto state = GaussianTrialState::fromPhi(phi);
    require(expectedParity == 0 || state.fermionParity() == expectedParity,
            "CDW Phi_0 parity does not match the explicit parity policy");
    return state;
}

inline DataType staggeredDensity(const MatType &green) {
    const int L = int(green.rows()/2);
    DataType result = 0;
    for (int j = 0; j < L; ++j)
        result += (j % 2 ? -1.0 : 1.0) * DataType(0,0.5) * green(j,L+j);
    return result / double(L);
}

inline PureFastConfiguration makeProtocolBranch(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        PureBranch branch, std::mt19937_64 &rng) {
    require(protocol.slices() > 0 &&
                std::abs(protocol.deltaTau() - p.dt) < 1e-12,
            "protocol Delta_tau does not match model dt");
    const MatType kinetic = kineticGenerator(p, p.t, p.delta, p.mu);
    const MatType halfKinetic = exponential(kinetic, -0.5 * p.dt);
    const auto counts = pureProjectorCheckerboardBondCounts(p.L, p.boundary);
    std::uniform_int_distribution<int> bit(0, 1);
    PureFastConfiguration result;
    int factor = 0;
    auto push = [&](const MatType &matrix, int field, int slice, int layer,
                    int auxiliary, const std::string &label) {
        result.slices.emplace_back(matrix, 1.0, label);
        result.hs_fields.push_back(field);
        result.locations.push_back(
            {branch, slice, factor++, layer, auxiliary});
    };
    // The validated Majorana-to-Fock representation is order reversing:
    // flattened one-particle factors B_0,...,B_n represent the dense Fock
    // product U_0...U_n.  Store late imaginary-time slices first so that the
    // represented ket is the time-ordered U(tau_f,0)=U_{n}...U_0.
    for (int slice = protocol.slices() - 1; slice >= 0; --slice) {
        const double interaction = protocol.midpointValue(slice);
        push(halfKinetic, 0, slice, -1, -1,
             "tau:" + std::to_string(slice) + ":K/2:first");
        for (int layer = 0; layer < 2; ++layer) {
            const int count = layer == 0 ? counts.first : counts.second;
            for (int auxiliary = 0; auxiliary < count; ++auxiliary) {
                const int sigma = bit(rng) ? 1 : -1;
                push(localHsFactor(p, interaction, layer, auxiliary, sigma),
                     sigma, slice, layer, auxiliary,
                     "tau:" + std::to_string(slice) + ":V" +
                         std::to_string(layer) + ":" +
                         std::to_string(auxiliary));
            }
        }
        push(halfKinetic, 0, slice, -1, -1,
             "tau:" + std::to_string(slice) + ":K/2:last");
    }
    return result;
}

inline PureFastConfiguration makeStaticBranch(
        const ModelParameters &p, double interaction, int slices,
        PureBranch branch, std::mt19937_64 &rng) {
    require(slices > 0, "static branch requires positive slice count");
    const MatType kinetic = kineticGenerator(p, p.t, p.delta, p.mu);
    const MatType halfKinetic = exponential(kinetic, -0.5 * p.dt);
    const auto counts = pureProjectorCheckerboardBondCounts(p.L, p.boundary);
    std::uniform_int_distribution<int> bit(0, 1);
    PureFastConfiguration result;
    int factor = 0;
    auto push = [&](const MatType &matrix, int field, int slice, int layer,
                    int auxiliary, const std::string &label) {
        result.slices.emplace_back(matrix, 1.0, label);
        result.hs_fields.push_back(field);
        result.locations.push_back(
            {branch, slice, factor++, layer, auxiliary});
    };
    for (int slice = 0; slice < slices; ++slice) {
        push(halfKinetic, 0, slice, -1, -1,
             "tau:" + std::to_string(slice) + ":K/2:first");
        for (int layer = 0; layer < 2; ++layer) {
            const int count = layer == 0 ? counts.first : counts.second;
            for (int auxiliary = 0; auxiliary < count; ++auxiliary) {
                const int sigma = bit(rng) ? 1 : -1;
                push(localHsFactor(p, interaction, layer, auxiliary, sigma),
                     sigma, slice, layer, auxiliary,
                     "tau:" + std::to_string(slice) + ":V" +
                         std::to_string(layer) + ":" +
                         std::to_string(auxiliary));
            }
        }
        push(halfKinetic, 0, slice, -1, -1,
             "tau:" + std::to_string(slice) + ":K/2:last");
    }
    return result;
}

inline PureDrivenContourResult makeDrivenContour(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        std::mt19937_64 &rng) {
    PureFastConfiguration ket =
        makeProtocolBranch(p, protocol, PureBranch::Ket, rng);
    PureFastConfiguration bra =
        makeProtocolBranch(p, protocol, PureBranch::Bra, rng);
    return pureProjectorIndependentBraKetContour(
        std::move(ket), std::move(bra));
}

inline PureDrivenContourResult makeMirroredDrivenContour(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        std::mt19937_64 &rng) {
    PureFastConfiguration ket =
        makeProtocolBranch(p, protocol, PureBranch::Ket, rng);
    PureFastConfiguration bra = ket;
    for (PureSliceLocation &location : bra.locations)
        location.branch = PureBranch::Bra;
    PureDrivenContourResult result = pureProjectorIndependentBraKetContour(
        std::move(ket), std::move(bra));
    result.bra_ket_hs_independent = false;
    result.contour_order =
        "bra_strict_reverse_adjoint_then_same_hs_ket_canonical_action";
    return result;
}

inline PureDrivenContourResult makeStaticIndependentContour(
        const ModelParameters &p, double interaction, int slices,
        std::mt19937_64 &rng) {
    PureFastConfiguration ket =
        makeStaticBranch(p, interaction, slices, PureBranch::Ket, rng);
    PureFastConfiguration bra =
        makeStaticBranch(p, interaction, slices, PureBranch::Bra, rng);
    // Keep the independently built static branches, but assemble their full
    // contour with the same corrected block order as the driven path.  This
    // helper is a regression oracle; duplicating the pre-fix concatenation
    // here would compare different contours rather than static versus driven
    // construction of the same contour.
    return pureProjectorIndependentBraKetContour(
        std::move(ket), std::move(bra));
}

inline std::vector<int> proposalIndices(const PureFastConfiguration &configuration) {
    std::vector<int> result;
    for (int index = 0; index < int(configuration.locations.size()); ++index)
        if (configuration.locations[index].aux >= 0) result.push_back(index);
    return result;
}

inline PureFastProposal flipProposal(
        const ModelParameters &p, const PureImaginaryTimeProtocol &protocol,
        const PureFastConfiguration &configuration, int index, double uniform) {
    const PureSliceLocation &location = configuration.locations[index];
    const int newHs = -configuration.hs_fields[index];
    MatType factor = localHsFactor(
        p, protocol.midpointValue(location.slice), location.bond,
        location.aux, newHs);
    if (location.branch == PureBranch::Bra) factor.adjointInPlace();
    PureFastProposal proposal;
    proposal.index = index;
    proposal.new_hs = newHs;
    proposal.new_factor = std::move(factor);
    proposal.new_eta = 1.0;
    proposal.uniform = uniform;
    if(p.hs==1)proposal.positive_ratio_log_bound=
        2*std::acosh(std::exp(.5*protocol.midpointValue(location.slice)*p.dt));
    return proposal;
}

// Run-scoped, compressed cache for fixed-interaction density proposals. The
// bound model/protocol cannot change underneath a key; ramps retain the
// original path so cache memory never grows with the number of time slices.
class PureDrivenProposalFactory {
    ModelParameters model_;
    PureImaginaryTimeProtocol protocol_;
    bool enabled_;
    std::map<std::tuple<int,int,int,int>,PureFastProposal> cache_;
public:
    PureDrivenProposalFactory(const ModelParameters &p,const PureImaginaryTimeProtocol &protocol,bool enabled=true)
        :model_(p),protocol_(protocol),enabled_(enabled&&p.hs==1&&p.analytic_hs_factors&&
            protocol.kind()!=PureImaginaryTimeProtocolKind::LinearRamp){}
    bool enabled()const{return enabled_;}
    std::size_t cachedFactorCount()const{return cache_.size();}
    PureFastProposal proposal(const PureFastConfiguration &cfg,int index,double uniform,bool materialize=true){
        if(!enabled_)return flipProposal(model_,protocol_,cfg,index,uniform);
        const auto &location=cfg.locations.at(index);
        const auto key=std::make_tuple(int(location.branch),location.bond,location.aux,-cfg.hs_fields.at(index));
        auto found=cache_.find(key);
        if(found==cache_.end()){
            auto value=flipProposal(model_,protocol_,cfg,index,uniform);
            value.cached_factor=std::make_shared<const PureFactorMatrix>(value.new_factor);
            value.new_factor.resize(0,0);
            found=cache_.emplace(key,std::move(value)).first;
        }
        auto result=found->second;result.index=index;result.uniform=uniform;
        if(materialize)result.new_factor=*result.cached_factor;
        else result.immutable_factor_only=true;
        return result;
    }
};

inline bool positiveHs1RatioContract(const ModelParameters &p,const GaussianTrialState &trial) {
    if(!(p.hs==1&&p.boundary==0&&p.L%2==0&&p.t==1&&p.delta==1&&p.mu==0&&trial.fermionParity()==-1))return false;
    const MatType expected=.5*kineticGenerator(p,1,1,0);
    return (trial.G_T-expected).norm()/std::max(1.,expected.norm())<1e-12;
}

inline DataType energy(const ModelParameters &p, double interaction,
                       const MatType &green) {
    auto value = [&](double V, double delta, double mu) {
        SpinlessTvChainUtils config(
            p.L, p.dt, V, 2, p.boundary, delta, mu, p.hs);
        return config.energyFromGreensFunc(green);
    };
    const DataType base = value(0.0, 0.0, 0.0);
    return p.t * base + p.delta * (value(0.0, 1.0, 0.0) - base) +
           p.mu * (value(0.0, 0.0, 1.0) - base) +
           (value(interaction, 0.0, 0.0) - base);
}

struct Observables {
    double energy = 0.0;
    double S_pi = 0.0;
    double S_pi_dq = 0.0;
    double R_CDW = 0.0;
    double parity = 0.0;
};

inline Observables observe(const ModelParameters &p, double interaction,
                           const MatType &green,int conservedParity=0) {
    Observables result;
    result.energy = energy(p, interaction, green).real();
    result.S_pi = pureProjectorStructureFactor(green, p.L, kPi).real();
    result.S_pi_dq = pureProjectorStructureFactor(
        green, p.L, kPi - 2.0 * kPi / p.L).real();
    result.R_CDW = 1.0 - result.S_pi_dq / result.S_pi;
    if(conservedParity!=0){
        require(conservedParity==1||conservedParity==-1,"invalid conserved fermion parity");
        // Every kinetic/HS factor is even in Majoranas. For a parity eigenstate
        // at both endpoints, <bra|P|ket>/<bra|ket> equals that initial parity.
        // This is distinct from the fluctuating physical Monte Carlo weight sign.
        result.parity=conservedParity;
        return result;
    }
    const PurePhysicalParityResult parity = pureProjectorPhysicalParity(green);
    require(parity.ok(), "fermion parity measurement is untrusted");
    result.parity = parity.physical_parity;
    return result;
}

inline std::uint64_t hashRng(const std::mt19937_64 &rng) {
    std::ostringstream state;
    state << rng;
    std::uint64_t hash = 1469598103934665603ULL;
    for (unsigned char byte : state.str()) {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    return hash;
}

inline double relativeError(const MatType &left, const MatType &right) {
    return (left - right).norm() / std::max(1.0, right.norm());
}

}  // namespace driven_stage1

#endif
