#ifndef KITAEV_SSH_PHONON_PROJECTOR_H
#define KITAEV_SSH_PHONON_PROJECTOR_H

#include "phonon_path.h"
#include "path_profile.h"
#include "free_projector.h"
#include "pure_projector_fast.h" // existing stabilized, signed Pfaffian reference
#include "bipartite_spectrum.h"
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <functional>
#include <memory>
#include <sstream>
#include <type_traits>
#include <utility>
namespace kitaev_ssh {
#ifdef KITAEV_SSH_TEST_FAILURE_INJECTION
// Test binaries only; no environment switch or injection API in production.
inline std::function<void(const char*,bool)> testNumericalFailure;
#endif

// P2 defaults explicitly differ from the frozen P1 baseline's Jz=0.6.
inline Couplings p2Couplings(double g) { return {1.,1.,.5,g}; }
// A generator keeps the original callable interface while allowing model
// implementations to fill a caller-owned matrix. Generic test lambdas use
// the value-returning fallback; the honeycomb path uses the output-buffer form.
struct MajoranaGenerator {
    using ValueFunction=std::function<MatType(const std::vector<double>&)>;
    using IntoFunction=std::function<void(const std::vector<double>&,MatType&)>;
    ValueFunction value;
    IntoFunction into;
    bool validated_dense=false;
    bool bipartite_real=false;
    MajoranaGenerator()=default;
    MajoranaGenerator(ValueFunction f):value(std::move(f)) {}
    MajoranaGenerator(ValueFunction f,IntoFunction output,bool dense=false,bool bipartite=false)
        :value(std::move(f)),into(std::move(output)),validated_dense(dense),bipartite_real(bipartite) {}
    template<class F,
             typename std::enable_if<!std::is_same<typename std::decay<F>::type,
                                                    MajoranaGenerator>::value,int>::type=0>
    MajoranaGenerator(F &&f):value(std::forward<F>(f)) {}
    MatType operator()(const std::vector<double> &x) const {
        if(!value) throw std::logic_error("empty Majorana generator");
        return value(x);
    }
    void generateInto(const std::vector<double> &x,MatType &out) const {
        if(into) into(x,out);
        else out=(*this)(x);
    }
    bool producesValidatedDense() const { return validated_dense; }
};
inline MajoranaGenerator honeycombGenerator(const KitaevSSHModel &model) {
    // A3: reuse the immutable model's lattice/couplings/gauge; only x varies.
    MajoranaGenerator::ValueFunction value=[model](const std::vector<double> &x) {
        return model.frameworkGeneratorWith(x);
    };
    MajoranaGenerator::IntoFunction into=[model](const std::vector<double> &x,MatType &out) {
        model.frameworkGeneratorWith(x,out);
    };
    return MajoranaGenerator(std::move(value),std::move(into),true,true);
}
inline PureProjectorSlice phononNodeSlice(const MatType &matrix,const MajoranaGenerator &generator,
                                          std::string label) {
    if(generator.producesValidatedDense())
        return PureProjectorSlice::fromValidatedDense(matrix,1.,std::move(label));
    return PureProjectorSlice(matrix,1.,std::move(label));
}
inline PureProjectorSlice phononNodeSlice(MatType &&matrix,const MajoranaGenerator &generator,
                                          std::string label) {
    if(generator.producesValidatedDense())
        return PureProjectorSlice::fromValidatedDense(std::move(matrix),1.,std::move(label));
    return PureProjectorSlice(matrix,1.,std::move(label));
}
// Stored in action order: F(X0), F(X1), F(X1), F(X2), ..., F(XM).
// Link l (1 based) uses factors 2*l-2,2*l-1. Center cut = M factors.
inline int factorNode(int factor,int intervals) {
    if (factor<0 || factor>=2*intervals) throw std::out_of_range("phonon factor index");
    return (factor+1)/2;
}
inline std::vector<int> nodeFactors(int node,int intervals) {
    if (node<0 || node>intervals) throw std::out_of_range("phonon node index");
    if (node==0) return {0};
    if (node==intervals) return {2*intervals-1};
    return {2*node-1,2*node};
}
struct MajoranaGenerationWorkspace {
    MatType generator;
    MatType half_step;
    BipartiteSpectrum bipartite_spectrum;
};
inline void nodeHalfStepInto(const PhononPath &path,int node,int dimension,
                             const MajoranaGenerator &generator,MatType &out,
                             MajoranaGenerationWorkspace &workspace) {
    bool half_step_invalid=false;
    profilePath("gen_generator",[&]{ generator.generateInto(path.node(node),workspace.generator); return true; });
    if(!generator.producesValidatedDense())
        profilePath("gen_validate",[&]{ KitaevSSHModel::validateGenerator(workspace.generator,dimension); return true; });
    profilePath("gen_expm",[&]{ workspace.half_step=expm(workspace.generator,-path.dt()/2); return true; });
    profilePath("gen_hermiticity",[&]{
        if (!workspace.half_step.allFinite() ||
            (workspace.half_step-workspace.half_step.adjoint()).norm()>
                1e-12*std::max(1.,workspace.half_step.norm()))
            half_step_invalid=true;
        return true;
    });
    if (half_step_invalid) throw std::runtime_error("invalid Majorana half-step");
    out=workspace.half_step;
}
inline MatType nodeHalfStep(const PhononPath &path,int node,int dimension,
                            const MajoranaGenerator &generator) {
    MajoranaGenerationWorkspace workspace;
    MatType out;
    nodeHalfStepInto(path,node,dimension,generator,out,workspace);
    return out;
}
inline std::vector<PureProjectorSlice> pathFactors(const PhononPath &path,int dimension,
                                                   const MajoranaGenerator &generator) {
    std::vector<PureProjectorSlice> factors(2*path.intervals());
    MajoranaGenerationWorkspace workspace;
    MatType half_step;
    for (int j=0;j<path.nodes();++j) {
        nodeHalfStepInto(path,j,dimension,generator,half_step,workspace);
        const PureProjectorSlice node=phononNodeSlice(half_step,generator,
                                                       "phonon_node:"+std::to_string(j));
        for (int f:nodeFactors(j,path.intervals())) factors[f]=node;
    }
    return factors;
}

// One spectrum per node (reduced real bipartite SVD for eligible large
// honeycomb generators, full Hermitian reference otherwise) yields the
// half-step h=exp(-dt K/2), the merged full node factor B=exp(-dt K) and the
// Cayley Green slices C=2(I+B)^{-1}-I=U tanh(a eps/2) U^dag (a=dt for merged
// interior nodes, a=dt/2 for the half-step endpoints). The canonical half
// factors still feed the endpoint center Green; the
// merged weight order plus an aligned Cayley array feed the signed weight.
struct PreparedOwnedPath {
    std::vector<PureProjectorSlice> factors;
    std::vector<PureProjectorSlice> weightOrder;
    std::vector<MatType> cayley;
};
struct PreparedOwnedNode {
    MatType half_step;
    MatType merged_factor;
};
inline void assembleOwnedHalfStep(PhononNodeSpectrum &spectral,double dt,
                                   MajoranaGenerationWorkspace &generation) {
    profilePath("gen_half_assemble",[&]{
        spectral.exponential(dt/2,generation.half_step); return true;
    });
    const bool valid=profilePath("gen_hermiticity",[&]{
        return generation.half_step.allFinite() &&
            (generation.half_step-generation.half_step.adjoint()).norm()<=
                1e-12*std::max(1.,generation.half_step.norm());
    });
    if(!valid) throw std::runtime_error("invalid Majorana half-step");
}
// Keep the deferred local path on the same spectral construction as the
// existing full prepared path. This is deliberately separate from the eager
// local h*h path below; the two branches have different floating-point
// arithmetic and are not interchangeable merely because they are equivalent
// mathematically.
inline PreparedOwnedNode prepareOwnedNode(const PhononPath &path,int node,int dimension,
                                          const MajoranaGenerator &generator,
                                          MajoranaGenerationWorkspace &generation,
                                          bool wantHalfSteps=true) {
    PhononNodeSpectrum spectral;
    const bool interior=node>0 && node<path.intervals();
    bool haveHalf=wantHalfSteps || !interior;
    profilePath("gen_generator",[&]{ generator.generateInto(path.node(node),generation.generator); return true; });
    if(!generator.producesValidatedDense())
        profilePath("gen_validate",[&]{ KitaevSSHModel::validateGenerator(generation.generator,dimension); return true; });
    profilePath("gen_expm",[&]{
        if(!spectral.compute(generation.generator,generator.bipartite_real,generation.bipartite_spectrum))
            throw std::runtime_error("Majorana generator eigendecomposition failed");
        if(haveHalf) assembleOwnedHalfStep(spectral,path.dt(),generation);
        return true;
    });
    PreparedOwnedNode result;
    if(interior) {
        spectral.exponential(path.dt(),result.merged_factor);
        if(!result.merged_factor.allFinite()) {
            if(!haveHalf) {assembleOwnedHalfStep(spectral,path.dt(),generation);haveHalf=true;}
            result.merged_factor=(generation.half_step*generation.half_step).eval();
        }
        if(!result.merged_factor.allFinite() ||
           (result.merged_factor-result.merged_factor.adjoint()).norm()>
               1e-12*std::max(1.,result.merged_factor.norm()))
            throw std::runtime_error("invalid merged Majorana node factor");
    } else result.merged_factor=generation.half_step;
    if(wantHalfSteps) result.half_step=generation.half_step;
    return result;
}
inline PreparedOwnedPath prepareOwnedPath(const PhononPath &path,int dimension,
                                          const MajoranaGenerator &generator,
                                          bool wantCayley=true,
                                          MajoranaGenerationWorkspace *generationWorkspace=nullptr,
                                          bool wantHalfSteps=true) {
    PreparedOwnedPath out;
    const int M=path.intervals();
    if(wantHalfSteps) out.factors.resize(2*M);
    std::vector<MatType> weightB(path.nodes()),weightC(path.nodes());
    MajoranaGenerationWorkspace localGeneration;
    MajoranaGenerationWorkspace &generation=generationWorkspace?*generationWorkspace:localGeneration;
    for(int j=0;j<path.nodes();++j){
        PhononNodeSpectrum spectral;
        const bool interior=(j>0&&j<M);
        const bool haveHalf=wantHalfSteps || !interior;
        profilePath("gen_generator",[&]{ generator.generateInto(path.node(j),generation.generator); return true; });
        if(!generator.producesValidatedDense())
            profilePath("gen_validate",[&]{ KitaevSSHModel::validateGenerator(generation.generator,dimension); return true; });
        profilePath("gen_expm",[&]{
            if(!spectral.compute(generation.generator,generator.bipartite_real,generation.bipartite_spectrum))
                throw std::runtime_error("Majorana generator eigendecomposition failed");
            if(haveHalf) assembleOwnedHalfStep(spectral,path.dt(),generation);
            return true;
        });
        if(wantHalfSteps) {
            const PureProjectorSlice node=phononNodeSlice(generation.half_step,generator,
                                                           "phonon_node:"+std::to_string(j));
            for(int f:nodeFactors(j,M)) out.factors[f]=node;
        }
        const double a=interior?path.dt():path.dt()/2;
        // Spectral condition of I+B(a); the FullPivLU reference rejects when
        // its pivot ratio is at/below ~D*eps. Use a conservative factor of 2:
        // only bypass the LU solve when the node is plainly nonsingular, else
        // emit an empty sentinel so the weight routine falls back to the
        // authoritative FullPivLU path (identical decision to the baseline).
        bool spectrallySafe=false;
        if(wantCayley){
            spectrallySafe=spectral.cayleySafe(a,dimension);
        }
        profilePath("gen_full_assemble",[&]{
            if(interior){
                MatType B;spectral.exponential(path.dt(),B);
                if(!B.allFinite()) {
                    // Preserve the original overflow fallback, materializing
                    // and validating this node's half-step only if needed.
                    if(!haveHalf) assembleOwnedHalfStep(spectral,path.dt(),generation);
                    B=(generation.half_step*generation.half_step).eval();
                }
                weightB[j]=std::move(B);
                if(wantCayley&&spectrallySafe){
                    spectral.cayley(path.dt(),weightC[j]);
                }
            } else {
                weightB[j]=generation.half_step;
                if(wantCayley&&spectrallySafe){
                    spectral.cayley(path.dt()/2,weightC[j]);
                }
            }
            return true;
        });
    }    out.weightOrder.reserve(path.nodes());out.cayley.reserve(path.nodes());    for(int j=M;j>=0;--j){
        const bool interior=(j>0&&j<M);
        if(interior && (!weightB[j].allFinite() ||
            (weightB[j]-weightB[j].adjoint()).norm()>1e-12*std::max(1.,weightB[j].norm())))
            throw std::runtime_error("invalid merged Majorana node factor");
        out.weightOrder.push_back(phononNodeSlice(std::move(weightB[j]),generator,
                                                  "owned_node_weight:"+std::to_string(j)));
        out.cayley.push_back(std::move(weightC[j]));
    }
    return out;
}
struct PathWeight {
    double log_phonon=0., log_abs_matter=0.;
    DataType matter_phase=0.;
    MatType center_green;
    // Task B validity flags. Eager snapshots have all three true; deferred
    // magnitude snapshots have only magnitude_valid until ensureMeasured()
    // recomputes the full signed weight and center Green. Readers of the
    // phase or Green MUST check the corresponding flag; no default phase=1.
    bool magnitude_valid=false, phase_valid=false, green_valid=false;
    double logAbsTotal() const { return log_phonon+log_abs_matter; }
};

// A deterministic identity for the immutable snapshot that owns a PathWeight.
// This is an audit/association hash, not a cryptographic digest: it binds all
// path metadata, trial shifts and coordinates so a valid phase/Green cache is
// visibly tied to the configuration that produced it.
inline std::uint64_t phononPathConfigurationHash(const PhononPath &path) {
    std::uint64_t hash=1469598103934665603ULL;
    auto mix=[&](std::uint64_t value) {
        hash^=value; hash*=1099511628211ULL;
    };
    auto mixInt=[&](int value) { mix(std::uint64_t(std::int64_t(value))); };
    auto mixDouble=[&](double value) {
        std::uint64_t bits=0;
        static_assert(sizeof(bits)==sizeof(value),"unexpected double size");
        std::memcpy(&bits,&value,sizeof(bits));
        mix(bits);
    };
    mixInt(path.bonds()); mixInt(path.intervals()); mixDouble(path.dt());
    mixDouble(path.parameters().mass); mixDouble(path.parameters().omega);
    for(double shift:path.trialShifts()) mixDouble(shift);
    for(int node=0;node<path.nodes();++node)
        for(double value:path.node(node)) mixDouble(value);
    return hash;
}

struct PhononWeightOptions {
    bool singular_values_only=false;
    bool merge_node_half_steps=false;
    // EXPERIMENTAL, default false: compute only the certified magnitude weight
    // on updates and defer the Pfaffian phase + center Green to measurement
    // strides (exact for the magnitude-target walk; see PHASE_ON_DEMAND.md).
    bool defer_phase_measurement=false;
    // Eager-local-only opt-in. Deferred and whole-bond proposals keep their
    // complete reference paths.
    bool reuse_signed_prefix_checkpoints=false;
    bool reuse_lu_cayley=false; // eager local merged factors only, default off
    // Signed phase is always measured. Energy/fermion consumers additionally
    // request the center Green; phonon-only runs may explicitly omit it.
    bool compute_measurement_green=true;
};
struct PhononCacheDiagnostics {
    std::uint64_t local_proposals=0;
    std::uint64_t local_node_rebuilds=0;
    std::uint64_t local_node_reuses=0;
    std::uint64_t local_merged_squares=0;
    std::uint64_t bond_refresh_full_rebuilds=0;
};
inline std::string weightFailure(const char *where,const PureProjectorWeightResult &r) {
    std::ostringstream out; out<<std::setprecision(17)<<where<<": "<<pureProjectorWeightStatusName(r.status)
        <<" slice="<<r.first_failing_slice<<" rank="<<r.overlap_rank<<" rcond="<<r.overlap_rcond
        <<" overlap_residual="<<r.overlap_residual<<" green_residual="<<r.green_residual;
    if(r.green_status!=PureProjectorStatus::invalid_dimension)
        out<<" green_status="<<pureProjectorStatusName(r.green_status)
           <<" green_rcond="<<r.green_overlap_rcond
           <<" green_solve_residual="<<r.green_solve_residual
           <<" green_skew_residual="<<r.green_skew_residual
           <<" green_diagonal_residual="<<r.green_diagonal_residual;
    return out.str();
}
inline PathWeight evaluatePath(const PhononPath &path,const GaussianTrialState &trial,
                              const std::vector<PureProjectorSlice> &factors,int block,
                              bool singularValuesOnly=false,
                              const std::vector<PureProjectorSlice> *ownedNodeWeightOrder=nullptr,
                              bool factorsPrevalidated=false,
                              const PreparedGreenTrial *preparedTrial=nullptr,
                              const std::vector<MatType> *preparedCayley=nullptr,
                              PureStableWeightWorkspace *weightWorkspace=nullptr,
                              const PureProjectorSignedPrefixCache *prefixCache=nullptr,
                              PureProjectorSignedPrefixCache *updatedPrefixCache=nullptr,
                              int requestedResumeIndex=-1,
                              PureProjectorLUCayleyCache *cayleyCache=nullptr,
                              bool computeGreen=true) {    if (block<=0 || factors.size()!=std::size_t(2*path.intervals())) throw std::invalid_argument("invalid path factors/block");    // Legacy weight order is the conjugate matter overlap for Hermitian
    // half-steps listed in physical action order. Reverse only the order for
    // this call (NOT scalar eta phases); independent complex Fock test fixes
    // the orientation. Keep the live factor list in chronological action order.
    if (!factorsPrevalidated) profilePath("factor_validation_and_temporary_copy",[&] {
        for (const auto &factor:factors) {
            const MatType b=factor.matrix;
            if (b.rows()!=trial.Phi.rows() || b.cols()!=trial.Phi.rows() || !b.allFinite() ||
                (b-b.adjoint()).norm()>1e-12*std::max(1.,b.norm()))
                throw std::invalid_argument("P2 expects Hermitian Majorana half-steps");
        }
        return true;
    });
    // The owned-node path is already in the required order. Keep it by
    // reference; copying every dense factor for every proposal is pure
    // overhead. The generic path still materializes its reversed order.
    std::vector<PureProjectorSlice> reversedOrder;
    profilePath("reverse_factor_copy",[&] {
        if (!ownedNodeWeightOrder) reversedOrder.assign(factors.rbegin(),factors.rend());
        return true;
    });
    const auto &weightOrder=ownedNodeWeightOrder ? *ownedNodeWeightOrder : reversedOrder;
    PureFastOptions weightOptions;
    weightOptions.reference_singular_values_only=singularValuesOnly;
    const auto matter=profilePath("physical_pfaffian",[&] {return pureProjectorStableReferenceWeight(trial,weightOrder,weightOptions,preparedTrial,preparedCayley,weightWorkspace,prefixCache,updatedPrefixCache,requestedResumeIndex,cayleyCache);});    if (!matter.ok()) throw std::runtime_error(weightFailure("phonon matter overlap",matter));    // General paths are NOT U^dagger U. Keep the Pfaffian phase above, even
    // though endpoint rebuilding below only computes |weight| and Green.
    if(!computeGreen) {
        PathWeight w{path.logWeight(),matter.log_abs_weight,matter.complex_phase,MatType{}};
        w.magnitude_valid=w.phase_valid=true; w.green_valid=false;
        return w;
    }
    const auto center=profilePath("center_green",[&] {return pureProjectorEndpointRebuild(trial,factors,path.intervals(),block,PureProjectorOptions{});});
    if (!center.ok()) {
        std::ostringstream out; out<<std::setprecision(17)<<"phonon center Green: "<<center.message
            <<" rank="<<center.overlap_rank<<" rcond="<<center.overlap_rcond
            <<" solve_residual="<<center.solve_residual<<" green_residual="<<center.green_residual;
        throw std::runtime_error(out.str());
    }
    if (std::abs(matter.log_abs_weight-center.log_abs_weight)>1e-9*std::max(1.,std::abs(matter.log_abs_weight)))
        throw std::runtime_error("signed overlap/endpoint magnitude mismatch");
    // beta_R=R^dagger c annihilates |R>, <L| beta_L^dagger=0. Therefore
    // G_phys = I-2 L(R^dagger L)^(-1)R^dagger, the adjoint of the core's
    // I-2 R(L^dagger R)^(-1)L^dagger. Equal bra/ket P1 tests cannot see this.
    PathWeight w{path.logWeight(),matter.log_abs_weight,matter.complex_phase,center.green.adjoint().eval()};
    w.magnitude_valid=w.phase_valid=w.green_valid=true;    return w;
}
// Task B: magnitude-only evaluation for deferred mode. Same factor gates, no
// Cayley/Green/Pfaffian/center-Green work; phase and Green are invalid until
// the full signed evaluation is run at measurement time.
inline PathWeight evaluatePathMagnitude(const PhononPath &path,const GaussianTrialState &trial,
                              const std::vector<PureProjectorSlice> &factors,int block,
                              const std::vector<PureProjectorSlice> *ownedNodeWeightOrder=nullptr,
                              PureStableWeightWorkspace *weightWorkspace=nullptr,
                              bool factorsPrevalidated=false) {
    // Empty canonical factors explicitly mean a certified merged-only path.
    // No generic/unmerged caller may bypass the canonical shape gate.
    const bool mergedOnly=factors.empty() && factorsPrevalidated && ownedNodeWeightOrder &&
        ownedNodeWeightOrder->size()==std::size_t(path.nodes());
    if (block<=0 || (!mergedOnly && factors.size()!=std::size_t(2*path.intervals())))
        throw std::invalid_argument("invalid path factors/block");
    // Owned factors were validated at construction, and the changed node is
    // validated before insertion. The propagation below still checks every
    // factor's shape and finite entries; generic callers keep this full gate.
    if(!factorsPrevalidated) {
        for (const auto &factor:factors) {
            const MatType b=factor.matrix;
            if (b.rows()!=trial.Phi.rows() || b.cols()!=trial.Phi.rows() || !b.allFinite() ||
                (b-b.adjoint()).norm()>1e-12*std::max(1.,b.norm()))
                throw std::invalid_argument("P2 expects Hermitian Majorana half-steps");
        }
    }
    std::vector<PureProjectorSlice> reversedOrder;
    if(!ownedNodeWeightOrder) reversedOrder.assign(factors.rbegin(),factors.rend());
    const auto &weightOrder=ownedNodeWeightOrder ? *ownedNodeWeightOrder : reversedOrder;
    PureFastOptions weightOptions;
    weightOptions.reference_singular_values_only=true;
    weightOptions.real_magnitude=factorsPrevalidated;
    const auto matter=profilePath("physical_pfaffian",[&] {
        return pureProjectorMagnitudeWeight(trial,weightOrder,weightOptions,weightWorkspace);});
    if (!matter.ok()) throw std::runtime_error(weightFailure("phonon matter magnitude",matter));
    PathWeight w{path.logWeight(),matter.log_abs_weight,DataType(0),MatType{}};
    w.magnitude_valid=true;w.phase_valid=w.green_valid=false;    return w;
}
inline std::vector<PureProjectorSlice> reverseFactors(
    const std::vector<PureProjectorSlice> &factors) {
    return std::vector<PureProjectorSlice>(factors.rbegin(),factors.rend());
}
class CoupledPhononPath {
public:
    struct Snapshot {
        PhononPath path;
        // Empty only for deferred merged proposals: weight_order owns all
        // certified propagators. ensureMeasured() rebuilds canonical factors
        // before signed/Green evaluation; eager and recovery paths retain them.
        std::vector<PureProjectorSlice> factors;
        // Immutable order used by the matter-weight routine. A merged eager
        // local cache uses h*h, while deferred/A2 uses the existing spectral
        // construction; unmerged paths use the ordinary reversed factors.
        std::vector<PureProjectorSlice> weight_order;
        bool weight_order_is_local=false;
        std::shared_ptr<const PureProjectorSignedPrefixCache> signed_prefix_cache;
        std::shared_ptr<const PureProjectorLUCayleyCache> lu_cayley_cache;
        PathWeight weight;
        std::uint64_t configuration_hash=0;
    };
    class Proposal {
    public:
        Proposal(const Proposal&)=delete;
        Proposal &operator=(const Proposal&)=delete;
        Proposal(Proposal&&)=default;
        Proposal &operator=(Proposal&&)=default;
        const Snapshot &candidate() const {
            if (!next_) throw std::logic_error("consumed phonon proposal");
            return *next_;
        }
        double localPhononLogRatio() const { candidate(); return local_phonon_; }
        double logAbsRatio() const { return candidate().weight.logAbsTotal()-base_->weight.logAbsTotal(); }
        DataType phaseRatio() const {
            const Snapshot &c=candidate();
            if(!c.weight.phase_valid||!base_->weight.phase_valid)
                throw std::logic_error("phase requested on a deferred magnitude-only snapshot");
            return c.weight.matter_phase/base_->weight.matter_phase;
        }    private:
        friend class CoupledPhononPath;
        std::shared_ptr<const Snapshot> base_,next_;
        double local_phonon_;
        Proposal(std::shared_ptr<const Snapshot> base,std::shared_ptr<const Snapshot> next,double local)
            : base_(std::move(base)),next_(std::move(next)),local_phonon_(local) {}
    };
    CoupledPhononPath(PhononPath path,GaussianTrialState trial,MajoranaGenerator generator,int block=4,
                     PhononWeightOptions options={},bool uncachedInitial=false)
        : trial_(std::move(trial)),generator_(std::move(generator)),block_(block),options_(options) {
        testFailure("initialization",uncachedInitial);
        state_=std::make_shared<const Snapshot>(uncachedInitial?uncachedReference(path):fullReference(path));
    }
    CoupledPhononPath(const CoupledPhononPath&)=delete;
    CoupledPhononPath &operator=(const CoupledPhononPath&)=delete;
    const Snapshot &state() const { return *state_; }
    std::uint64_t configurationHash() const { return state_->configuration_hash; }
    const PhononCacheDiagnostics &cacheDiagnostics() const { return cacheDiagnostics_; }
    // Task B: full signed weight + center Green for the current configuration,
    // computed at a measurement stride. Consumes no RNG; an immutable snapshot
    // is replaced atomically only after the reference evaluation succeeds.
    const Snapshot &ensureMeasured(bool reference=false) {
        if(state_->weight.phase_valid &&
           (!options_.compute_measurement_green || state_->weight.green_valid)) return *state_;
        testFailure("measurement",reference);
        if(reference) {
            auto next=std::make_shared<const Snapshot>(uncachedReference(state_->path,true,options_.compute_measurement_green));
            state_=std::move(next);
            return *state_;
        }
        if(!options_.defer_phase_measurement) {
            const auto center=profilePath("center_green",[&] {
                return pureProjectorEndpointRebuild(trial_,state_->factors,
                    state_->path.intervals(),block_,PureProjectorOptions{});
            });
            if(!center.ok()) {
                std::ostringstream out;
                out<<std::setprecision(17)<<"phonon center Green: "<<center.message
                   <<" rank="<<center.overlap_rank<<" rcond="<<center.overlap_rcond
                   <<" solve_residual="<<center.solve_residual
                   <<" green_residual="<<center.green_residual;
                throw std::runtime_error(out.str());
            }
            if(std::abs(state_->weight.log_abs_matter-center.log_abs_weight)>
               1e-9*std::max(1.,std::abs(state_->weight.log_abs_matter)))
                throw std::runtime_error("signed overlap/endpoint magnitude mismatch");
            auto weight=state_->weight;
            weight.center_green=center.green.adjoint().eval();
            weight.green_valid=true;
            state_=std::make_shared<const Snapshot>(makeSnapshot(
                state_->path,state_->factors,std::move(weight),state_->weight_order));
            return *state_;
        }
        if(options_.merge_node_half_steps && !a2Disabled()) {
            auto prepared=prepareOwnedPath(state_->path,trial_.Phi.rows(),generator_,true,
                                           &generationWorkspace_);
            auto factors=std::move(prepared.factors);
            auto weight=evaluatePath(state_->path,trial_,factors,block_,
                options_.singular_values_only,&prepared.weightOrder,true,preparedTrial(),&prepared.cayley,
                &weightWorkspace_,nullptr,nullptr,-1,nullptr,
                options_.compute_measurement_green);
            state_=std::make_shared<const Snapshot>(makeSnapshot(
                state_->path,std::move(factors),std::move(weight),
                std::move(prepared.weightOrder)));
        } else {
            // The deferred candidate already stores the canonical half-step
            // factors. Reuse them for the full signed measurement path.
            auto weight=evaluatePath(state_->path,trial_,state_->factors,block_,
                options_.singular_values_only,&state_->weight_order,true,preparedTrial(),nullptr,
                &weightWorkspace_,nullptr,nullptr,-1,nullptr,
                options_.compute_measurement_green);
            state_=std::make_shared<const Snapshot>(makeSnapshot(
                state_->path,state_->factors,std::move(weight),state_->weight_order));
        }
        return *state_;
    }
    bool deferredPhase() const { return options_.defer_phase_measurement; }
    Snapshot fullReference(const PhononPath &path,bool wantPrefixCache=true,bool computeGreen=true) const {
        if(options_.defer_phase_measurement) {
            if(options_.merge_node_half_steps && !a2Disabled()) {
                auto prepared=profilePath("factor_generation_and_storage",[&] {
                    return prepareOwnedPath(path,trial_.Phi.rows(),generator_,/*wantCayley=*/false,
                                            &generationWorkspace_,/*wantHalfSteps=*/false);});
                auto factors=std::move(prepared.factors);
                auto weight=evaluatePathMagnitude(path,trial_,factors,block_,&prepared.weightOrder,
                                                  &weightWorkspace_,true);
                return profilePath("snapshot_path_copy",[&] {
                    return makeSnapshot(path,std::move(factors),std::move(weight),
                                        std::move(prepared.weightOrder));});
            }
            auto factors=profilePath("factor_generation_and_storage",[&] {
                return pathFactors(path,trial_.Phi.rows(),generator_);});
            auto weightOrder=reverseFactors(factors);
            auto weight=evaluatePathMagnitude(path,trial_,factors,block_,&weightOrder,
                                              &weightWorkspace_,true);
            return profilePath("snapshot_path_copy",[&] {
                return makeSnapshot(path,std::move(factors),std::move(weight),std::move(weightOrder));});
        }
        if(options_.merge_node_half_steps && !a2Disabled()) {
            auto prepared=profilePath("factor_generation_and_storage",[&] {
                return prepareOwnedPath(path,trial_.Phi.rows(),generator_,true,
                                        &generationWorkspace_);});
            auto factors=std::move(prepared.factors);
            auto prefixCache=wantPrefixCache&&options_.reuse_signed_prefix_checkpoints
                ?std::make_shared<PureProjectorSignedPrefixCache>():nullptr;
            auto weight=evaluatePath(path,trial_,factors,block_,options_.singular_values_only,
                                     &prepared.weightOrder,true,preparedTrial(),&prepared.cayley,
                                     &weightWorkspace_,nullptr,prefixCache.get(),-1,nullptr,computeGreen);
            return profilePath("snapshot_path_copy",[&] {
                return makeSnapshot(path,std::move(factors),std::move(weight),
                                    std::move(prepared.weightOrder),false,
                                    std::move(prefixCache));});
        }
        auto factors=profilePath("factor_generation_and_storage",[&] {
            return pathFactors(path,trial_.Phi.rows(),generator_);});
        auto weightOrder=reverseFactors(factors);
        auto weight=evaluatePath(path,trial_,factors,block_,options_.singular_values_only,
                                 &weightOrder,true,preparedTrial(),nullptr,&weightWorkspace_,nullptr,nullptr,-1,nullptr,computeGreen);
        return profilePath("snapshot_path_copy",[&] {
            return makeSnapshot(path,std::move(factors),std::move(weight),std::move(weightOrder));});
    }
    // Fresh canonical half-steps and full propagation: no accepted factors,
    // merged weight order, prepared trial/Cayley data or prefix checkpoints.
    // Deferred transitions retain magnitude-only semantics; measurements must
    // certify the signed weight and requested Green. No stability thresholds change.
    Snapshot uncachedReference(const PhononPath &path,bool measurement=false,
                               bool computeGreen=true) const {
        auto factors=pathFactors(path,trial_.Phi.rows(),generator_);
        auto weight=options_.defer_phase_measurement&&!measurement
            ?evaluatePathMagnitude(path,trial_,factors,block_)
            :evaluatePath(path,trial_,factors,block_,false,nullptr,false,
                          nullptr,nullptr,nullptr,nullptr,nullptr,-1,nullptr,computeGreen);
        // Reconstruct future proposal bookkeeping from these fresh factors,
        // only after certification. It is not used to evaluate this weight.
        std::vector<PureProjectorSlice> order;
        if(options_.merge_node_half_steps&&!a2Disabled()) {
            for(int j=path.intervals();j>=0;--j) {
                MatType b=factors[nodeFactors(j,path.intervals()).front()].matrix;
                if(j>0&&j<path.intervals()) b=(b*b).eval();
                if(!b.allFinite()||(b-b.adjoint()).norm()>1e-12*std::max(1.,b.norm()))
                    throw std::runtime_error("invalid rebuilt merged Majorana node factor");
                order.push_back(phononNodeSlice(
                    b,generator_,"owned_node_weight:"+std::to_string(j)));
            }
        } else {
            // The accepted snapshot must remain ready for the next local
            // proposal even though this recovery evaluation deliberately used
            // no cached order. Unmerged proposals update slots in this reversed
            // factor list in place.
            order=reverseFactors(factors);
        }
        const bool localOrder=options_.merge_node_half_steps&&!a2Disabled();
        return makeSnapshot(path,std::move(factors),std::move(weight),
                            std::move(order),localOrder);
    }
#ifdef KITAEV_SSH_TEST_FAILURE_INJECTION
    std::function<void(const char*,bool)> failureInjection;
#endif
    Proposal propose(int node,int bond,double value,bool reference=false) const {
        const auto &old=*state_;
        const double local=old.path.localLogRatio(node,bond,value);
        auto path=old.path.withCoordinate(node,bond,value);
        ++cacheDiagnostics_.local_proposals;
        testFailure("local",reference);
        if(reference) return Proposal(state_,std::make_shared<const Snapshot>(uncachedReference(path)),local);
        auto weightOrder=old.weight_order;
        const int intervals=path.intervals();
        if(options_.defer_phase_measurement && options_.merge_node_half_steps && !a2Disabled()) {
            // Deferred updates retain the existing spectral node construction.
            // Only the changed node is regenerated; all other immutable records
            // are shared through PureFactorMatrix's storage.
            auto prepared=prepareOwnedNode(path,node,trial_.Phi.rows(),generator_,
                                           generationWorkspace_,/*wantHalfSteps=*/false);
            weightOrder[intervals-node]=phononNodeSlice(
                std::move(prepared.merged_factor),generator_,"owned_node_weight:"+std::to_string(node));
            ++cacheDiagnostics_.local_node_rebuilds;
            cacheDiagnostics_.local_node_reuses+=std::uint64_t(path.nodes()-1);
            auto weight=evaluatePathMagnitude(path,trial_,{},block_,&weightOrder,
                                              &weightWorkspace_,true);
            auto next=std::make_shared<const Snapshot>(makeSnapshot(
                std::move(path),{},std::move(weight),std::move(weightOrder),false));
            return Proposal(state_,std::move(next),local);
        }
        auto factors=old.factors;
        MatType b;
        nodeHalfStepInto(path,node,trial_.Phi.rows(),generator_,b,generationWorkspace_);
        for(int f:nodeFactors(node,intervals))
            factors[f]=phononNodeSlice(b,generator_,"phonon_node:"+std::to_string(node));
        if(options_.merge_node_half_steps && !a2Disabled()) {
            const bool canResume=old.weight_order_is_local;
            if(old.weight_order_is_local) {
                // Eager local proposals keep their original h*h arithmetic.
                // Reuse every accepted merged node except the changed one.
                weightOrder=old.weight_order;
                weightOrder[intervals-node]=localMergedNode(factors,node,intervals);
                cacheDiagnostics_.local_node_reuses+=std::uint64_t(path.nodes()-1);
            } else {
                // The first eager local proposal establishes the h*h cache.
                // This preserves the pre-cache arithmetic exactly, then makes
                // later proposals incremental.
                weightOrder=buildLocalMergedOrder(factors,intervals);
            }
            std::shared_ptr<PureProjectorSignedPrefixCache> prefixCache;
            if(options_.reuse_signed_prefix_checkpoints)
                prefixCache=std::make_shared<PureProjectorSignedPrefixCache>();
            const int changedWeightIndex=ownedWeightIndexForNode(weightOrder,node);
            std::shared_ptr<PureProjectorLUCayleyCache> cayleyCache;
            if(options_.reuse_lu_cayley) {
                cayleyCache=canResume&&old.lu_cayley_cache
                    ?std::make_shared<PureProjectorLUCayleyCache>(*old.lu_cayley_cache)
                    :std::make_shared<PureProjectorLUCayleyCache>(weightOrder.size());
                // The changed factor must be solved again even if it happens to
                // have the same node label. Storage identity is the authority.
                (*cayleyCache)[changedWeightIndex].reset();
            }
            auto weight=evaluatePath(path,trial_,factors,block_,options_.singular_values_only,
                                     &weightOrder,true,preparedTrial(),nullptr,&weightWorkspace_,
                                     canResume?old.signed_prefix_cache.get():nullptr,
                                     prefixCache.get(),canResume?changedWeightIndex:-1,
                                     cayleyCache.get());
            auto next=std::make_shared<const Snapshot>(makeSnapshot(
                std::move(path),std::move(factors),std::move(weight),
                std::move(weightOrder),true,std::move(prefixCache),
                std::move(cayleyCache)));
            return Proposal(state_,std::move(next),local);
        }
        // Unmerged and A2-disabled paths retain their ordinary reversed order,
        // but update only the two factor slots owned by this node.
        for(int f:nodeFactors(node,intervals))
            weightOrder[2*intervals-1-f]=factors[f];
        auto weight=options_.defer_phase_measurement
            ? evaluatePathMagnitude(path,trial_,factors,block_,&weightOrder,&weightWorkspace_,true)
            : evaluatePath(path,trial_,factors,block_,options_.singular_values_only,
                           &weightOrder,true,preparedTrial(),nullptr,&weightWorkspace_);
        auto next=std::make_shared<const Snapshot>(makeSnapshot(
            std::move(path),std::move(factors),std::move(weight),
            std::move(weightOrder),false));
        return Proposal(state_,std::move(next),local);
    }
    // Whole single-bond path refresh, still ONE atomic snapshot transaction.
    // Every time node changes; the full reference is the minimal safe rebuild.
    Proposal proposeBondPath(int bond,const std::vector<double> &values,bool reference=false) const {
        const auto path=profilePath("proposal_path_copy",[&] {return state_->path.withBondCoordinates(bond,values);});
        ++cacheDiagnostics_.bond_refresh_full_rebuilds;
        testFailure("bond-refresh",reference);
        auto next=std::make_shared<const Snapshot>(reference?uncachedReference(path):fullReference(path,false,false));
        const double delta=next->weight.log_phonon-state_->weight.log_phonon;
        return Proposal(state_,std::move(next),delta);
    }
    Proposal proposeQModePath(const PhononPath &path,bool reference=false) const {
        testFailure("q-mode",reference);
        auto next=std::make_shared<const Snapshot>(reference?uncachedReference(path):fullReference(path,false,false));
        const double delta=next->weight.log_phonon-state_->weight.log_phonon;
        return Proposal(state_,std::move(next),delta);
    }
    void accept(Proposal &proposal) {
        checkProposal(proposal);
        state_=std::move(proposal.next_); // single nonthrowing snapshot commit
        proposal.base_.reset();
    }
    void reject(Proposal &proposal) {
        checkProposal(proposal); proposal.next_.reset(); proposal.base_.reset();
    }
    // A1: lazily prepared, immutable thin QR of the fixed trial state. The
    // trial is a const member set at construction and never changes, so the
    // preparation is valid for the engine's whole lifetime and is independent
    // of candidate/accepted configurations (no stale-state risk).
    static bool a1Disabled() {
        static const bool v=[] { const char *e=std::getenv("PFQMC_NO_A1"); return e&&*e&&*e!='0'; }();
        return v;
    }
    static bool a2Disabled() {
        static const bool v=[] { const char *e=std::getenv("PFQMC_NO_A2"); return e&&*e&&*e!='0'; }();
        return v;
    }
    const PreparedGreenTrial *preparedTrial() const {
        if(a1Disabled()) return nullptr;
        if(!preparedTrial_) {            auto prepared=std::make_shared<PreparedGreenTrial>(
                prepareGreenTrial(trial_,pureFastGreenOptions(PureFastOptions{})));
            preparedTrial_=std::move(prepared);
        }
        return preparedTrial_->ok()?preparedTrial_.get():nullptr;
    }
private:
    void testFailure(const char *operation,bool reference) const {
#ifdef KITAEV_SSH_TEST_FAILURE_INJECTION
        if(failureInjection) failureInjection(operation,reference);
        if(testNumericalFailure) testNumericalFailure(operation,reference);
#else
        (void)operation; (void)reference;
#endif
    }
    static Snapshot makeSnapshot(PhononPath path,std::vector<PureProjectorSlice> factors,
                                 PathWeight weight,
                                 std::vector<PureProjectorSlice> weightOrder={},
                                 bool weightOrderIsLocal=false,
                                 std::shared_ptr<const PureProjectorSignedPrefixCache> prefixCache=nullptr,
                                 std::shared_ptr<const PureProjectorLUCayleyCache> cayleyCache=nullptr) {
        const auto hash=phononPathConfigurationHash(path);
        return Snapshot{std::move(path),std::move(factors),std::move(weightOrder),
                        weightOrderIsLocal,std::move(prefixCache),std::move(cayleyCache),
                        std::move(weight),hash};
    }
    PureProjectorSlice localMergedNode(const std::vector<PureProjectorSlice> &factors,
                                       int node,int intervals) const {
        MatType b=factors[nodeFactors(node,intervals).front()].matrix;
        if(node>0 && node<intervals) {
            profilePath("merge_square_gemm",[&]{ b=(b*b).eval(); return true; });
            ++cacheDiagnostics_.local_merged_squares;
            profilePath("merge_hermiticity",[&]{
                if(!b.allFinite() || (b-b.adjoint()).norm()>1e-12*std::max(1.,b.norm()))
                    throw std::runtime_error("invalid merged Majorana node factor");
                return true;
            });
        }
        ++cacheDiagnostics_.local_node_rebuilds;
        return phononNodeSlice(b,generator_,"owned_node_weight:"+std::to_string(node));
    }
    std::vector<PureProjectorSlice> buildLocalMergedOrder(
        const std::vector<PureProjectorSlice> &factors,int intervals) const {
        std::vector<PureProjectorSlice> result;
        result.reserve(intervals+1);
        profilePath("merge_owned_node_factors",[&]{
            for(int node=intervals;node>=0;--node)
                result.push_back(localMergedNode(factors,node,intervals));
            return true;
        });
        return result;
    }
    static int ownedWeightIndexForNode(
        const std::vector<PureProjectorSlice> &weightOrder,int node) {
        const std::string label="owned_node_weight:"+std::to_string(node);
        for(int index=0;index<int(weightOrder.size());++index)
            if(weightOrder[index].label==label)return index;
        throw std::logic_error("missing owned node in weight order: "+label);
    }
    GaussianTrialState trial_;
    MajoranaGenerator generator_;
    mutable MajoranaGenerationWorkspace generationWorkspace_;
    mutable PureStableWeightWorkspace weightWorkspace_;
    mutable std::shared_ptr<const PreparedGreenTrial> preparedTrial_;    int block_;
    PhononWeightOptions options_;
    std::shared_ptr<const Snapshot> state_;
    mutable PhononCacheDiagnostics cacheDiagnostics_;
    PathWeight evaluateOwned(const PhononPath &path,const std::vector<PureProjectorSlice> &factors) const {
        if(!options_.merge_node_half_steps || a2Disabled())
            return evaluatePath(path,trial_,factors,block_,options_.singular_values_only,nullptr,true,
                                preparedTrial(),nullptr,&weightWorkspace_);        // Only this class's own canonical exp(-h K(X_j)/2), eta=1 factors
        // reach here. No generic-factor or approximate-matrix matching.
        // Scalar phases in arbitrary factors still use evaluatePath reference.
        auto nodes=profilePath("merge_owned_node_factors",[&] {
            std::vector<PureProjectorSlice> result;
            result.reserve(path.nodes());
            for(int j=path.intervals();j>=0;--j) {
                MatType b=factors[nodeFactors(j,path.intervals()).front()].matrix;
                if(j>0 && j<path.intervals()) {
                    profilePath("merge_square_gemm",[&]{ b=(b*b).eval(); return true; });
                    profilePath("merge_hermiticity",[&]{
                        if(!b.allFinite() || (b-b.adjoint()).norm()>1e-12*std::max(1.,b.norm()))
                            throw std::runtime_error("invalid merged Majorana node factor");
                        return true;
                    });
                }
                result.emplace_back(b,1.,"owned_node_weight:"+std::to_string(j));
            }            return result;
        });
        return evaluatePath(path,trial_,factors,block_,options_.singular_values_only,&nodes,true,
                            preparedTrial(),nullptr,&weightWorkspace_);
    }    void checkProposal(const Proposal &proposal) const {
        if (!proposal.next_ || proposal.base_!=state_) throw std::logic_error("foreign, stale or consumed phonon proposal");
    }
};

} // namespace kitaev_ssh
#endif
