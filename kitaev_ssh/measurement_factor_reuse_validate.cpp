// Deterministic measurement-factor provenance checks. No sampler, MC or RNG.
#include "phonon_projector.h"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>

using namespace kitaev_ssh;
namespace {
int checks=0,measurements=0,reusedMeasurements=0,rebuiltMeasurements=0;
double maxLog=0.,maxPhase=0.,maxGreen=0.;

void require(bool ok,const std::string &what) {
    ++checks;if(!ok)throw std::runtime_error(what);
}
bool same(const MatType &a,const MatType &b) {
    return a.rows()==b.rows()&&a.cols()==b.cols()&&
        (a.size()==0||std::memcmp(a.data(),b.data(),sizeof(DataType)*a.size())==0);
}
bool sameDouble(double a,double b) {return std::memcmp(&a,&b,sizeof(a))==0;}
bool sameComplex(DataType a,DataType b) {return std::memcmp(&a,&b,sizeof(a))==0;}
std::uint64_t hash(const PureFactorMatrix &a) {
    std::uint64_t result=1469598103934665603ULL;a.mixHash(result);return result;
}
struct FrozenSlice {
    MatType matrix;
    std::uint64_t coefficientHash;
    DataType eta;
    std::string label;
    bool validated;
    double identity;
};
std::vector<FrozenSlice> freeze(const std::vector<PureProjectorSlice> &slices) {
    std::vector<FrozenSlice> out;
    for(const auto &s:slices)
        out.push_back({MatType(s.matrix),hash(s.matrix),s.eta,s.label,
                       s.validated_dense,s.identity_residual});
    return out;
}
void frozen(const std::vector<PureProjectorSlice> &slices,
            const std::vector<FrozenSlice> &saved,const std::string &where) {
    require(slices.size()==saved.size(),where+": slice count changed");
    for(std::size_t i=0;i<slices.size();++i) {
        const auto &s=slices[i];const auto &old=saved[i];
        require(same(MatType(s.matrix),old.matrix),where+": complex coefficients/signed zeros changed");
        require(hash(s.matrix)==old.coefficientHash,where+": coefficient hash changed");
        require(sameComplex(s.eta,old.eta)&&s.label==old.label&&
                s.validated_dense==old.validated&&sameDouble(s.identity_residual,old.identity),
                where+": slice metadata changed");
    }
}
void equalSlices(const std::vector<PureProjectorSlice> &a,
                 const std::vector<PureProjectorSlice> &b,const std::string &where) {
    frozen(a,freeze(b),where);
}
void compareWeight(const PathWeight &a,const PathWeight &b,bool green,
                   const std::string &where) {
    require(a.magnitude_valid&&a.phase_valid&&b.magnitude_valid&&b.phase_valid,
            where+": signed weight invalid");
    require(sameDouble(a.log_phonon,b.log_phonon),where+": phonon action changed");
    const double log=std::abs(a.log_abs_matter-b.log_abs_matter);
    const double phase=std::abs(a.matter_phase-b.matter_phase);
    maxLog=std::max(maxLog,log);maxPhase=std::max(maxPhase,phase);
    require(std::isfinite(log)&&log<2e-10,where+": signed log weight differs");
    require(std::isfinite(phase)&&phase<2e-10,where+": signed phase differs");
    require(a.green_valid==green&&b.green_valid==green,where+": Green validity differs");
    if(green) {
        require(a.center_green.rows()==b.center_green.rows()&&
                a.center_green.cols()==b.center_green.cols()&&a.center_green.size()>0,
                where+": full Green dimensions differ");
        const double err=(a.center_green-b.center_green).cwiseAbs().maxCoeff();
        maxGreen=std::max(maxGreen,err);
        require(std::isfinite(err)&&err<2e-10,where+": full complex Green differs");
    } else require(a.center_green.size()==0&&b.center_green.size()==0,
                   where+": omitted Green was populated");
}
PhononPath pathFor(int bonds,int intervals=4,double dt=.08) {
    PhononPath path(bonds,intervals,dt,{7.,.6});
    std::vector<std::vector<double>> x(path.nodes(),std::vector<double>(path.bonds()));
    for(int j=0;j<path.nodes();++j)for(int b=0;b<path.bonds();++b)
        x[j][b]=.07*std::sin(.21*j+.13*b)+.01*std::cos(.37*j-.09*b);
    return path.withCoordinates(x);
}
std::vector<double> bondCoordinates(const PhononPath &path,int bond) {
    std::vector<double> values(path.nodes());
    for(int j=0;j<path.nodes();++j)values[j]=path.coordinate(j,bond);
    return values;
}
PhononWeightOptions options(bool green=true,bool fast=true,bool merged=true,bool deferred=true) {
    PhononWeightOptions out;out.singular_values_only=true;
    out.merge_node_half_steps=merged;out.defer_phase_measurement=deferred;
    out.reuse_signed_prefix_checkpoints=true;out.reuse_lu_cayley=true;
    out.compute_measurement_green=green;out.center_green_fast=fast;return out;
}
GaussianTrialState pairTrial(int pairs) {
    MatType phi=MatType::Zero(2*pairs,pairs);
    const double scale=1./std::sqrt(2.);
    for(int i=0;i<pairs;++i){phi(i,i)=scale;phi(pairs+i,i)=DataType(0.,scale);}
    return GaussianTrialState::fromPhi(phi);
}

// Rebuild with the original no-source prepare path. This oracle has no
// accepted factor owners or private snapshot provenance to consult.
void measureAgainstOriginal(CoupledPhononPath &state,const GaussianTrialState &trial,
                            const MajoranaGenerator &generator,bool expectReuse,
                            bool green=true,bool fast=true,const std::string &where="measurement") {
    const auto before=state.state();
    const auto savedOrder=freeze(before.weight_order),savedFactors=freeze(before.factors);
    const auto configuration=state.configurationHash();
    const auto counters=state.cacheDiagnostics();
    const auto original=prepareOwnedPath(before.path,trial.Phi.rows(),generator,true);
    const auto reference=evaluatePath(before.path,trial,original.factors,2,true,
        &original.weightOrder,true,state.preparedTrial(),&original.cayley,
        nullptr,nullptr,nullptr,-1,nullptr,green,fast);
    const auto &measured=state.ensureMeasured();++measurements;
    const auto after=state.cacheDiagnostics();
    require(after.measurement_spectral_reuses==counters.measurement_spectral_reuses+(expectReuse?1:0)&&
            after.measurement_full_rebuilds==counters.measurement_full_rebuilds+(expectReuse?0:1)&&
            after.measurement_merged_node_reuses==counters.measurement_merged_node_reuses+
                (expectReuse?std::uint64_t(before.path.nodes()):0),
            where+": reuse/rebuild counters disagree with selected source");
    require(measured.configuration_hash==configuration&&
            phononPathConfigurationHash(measured.path)==configuration,
            where+": measurement changed configuration identity");
    equalSlices(measured.factors,original.factors,where+": canonical halves");
    equalSlices(measured.weight_order,original.weightOrder,where+": merged spectral B");
    compareWeight(measured.weight,reference,green,where+": original prepare");
    const auto independent=evaluatePath(before.path,trial,pathFactors(before.path,trial.Phi.rows(),generator),2,
        false,nullptr,false,nullptr,nullptr,nullptr,nullptr,nullptr,-1,nullptr,green,false);
    compareWeight(measured.weight,independent,green,where+": independent canonical reference");
    require(before.weight_order.size()==measured.weight_order.size(),where+": merged order size changed");
    for(std::size_t i=0;i<measured.weight_order.size();++i)
        require(measured.weight_order[i].matrix.sharesStorageWith(before.weight_order[i].matrix)==expectReuse,
                where+(expectReuse?": certified factor was rebuilt":": uncertified factor was reused"));
    if(expectReuse)++reusedMeasurements;else ++rebuiltMeasurements;
    frozen(before.weight_order,savedOrder,where+": held snapshot order");
    frozen(before.factors,savedFactors,where+": held snapshot halves");
    const auto *measuredAddress=&measured;
    require(&state.ensureMeasured()==measuredAddress,where+": repeated measurement replaced valid snapshot");
    require(state.cacheDiagnostics().measurement_spectral_reuses==after.measurement_spectral_reuses&&
            state.cacheDiagnostics().measurement_full_rebuilds==after.measurement_full_rebuilds&&
            state.cacheDiagnostics().measurement_merged_node_reuses==after.measurement_merged_node_reuses,
            where+": repeated measurement changed reuse counters");
}

void honeycombTransitions(int L,bool green,bool fast) {
    KitaevSSHModel model(HoneycombLattice(L,L,0),p2Couplings(.8));
    const auto trial=freeGroundState(model);const auto generator=honeycombGenerator(model);
    const auto path=pathFor(model.lattice().bondCount());
    CoupledPhononPath state(path,trial,generator,2,options(green,fast));
    const auto oldest=state.state();const auto oldestOrder=freeze(oldest.weight_order);
    const auto initial=state.configurationHash();

    // Public copies are editable, but neither Proposal's commit pointers nor
    // the accepted immutable Snapshot can be created/adopted by a caller.
    static_assert(!std::is_constructible<CoupledPhononPath::Proposal,
                  std::shared_ptr<const CoupledPhononPath::Snapshot>,
                  std::shared_ptr<const CoupledPhononPath::Snapshot>,double>::value,
                  "proposal construction must remain private");
    static_assert(!std::is_constructible<OwnedMeasurementReuse,
                  const std::vector<PureProjectorSlice>&>::value,
                  "measurement source token construction must remain private");
    static_assert(!std::is_default_constructible<OwnedMeasurementReuse>::value,
                  "empty caller-created measurement source token must be impossible");
    auto edited=oldest;
    edited.weight_order[0].matrix=MatType::Identity(trial.Phi.rows(),trial.Phi.rows());
    edited.weight_order[1].eta=DataType(-1.,-0.);
    edited.weight_order[2].label="owned_node_weight:forged";
    edited.weight_order[3].identity_residual=0.;
    edited.weight_order[4].validated_dense=false;
    edited.weight_order_is_local=true;edited.configuration_hash=0;
    frozen(state.state().weight_order,oldestOrder,"external Snapshot edits");
    require(state.configurationHash()==initial,"external Snapshot edits changed accepted hash");
    measureAgainstOriginal(state,trial,generator,true,green,fast,"initial honeycomb");

    for(int node:{0,path.intervals()/2,path.intervals()}) {
        const auto accepted=state.state();const auto acceptedOrder=freeze(accepted.weight_order);
        const double value=accepted.path.coordinate(node,0)+.013*(node+1);
        auto rejected=state.propose(node,0,value);
        const auto heldCandidate=rejected.candidate();const auto heldOrder=freeze(heldCandidate.weight_order);
        require(!heldCandidate.weight.phase_valid&&!heldCandidate.weight.green_valid,
                "local deferred proposal populated signed caches");
        state.reject(rejected);
        require(state.configurationHash()==accepted.configuration_hash,"reject changed accepted identity");
        frozen(state.state().weight_order,acceptedOrder,"rejected local proposal");
        frozen(heldCandidate.weight_order,heldOrder,"held rejected candidate");
        auto next=state.propose(node,0,value);
        for(int j=0;j<path.nodes();++j)
            require(next.candidate().weight_order[path.intervals()-j].matrix.sharesStorageWith(
                        accepted.weight_order[path.intervals()-j].matrix)==(j!=node),
                    "local replacement changed unrelated owner");
        state.accept(next);
        measureAgainstOriginal(state,trial,generator,true,green,fast,"accepted local honeycomb");
        frozen(accepted.weight_order,acceptedOrder,"held previous accepted snapshot");
    }
    auto bondValues=bondCoordinates(state.state().path,1);
    for(std::size_t j=0;j<bondValues.size();++j)bondValues[j]+=.012*std::cos(.31*j);
    auto bond=state.proposeBondPath(1,bondValues);
    const auto heldBond=bond.candidate();const auto heldBondOrder=freeze(heldBond.weight_order);
    // A second independent candidate must not overwrite the first one's pool owners.
    auto second=state.proposeBondPath(2,bondValues);
    state.reject(second);frozen(heldBond.weight_order,heldBondOrder,"second proposal/reject");
    state.accept(bond);
    measureAgainstOriginal(state,trial,generator,true,green,fast,"accepted whole-bond honeycomb");
    frozen(heldBond.weight_order,heldBondOrder,"held whole-bond candidate after measurement");

    auto collectivePath=state.state().path;
    std::vector<std::vector<double>> collective(collectivePath.nodes());
    for(int j=0;j<collectivePath.nodes();++j) {
        collective[j]=collectivePath.node(j);
        for(int b=0;b<collectivePath.bonds();++b)
            collective[j][b]+=.004*std::sin(.17*j+.29*b);
    }
    auto collectiveProposal=state.proposeQModePath(collectivePath.withCoordinates(collective));
    state.accept(collectiveProposal);
    measureAgainstOriginal(state,trial,generator,true,green,fast,"accepted collective honeycomb");
    frozen(oldest.weight_order,oldestOrder,"oldest snapshot after transitions");
}

void recoveryAndMixed() {
    KitaevSSHModel model(HoneycombLattice(2,2,0),p2Couplings(.8));
    const auto trial=freeGroundState(model);const auto generator=honeycombGenerator(model);
    const auto path=pathFor(model.lattice().bondCount());
    CoupledPhononPath recovered(path,trial,generator,2,options(),true);
    const auto old=recovered.state();const auto oldOrder=freeze(old.weight_order);
    require(old.weight_order_is_local,"uncached recovery fixture has no half-product order");
    auto mixed=recovered.propose(2,0,path.coordinate(2,0)+.019);
    const auto mixedSnapshot=mixed.candidate();const auto mixedOrder=freeze(mixedSnapshot.weight_order);
    recovered.accept(mixed);
    measureAgainstOriginal(recovered,trial,generator,false,true,true,"uncached plus spectral local mixed path");
    frozen(mixedSnapshot.weight_order,mixedOrder,"held mixed snapshot");
    frozen(old.weight_order,oldOrder,"held uncached snapshot");
    auto nowSpectral=recovered.propose(1,0,recovered.state().path.coordinate(1,0)-.011);
    recovered.accept(nowSpectral);
    measureAgainstOriginal(recovered,trial,generator,true,true,true,"post-rebuild spectral local path");

    auto referenceProposal=recovered.propose(2,1,recovered.state().path.coordinate(2,1)+.014,true);
    recovered.accept(referenceProposal);
    measureAgainstOriginal(recovered,trial,generator,false,true,true,"accepted uncached reference proposal");
    auto referenceBond=recovered.proposeBondPath(0,bondCoordinates(recovered.state().path,0),true);
    recovered.accept(referenceBond);
    const auto beforeReference=recovered.state();const auto beforeReferenceOrder=freeze(beforeReference.weight_order);
    const auto referenceWeight=recovered.uncachedReference(beforeReference.path,true,true).weight;
    const auto &referenceMeasured=recovered.ensureMeasured(true);
    compareWeight(referenceMeasured.weight,referenceWeight,true,"explicit reference measurement");
    for(std::size_t i=0;i<beforeReference.weight_order.size();++i)
        require(!referenceMeasured.weight_order[i].matrix.sharesStorageWith(beforeReference.weight_order[i].matrix),
                "explicit reference measurement reused source owner");
    frozen(beforeReference.weight_order,beforeReferenceOrder,"explicit reference retained snapshot");
    auto afterReference=recovered.propose(1,0,recovered.state().path.coordinate(1,0)+.009);
    recovered.accept(afterReference);
    measureAgainstOriginal(recovered,trial,generator,false,true,true,"reference measurement then mixed local");
}

void genericAndIdentity() {
    // A stateful callable can change K while every public generator flag and
    // X stays fixed. It cannot receive the immutable honeycomb source proof.
    auto scale=std::make_shared<double>(1.);
    MajoranaGenerator::IntoFunction into=[scale](const std::vector<double> &x,MatType &k) {
        k.setZero(8,8);
        for(int i=0;i<4;++i)for(int j=0;j<4;++j) {
            const double t=(*scale)*((i==j?.6:.04*std::cos(.3*i+.2*j))+.03*x[0]);
            k(i,4+j)=DataType(0.,-t);k(4+j,i)=DataType(0.,t);
        }
        // Exact intra-sublattice coupling requires the full Hermitian path.
        k(0,1)=DataType(0.,-.031);k(1,0)=DataType(0.,.031);
    };
    MajoranaGenerator generic([into](const std::vector<double> &x){MatType k;into(x,k);return k;},into,true,true);
    require(!generic.hasImmutableSpectralSource(),"public generator flags forged immutable source proof");
    const auto trial=pairTrial(4);const auto path=pathFor(1);
    CoupledPhononPath state(path,trial,generic,2,options());
    measureAgainstOriginal(state,trial,generic,false,true,true,"generic nonreduced source");
    auto p=state.propose(2,0,path.coordinate(2,0)+.017);state.accept(p);
    const auto stale=freeze(state.state().weight_order);*scale=1.09;
    measureAgainstOriginal(state,trial,generic,false,true,true,"stateful generic changed source");
    bool changed=false;for(std::size_t i=0;i<stale.size();++i)
        changed=changed||!same(MatType(state.state().weight_order[i].matrix),stale[i].matrix);
    require(changed,"stateful generic fixture did not change merged coefficients");

    MajoranaGenerator zero([](const std::vector<double> &){return MatType::Zero(8,8).eval();});
    CoupledPhononPath identity(path,trial,zero,2,options());
    measureAgainstOriginal(identity,trial,zero,false,true,true,"identity general source");
    for(const auto &s:identity.state().weight_order) {
        require(sameDouble(s.identity_residual,-1.),"generic identity residual metadata changed");
        require(same(MatType(s.matrix),MatType::Identity(8,8)),"identity coefficients changed");
    }
}

void immutableSourceMutations() {
    KitaevSSHModel model(HoneycombLattice(2,2,0),p2Couplings(.8));
    const auto trial=freeGroundState(model);const auto canonical=honeycombGenerator(model);
    const auto path=pathFor(model.lattice().bondCount());
    require(canonical.hasImmutableSpectralSource(),"factory source identity is absent");
    const auto copy=canonical;require(copy.hasImmutableSpectralSource(),"copy lost immutable source identity");
    for(int edit=0;edit<6;++edit) {
        auto generator=canonical;
        if(edit==0)generator.value=[canonical](const std::vector<double> &x){return canonical(x);};
        if(edit==1)generator.into=[canonical](const std::vector<double> &x,MatType &k){canonical.generateInto(x,k);};
        if(edit==2)generator.validated_dense=false;
        if(edit==3)generator.bipartite_real=false;
        if(edit==4)generator.value=honeycombGenerator(model).value;
        if(edit==5)generator.into=honeycombGenerator(model).into;
        require(!generator.hasImmutableSpectralSource(),"edited public generator retained immutable source proof");
        CoupledPhononPath state(path,trial,generator,2,options());
        measureAgainstOriginal(state,trial,generator,false,true,true,"edited canonical generator source");
    }
}

void eagerAndUnmergedControls() {
    KitaevSSHModel model(HoneycombLattice(2,2,0),p2Couplings(.8));
    const auto trial=freeGroundState(model);const auto generator=honeycombGenerator(model);
    const auto path=pathFor(model.lattice().bondCount());
    for(bool deferred:{false,true})for(bool merged:{false,true}) {
        if(deferred&&merged)continue;
        CoupledPhononPath state(path,trial,generator,2,options(true,true,merged,deferred));
        auto p=state.proposeBondPath(0,bondCoordinates(path,0));state.accept(p);
        const auto before=state.state();const auto savedOrder=freeze(before.weight_order);
        const auto expected=deferred?evaluatePath(before.path,trial,before.factors,2,true,
            &before.weight_order,true,state.preparedTrial(),nullptr,nullptr,nullptr,nullptr,-1,nullptr,true,true)
            :state.fullReference(before.path,false,true).weight;
        const auto &measured=state.ensureMeasured();compareWeight(measured.weight,expected,true,"eager/unmerged control");
        frozen(measured.weight_order,savedOrder,"eager/unmerged retained weight order");
        for(std::size_t i=0;i<before.weight_order.size();++i)
            require(measured.weight_order[i].matrix.sharesStorageWith(before.weight_order[i].matrix),
                    "eager/unmerged existing owner behavior changed");
    }
}

template<class F> void overflowFailure(F &&operation) {
    bool runtime=false;std::string message;
    try{operation();}catch(const std::invalid_argument &){throw std::runtime_error("overflow bypassed runtime recovery gate");}
    catch(const std::runtime_error &e){runtime=true;message=e.what();}
    require(runtime&&message=="invalid merged Majorana node factor","overflow recovery failure semantics changed");
}
void overflowRecovery() {
    MajoranaGenerator::IntoFunction into=[](const std::vector<double> &x,MatType &k) {
        k.setZero(4,4);for(int i=0;i<2;++i){k(i,2+i)=DataType(0.,-x[0]);k(2+i,i)=DataType(0.,x[0]);}
    };
    MajoranaGenerator generator([into](const std::vector<double> &x){MatType k;into(x,k);return k;},into,true,true);
    const auto path=PhononPath(1,2,1.,{1.,1.}).withCoordinates({{.02},{.1},{.02}});
    const auto extreme=path.withCoordinates({{.02},{800.},{.02}});
    BipartiteSpectrum spectral;require(spectral.compute(generator(extreme.node(1))),"overflow fixture spectral setup failed");
    MatType half,full;spectral.exponential(.5,half);spectral.exponential(1.,full);
    require(half.allFinite()&&!full.allFinite()&&!(half*half).eval().allFinite(),"overflow fixture misses fallback");
    CoupledPhononPath state(path,pairTrial(2),generator,2,options());
    const auto saved=state.state();const auto savedOrder=freeze(saved.weight_order);
    const auto savedHash=state.configurationHash();
    overflowFailure([&]{(void)prepareOwnedPath(extreme,4,generator,false,nullptr,false);});
    overflowFailure([&]{(void)state.proposeBondPath(0,{.02,800.,.02});});
    overflowFailure([&]{(void)state.propose(1,0,800.);});
    require(state.configurationHash()==savedHash,"overflow changed accepted configuration");
    frozen(state.state().weight_order,savedOrder,"overflow accepted snapshot");
    frozen(saved.weight_order,savedOrder,"overflow held snapshot");
    auto normal=state.proposeBondPath(0,{.02,.11,.02});state.accept(normal);
    measureAgainstOriginal(state,pairTrial(2),generator,false,true,true,"normal recovery after overflow");
    frozen(saved.weight_order,savedOrder,"successful recovery held snapshot");
}

void finiteOverflowFallback() {
    // A delocalized positive mode has exp(711)=inf, but division among the
    // matrix entries leaves h*h finite. This is a successful fallback, so
    // its explicit product provenance must survive publication. No weight
    // is evaluated for these deliberately extreme, ill-conditioned factors.
    for(int dimension:{4,32}) {
        const int n=dimension/2;const bool reduced=dimension==32;
        MajoranaGenerator::IntoFunction into=[n](const std::vector<double> &x,MatType &k) {
            k.setZero(2*n,2*n);
            for(int i=0;i<n;++i)for(int j=0;j<n;++j) {
                k(i,n+j)=DataType(0.,-x[0]);k(n+j,i)=DataType(0.,x[0]);
            }
        };
        MajoranaGenerator generator([into](const std::vector<double> &x){MatType k;into(x,k);return k;},
                                    into,true,reduced);
        const auto path=PhononPath(1,2,1.,{1.,1.}).withCoordinates({{.02},{711./n},{.02}});
        BipartiteSpectrum bipartite;PhononNodeSpectrum spectral;
        require(spectral.compute(generator(path.node(1)),reduced,bipartite),"finite fallback fixture eigendecomposition failed");
        require(spectral.reduced()==reduced,"finite fallback fixture selected wrong spectrum");
        MatType half,full;spectral.exponential(.5,half);spectral.exponential(1.,full);
        const MatType expected=(half*half).eval();
        require(half.allFinite()&&!full.allFinite()&&expected.allFinite(),
                "finite fallback fixture misses finite half-product/nonfinite full distinction");
        MajoranaGenerationWorkspace generation;
        const auto node=prepareOwnedNode(path,1,dimension,generator,generation,true);
        require(!node.spectral,"finite full-assembly fallback certified a spectral node");
        require(same(node.half_step,half)&&same(node.merged_factor,expected),
                "finite node fallback changed original half/product complex coefficients");
        const auto original=prepareOwnedPath(path,dimension,generator,false);
        require(!original.all_spectral,"finite unpooled fallback certified an all-spectral path");
        require(same(MatType(original.weightOrder[1].matrix),expected),
                "finite unpooled fallback changed half-product coefficients");
        PureFactorMatrix::DenseStoragePool pool;
        const auto pooled=prepareOwnedPath(path,dimension,generator,false,nullptr,true,&pool);
        require(!pooled.all_spectral,"finite pooled fallback certified an all-spectral path");
        equalSlices(pooled.factors,original.factors,"finite fallback pooled canonical halves");
        equalSlices(pooled.weightOrder,original.weightOrder,"finite fallback pooled full factors");
        require(hash(pooled.weightOrder[1].matrix)==hash(PureFactorMatrix::fromValidatedDense(expected)),
                "finite fallback coefficient hash differs from original product");
    }
}
}
int main() {
    try {
        honeycombTransitions(2,true,true);
        honeycombTransitions(4,false,false);
        recoveryAndMixed();genericAndIdentity();immutableSourceMutations();
        eagerAndUnmergedControls();overflowRecovery();finiteOverflowFallback();
        std::cout<<std::setprecision(17)<<"{\"checks\":"<<checks<<",\"measurements\":"<<measurements
                 <<",\"reused_measurements\":"<<reusedMeasurements<<",\"rebuilt_measurements\":"<<rebuiltMeasurements
                 <<",\"max_signed_log_difference\":"<<maxLog<<",\"max_phase_difference\":"<<maxPhase
                 <<",\"max_full_green_difference\":"<<maxGreen<<"}\nMEASUREMENT_FACTOR_REUSE_VALIDATE_OK\n";
    } catch(const std::exception &e) {
        std::cerr<<"MEASUREMENT_FACTOR_REUSE_VALIDATE_FAIL: "<<e.what()<<'\n';return 1;
    }
}
