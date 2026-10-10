// Same physical proposals against freshly generated canonical half-steps.
// This tests the cache/update arithmetic, not an independent physical model;
// pfqmc_fock_validate supplies the independent many-body reference.
#include "phonon_mc.h"
#include <iostream>
#include <map>
#include <sstream>
using namespace kitaev_ssh;
namespace {
struct Metric { double error=0., tolerance=0.; };
std::map<std::string,Metric> metrics;
void check(const std::string &name,double error,double tolerance) {
    auto &m=metrics[name];m.error=std::max(m.error,error);m.tolerance=tolerance;
    if(!std::isfinite(error)||error>tolerance)
        throw std::runtime_error(name+" error="+std::to_string(error));
}
void demand(bool ok,const std::string &name) {check(name,ok?0.:1.,0.);}
void weights(const PathWeight &actual,const PathWeight &reference) {
    demand(actual.magnitude_valid&&reference.magnitude_valid,"certified_magnitude");
    check("log_matter_weight",std::abs(actual.log_abs_matter-reference.log_abs_matter),2e-10);
    check("log_phonon_weight",std::abs(actual.log_phonon-reference.log_phonon),2e-10);
    if(actual.phase_valid) {
        demand(reference.phase_valid,"reference_signed");
        check("signed_phase",std::abs(actual.matter_phase-reference.matter_phase),2e-10);
    }
    if(actual.green_valid) {
        demand(reference.green_valid,"reference_green");
        check("center_green",(actual.center_green-reference.center_green).cwiseAbs().maxCoeff(),2e-11);
    }
}
// Serialize values now, so later checks cannot be fooled by shallow copies
// sharing an accidentally mutated matrix or checkpoint backing allocation.
std::string fingerprint(const CoupledPhononPath::Snapshot &s) {
    std::ostringstream out;out<<std::setprecision(17);
    auto matrix=[&](const MatType &m){out<<m.rows()<<","<<m.cols()<<":";for(int j=0;j<m.size();++j)out<<m.data()[j]<<";";};
    for(int j=0;j<s.path.nodes();++j)for(int b=0;b<s.path.bonds();++b)out<<s.path.coordinate(j,b)<<";";
    for(const auto &f:s.factors){matrix(MatType(f.matrix));out<<f.eta;}
    for(const auto &f:s.weight_order){matrix(MatType(f.matrix));out<<f.eta;}
    out<<s.configuration_hash<<s.weight_order_is_local<<s.weight.logAbsTotal()<<s.weight.matter_phase
       <<s.weight.magnitude_valid<<s.weight.phase_valid<<s.weight.green_valid;matrix(s.weight.center_green);
    if(s.signed_prefix_cache)for(const auto &p:s.signed_prefix_cache->checkpoints) {
        out<<bool(p);if(!p)continue;matrix(p->q);
        out<<p->log_det_r<<p->previous_log_overlap<<p->det_r_phase<<p->previous_overlap_phase
           <<p->log_abs_weight<<p->complex_phase<<p->z2_sign<<p->overlap_rank<<p->overlap_rcond
           <<p->overlap_residual<<int(p->pfaffian_status)<<p->prefix_hash;
    }
    if(s.lu_cayley_cache)for(const auto &p:*s.lu_cayley_cache){out<<bool(p);if(p)matrix(p->value());}
    return out.str();
}
void unchanged(const CoupledPhononPath::Snapshot &a,const CoupledPhononPath::Snapshot &b) {
    demand(a.configuration_hash==b.configuration_hash,"rejection_configuration_hash");
    demand(a.signed_prefix_cache==b.signed_prefix_cache,"rejection_prefix_ownership");
    demand(a.lu_cayley_cache==b.lu_cayley_cache,"rejection_cayley_ownership");
    demand(a.weight_order_is_local==b.weight_order_is_local,"rejection_order_type");
    demand(a.weight.phase_valid==b.weight.phase_valid&&a.weight.green_valid==b.weight.green_valid,"rejection_validity");
    check("rejection_weight",std::abs(a.weight.logAbsTotal()-b.weight.logAbsTotal()),0.);
    check("rejection_phase",std::abs(a.weight.matter_phase-b.weight.matter_phase),0.);
    if(a.weight.green_valid)check("rejection_green",(a.weight.center_green-b.weight.center_green).norm(),0.);
    for(int j=0;j<a.path.nodes();++j)for(int bond=0;bond<a.path.bonds();++bond)
        check("rejection_coordinates",std::abs(a.path.coordinate(j,bond)-b.path.coordinate(j,bond)),0.);
    demand(a.factors.size()==b.factors.size()&&a.weight_order.size()==b.weight_order.size(),"rejection_factor_count");
    for(std::size_t j=0;j<a.factors.size();++j)
        check("rejection_factors",(MatType(a.factors[j].matrix)-MatType(b.factors[j].matrix)).norm(),0.);
    for(std::size_t j=0;j<a.weight_order.size();++j)
        check("rejection_weight_order",(MatType(a.weight_order[j].matrix)-MatType(b.weight_order[j].matrix)).norm(),0.);
}
void lazyHalfSteps() {
    const KitaevSSHModel model(HoneycombLattice(6,6,1),p2Couplings(.3));
    std::mt19937_64 rng(106109250);
    const auto path=PhononPath(model.lattice().bondCount(),4,.1,{1.7,.8}).freeSample(rng);
    const auto trial=freeGroundState(model);
    const auto generator=honeycombGenerator(model);
    for(const auto &g:{generator,MajoranaGenerator(generator.value,generator.into)}) {
        const auto canonical=prepareOwnedPath(path,trial.Phi.rows(),g,false);
        const auto lazy=prepareOwnedPath(path,trial.Phi.rows(),g,false,nullptr,false);
        demand(lazy.factors.empty()&&canonical.factors.size()==8,"lazy_preparation_shape");
        demand(lazy.weightOrder.size()==canonical.weightOrder.size(),"lazy_merged_shape");
        for(std::size_t j=0;j<lazy.weightOrder.size();++j)
            check("lazy_merged_arithmetic",(MatType(lazy.weightOrder[j].matrix)-
                MatType(canonical.weightOrder[j].matrix)).norm(),0.);
    }
    CoupledPhononPath state(path,trial,generator,4,{true,true,true});
    demand(state.state().factors.empty(),"initial_lazy_factors");
    const auto initial=state.state();
    const auto initialFingerprint=fingerprint(initial);
    auto rejected=state.propose(2,0,path.coordinate(2,0)+.1);
    demand(rejected.candidate().factors.empty(),"local_lazy_factors");
    state.reject(rejected);
    demand(fingerprint(state.state())==initialFingerprint,"lazy_rejection_deep_values");
    weights(state.ensureMeasured().weight,state.uncachedReference(path,true).weight);
    demand(state.state().factors.size()==8,"measurement_materializes_factors");
    auto measuredLocal=state.propose(2,0,path.coordinate(2,0)+.05);
    state.accept(measuredLocal);
    demand(state.state().factors.empty(),"measured_to_lazy_local");
    weights(state.ensureMeasured().weight,state.uncachedReference(state.state().path,true).weight);
    // A successfully recovered candidate owns canonical h*h bookkeeping.
    // Its next ordinary local update must safely return to merged-only form.
    auto recovered=state.propose(0,1,path.coordinate(0,1)+.03,true);
    state.accept(recovered);
    demand(state.state().factors.size()==8,"recovery_canonical_factors");
    auto afterRecovery=state.propose(4,2,path.coordinate(4,2)+.02);
    state.accept(afterRecovery);
    demand(state.state().factors.empty(),"recovered_to_lazy_local");
    weights(state.ensureMeasured().weight,state.uncachedReference(state.state().path,true).weight);
    const auto changed=state.state().path.withCoordinate(2,3,path.coordinate(2,3)+.04);
    auto collective=state.proposeQModePath(changed);
    state.accept(collective);
    demand(state.state().factors.empty(),"collective_lazy_factors");
    weights(state.ensureMeasured().weight,state.uncachedReference(changed,true).weight);
    // Empty factors are accepted only with a certified complete merged order.
    bool failed=false;
    try {evaluatePathMagnitude(path,trial,{},4);}catch(const std::invalid_argument&){failed=true;}
    demand(failed,"missing_canonical_factors_rejected");
    const std::vector<PureProjectorSlice> incomplete;
    failed=false;
    try {evaluatePathMagnitude(path,trial,{},4,&incomplete,nullptr,true);}
    catch(const std::invalid_argument&){failed=true;}
    demand(failed,"incomplete_merged_factors_rejected");
    std::cout<<"PASS L6 lazy half-step measurement/recovery/local/collective transitions\n";
}
void run(int shear,int block,PhononWeightOptions options,int &accepted,int &rejected) {
    const int Ntau=8;const double dtau=.125,tau=Ntau*dtau,tau_mid=tau/2;
    KitaevSSHModel model(HoneycombLattice(3,2,shear),p2Couplings(.3));
    std::mt19937_64 rng(802715+shear);
    auto path=PhononPath(model.lattice().bondCount(),Ntau,dtau,{1.7,.8}).freeSample(rng);
    CoupledPhononPath state(path,freeGroundState(model),honeycombGenerator(model),block,options);
    PhononSampler sampler(state,rng);
    weights(state.state().weight,state.uncachedReference(path,true).weight);
    for(int step=0;step<18;++step) {
        const auto before=state.state();
        const auto beforeFingerprint=fingerprint(before);
        const int node=step%3==0?0:(step%3==1?Ntau/2:Ntau),bond=(step*7)%path.bonds();
        std::vector<double> values(path.nodes());
        for(int j=0;j<path.nodes();++j)
            values[j]=before.path.coordinate(j,bond)+.17*std::sin(.7*(step+1)+.3*j);
        const bool whole=step%4==3;
        const auto changed=whole?before.path.withBondCoordinates(bond,values):before.path.withCoordinate(node,bond,values[node]);
        const auto reference=state.uncachedReference(changed,true);
        auto proposal=whole?state.proposeBondPath(bond,values):state.propose(node,bond,values[node]);
        demand(fingerprint(state.state())==beforeFingerprint,"proposal_does_not_mutate_accepted_cache");
        weights(proposal.candidate().weight,reference.weight);
        const auto oldReference=state.uncachedReference(before.path,true);
        const double expected=reference.weight.log_abs_matter-oldReference.weight.log_abs_matter;
        const double actual=proposal.candidate().weight.log_abs_matter-before.weight.log_abs_matter;
        check("proposal_log_ratio",std::abs(actual-expected),4e-10);
        check("Metropolis_probability",std::abs(std::exp(std::min(0.,actual))-std::exp(std::min(0.,expected))),4e-10);
        if(proposal.candidate().weight.phase_valid&&before.weight.phase_valid)
            check("proposal_signed_phase_ratio",std::abs(proposal.phaseRatio()-reference.weight.matter_phase/oldReference.weight.matter_phase),4e-10);
        const auto candidateHash=proposal.candidate().configuration_hash;
        const auto candidatePrefix=proposal.candidate().signed_prefix_cache;
        const auto candidateCayley=proposal.candidate().lu_cayley_cache;
        if(step%3==0) {
            state.reject(proposal);++rejected;unchanged(before,state.state());
            demand(fingerprint(state.state())==beforeFingerprint,"rejected_cache_deep_values");
        } else {
            state.accept(proposal);++accepted;
            demand(state.configurationHash()==candidateHash,"accept_configuration_hash");
            demand(state.state().signed_prefix_cache==candidatePrefix,"accept_prefix_ownership");
            demand(state.state().lu_cayley_cache==candidateCayley,"accept_cayley_ownership");
            weights(state.ensureMeasured().weight,reference.weight);
        }
        weights(state.state().weight,state.uncachedReference(state.state().path,true).weight);
        // Run real MH too: both acceptance and general proposal-density ratio
        // must agree with a complete reconstruction of this same candidate.
        const auto prior=state.state();
        const auto priorFingerprint=fingerprint(prior);
        std::vector<double> mhValues(prior.path.nodes());
        for(int j=0;j<prior.path.nodes();++j)mhValues[j]=prior.path.coordinate(j,bond)+.09*std::cos(step+.2*j);
        const auto mhPath=whole?prior.path.withBondCoordinates(bond,mhValues):prior.path.withCoordinate(node,bond,mhValues[node]);
        const auto mhRef=state.uncachedReference(mhPath,true);
        const double ratio=mhRef.weight.log_abs_matter-state.uncachedReference(prior.path,true).weight.log_abs_matter;
        const double u=step%2?.999999:1e-12;
        const auto move=whole?sampler.attemptBond(bond,mhValues,u):sampler.attempt(node,bond,mhPath.coordinate(node,bond),u);
        check("sampler_log_ratio",std::abs(move.log_ratio-ratio),4e-10);
        check("sampler_general_MH",std::abs(move.general_log_ratio-ratio),4e-10);
        demand(move.accepted==(std::log(u)<std::min(0.,ratio)),"sampler_decision");
        if(move.accepted){++accepted;weights(state.ensureMeasured().weight,mhRef.weight);}
        else {++rejected;unchanged(prior,state.state());
            demand(fingerprint(state.state())==priorFingerprint,"MH_rejected_cache_deep_values");}
    }
    std::cout<<"PASS cache sequence shear="<<shear<<" block="<<block<<" tau="<<tau<<" dtau="<<dtau<<" Ntau="<<Ntau<<" tau_mid="<<tau_mid
        <<" svd="<<options.singular_values_only<<" merge="<<options.merge_node_half_steps<<" defer="<<options.defer_phase_measurement
        <<" prefix="<<options.reuse_signed_prefix_checkpoints<<" cayley="<<options.reuse_lu_cayley<<'\n';
}
}
int main() {
    try {
        int accepted=0,rejected=0;
        lazyHalfSteps();
        const std::vector<PhononWeightOptions> settings{{false,false},{true,false},{true,true},{true,true,false,true},
            {true,true,false,true,true},{true,true,true,true,true}};
        for(int shear:{0,1})for(int block:{1,3})for(auto options:settings)run(shear,block,options,accepted,rejected);
        demand(accepted>0&&rejected>0,"both_transaction_outcomes");
        for(const auto &m:metrics)std::cout<<"PASS "<<m.first<<" reference_error=0 max_PfQMC_error="<<m.second.error<<" tolerance="<<m.second.tolerance<<'\n';
        std::cout<<"PASS cache transactions accepted="<<accepted<<" rejected="<<rejected<<'\n';
        return 0;
    }catch(const std::exception &e){std::cerr<<"FAIL cache validation: "<<e.what()<<'\n';return 1;}
}
