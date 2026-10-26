#include "phonon_mc.h"
#include "numerical_policy.h"
#include "physical_matter_trial.h"
#include "vbs_observables.h"
#include "vbs_susceptibility.h"
#include "vbs_amplitude.h"
#include "vbs_field.h"
#include "vbs_candidates.h"
#include "fermion_observables.h"
#include "bond_diagnostics.h"
#include "total_energy.h"
#include <chrono>
#include <ctime>
#include "output_path.h"
#include "mc_checkpoint.h"
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>

using namespace kitaev_ssh;
namespace {
// JSON diagnostics contain replay inputs and scalar certificates, never matrices.
void jsonString(std::ostream &out,const std::string &s) {
    out<<'"';
    for(unsigned char c:s) {
        if(c=='"'||c=='\\') out<<'\\'<<char(c);
        else if(c<32) { const char *hex="0123456789abcdef";out<<"\\u00"<<hex[c>>4]<<hex[c&15]; }
        else out<<char(c);
    }
    out<<'"';
}
void writeResult(const std::string &prefix,const std::string &json) {
    const std::string temp=prefix+"_results.json.tmp";
    std::ofstream out(temp);out<<json;out.close();
    if(!out || std::rename(temp.c_str(),(prefix+"_results.json").c_str())!=0)
        throw std::runtime_error("atomic results write failed");
}
struct RunNumerics {
    std::map<std::string,MoveCounts> stages;
    std::map<std::string,int> saved;
    std::string argv_json;
    bool owns_output=false;
    void metadata(std::ostream &out) const {
        std::uint64_t total=0;for(const auto &entry:stages)total+=entry.second.first_failures;
        out<<"\"numerical_failure_policy\":\""<<numericalFailurePolicy<<"\",\"numerical_failures\":"<<total
           <<",\"numerical_failures_definition\":\"first failed numerical evaluations, including recovered and initialization attempts; retries counted separately\""
           <<",\"numerical_diagnostics\":{";
        bool comma=false;
        for(const std::string stage:{"initialization","warmup","measurement"}) {
            if(comma)out<<',';comma=true;jsonString(out,stage);
            const auto it=stages.find(stage);const MoveCounts c=it==stages.end()?MoveCounts{}:it->second;
            out<<":{\"first_failures\":"<<c.first_failures<<",\"reference_recovered\":"<<c.reference_recovered
               <<",\"reference_failed\":"<<c.reference_failed<<'}';
        }
        out<<"},\"argv\":"<<argv_json;
    }
    void status(const std::string &prefix,const char *status,const std::string &error="") const {
        std::ostringstream out;out<<"{\"status\":\""<<status<<"\",\"error\":";jsonString(out,error);out<<',';
        metadata(out);out<<",\"partial_blocks_are_not_completed_results\":true}\n";writeResult(prefix,out.str());
    }
    void record(const std::string &prefix,const NumericalEvent &e,const PhononPath &path) {
        auto &c=stages[e.stage];++c.first_failures;
        if(e.recovered)++c.reference_recovered;else ++c.reference_failed;
        const bool fatal=!e.recovered&&e.stage!="initialization";
        if(saved[e.stage]>=3&&!fatal)return;
        ++saved[e.stage];
        std::ofstream out(prefix+"_numerical.jsonl",std::ios::app);
        out<<std::setprecision(17)<<"{\"stage\":";jsonString(out,e.stage);
        out<<",\"operation\":";jsonString(out,e.operation);
        out<<",\"reference_recovered\":"<<(e.recovered?"true":"false")<<",\"proposal_index\":"<<e.proposal_index
           <<",\"node\":"<<e.node<<",\"bond\":"<<e.bond<<",\"acceptance_u\":"<<e.u<<",\"candidate_values\":[";
        for(std::size_t i=0;i<e.values.size();++i){if(i)out<<',';out<<e.values[i];}
        out<<"],\"first_error\":";jsonString(out,e.first_error);
        out<<",\"reference_error\":";jsonString(out,e.reference_error);
        out<<",\"path_before_proposal_or_initial_draw\":[";
        for(int n=0;n<path.nodes();++n){if(n)out<<',';out<<'[';for(int b=0;b<path.bonds();++b){if(b)out<<',';out<<path.coordinate(n,b);}out<<']';}
        out<<"],\"intervals\":"<<path.intervals()<<",\"dt\":"<<path.dt()<<",\"mass\":"<<path.parameters().mass
           <<",\"omega_ph\":"<<path.parameters().omega<<",\"trial_shifts\":[";
        for(int b=0;b<path.bonds();++b){if(b)out<<',';out<<path.trialShift(b);}out<<"],";
        metadata(out);out<<"}\n";out.close();
        if(!out)throw std::runtime_error("numerical diagnostic write failed");
        std::cerr<<"NUMERICAL "<<e.stage<<' '<<e.operation<<" reference="<<(e.recovered?"recovered":"failed")
                 <<" proposal="<<e.proposal_index<<" diagnostics="<<prefix<<"_numerical.jsonl\n";
    }
};
void number(std::ostream &out,double x) { if(std::isfinite(x)) out<<x; else out<<"null"; }

// Process totals include initialization and numerical retries. These counters
// describe arithmetic recovery, not Monte Carlo equilibration or acceptance.
void writeRoundoffDiagnostics(std::ostream &out) {
    const auto &q=pureCanonicalCorrectionStats;
    out<<",\"canonical_frame_roundoff\":{\"scope\":\"process_lifetime\""
       <<",\"complex_qr_calls\":"<<q.complex_calls<<",\"real_qr_calls\":"<<q.real_calls
       <<",\"extended_frame_calls\":"<<q.extended_calls<<",\"failures\":"<<q.failures;
    for(const auto &item:std::vector<std::pair<const char*,double>>{
        {"maximum_relative_correction",q.maximum_relative_correction},
        {"maximum_input_isotropy",q.maximum_input_isotropy},
        {"maximum_output_isotropy",q.maximum_output_isotropy},
        {"maximum_output_orthonormality",q.maximum_output_orthonormality}}) {
        out<<",\""<<item.first<<"\":";number(out,item.second);
    }
    const auto &g=pureProjectorGreenPrecisionStats;
    out<<"},\"green_precision_recovery\":{\"scope\":\"process_lifetime\""
       <<",\"stored_matrix_ld_checks\":"<<g.long_double_residual_checks
       <<",\"stored_matrix_ld_passes\":"<<g.long_double_residual_passes
       <<",\"canonical_ld_rebuild_attempts\":"<<g.canonical_ld_rebuild_attempts
       <<",\"canonical_ld_rebuild_assemblies\":"<<g.canonical_ld_rebuild_assemblies
       <<",\"canonical_ld_rebuild_certified\":"<<g.canonical_ld_rebuild_certified
       <<",\"maximum_relative_frame_correction\":";
    number(out,g.maximum_canonical_relative_correction);
    const auto &p=pure_projector_augmented_pf::statistics;
    out<<"},\"augmented_phase_recovery\":{\"scope\":\"process_lifetime\""
       <<",\"attempts\":"<<p.attempts<<",\"accepted\":"<<p.accepted
       <<",\"failed\":"<<p.failed<<",\"raw_b_identity_failures\":"<<p.raw_b_identity_failures
       <<",\"unbounded_phase\":"<<p.unbounded_phase;
    for(const auto &item:std::vector<std::pair<const char*,double>>{
        {"maximum_phase_bound_radians",p.max_phase_bound},
        {"maximum_bound_to_branch_margin",p.max_branch_fraction},
        {"maximum_raw_b_log_error",p.max_raw_b_log_error},
        {"maximum_raw_b_phase_error",p.max_raw_b_phase_error},
        {"maximum_same_c_log_error",p.max_same_c_log_error},
        {"maximum_same_c_phase_error",p.max_same_c_phase_error},
        {"maximum_relative_frame_correction",p.max_relative_frame_correction},
        {"maximum_augmented_condition_estimate",p.max_k_condition},
        {"maximum_overlap_condition_estimate",p.max_m_condition},
        {"maximum_cayley_factor_difference",p.max_cayley_factor_difference},
        {"seconds",p.total_seconds}}) {
        out<<",\""<<item.first<<"\":";number(out,item.second);
    }
    out<<'}';
}
void estimate(std::ostream &out,const ComplexEstimate &e,bool defined=true) {
    out<<"{\"value\":["; if(defined) out<<e.value.real()<<','<<e.value.imag(); else out<<"null,null";
    out<<"],\"stderr_real\":"; number(out,e.error_real); out<<",\"stderr_imag\":"; number(out,e.error_imag);
    out<<",\"cov_real_imag\":"; number(out,e.covariance); out<<'}';
}
// SE is the joint complex mean uncertainty sqrt(Var(Re)+Var(Im)).
// The modulus SE separately reuses the same delete-one-block replicas.
void phaseDiagnostics(std::ostream &out,const PhaseSummary &s) {
    const double magnitude=std::abs(s.phase.value);
    const double se=std::hypot(s.phase.error_real,s.phase.error_imag);
    out<<"{\"positive_count\":"<<s.positive<<",\"negative_count\":"<<s.negative
       <<",\"complex_count\":"<<s.nonreal<<",\"reality_tolerance\":"<<phaseRealityTolerance
       <<",\"magnitude\":"; number(out,magnitude);
    out<<",\"stderr_complex\":"; number(out,se);
    out<<",\"stderr_magnitude\":"; number(out,s.phase_magnitude_error);
    out<<",\"denominator_snr\":";
    number(out,se>0?magnitude/se:std::numeric_limits<double>::quiet_NaN());
    out<<",\"snr_status\":\""<<(std::isfinite(se)?(se>0?"finite":(magnitude>0?"zero_observed_variance":"undefined_zero_mean")):"insufficient_blocks")<<"\"}";
}
void estimates(std::ostream &out,const std::vector<std::string> &names,const PhaseSummary &summary) {
    out<<'{';
    for(std::size_t n=0;n<names.size();++n) {
        if(n)out<<',';out<<'"'<<names[n]<<"\":";
        estimate(out,n<summary.observables.size()?summary.observables[n]:ComplexEstimate{},summary.ratio_defined);
    }
    out<<'}';
}
// Categorical/count entries (winner index, participant count) are labels, not
// continuous observables: their values are written but their jackknife
// "errors" are suppressed (null), because reordering the q-file changes the
// numbers without changing any physics. The selection-stability diagnostics
// live in the structured "candidates" JSON section instead.
bool isCategoricalDerivedName(const std::string &name) {
    return name=="candidate_max_winner_index" || name=="candidate_participant_count";
}
void estimatesDerived(std::ostream &out,const std::vector<std::string> &names,const PhaseSummary &summary) {
    out<<'{';
    for(std::size_t n=0;n<names.size();++n) {
        if(n)out<<',';out<<'"'<<names[n]<<"\":";
        if(n<summary.observables.size() && isCategoricalDerivedName(names[n])) {
            auto e=summary.observables[n];
            e.error_real=e.error_imag=e.covariance=std::numeric_limits<double>::quiet_NaN();
            estimate(out,e,summary.ratio_defined);
        } else {
            estimate(out,n<summary.observables.size()?summary.observables[n]:ComplexEstimate{},summary.ratio_defined);
        }
    }
    out<<'}';
}
void qPointJson(std::ostream &json,const VBSQPoint &q,int j) {
    json<<"{\"index\":"<<j<<",\"label\":\""<<q.label<<"\""
        <<",\"input_fraction\":["<<q.u_input<<','<<q.v_input<<"]"
        <<",\"reciprocal_fraction\":["<<q.u<<','<<q.v<<"]"
        <<",\"phase_fraction\":["<<q.phase_u<<','<<q.phase_v<<"]"
        <<",\"cartesian\":["<<q.cartesian.x<<','<<q.cartesian.y<<"]"
        <<",\"on_grid\":"<<(q.on_grid?"true":"false");
    if(q.on_grid) json<<",\"n\":["<<q.n1<<','<<q.n2<<"]"; else json<<",\"n\":null";
    json<<",\"grid_snap_residual\":"<<q.snap_residual
        <<",\"W_L\":["<<q.window.real()<<','<<q.window.imag()<<"]"
        <<",\"W_L_abs2\":"<<std::norm(q.window)<<"}";
}
void writeDoubleOrNull(std::ostream &json,double x) {
    if(std::isfinite(x)) json<<x; else json<<"null";
}
void writeCandidatesJson(std::ostream &json,const VBSObservables &vbs,
                         const CandidateDiagnostics &diag) {
    json<<",\"candidates\":{\"schema\":\"ivbs-candidates-v2\",\"schema_revision\":2,\"level\":\"finest_blocks\""
        <<",\"ratio_defined\":"<<(diag.ratio_defined?"true":"false")
        <<",\"center_defined\":"<<(diag.center_defined?"true":"false")
        <<",\"center_reason\":\""<<diag.center_reason<<"\""
        <<",\"jackknife_defined\":"<<(diag.jackknife_defined?"true":"false")
        <<",\"jackknife_reason\":\""<<diag.jackknife_reason<<"\""
        <<",\"stability_defined\":"<<(diag.stability_defined?"true":"false")
        <<",\"stability_reason\":\""<<diag.stability_reason<<"\""
        <<",\"total_replicas\":"<<diag.total_replicas
        <<",\"valid_replicas\":"<<diag.valid_replicas
        <<",\"invalid_replicas\":"<<diag.invalid_replicas
        <<",\"invalid_replica_reasons\":[";
    for(std::size_t r=0;r<diag.invalid_replica_reasons.size();++r) {
        if(r) json<<',';
        json<<'\"'<<diag.invalid_replica_reasons[r]<<'\"';
    }
    json<<"]"
        <<",\"excluded_Q0_rule\":\"mod reciprocal-lattice normalization + on_grid test; (n1,n2)==(0,0); "
        "a momentum congruent to 0 is never off-grid\""
        <<",\"scoring\":\"real part of the reweighted ensemble mz2=Szz/Nc; imaginary parts kept as "
        "diagnostics; no |value|, no clipping of negatives\""
        <<",\"tie_tolerance\":{\"abs\":"<<candidateTieAbs<<",\"rel\":"<<candidateTieRel<<"}"
        <<",\"winner_set_rule\":\"all candidates within tie tolerance of the group maximum; "
        "no random tie-breaking; winner_display_index is the input-order first winner, a display "
        "number without physical error\""
        <<",\"stability_definition\":\"delete-one-block replicas re-select peaks from their own "
        "reweighted ratios; frequencies are resampling selection-stability diagnostics, not "
        "probabilities of a phase or confidence intervals; tie-membership frequencies may sum to "
        "more than 1; invalid replicas are counted, not silently dropped\""
        <<",\"legacy\":{\"candidate_max_mz2\":\"mixed grid+off-grid raw maximum over non-Q=0 "
        "candidates, definition unchanged from the previous release; kept for backward "
        "compatibility, prefer candidates.groups\"}"
        <<",\"groups\":[";
    for(std::size_t g=0;g<diag.groups.size();++g) {
        const auto &G=diag.groups[g];
        if(g)json<<',';
        json<<"{\"name\":\""<<G.name<<"\",\"estimator\":\""<<G.estimator<<"\""
            <<",\"available\":"<<(G.available?"true":"false")
            <<",\"participants\":[";
        for(std::size_t t=0;t<G.participants.size();++t){if(t)json<<',';json<<G.participants[t];}
        json<<"],\"participant_labels\":[";
        for(std::size_t t=0;t<G.participant_labels.size();++t){if(t)json<<',';json<<'"'<<G.participant_labels[t]<<'"';}
        json<<"],\"center_defined\":"<<(G.center_defined?"true":"false")
            <<",\"max_mz2_real\":";writeDoubleOrNull(json,G.max_mz2_real);
        json<<",\"winner_display_index\":";
        if(G.winner_display_index>=0) json<<G.winner_display_index; else json<<"null";
        json<<",\"winner_set\":[";
        for(std::size_t t=0;t<G.winner_set.size();++t){if(t)json<<',';json<<G.winner_set[t];}
        json<<"],\"winner_labels\":[";
        for(std::size_t t=0;t<G.winner_labels.size();++t){if(t)json<<',';json<<'"'<<G.winner_labels[t]<<'"';}
        json<<"],\"center_reason\":\""<<G.center_reason<<"\""
            <<",\"total_replicas\":"<<G.total_replicas
            <<",\"valid_replicas\":"<<G.valid_replicas<<",\"invalid_replicas\":"<<G.invalid_replicas
            <<",\"jackknife_defined\":"<<(G.jackknife_defined?"true":"false")
            <<",\"jackknife_reason\":\""<<G.jackknife_reason<<"\""
            <<",\"stability_defined\":"<<(G.stability_defined?"true":"false")
            <<",\"stability_reason\":\""<<G.stability_reason<<"\""
            <<",\"conditional_on_valid_replicas\":"<<(G.conditional_on_valid_replicas?"true":"false")
            <<",\"unique_win_counts\":[";
        for(std::size_t t=0;t<G.unique_win_counts.size();++t){if(t)json<<',';json<<G.unique_win_counts[t];}
        json<<"],\"tie_membership_counts\":[";
        for(std::size_t t=0;t<G.tie_membership_counts.size();++t){if(t)json<<',';json<<G.tie_membership_counts[t];}
        json<<"],\"unique_win_frequency\":[";
        for(std::size_t t=0;t<G.unique_win_frequency.size();++t){if(t)json<<',';writeDoubleOrNull(json,G.unique_win_frequency[t]);}
        json<<"],\"tie_membership_frequency\":[";
        for(std::size_t t=0;t<G.tie_membership_frequency.size();++t){if(t)json<<',';writeDoubleOrNull(json,G.tie_membership_frequency[t]);}
        json<<"],\"class_tie_counts\":[";
        for(std::size_t t=0;t<G.class_tie_counts.size();++t){if(t)json<<',';json<<G.class_tie_counts[t];}
        json<<"],\"class_full_set_counts\":[";
        for(std::size_t t=0;t<G.class_full_set_counts.size();++t){if(t)json<<',';json<<G.class_full_set_counts[t];}
        json<<"],\"class_tie_frequency\":[";
        for(std::size_t t=0;t<G.class_tie_frequency.size();++t){if(t)json<<',';writeDoubleOrNull(json,G.class_tie_frequency[t]);}
        json<<"],\"class_full_set_frequency\":[";
        for(std::size_t t=0;t<G.class_full_set_frequency.size();++t){if(t)json<<',';writeDoubleOrNull(json,G.class_full_set_frequency[t]);}
        json<<"]}";
    }
    json<<"],\"points\":[";
    for(std::size_t j=0;j<diag.points.size();++j) {
        const auto &p=diag.points[j];
        if(j)json<<',';
        json<<"{\"index\":"<<p.index<<",\"label\":\""<<p.label<<"\",\"group\":\""<<p.group<<"\""
            <<",\"u\":"<<p.u<<",\"v\":"<<p.v<<",\"on_grid\":"<<(p.on_grid?"true":"false")
            <<",\"W_L_abs2\":"<<p.W_abs2
            <<",\"mz2_raw\":{\"defined\":"<<(p.raw_center_defined?"true":"false")
            <<",\"reason\":\""<<p.raw_center_reason<<"\",\"real\":";writeDoubleOrNull(json,p.mz2_raw_real);
        json<<",\"imag_diag\":";writeDoubleOrNull(json,p.mz2_raw_imag);
        json<<",\"stderr_real\":";writeDoubleOrNull(json,p.mz2_raw_error_real);
        json<<",\"stderr_imag\":";writeDoubleOrNull(json,p.mz2_raw_error_imag);
        json<<",\"cov_real_imag\":";writeDoubleOrNull(json,p.mz2_raw_covariance);
        json<<",\"jackknife_defined\":"<<(p.raw_jackknife_defined?"true":"false")
            <<",\"jackknife_reason\":\""<<p.raw_jackknife_reason<<"\""
            <<"},\"mz2_projected\":{\"defined\":"<<(p.projected_center_defined?"true":"false")
            <<",\"reason\":\""<<p.projected_center_reason<<"\",\"real\":";writeDoubleOrNull(json,p.mz2_proj_real);
        json<<",\"imag_diag\":";writeDoubleOrNull(json,p.mz2_proj_imag);
        json<<",\"stderr_real\":";writeDoubleOrNull(json,p.mz2_proj_error_real);
        json<<",\"stderr_imag\":";writeDoubleOrNull(json,p.mz2_proj_error_imag);
        json<<",\"cov_real_imag\":";writeDoubleOrNull(json,p.mz2_proj_covariance);
        json<<",\"jackknife_defined\":"<<(p.projected_jackknife_defined?"true":"false")
            <<",\"jackknife_reason\":\""<<p.projected_jackknife_reason<<"\""
            <<"},\"plus_minus_class\":";
        if(p.plus_minus_class>=0) json<<p.plus_minus_class; else json<<"null";
        json<<",\"class_self_conjugate\":"<<(p.class_self_conjugate?"true":"false")<<'}';
    }
    json<<"],\"plus_minus_classes\":[";
    for(std::size_t c=0;c<diag.classes.size();++c) {
        const auto &k=diag.classes[c];
        if(c)json<<',';
        json<<"{\"id\":"<<c<<",\"members\":[";
        for(std::size_t t=0;t<k.members.size();++t){if(t)json<<',';json<<k.members[t];}
        json<<"],\"member_labels\":[";
        for(std::size_t t=0;t<k.members.size();++t){if(t)json<<',';json<<'"'<<diag.points[k.members[t]].label<<'"';}
        json<<"],\"self_conjugate\":"<<(k.self_conjugate?"true":"false")
            <<",\"canonical_fraction\":["<<k.u<<','<<k.v<<"]}";
    }
    json<<"]}";
}
}
int main(int argc,char **argv) {
    std::string prefix="kitaev_ssh_mc";
    RunNumerics numerics;
    {std::ostringstream out;out<<'[';for(int i=0;i<argc;++i){if(i)out<<',';jsonString(out,argv[i]);}out<<']';numerics.argv_json=out.str();}
    try {
        const auto start=std::chrono::steady_clock::now();
        const auto cpuStart=std::clock();
        std::map<std::string,std::string> args;
        const std::set<std::string> allowed{"Lx","Ly","shear","Jx","Jy","Jz","intervals","dt","mass","omega-ph","g","seed",
            "warmup","sweeps","measure-every","block-samples","qr-block","output","vbs","vbs-susceptibility","vbs-amplitude","vbs-amplitude-response","vbs-mode","vbs-q-file",
            "proposal","measure-energy","trial-dx","trial-dy","trial-dz","weight-engine","defer-phase","prefix-cache","cayley-cache","fermions","warmup-observable-every","acceptance-every","q-mode-every","q-mode-n1","q-mode-n2","q-mode-beta","q-diagnostic-every",
            "vbs-field-h","vbs-field-n1","vbs-field-n2","vbs-field-phase",
            "bond-beta","bond-guide","bond-guide-x","bond-guide-y","bond-guide-z","q-mode-pool","q-mode-bonds","center-green-fast","numerics","phonon-diagnostics","bond-diagnostics","validation-series","path-series-bond","matter-trial","gamma-every","gamma-beta","q-mode-extra-pool","initial-path","mixing-diagnostics","checkpoint-file","checkpoint-every","resume","measure-total-energy"};
        for(int i=1;i<argc;++i) {
            const std::string key=argv[i];
            if(key=="--help") {
                std::cout<<"Static coupled |W| MC with complex phase reweighting; fixed-gauge matter sector.\n"
                    <<"--matter-trial physical|unrestricted (default physical): fixed free trial in the geometry/gauge-compatible matter parity; unrestricted is a legacy control.\n"
                    <<"--Lx 2 --Ly 2 --shear 0 --Jx 1 --Jy 1 --Jz 0.5 --intervals 4 --dt 0.2\n"
                    <<"--mass 4 --omega-ph 0.5 --g 0.3 --seed 42 --warmup 100 --sweeps 400\n"
                    <<"--measure-every 2 --block-samples 10 --qr-block 4 --vbs 1 --output PREFIX\n"
                    <<"--bond-diagnostics 0|1 appends center K_b=i*c_A*c_B and X_b*K_b per physical bond and direction; no gauge factor in K_b (default 1; 0 disables).\n"
                    <<"--phonon-diagnostics 1 appends node X/X2, center-node same/neighbor-bond correlations and matter log-weight to raw phase blocks (default 0).\n"
                    <<"--vbs-susceptibility 1 saves finite-projection static response moments (requires --vbs 1 --vbs-mode all).\n"
                    <<"--vbs-amplitude 1 saves A=|D|, A2 and per-configuration A4 for all allowed q and bonds; requires --vbs 1 --vbs-mode all.\n"
                    <<"--vbs-amplitude-response 1 additionally saves IA=sum_j t_j |D(node_j)| and A(center)*IA, with trapezoidal weights; requires --vbs-amplitude 1. No physical source is applied.\n"
                    <<"--vbs-field-h h --vbs-field-n1 n1 --vbs-field-n2 n2 --vbs-field-phase phi: H_h=H-h sum_r cos(q.r+phi) X_z(r), fixed trial; requires VBS susceptibility. Explicit h=0 measures the same selected channel.\n"
                    <<"--vbs 1 is the master VBS switch (default 0). With --vbs 1:\n"
                    <<"  --vbs-mode all      (default) full PBC momentum grid via FFT and full real-space C;\n"
                    <<"  --vbs-mode selected --vbs-q-file q.txt: only the listed momenta (finite-window\n"
                    <<"                       Fourier estimator; no full FFT/real-space C; phonon stays PBC);\n"
                    <<"  --vbs-mode both      full grid block followed by the selected block.\n"
                    <<"q-file lines: 'label u v' with Q=u*b1+v*b2 (fractional reciprocal coords);\n"
                    <<"blank lines and '#' comments are allowed. Off-grid points are reported, never snapped\n"
                    <<"to fabricated integer indices. --shear is an integer tilted torus, not a U(1) twist.\n"
                    <<"--fermions 1 adds all cell-displacement AA/AB/BA/BB bare and directed gauge-string Majorana correlations (default 0).\n"
                    <<"--measure-energy 0 omits Majorana energy; center Green is still required when --fermions 1 or --bond-diagnostics 1 is enabled.\n"
                    <<"--acceptance-every N flushes stage-cumulative acceptance to PREFIX_acceptance.csv every N completed sweeps (0 disables).\n"
                    <<"--gamma-every N --gamma-beta beta: optional separate xyz real Gamma path updates (default off). --q-mode-extra-pool file preserves primary cadence. --initial-path file loads matching coordinates only. --mixing-diagnostics 1 writes all-q raw traces and formal counts.\n"
                    <<"--q-mode-every N --q-mode-n1 n1 --q-mode-n2 n2 --q-mode-beta beta: optional nonzero allowed z-bond spatial pCN mode (default interval 0, beta 1).\n"
                    <<"--q-diagnostic-every N writes raw measurement diagnostics to PREFIX_q_diagnostic.csv; N must be a multiple of measure-every (default 0).\n"
                    <<"--warmup-observable-every N records raw per-configuration all-q-max every N warmup sweeps (diagnostic only).\n"
                    <<"--proposal local (default), bond-refresh (legacy), or bond-guided (fixed-reference pCN).\n"
                    <<"--validation-series 1 writes compact unweighted values plus common phase each measurement; --path-series-bond b appends its full time path. Default off.\n"
                    <<"--bond-beta beta (default 1); --bond-guide bare|zero (default bare); --bond-guide-x/y/z f overrides. Proposal-only forces, NOT external fields or changed trial states.\n"
                    <<"--q-mode-bonds z|xyz|... (default z); --q-mode-pool file: fixed integer n1 n2 pairs, one uniformly chosen per q update; no Gamma or conjugate duplicates. With a pool, q diagnostics use its first entry.\n"
                    <<"Numerics are built in: combined weight, deferred signed phase, applicable prefix/Cayley caches and fast center Green, with unchanged strict gates.\n"
                    <<"--numerics reference: diagnostic unmerged/eager configuration with caches and fast center Green off (default production).\n"
                    <<"Legacy per-feature numerical options remain accepted for diagnostic/compatibility overrides; ordinary production scripts need none.\n"
                    <<"--trial-dx 0 --trial-dy 0 --trial-dz 0: equal bra/ket endpoint shifts, uniform by bond direction; not static X or a Hamiltonian shift.\n"
                    <<"Numerical policy: same candidate, one uncached reference retry, then abort; never silently reject.\n"
                    <<"Writes PREFIX_blocks.csv, PREFIX_results.json and representative PREFIX_numerical.jsonl (refuses overwrite).\n";
                std::cout<<"SIGUSR1 saves and continues; SIGUSR2/SIGTERM save and stop (exit75), at a completed sweep boundary.\n"
                    <<"--checkpoint-file FILE.dat (default PREFIX_checkpoint.dat); --checkpoint-every N (default0, signals only).\n"
                    <<"--measure-total-energy 0|1 (default 0): six extensive energy components at symmetric adjacent phonon-link midpoints; legacy majorana_energy stays at the center node.\n"
                    <<"--resume FILE --output NEW_PREFIX restores path, RNG/distribution, sweep progress and raw statistics. Original options are inherited and locked; --sweeps sets the TOTAL target.\n"
                    <<"Resume requires the identical executable/runtime; caches rebuild with strict numerical gates. Bitwise interacting trajectories across cache rebuilds are not guaranteed.\n";
                return 0;
            }
            if(key.rfind("--",0)!=0 || !allowed.count(key.substr(2)) || args.count(key.substr(2)) || i+1>=argc)
                throw std::invalid_argument("unknown/repeated/incomplete argument: "+key);
            args[key.substr(2)]=argv[++i];
        }
        std::unique_ptr<MCCheckpoint> resumed;
        if(args.count("resume")) {
            resumed=std::make_unique<MCCheckpoint>(MCCheckpoint::load(args.at("resume")));
            checkpointRequire(resumed->executable==checkpointExecutableIdentity(),"resume executable differs from saved binary");
            checkpointMergeArguments(args,*resumed);
        }
        if(args.count("output")) prefix=args.at("output");
        const std::string checkpointFile=args.count("checkpoint-file")?args.at("checkpoint-file"):prefix+"_checkpoint.dat";
        // All ordinary output streams end in .csv/.json/.jsonl/.path. Reserve
        // .dat so a mistaken --checkpoint-file cannot overwrite any of them.
        checkpointRequire(checkpointFile.size()>=4 && checkpointFile.compare(checkpointFile.size()-4,4,".dat")==0,"checkpoint file must end in .dat (separate from ordinary outputs)");
        std::map<std::string,std::string> checkpointInputs;
        auto inputText=[&](const std::string &key) {
            const auto found=resumed?resumed->inputs.find(key):checkpointInputs.end();
            auto text=resumed && found!=resumed->inputs.end()?found->second:checkpointReadFile(args.at(key));
            checkpointInputs[key]=text;return text;
        };
        auto real=[&](const std::string &key,double fallback) {
            if(!args.count(key)) return fallback;
            std::size_t end; const double v=std::stod(args.at(key),&end);
            if(end!=args.at(key).size() || !std::isfinite(v)) throw std::invalid_argument("invalid --"+key); return v;
        };
        auto integer=[&](const std::string &key,int fallback) {
            const double v=real(key,fallback);
            if(v!=std::floor(v) || v<0 || v>std::numeric_limits<int>::max()) throw std::invalid_argument("invalid --"+key); return int(v);
        };
        std::uint64_t seed=42;
        if(args.count("seed")) {
            std::size_t end; const auto &s=args.at("seed");
            if(s.empty() || s[0]=='-') throw std::invalid_argument("invalid seed");
            seed=std::stoull(s,&end); if(end!=s.size()) throw std::invalid_argument("invalid seed");
        }
        const SamplingSchedule schedule{integer("warmup",100),integer("sweeps",400),integer("measure-every",2),integer("block-samples",10)};
        schedule.validate();
        const int checkpointEvery=integer("checkpoint-every",0);
        SamplingProgress samplingProgress=resumed?resumed->progress:SamplingProgress{};
        checkpointRequire(samplingProgress.warmup_done>=0 && samplingProgress.warmup_done<=schedule.warmup &&
            samplingProgress.measurement_done>=0 && samplingProgress.measurement_done<=schedule.measurement_sweeps &&
            (samplingProgress.measurement_done==0 || samplingProgress.warmup_done==schedule.warmup),"saved progress outside requested schedule");
        const int validationSeries=integer("validation-series",0);
        if(validationSeries>1)throw std::invalid_argument("--validation-series must be 0 or 1");
        const bool pathSeries=args.count("path-series-bond")!=0;
        const int pathSeriesBond=integer("path-series-bond",0);
        const int bondDiagnosticsEnabled=integer("bond-diagnostics",1);
        if(bondDiagnosticsEnabled>1)throw std::invalid_argument("--bond-diagnostics must be 0 or 1");
        const int phononDiagnostics=integer("phonon-diagnostics",0);
        if(phononDiagnostics>1)throw std::invalid_argument("--phonon-diagnostics must be 0 or 1");
        const int vbsEnabled=integer("vbs",0);
        if(vbsEnabled>1)throw std::invalid_argument("--vbs must be 0 or 1");
        const int susceptibilityEnabled=integer("vbs-susceptibility",0);
        if(susceptibilityEnabled>1)throw std::invalid_argument("--vbs-susceptibility must be 0 or 1");
        if(susceptibilityEnabled && !vbsEnabled)throw std::invalid_argument("--vbs-susceptibility requires --vbs 1");
        const int amplitudeEnabled=integer("vbs-amplitude",0);
        const int amplitudeResponse=integer("vbs-amplitude-response",0);
        if(amplitudeEnabled>1 || amplitudeResponse>1)
            throw std::invalid_argument("amplitude switches must be 0 or 1");
        if(amplitudeEnabled && !vbsEnabled)
            throw std::invalid_argument("--vbs-amplitude requires --vbs 1");
        if(amplitudeResponse && !amplitudeEnabled)
            throw std::invalid_argument("--vbs-amplitude-response requires --vbs-amplitude 1");
        const int fermionsEnabled=integer("fermions",0);
        if(fermionsEnabled>1)throw std::invalid_argument("--fermions must be 0 or 1");
        const int measureEnergy=integer("measure-energy",1);
        if(measureEnergy>1)throw std::invalid_argument("--measure-energy must be 0 or 1");
        const int measureTotalEnergy=integer("measure-total-energy",0);
        if(measureTotalEnergy>1)throw std::invalid_argument("--measure-total-energy must be 0 or 1");
        const int acceptanceEvery=integer("acceptance-every",0);
        const int qDiagnosticEvery=integer("q-diagnostic-every",0);
        const int warmupObservableEvery=integer("warmup-observable-every",0);
        if(warmupObservableEvery>0 && vbsEnabled==0)
            throw std::invalid_argument("--warmup-observable-every requires --vbs 1");
        const auto proposalName=args.count("proposal")?args.at("proposal"):"local";
        if(proposalName!="local" && proposalName!="bond-refresh" && proposalName!="bond-guided")throw std::invalid_argument("unknown --proposal");
        const auto proposalMode=proposalName=="local"?PhononProposal::LocalConditional:
            (proposalName=="bond-refresh"?PhononProposal::BondRefresh:PhononProposal::GuidedBond);
        if(proposalMode!=PhononProposal::GuidedBond)
            for(const std::string name:{"bond-beta","bond-guide","bond-guide-x","bond-guide-y","bond-guide-z"})
                if(args.count(name))throw std::invalid_argument("--"+name+" requires --proposal bond-guided");
        Couplings j=p2Couplings(real("g",.3)); j.Jx=real("Jx",1.); j.Jy=real("Jy",1.); j.Jz=real("Jz",.5);
        const int lx=integer("Lx",2),ly=integer("Ly",2),shear=integer("shear",0);
        const KitaevSSHModel model(HoneycombLattice(lx,ly,shear),j);
        const auto trialSelection=selectFreeMatterTrial(model,args.count("matter-trial")?args.at("matter-trial"):"physical");
        if(pathSeries && pathSeriesBond>=model.lattice().bondCount())throw std::invalid_argument("--path-series-bond outside bond table");
        const std::string qBonds=args.count("q-mode-bonds")?args.at("q-mode-bonds"):"z";
        std::vector<QModeOptions> qPool;
        if(args.count("q-mode-pool")) {
            if(args.count("q-mode-n1")||args.count("q-mode-n2"))throw std::invalid_argument("pool and single q-mode indices are mutually exclusive");
            std::istringstream in(inputText("q-mode-pool"));
            qPool=readQModePool(in,model.lattice(),integer("q-mode-every",0),real("q-mode-beta",1.),qBonds);
        }
        const auto qOptions=qPool.empty()?makeQModeOptions(model.lattice(),integer("q-mode-n1",0),integer("q-mode-n2",0),integer("q-mode-every",0),real("q-mode-beta",1.),qBonds):qPool.front();
        const auto gammaOptions=makeGammaOptions(model.lattice(),integer("gamma-every",0),real("gamma-beta",.5));
        const int mixingDiagnostics=integer("mixing-diagnostics",0);
        if(mixingDiagnostics>1)throw std::invalid_argument("mixing-diagnostics must be 0 or 1");
        std::vector<QModeOptions> extraPool;
        if(args.count("q-mode-extra-pool")) {
            if(qPool.empty())throw std::invalid_argument("extra pool requires primary pool");
            std::istringstream extra(inputText("q-mode-extra-pool"));
            extraPool=readQModePool(extra,model.lattice(),qOptions.every,qOptions.beta,qBonds);
            std::stringstream all;all<<inputText("q-mode-pool")<<'\n'<<inputText("q-mode-extra-pool");
            readQModePool(all,model.lattice(),qOptions.every,qOptions.beta,qBonds); // cross-pool duplicate gate
        }
        GuidedBondOptions guide;
        std::string guideSource="disabled";
        std::array<double,3> guideByType{};
        if(proposalMode==PhononProposal::GuidedBond) {
            guide.beta=real("bond-beta",1.);
            guideSource=args.count("bond-guide")?args.at("bond-guide"):"bare";
            if(guideSource!="bare" && guideSource!="zero")throw std::invalid_argument("--bond-guide must be bare or zero");
            if(!std::isfinite(guide.beta)||guide.beta<=0||guide.beta>1)throw std::invalid_argument("bond-beta must be in (0,1]");
            // Approximate ONLY the proposal center. The exact MH ratio corrects
            // all deviations; the real matter trial below remains unchanged.
            if(guideSource=="bare") {
                Eigen::SelfAdjointEigenSolver<MatType> eig(model.frameworkGenerator());
                if(eig.info()!=Eigen::Success)throw std::runtime_error("bare proposal guide eigensolver failed");
                MatType occupied=MatType::Zero(model.majoranas(),model.majoranas());
                for(int n=0;n<model.majoranas();++n)if(eig.eigenvalues()[n]<-1e-10)
                    occupied.noalias()+=eig.eigenvectors().col(n)*eig.eigenvectors().col(n).adjoint();
                for(const auto &b:model.lattice().bonds())
                    guideByType[static_cast<int>(b.type)]+=2*j.g*model.gauge()[b.index]*std::imag(occupied(b.a,b.b))/model.lattice().cells();
            }
            for(int mu=0;mu<3;++mu)guideByType[mu]=real(std::string("bond-guide-")+"xyz"[mu],guideByType[mu]);
            guide.forces.resize(model.lattice().bondCount());
            for(const auto &b:model.lattice().bonds())guide.forces[b.index]=guideByType[static_cast<int>(b.type)];
        }
        if(qDiagnosticEvery && qOptions.basis.empty()) throw std::invalid_argument("q diagnostic requires nonzero momentum");
        if(qDiagnosticEvery && qDiagnosticEvery%schedule.stride) throw std::invalid_argument("q diagnostic interval must be a multiple of measure-every");
        const bool qOutputEnabled=qOptions.every>0 || qDiagnosticEvery>0;
        std::vector<std::pair<int,Complex>> diagnosticFourier;
        if(qDiagnosticEvery>0) {
            const double u=double(qOptions.n1)/lx,v=double(qOptions.n2)/ly-double(shear)*qOptions.n1/(double(lx)*ly);
            for(const auto &b:model.lattice().bonds()) if(b.type==BondType::Z)
                diagnosticFourier.emplace_back(b.index,std::exp(Complex(0,-2*std::acos(-1.)*(u*b.cell.x+v*b.cell.y)))/double(model.lattice().cells()));
        }
        const PhononParameters phonon(real("mass",4.),real("omega-ph",.5));
        const bool fieldEnabled=args.count("vbs-field-h")!=0;
        if(!fieldEnabled && (args.count("vbs-field-n1") || args.count("vbs-field-n2") || args.count("vbs-field-phase")))
            throw std::invalid_argument("VBS source momentum/phase requires --vbs-field-h");
        if(fieldEnabled && (!susceptibilityEnabled || !args.count("vbs-field-n1") || !args.count("vbs-field-n2")))
            throw std::invalid_argument("--vbs-field-h requires susceptibility and explicit allowed momentum n1,n2");
        const VBSField field(model.lattice(),integer("vbs-field-n1",0),integer("vbs-field-n2",0),real("vbs-field-phase",0),real("vbs-field-h",0));
        const double lambda=phonon.lambda(j.g);
        std::mt19937_64 rng(seed);
        const double shifts[3]={real("trial-dx",0),real("trial-dy",0),real("trial-dz",0)};
        std::vector<double> bondShifts(model.lattice().bondCount());
        for(const auto &b:model.lattice().bonds()) bondShifts[b.index]=shifts[static_cast<int>(b.type)];
        const std::string numericsRequested=args.count("numerics")?args.at("numerics"):"production";
        if(numericsRequested!="production" && numericsRequested!="reference")
            throw std::invalid_argument("--numerics must be production or reference");
        const auto &policy=numericsRequested=="production"?productionNumerics:referenceNumerics;
        bool numericalOverrides=false;
        for(const std::string name:{"weight-engine","defer-phase","prefix-cache","cayley-cache","center-green-fast"})
            numericalOverrides|=args.count(name)!=0;
        const std::string numericsProfile=numericalOverrides?"compatibility-overrides":policy.name;
        const std::string engine=args.count("weight-engine")?args.at("weight-engine"):policy.weight_engine;
        if(engine!="reference" && engine!="svd-values" && engine!="node-merge" && engine!="combined")
            throw std::invalid_argument("unknown --weight-engine");
        const int deferPhase=integer("defer-phase",policy.defer_phase);
        if(deferPhase>1) throw std::invalid_argument("--defer-phase must be 0 or 1");
        const int prefixCache=integer("prefix-cache",policy.prefix_cache);
        const int cayleyCache=integer("cayley-cache",policy.cayley_cache);
        const int centerGreenFast=integer("center-green-fast",policy.center_green_fast);
        if(centerGreenFast>1)throw std::invalid_argument("--center-green-fast must be 0 or 1");
        if(cayleyCache>1)throw std::invalid_argument("--cayley-cache must be 0 or 1");
        if(prefixCache>1) throw std::invalid_argument("--prefix-cache must be 0 or 1");
        if(deferPhase)
            std::cerr<<"[deferred phase-on-demand] magnitude-only updates, full Pfaffian phase/Green at measurement\n";
        const PhononWeightOptions weightOptions{
            engine=="svd-values" || engine=="combined",
            engine=="node-merge" || engine=="combined",
            deferPhase==1,
            prefixCache==1,cayleyCache==1,
            fermionsEnabled!=0 || measureEnergy!=0 || bondDiagnosticsEnabled!=0,centerGreenFast==1};
        const int qrBlock=integer("qr-block",4);
        if(outputPathExists(prefix+"_checkpoint.control.json") || outputPathExists(prefix+"_checkpoint.receipt.json") ||
           outputPathExists(checkpointFile))throw std::invalid_argument("checkpoint output exists; choose fresh output and checkpoint paths");
        if(outputPathExists(prefix+"_blocks.csv") || outputPathExists(prefix+"_results.json") ||
           outputPathExists(prefix+"_acceptance.csv") ||
           (qOutputEnabled && outputPathExists(prefix+"_q_acceptance.csv")) ||
           (qDiagnosticEvery>0 && outputPathExists(prefix+"_q_diagnostic.csv")) || outputPathExists(prefix+"_numerical.jsonl") || outputPathExists(prefix+"_results.json.tmp") || outputPathExists(prefix+"_matter_trial.json"))
            throw std::invalid_argument("output exists; choose another --output prefix");
        {
            std::ofstream trialRecord(prefix+"_matter_trial.json");
            trialRecord<<std::setprecision(17);writeMatterTrialSelection(trialRecord,trialSelection);trialRecord<<'\n';
            if(!trialRecord)throw std::runtime_error("cannot write matter trial identity");
        }
        numerics.status(prefix,"running");numerics.owns_output=true;
        std::cerr<<"Numerical failure policy: "<<numericalFailurePolicy<<'\n';
        // Initial-state search is separate from an MC transition. Retry the
        // same initial draw first; only initialization may draw another path.
        std::unique_ptr<CoupledPhononPath> ownedState;
        CheckpointSignals checkpointSignals;
        long long untrustedInitial=0; const int maxInitialAttempts=100000;
        for(int initAttempt=0; initAttempt<=maxInitialAttempts; ++initAttempt) {
            auto initial=PhononPath(model.lattice().bondCount(),integer("intervals",4),real("dt",.2),phonon,bondShifts,fieldEnabled?field.forces():std::vector<double>{});
            if(resumed) {
                std::istringstream in(resumed->path);auto loaded=PhononPath::load(in);checkpointEnd(in);
                if(loaded.bonds()!=initial.bonds() || loaded.intervals()!=initial.intervals() || loaded.dt()!=initial.dt() ||
                   loaded.parameters().mass!=phonon.mass || loaded.parameters().omega!=phonon.omega ||
                   loaded.trialShifts()!=initial.trialShifts() || loaded.forces()!=initial.forces())
                    throw std::invalid_argument("checkpoint path physical/trial metadata mismatch");
                initial=std::move(loaded);
            } else initial=initial.freeSample(rng);
            if(args.count("initial-path") && !resumed) {
                std::ifstream f(args.at("initial-path"));if(!f)throw std::invalid_argument("cannot open initial path");
                const auto loaded=PhononPath::load(f);f>>std::ws;if(!f.eof())throw std::invalid_argument("trailing initial-path data");
                if(loaded.bonds()!=initial.bonds() || loaded.intervals()!=initial.intervals() || loaded.dt()!=initial.dt() ||
                   loaded.parameters().mass!=phonon.mass || loaded.parameters().omega!=phonon.omega ||
                   loaded.trialShifts()!=initial.trialShifts() || loaded.forces()!=initial.forces())
                    throw std::invalid_argument("initial path physical/trial metadata mismatch");
                initial=loaded;
            }
            NumericalEvent event{"initialization","initial-state","","",false,std::uint64_t(initAttempt+1)};
            try {
                ownedState=std::make_unique<CoupledPhononPath>(
                    initial,trialSelection.trial,honeycombGenerator(model),qrBlock,weightOptions);
            } catch(const std::runtime_error &e) {event.first_error=e.what();}
            if(ownedState)break;
            try {
                ownedState=std::make_unique<CoupledPhononPath>(
                    initial,trialSelection.trial,honeycombGenerator(model),qrBlock,weightOptions,true);
            } catch(const std::exception &e) {event.reference_error=e.what();}
            event.recovered=bool(ownedState);numerics.record(prefix,event,initial);
            if(ownedState)break;
            if(args.count("initial-path") || resumed)throw std::runtime_error("supplied initial/checkpoint path failed certification; no redraw permitted");
            ++untrustedInitial;
            if(initAttempt==maxInitialAttempts)
                throw std::runtime_error("initial trusted state unavailable: "+event.reference_error);
        }
        CoupledPhononPath &state=*ownedState;
        if(args.count("initial-path") || resumed){std::ofstream f(prefix+"_initial.path");state.state().path.save(f);}
        // VBS measurement selection. --vbs is the master switch; measurement mode
        // never feeds back into the sampler (no RNG draws, no weight changes).
        const auto vbsModeNameStr=args.count("vbs-mode")?args.at("vbs-mode"):"all";
        const VBSMode vbsMode=vbsEnabled?parseVBSMode(vbsModeNameStr):VBSMode::All;
        if(!vbsEnabled && (args.count("vbs-mode") || args.count("vbs-q-file")))
            throw std::invalid_argument("--vbs-mode/--vbs-q-file require --vbs 1");
        if(vbsEnabled && vbsMode==VBSMode::All && args.count("vbs-q-file"))
            throw std::invalid_argument("--vbs-q-file is only used with --vbs-mode selected|both");
        std::vector<VBSQPoint> vbsQ;
        if(vbsEnabled && (vbsMode==VBSMode::Selected || vbsMode==VBSMode::Both)) {
            if(!args.count("vbs-q-file"))
                throw std::invalid_argument("--vbs-mode selected|both requires --vbs-q-file");
            std::istringstream in(inputText("vbs-q-file"));vbsQ=parseVBSQPoints(in,model.lattice());
        }
        if(susceptibilityEnabled && vbsMode!=VBSMode::All)
            throw std::invalid_argument("--vbs-susceptibility requires --vbs-mode all");
        if(amplitudeEnabled && vbsMode!=VBSMode::All)
            throw std::invalid_argument("--vbs-amplitude requires --vbs-mode all");
        const VBSObservables vbs(model.lattice(),vbsMode,vbsQ);
        const FermionObservables fermions(model);
        auto names=basicObservableNames(measureEnergy!=0);
        if(fermionsEnabled) {const auto extra=fermions.names();names.insert(names.end(),extra.begin(),extra.end());}
        const auto vbsOffset=names.size();
        if(vbsEnabled) { const auto extra=vbs.names(); names.insert(names.end(),extra.begin(),extra.end()); }
        const auto susceptibilityOffset=names.size();
        if(susceptibilityEnabled) {const auto extra=VBSSusceptibility(vbs).names();names.insert(names.end(),extra.begin(),extra.end());}
        const auto fieldOffset=names.size();
        if(fieldEnabled) {const auto extra=field.names();names.insert(names.end(),extra.begin(),extra.end());}
        const auto amplitudeOffset=names.size();
        if(amplitudeEnabled) {const auto extra=VBSAmplitude(vbs,amplitudeResponse).names();names.insert(names.end(),extra.begin(),extra.end());}
        if(phononDiagnostics) {
            for(int j=0;j<state.state().path.nodes();++j)
                for(const auto *stem:{"PhDiagX", "PhDiagX2", "PhDiagCenterC", "PhDiagOtherBondC"})
                    names.push_back(std::string(stem)+"_j"+std::to_string(j));
            names.push_back("PhDiagLogMatter");
        }
        if(bondDiagnosticsEnabled) {const auto extra=bondDiagnosticNames(model);names.insert(names.end(),extra.begin(),extra.end());}
        const auto totalEnergyOffset=names.size();
        if(measureTotalEnergy) {const auto extra=totalEnergyNames();names.insert(names.end(),extra.begin(),extra.end());}
        // Optional raw per-measurement stream. No RNG, state, or estimator changes.
        // Values are UNWEIGHTED; the adjacent common complex phase must be applied.
        std::ofstream validationOutput;
        std::vector<std::size_t> validationColumns;
        if(validationSeries || pathSeries) {
            validationOutput.open(prefix+"_validation_series.csv");
            if(!validationOutput)throw std::runtime_error("cannot open validation series");
            validationOutput<<"sweep,phase_real,phase_imag"<<std::setprecision(17);
            const std::vector<std::string> selected={"X_x","X_y","X_z","X2_x","X2_y","X2_z","majorana_energy",
                "BondK_b2","BondXK_b2","BondK_x","BondK_y","BondK_z","BondXK_x","BondXK_y","BondXK_z",
                "S_q0_zz","S_q3_zz","AmpA2_q0_z","AmpA4_q0_z","AmpA2_q3_z","AmpA4_q3_z"};
            if(validationSeries)for(const auto &name:selected) {
                const auto it=std::find(names.begin(),names.end(),name);
                if(it!=names.end()) {validationColumns.push_back(std::distance(names.begin(),it));validationOutput<<','<<name<<"_real,"<<name<<"_imag";}
            }
            if(pathSeries)for(int j=0;j<=integer("intervals",4);++j)validationOutput<<",PathX_j"<<j;
            validationOutput<<'\n';
        }
        PhononSampler sampler(state,rng,proposalMode,
            [&](const NumericalEvent &e,const PhononPath &p){numerics.record(prefix,e,p);},qOptions,guide,qPool,gammaOptions,extraPool);
        PhaseStatistics statistics(names,schedule.block_samples);
        if(resumed) {
            std::istringstream samplerIn(resumed->sampler),statsIn(resumed->statistics);
            sampler.loadCheckpoint(samplerIn);statistics.loadCheckpoint(statsIn);
            checkpointRequire((samplingProgress.measurement_done==0 || sampler.stage()=="measurement") &&
                (samplingProgress.warmup_done==schedule.warmup || sampler.stage()=="warmup"),"sampler stage/progress mismatch");
            checkpointRequire(statistics.samplesIncludingPending()==std::uint64_t(samplingProgress.measurement_done/schedule.stride),"progress/statistics count mismatch");
        }
        double blockFileOutputSeconds=0.,totalEnergySeconds=0.;std::uint64_t totalEnergyMeasurements=0;
        const auto blockHeaderBegin=std::chrono::steady_clock::now();
        std::ofstream blocks(prefix+"_blocks.csv");
        if(!blocks) throw std::runtime_error("cannot open block output");
        blocks<<std::setprecision(17); statistics.writeHeader(blocks);
        if(resumed)statistics.writeAllBlocks(blocks);
        blocks.flush();
        if(!blocks) throw std::runtime_error("block header write failed");
        std::ofstream acceptanceOutput;
        if(acceptanceEvery>0) {
            acceptanceOutput.open(prefix+"_acceptance.csv");
            if(!acceptanceOutput) throw std::runtime_error("cannot open acceptance output");
            acceptanceOutput<<"stage,sweep,proposed,accepted,metropolis_rejected,untrusted_rejected,first_failures,reference_recovered,reference_failed,acceptance,elapsed_seconds\n"<<std::setprecision(17);
            acceptanceOutput.flush();
            if(!acceptanceOutput) throw std::runtime_error("acceptance header write failed");
        }
        std::ofstream qAcceptanceOutput,qDiagnosticOutput;
        if(qOutputEnabled) {
            qAcceptanceOutput.open(prefix+"_q_acceptance.csv");
            if(!qAcceptanceOutput) throw std::runtime_error("cannot open q acceptance output");
            qAcceptanceOutput<<"stage,sweep,base_proposed,base_accepted,base_metropolis_rejected,base_untrusted_rejected,base_first_failures,base_reference_recovered,base_reference_failed,q_proposed,q_accepted,q_metropolis_rejected,q_untrusted_rejected,q_first_failures,q_reference_recovered,q_reference_failed\n"<<std::setprecision(17);
        }
        if(qDiagnosticEvery>0) {
            qDiagnosticOutput.open(prefix+"_q_diagnostic.csv");
            if(!qDiagnosticOutput) throw std::runtime_error("cannot open q diagnostic output");
            qDiagnosticOutput<<"stage,sweep,cpu_seconds,elapsed_seconds,phase_real,phase_imag,m2_q,D_real,D_imag";
            if(measureEnergy) qDiagnosticOutput<<",majorana_energy_real,majorana_energy_imag";
            qDiagnosticOutput<<",X_z_real,X_z_imag,proposed,accepted,q_proposed,q_accepted\n"<<std::setprecision(17);
        }
        std::ofstream mixingOutput,modeOutput;
        if(mixingDiagnostics) {
            mixingOutput.open(prefix+"_mixing.csv");modeOutput.open(prefix+"_mode_counts.csv");
            if(!mixingOutput||!modeOutput)throw std::runtime_error("cannot open mixing diagnostics");
            mixingOutput<<std::setprecision(17)<<"sweep,phase_real,phase_imag,cpu_seconds,wall_seconds,X_x,X_y,X_z";
            for(int n1=0;n1<lx;++n1)for(int n2=0;n2<ly;++n2)for(char mu:std::string("xyz"))
                mixingOutput<<",D_"<<n1<<'_'<<n2<<'_'<<mu<<"_re,D_"<<n1<<'_'<<n2<<'_'<<mu<<"_im,M2_"<<n1<<'_'<<n2<<'_'<<mu;
            mixingOutput<<'\n';
            modeOutput<<std::setprecision(17)<<"stage,sweep,kernel,n1,n2,bonds,rank,proposed,accepted,metropolis_rejected,first_failures,reference_recovered,reference_failed,seconds\n";
        }
        MoveCounts qWarm,totalWarm,gammaWarm;
        int observationIndex=samplingProgress.measurement_done/schedule.stride;
        std::ofstream warmupObservables;
        if(warmupObservableEvery>0) {
            warmupObservables.open(prefix+"_warmup.csv");
            if(!warmupObservables) throw std::runtime_error("cannot open warmup observable output");
            warmupObservables<<"sweep,raw_all_q_max,q_index\n"<<std::setprecision(17);
        }
        blockFileOutputSeconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-blockHeaderBegin).count();
        const auto samplingBegin=std::chrono::steady_clock::now();
        SamplingTimings samplingTimings;
        std::uint64_t checkpointGeneration=resumed?resumed->generation:0;
        if(resumed) {
            std::istringstream in(resumed->runner);
            qWarm=PhononSampler::readCounts(in);totalWarm=PhononSampler::readCounts(in);gammaWarm=PhononSampler::readCounts(in);
            std::size_t n;checkpointRequire(bool(in>>untrustedInitial>>n) && untrustedInitial>=0 && n<=3,"invalid runner counters");
            for(std::size_t i=0;i<n;++i) {
                std::string stage;checkpointRequire(bool(in>>std::quoted(stage)) && (stage=="initialization"||stage=="warmup"||stage=="measurement"),"invalid numerical stage");
                const auto inherited=PhononSampler::readCounts(in);auto &current=numerics.stages[stage];
                // Do not hide a numerical recovery during checkpoint-path reconstruction.
                current.first_failures+=inherited.first_failures;current.reference_recovered+=inherited.reference_recovered;current.reference_failed+=inherited.reference_failed;
            }
            checkpointEnd(in);
        }
        const auto executableIdentity=checkpointExecutableIdentity();
        char checkpointHost[256]{};
        if(gethostname(checkpointHost,sizeof(checkpointHost)-1)!=0)throw std::runtime_error("checkpoint hostname lookup failed");
        const auto processStat=checkpointReadFile("/proc/self/stat");
        std::istringstream statFields(processStat.substr(processStat.rfind(')')+2));
        std::string startTicks,statWord;for(int field=3;field<=22;++field){statFields>>statWord;if(field==22)startTicks=statWord;}
        auto receipt=[&](bool stopped) {
            std::ostringstream out;out<<"{\"pid\":"<<getpid()<<",\"hostname\":";jsonString(out,checkpointHost);
            out<<",\"process_start_ticks\":";jsonString(out,startTicks);
            out<<",\"pbs_jobid\":";jsonString(out,std::getenv("PBS_JOBID")?std::getenv("PBS_JOBID"):"");
            out<<",\"checkpoint_file\":";jsonString(out,checkpointFile);
            out<<",\"generation\":"<<checkpointGeneration<<",\"warmup_done\":"<<samplingProgress.warmup_done
               <<",\"measurement_done\":"<<samplingProgress.measurement_done<<",\"stopped\":"<<(stopped?"true":"false")
               <<",\"exit_code_if_stopped\":75}\n";return out.str();
        };
        checkpointAtomicWrite(prefix+"_checkpoint.control.json",receipt(false));
        auto saveCheckpoint=[&](bool stop) {
            for(auto *stream:{&blocks,&validationOutput,&acceptanceOutput,&qAcceptanceOutput,&qDiagnosticOutput,&mixingOutput,&modeOutput,&warmupObservables})
                if(stream->is_open()){stream->flush();if(!*stream)throw std::runtime_error("flush before checkpoint failed");}
            MCCheckpoint c;c.executable=executableIdentity;c.origin=prefix;c.args=args;c.inputs=checkpointInputs;c.progress=samplingProgress;c.generation=++checkpointGeneration;
            std::ostringstream pathOut,samplerOut,statsOut,runner;
            state.state().path.save(pathOut);sampler.saveCheckpoint(samplerOut);statistics.saveCheckpoint(statsOut);
            PhononSampler::writeCounts(runner,qWarm);PhononSampler::writeCounts(runner,totalWarm);PhononSampler::writeCounts(runner,gammaWarm);
            runner<<untrustedInitial<<' '<<numerics.stages.size()<<'\n';
            for(const auto &p:numerics.stages){runner<<std::quoted(p.first)<<' ';PhononSampler::writeCounts(runner,p.second);}
            c.path=pathOut.str();c.sampler=samplerOut.str();c.statistics=statsOut.str();c.runner=runner.str();c.save(checkpointFile);
            checkpointAtomicWrite(prefix+"_checkpoint.receipt.json",receipt(stop));
            std::cerr<<"CHECKPOINT_SAVED generation="<<checkpointGeneration<<" warmup="<<samplingProgress.warmup_done<<" measurement="<<samplingProgress.measurement_done<<" file="<<checkpointFile<<" stop="<<stop<<'\n';
        };
        const auto warm=sample(sampler,schedule,statistics,[&](const auto &s){
            auto values=basicObservables(model,s,measureEnergy!=0);
            const int sweep=++observationIndex*schedule.stride;
            if(qDiagnosticEvery>0 && sweep%qDiagnosticEvery==0) {
                Complex d=0.;
                for(const auto &term:diagnosticFourier)
                    d+=s.path.coordinate(s.path.center(),term.first)*term.second;
                const auto c=sampler.counts()-totalWarm,q=sampler.qCounts()-qWarm;
                qDiagnosticOutput<<"measurement,"<<sweep<<','<<double(std::clock()-cpuStart)/CLOCKS_PER_SEC<<','
                    <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<','
                    <<s.weight.matter_phase.real()<<','<<s.weight.matter_phase.imag()<<','<<std::norm(d)<<','<<d.real()<<','<<d.imag();
                if(measureEnergy) qDiagnosticOutput<<','<<values[6].real()<<','<<values[6].imag();
                qDiagnosticOutput<<','<<values[2].real()<<','<<values[2].imag()<<','<<c.proposed<<','<<c.accepted<<','<<q.proposed<<','<<q.accepted<<'\n';
                qDiagnosticOutput.flush();
                if(!qDiagnosticOutput) throw std::runtime_error("q diagnostic output write failed");
            }
            if(mixingDiagnostics) {
                mixingOutput<<sweep<<','<<s.weight.matter_phase.real()<<','<<s.weight.matter_phase.imag()<<','
                    <<double(std::clock()-cpuStart)/CLOCKS_PER_SEC<<','<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                for(int mu=0;mu<3;++mu)mixingOutput<<','<<values[mu].real();
                for(int n1=0;n1<lx;++n1)for(int n2=0;n2<ly;++n2) {
                    std::array<Complex,3> d{};
                    const double u=double(n1)/lx,v=double(n2)/ly-double(shear)*n1/(double(lx)*ly);
                    for(const auto &b:model.lattice().bonds())d[int(b.type)]+=s.path.coordinate(s.path.center(),b.index)*std::exp(Complex(0,-2*pi*(u*b.cell.x+v*b.cell.y)))/double(model.lattice().cells());
                    for(const auto &z:d)mixingOutput<<','<<z.real()<<','<<z.imag()<<','<<std::norm(z);
                }
                mixingOutput<<'\n';mixingOutput.flush();if(!mixingOutput)throw std::runtime_error("mixing stream write failed");
            }
            if(fermionsEnabled) {const auto extra=fermions.measure(s.weight.center_green);values.insert(values.end(),extra.begin(),extra.end());}
            if(vbsEnabled) {
                const auto extra=vbs.measure(s.path.node(s.path.center()));
                values.insert(values.end(),extra.begin(),extra.end());
            }
            if(susceptibilityEnabled) {const auto extra=VBSSusceptibility(vbs).measure(s.path);values.insert(values.end(),extra.begin(),extra.end());}
            if(fieldEnabled) {const auto extra=field.measure(s.path);values.insert(values.end(),extra.begin(),extra.end());}
            if(amplitudeEnabled) {const auto extra=VBSAmplitude(vbs,amplitudeResponse).measure(s.path);values.insert(values.end(),extra.begin(),extra.end());}
            if(phononDiagnostics) {
                const auto &p=s.path;
                for(int j=0;j<p.nodes();++j) {
                    double x=0,x2=0,c=0,other=0;
                    for(int b=0;b<p.bonds();++b) {
                        const double y=p.coordinate(j,b)/p.bonds();
                        x+=y;x2+=y*p.coordinate(j,b);
                        c+=y*p.coordinate(p.center(),b);
                        other+=y*p.coordinate(p.center(),(b+1)%p.bonds());
                    }
                    values.insert(values.end(),{x,x2,c,other});
                }
                values.push_back(s.weight.log_abs_matter);
            }
            if(bondDiagnosticsEnabled) {const auto extra=bondDiagnostics(model,s);values.insert(values.end(),extra.begin(),extra.end());}
            if(measureTotalEnergy) {
                const auto energyBegin=std::chrono::steady_clock::now();
                const auto extra=sampler.totalEnergy(model);values.insert(values.end(),extra.begin(),extra.end());
                totalEnergySeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-energyBegin).count();
                ++totalEnergyMeasurements;
            }
            if(validationOutput.is_open()) {
                validationOutput<<sweep<<','<<s.weight.matter_phase.real()<<','<<s.weight.matter_phase.imag();
                for(const auto index:validationColumns)validationOutput<<','<<values.at(index).real()<<','<<values.at(index).imag();
                if(pathSeries)for(int j=0;j<s.path.nodes();++j)validationOutput<<','<<s.path.coordinate(j,pathSeriesBond);
                validationOutput<<'\n';
                if(!validationOutput)throw std::runtime_error("validation series write failed");
            }
            return values;
        },
            [&](const auto &stats){
                const auto blockOutputBegin=std::chrono::steady_clock::now();
                stats.writeLastBlock(blocks);
                blockFileOutputSeconds+=std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-blockOutputBegin).count();
            },&samplingTimings,{},[&](const char *stage,int sweep,const MoveCounts &counts) {
                const bool inWarmup=std::string(stage)=="warmup";
                const int endSweep=inWarmup?schedule.warmup:schedule.measurement_sweeps;
                if(mixingDiagnostics && !inWarmup && (sweep%100==0 || sweep==endSweep)) {
                    auto write=[&](const char *kernel,const QModeOptions &q,const MoveCounts &c,double seconds) {
                        modeOutput<<"measurement,"<<sweep<<','<<kernel<<','<<q.n1<<','<<q.n2<<','<<q.bond_types<<','<<q.basis.size()<<','<<c.proposed<<','<<c.accepted<<','<<c.metropolis_rejected<<','<<c.first_failures<<','<<c.reference_recovered<<','<<c.reference_failed<<','<<seconds<<'\n';
                    };
                    for(std::size_t k=0;k<qPool.size()+extraPool.size();++k)write(k<qPool.size()?"primary":"extra",k<qPool.size()?qPool[k]:extraPool[k-qPool.size()],sampler.qMeasuredCounts()[k],sampler.qMeasuredSeconds()[k]);
                    if(gammaOptions.every)for(int mu=0;mu<3;++mu)write("gamma",gammaOptions.modes[mu],sampler.gammaMeasuredCounts()[mu],sampler.gammaSeconds()[mu]);
                    modeOutput.flush();if(!modeOutput)throw std::runtime_error("mode counts write failed");
                }
                if(inWarmup) {
                    qWarm=sampler.qCounts();totalWarm=sampler.counts();gammaWarm=sampler.gammaCounts();
                    // This is a raw path diagnostic, with no phase reweighting
                    // or Green observable. It must not materialize either.
                    if(warmupObservableEvery>0 && sweep%warmupObservableEvery==0) {
                        const auto &p=sampler.state().path;
                        const auto value=vbs.rawAllQMax(p.node(p.center()));
                        warmupObservables<<sweep<<','<<value.first<<','<<value.second<<'\n';
                        if(!warmupObservables) throw std::runtime_error("warmup observable output write failed");
                    }
                }
                const int qEvery=acceptanceEvery>0?acceptanceEvery:(qDiagnosticEvery>0?qDiagnosticEvery:qOptions.every);
                if(qOutputEnabled && (sweep%qEvery==0 || sweep==endSweep)) {
                    const auto q=inWarmup?sampler.qCounts():sampler.qCounts()-qWarm;
                    const auto gamma=inWarmup?sampler.gammaCounts():sampler.gammaCounts()-gammaWarm;
                    const auto base=counts-q-gamma;
                    qAcceptanceOutput<<stage<<','<<sweep;
                    for(const auto &c:{base,q}) qAcceptanceOutput<<','<<c.proposed<<','<<c.accepted<<','<<c.metropolis_rejected<<','<<c.untrusted_rejected<<','<<c.first_failures<<','<<c.reference_recovered<<','<<c.reference_failed;
                    qAcceptanceOutput<<'\n';qAcceptanceOutput.flush();
                    if(!qAcceptanceOutput) throw std::runtime_error("q acceptance output write failed");
                }
                if(acceptanceEvery>0 && (sweep%acceptanceEvery==0 || sweep==endSweep)) {
                    acceptanceOutput<<stage<<','<<sweep<<','<<counts.proposed<<','<<counts.accepted
                        <<','<<counts.metropolis_rejected<<','<<counts.untrusted_rejected
                        <<','<<counts.first_failures<<','<<counts.reference_recovered<<','<<counts.reference_failed
                        <<','<<counts.acceptance()<<','
                        <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
                    acceptanceOutput.flush();
                    if(!acceptanceOutput) throw std::runtime_error("acceptance output write failed");
                }
                const auto request=checkpointSignals.consume();
                if(request.first || (checkpointEvery>0 && sweep%checkpointEvery==0))saveCheckpoint(request.second);
                if(request.second)throw CheckpointStopped{};
            },&samplingProgress);
        if(warmupObservables) warmupObservables.close();
        const double initializationSeconds=std::chrono::duration<double>(samplingBegin-start).count();
        const double samplingWallSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-samplingBegin).count();
        const auto statisticsFinalizeBegin=std::chrono::steady_clock::now();
        const auto summary=statistics.summary(), doubled=statistics.doubledBlockSummary();
        const double statisticsFinalizeSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-statisticsFinalizeBegin).count();
        const auto all=sampler.counts(); const MoveCounts production=all-warm;
        // Keep elapsed_seconds at its Stage 1 boundary: sampling plus the
        // primary and doubled summaries, before VBS/JSON/output work.
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        const auto blockCloseBegin=std::chrono::steady_clock::now();
        blocks.flush();
        blocks.close();
        if(!blocks) throw std::runtime_error("block output close failed");
        blockFileOutputSeconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-blockCloseBegin).count();
        const auto &finalPath=state.state().path;
        double vbsBlockingSeconds=0.,candidateDiagnosticsSeconds=0.;
        const auto jsonBuildBegin=std::chrono::steady_clock::now();
        std::ostringstream json; json<<std::setprecision(17)<<std::boolalpha;
        json<<"{\"status\":\"complete\",\"scope\":\"fixed_gauge_matter_sector\",\"coupled_mc\":true,\"seed\":"<<seed
            <<",\"Lx\":"<<model.lattice().Lx()<<",\"Ly\":"<<model.lattice().Ly()<<",\"Jx\":"<<j.Jx<<",\"Jy\":"<<j.Jy<<",\"Jz\":"<<j.Jz
            <<",\"g\":"<<j.g<<",\"mass\":"<<phonon.mass<<",\"omega_ph\":"<<phonon.omega<<",\"k\":"<<phonon.spring()<<",\"lambda\":"<<lambda
            <<",\"phonon_trial\":\"equal_shifted_vacuum_endpoints\",\"trial_shift\":["<<shifts[0]<<','<<shifts[1]<<','<<shifts[2]<<']'
            <<",\"matter_trial\":\""<<(trialSelection.policy=="physical"?"physical_parity_free_state_static_X_zero":"unrestricted_free_ground_state_static_X_zero")<<"\",\"matter_trial_sector\":";
        writeMatterTrialSelection(json,trialSelection);
        if(resumed) {
            json<<",\"resume\":{\"checkpoint\":";jsonString(json,args.at("resume"));
            json<<",\"ancestor_output_prefix\":";jsonString(json,resumed->origin);
            json<<",\"inherited_samples\":"<<resumed->progress.measurement_done/schedule.stride
                <<",\"raw_blocks_include_ancestor\":true,\"diagnostic_streams_only_this_segment\":true,\"bitwise_interacting_replay_guaranteed\":false,\"effective_options\":{";
            bool comma=false;for(const auto &p:args){if(comma)json<<',';comma=true;jsonString(json,p.first);json<<':';jsonString(json,p.second);}json<<"}}";
        }
        json<<",\"spatial_boundary\":\"pbc\",\"pbc_geometry\":\""<<(shear==0?"axis_aligned":"sheared")<<"\""
            <<",\"shear\":"<<shear<<",\"geometry_id\":\""<<model.lattice().geometryId()<<"\""
            <<",\"torus_matrix\":[["<<lx<<','<<shear<<"],[0,"<<ly<<"]]"
            <<",\"period_vectors_cartesian\":[["<<model.lattice().periodVector(0).x<<','<<model.lattice().periodVector(0).y
            <<"],["<<model.lattice().periodVector(1).x<<','<<model.lattice().periodVector(1).y<<"]]"
            <<",\"momentum_grid_convention\":\"n1/Lx, n2/Ly-shear*n1/(Lx*Ly), index=n1*Ly+n2\""
            <<",\"wilson_loops\":[1,1],\"time_contour\":\"open\",\"intervals\":"<<state.state().path.intervals()
            <<",\"dt\":"<<state.state().path.dt()<<",\"tau\":"<<state.state().path.totalTau()<<",\"acceptance_every\":"<<acceptanceEvery<<",\"warmup_sweeps\":"<<schedule.warmup
            <<",\"measurement_sweeps\":"<<schedule.measurement_sweeps<<",\"measure_every\":"<<schedule.stride<<",\"block_samples\":"<<schedule.block_samples
            <<",\"bond_diagnostics\":"<<bondDiagnosticsEnabled
            <<",\"phonon_diagnostics\":"<<phononDiagnostics
            <<",\"qr_block\":"<<integer("qr-block",4)
            <<",\"numerics_requested\":\""<<numericsRequested<<"\""
            <<",\"numerics_profile\":\""<<numericsProfile<<"\""
            <<",\"weight_engine\":\""<<engine<<"\""
            <<",\"defer_phase\":"<<deferPhase
            <<",\"center_green_fast\":"<<centerGreenFast
            <<",\"center_green_qr_fast_calls\":"<<pure_projector_endpoint_detail::endpointAcceleratedQrCalls
            <<",\"center_green_qr_reference_calls\":"<<pure_projector_endpoint_detail::endpointReferenceQrCalls;
        json<<",\"gauge_u\":[";
        for(std::size_t b=0;b<model.gauge().size();++b){if(b)json<<',';json<<model.gauge()[b];}json<<']';
        writeRoundoffDiagnostics(json);
        json<<",\"cayley_cache\":"<<cayleyCache
            <<",\"cayley_lu_calls\":"<<cayleyLUSolveCount
            <<",\"cayley_cache_hits\":"<<cayleyCacheHitCount
            <<",\"cayley_cache_fallbacks\":"<<cayleyCacheFallbackCount
            <<",\"cayley_bytes_estimate\":"<<(state.state().lu_cayley_cache?pureProjectorCayleyBytes(*state.state().lu_cayley_cache):0)
            <<",\"prefix_resumes\":"<<signedPrefixResumeCount
            <<",\"signed_factor_steps\":"<<signedPropagationSliceCount
            <<",\"prefix_cache\":"<<prefixCache
            <<",\"validation_series\":"<<validationSeries
            <<",\"path_series_bond\":"<<(pathSeries?pathSeriesBond:-1)
            <<",\"proposal\":\""<<(proposalMode==PhononProposal::LocalConditional?"free_local_conditional_gaussian":(proposalMode==PhononProposal::BondRefresh?"free_whole_bond_refresh":"fixed_guided_bond_pcn"))
            <<"\",\"sweep\":\""<<(proposalMode==PhononProposal::LocalConditional?"random_permutation_all_coordinates":"random_permutation_all_bond_paths")
            <<"\",\"proposal_tuning\":false"
            <<",\"metropolis_rejected_warmup\":"<<warm.metropolis_rejected<<",\"metropolis_rejected_measurement\":"<<production.metropolis_rejected
            <<",\"proposals_warmup\":"<<warm.proposed<<",\"accepted_warmup\":"<<warm.accepted
            <<",\"proposals_measurement\":"<<production.proposed<<",\"accepted_measurement\":"<<production.accepted
            <<",\"untrusted_rejected_warmup\":"<<warm.untrusted_rejected<<",\"untrusted_rejected_measurement\":"<<production.untrusted_rejected<<",\"untrusted_rejected_initial\":"<<untrustedInitial<<",\"acceptance\":"<<production.acceptance()<<",\"samples\":"<<summary.samples<<",\"blocks\":"<<summary.blocks
            <<",\"final_center_node\":"<<finalPath.center()<<",\"final_phonon_path\":[";
        for(int node=0;node<finalPath.nodes();++node) {
            if(node) json<<',';
            json<<'[';
            for(int bond=0;bond<finalPath.bonds();++bond) {
                if(bond) json<<',';
                json<<finalPath.coordinate(node,bond);
            }
            json<<']';
        }
        json<<"],\"final_center_coordinates\":[";
        for(int bond=0;bond<finalPath.bonds();++bond) {
            if(bond) json<<',';
            json<<finalPath.coordinate(finalPath.center(),bond);
        }
        json<<']';
        json<<",\"guided_bond\":{\"enabled\":"<<(proposalMode==PhononProposal::GuidedBond)
            <<",\"beta\":"<<guide.beta<<",\"reference\":\""<<guideSource<<"\",\"force_by_type\":["
            <<guideByType[0]<<','<<guideByType[1]<<','<<guideByType[2]
            <<"],\"physical_source_changed\":false,\"trial_changed\":false,\"adaptation\":false}";
        json<<",\"gamma\":{\"every\":"<<gammaOptions.every<<",\"beta\":"<<gammaOptions.beta<<",\"proposal\":\"real_normalized_open_path_pCN_full_MH\",\"proposed_measurement\":"<<(sampler.gammaCounts()-gammaWarm).proposed<<",\"accepted_measurement\":"<<(sampler.gammaCounts()-gammaWarm).accepted<<"}";
        json<<",\"extra_q_pool_size\":"<<extraPool.size();
        json<<",\"q_mode_pool\":{\"bond_types\":\""<<qBonds<<"\",\"selection\":\"fixed_uniform\",\"counts_scope\":\"warmup_plus_measurement\",\"entries\":[";
        for(std::size_t k=0;k<qPool.size();++k) {
            if(k)json<<',';const auto &c=sampler.qPoolCounts()[k];
            json<<"{\"n1\":"<<qPool[k].n1<<",\"n2\":"<<qPool[k].n2<<",\"rank\":"<<qPool[k].basis.size()
                <<",\"proposed\":"<<c.proposed<<",\"accepted\":"<<c.accepted<<'}';
        }
        json<<"]},\"elapsed_seconds\":"<<elapsed<<',';
        json<<"\"q_mode\":{\"every\":"<<qOptions.every<<",\"n1\":"<<qOptions.n1<<",\"n2\":"<<qOptions.n2<<",\"beta\":"<<qOptions.beta<<",\"diagnostic_every\":"<<qDiagnosticEvery<<",\"basis_rank\":"<<qOptions.basis.size()<<",\"diagnostic_observables\":\"raw_unreweighted\",\"diagnostic_counts\":\"measurement_stage_cumulative\"},\"move_counts\":{";
        bool firstCounts=true;
        for(const auto &entry:std::vector<std::pair<std::string,MoveCounts>>{{"base_warmup",warm-qWarm-gammaWarm},{"q_warmup",qWarm},{"base_measurement",production-(sampler.qCounts()-qWarm)-(sampler.gammaCounts()-gammaWarm)},{"q_measurement",sampler.qCounts()-qWarm},{"gamma_warmup",gammaWarm},{"gamma_measurement",sampler.gammaCounts()-gammaWarm}}) {
            if(!firstCounts)json<<',';firstCounts=false;
            const auto &c=entry.second;jsonString(json,entry.first);
            json<<":{\"proposed\":"<<c.proposed<<",\"accepted\":"<<c.accepted<<",\"metropolis_rejected\":"<<c.metropolis_rejected<<",\"untrusted_rejected\":"<<c.untrusted_rejected<<",\"first_failures\":"<<c.first_failures<<",\"reference_recovered\":"<<c.reference_recovered<<",\"reference_failed\":"<<c.reference_failed<<'}';
        }
        json<<"},";
        numerics.metadata(json);
        json<<",\"ratio_defined\":"<<summary.ratio_defined
            <<",\"cache_diagnostics\":{\"local_proposals\":"<<state.cacheDiagnostics().local_proposals
            <<",\"local_node_rebuilds\":"<<state.cacheDiagnostics().local_node_rebuilds
            <<",\"local_node_reuses\":"<<state.cacheDiagnostics().local_node_reuses
            <<",\"local_merged_squares\":"<<state.cacheDiagnostics().local_merged_squares
            <<",\"bond_refresh_full_rebuilds\":"<<state.cacheDiagnostics().bond_refresh_full_rebuilds<<"}"
            <<",\"jackknife_defined\":"<<summary.jackknife_defined<<",\"phase_resolved_5sigma\":"<<summary.phase_resolved
            <<",\"too_few_blocks\":"<<summary.too_few_blocks<<",\"block_independence_verified\":false,\"average_phase\":";
        estimate(json,summary.phase); json<<",\"phase_diagnostics\":"; phaseDiagnostics(json,summary);
        json<<",\"observables\":{";
        for(std::size_t n=0;n<statistics.names().size();++n) {
            if(n) json<<','; json<<'"'<<statistics.names()[n]<<"\":"; estimate(json,summary.observables[n],summary.ratio_defined);
        }
        json<<"},\"doubled_block_count\":"<<doubled.blocks<<",\"doubled_block_observables\":{";
        for(std::size_t n=0;n<doubled.observables.size();++n) {
            if(n) json<<','; json<<'"'<<statistics.names()[n]<<"\":"; estimate(json,doubled.observables[n],doubled.ratio_defined);
        }
        json<<"},\"vbs_enabled\":"<<bool(vbsEnabled);
        if(fieldEnabled) {
            json<<",\"vbs_field\":{\"schema\":\"vbs-finite-field-v1\",\"h\":"<<field.h()
                <<",\"n\":["<<field.n1()<<','<<field.n2()<<"],\"q_fraction\":["<<field.u()<<','<<field.v()<<"]"
                <<",\"phase\":"<<field.phase()<<",\"profile_mean_square\":"<<field.profileMeanSquare()
                <<",\"raw_offset\":"<<fieldOffset
                <<",\"source\":\"H_h=H-h sum_r cos(q.r+phase) X_z(r); unrescaled X\""
                <<",\"trial_endpoints\":\"fixed under h differentiation\",\"translation_invariant\":"<<(field.h()==0.)
                <<",\"response\":\"joint phase ratios: <m*score>-<m><score>; score=d log W_h/dh with exact forced Mehler kernel\""
                <<",\"order\":\"absD2 and absD4 are raw fixed allowed q moments; B4_complex=absD4/absD2^2; finite field breaks symmetry\""
                <<",\"blocking\":[";
            auto fb=statistics.blocks();bool comma=false;
            while(!fb.empty()) {
                if(comma)json<<',';comma=true;
                const auto ds=summarizeDerivedBlocks(fb,[&](const auto &r){return VBSField::derive(r,fieldOffset);});
                json<<"{\"blocks\":"<<fb.size()<<",\"block_samples\":"<<fb.front().count<<",\"derived\":";
                estimatesDerived(json,VBSField::derivedNames(),ds);json<<'}';
                if(fb.size()%2 || fb.size()<4)break;fb=mergeAdjacentBlocks(fb);
            }
            json<<"]}";
        }
        if(susceptibilityEnabled) {
            const auto &p=finalPath;
            json<<",\"vbs_susceptibility\":{\"schema\":\"vbs-static-response-v1\",\"raw_offset\":"<<susceptibilityOffset
                <<",\"stride_per_q\":24,\"interpretation\":\"finite-projection central observable response; ground-state limit not assumed\""
                <<",\"source\":\"H_h=H-h sum_r f(r) X_z(r); cos f=cos(q.r); sin f=-sin(q.r)\""
                <<",\"trial_endpoints\":\"fixed under source differentiation\",\"tau_total\":"<<p.totalTau()
                <<",\"center_node\":"<<p.center()<<",\"center_tau\":"<<p.center()*p.dt()
                <<",\"endpoint_score_weight\":"<<VBSSusceptibility::nodeWeight(p,0)
                <<",\"interior_score_weight\":"<<2*VBSSusceptibility::nodeWeight(p,0)
                <<",\"score_weights\":\"c=tanh(omega*dt/2)/omega; endpoints c, interiors 2c; exact Mehler source score\""
                <<",\"normalization\":\"D=F(center)/Nc; I=sum_j w_j F(node_j)/Nc; Raw=Nc D_mu conj(I_nu)\""
                <<",\"connected\":\"joint phase ratios then Raw-Nc <D_mu><conj(I_nu)>; real sources use saved component means; joint delete-one-block errors; no Hermitization\"}";
        }
        if(amplitudeEnabled) {
            json<<",\"vbs_amplitude\":{\"schema\":\"vbs-amplitude-v1\",\"raw_offset\":"<<amplitudeOffset
                <<",\"stride_per_q\":"<<(amplitudeResponse?15:9)<<",\"response_enabled\":"<<bool(amplitudeResponse)
                <<",\"normalization\":\"D_mu=F_mu/Nc; A_mu=abs(D_mu); raw A2 and A4 are per configuration\""
                <<",\"source_epsilon\":0,\"source\":\"H_epsilon=H_0-epsilon*Nc*abs(D_mu(q))\""
                <<",\"trial_endpoints\":\"fixed under source differentiation\",\"tau_total\":"<<finalPath.totalTau()
                <<",\"center_node\":"<<finalPath.center()
                <<",\"integral_weights\":\"trapezoid: dt/2 endpoints, dt interiors; modulus before integration\""
                <<",\"finite_dt_source_definition\":\"W_epsilon=W_0*exp(epsilon*Nc*IA); symmetric diagonal insertions\""
                <<",\"interpretation\":\"finite-dt finite-projection modulus-source response; continuum and ground-state limits require convergence\""
                <<",\"derived\":\"common phase ratios first; VA=Nc*(<A2>-<A>*<A>), B4=<A4>/<A2>^2, chiA=Nc*(<AIA>-<A>*<IA>); joint errors in analyze_vbs_amplitude.py\"}";
        }
        if(vbsEnabled) {
            const bool gridBlock=(vbsMode==VBSMode::All || vbsMode==VBSMode::Both);
            const bool selBlock=(vbsMode==VBSMode::Selected || vbsMode==VBSMode::Both);
            json<<",\"vbs\":{\"normalization\":\"F=sum_R exp(-iqR)X; Delta=F/Nc; S=<F F^dagger>/Nc; mz2=Szz/Nc\","
                <<"\"eigenvalues_from\":\"Hermitian part of ensemble S, no negative clipping\",";
            if(vbsMode==VBSMode::All) {
                json<<"\"connected_background\":\"mean_Delta_mu(0)*mean_Delta_nu(0), joint jackknife\",";
            } else {
                json<<"\"connected_background\":\"ensemble connected C not computed in selected mode; "
                    <<"'proj' outputs project F_mu(Q)-W_L(Q)F_mu(0) per configuration, which is not the ensemble connected correlator\",";
            }
            json<<"\"Nc\":"<<model.lattice().cells();
            json<<",\"mode\":\""<<vbsModeName(vbsMode)<<"\""
                <<",\"estimator\":\""<<(vbsMode==VBSMode::All?"pbc_grid_fft":"finite_window_fourier")<<"\""
                <<",\"Nq\":"<<(vbsMode==VBSMode::All?model.lattice().cells():vbs.nq())
                <<",\"spatial_boundary_remains_pbc\":true,\"twist_applied\":false"
                <<",\"pbc_geometry\":\""<<(shear==0?"axis_aligned":"sheared")<<"\",\"shear\":"<<shear
                <<",\"geometry_id\":\""<<model.lattice().geometryId()<<"\""
                <<",\"real_space_correlation\":\""<<(vbsMode==VBSMode::Selected?"unavailable_selected_mode":"full_pbc_grid")<<"\""
                <<",\"window_definition\":\"W_L(Q)=(1/Nc)sum_R exp(-iQ.R) over rectangular representatives; W=0 is not an on-grid test\""
                <<",\"projected_definition\":\"F_proj_mu(Q)=F_mu(Q)-W_L(Q)F_mu(0), applied per configuration; not the ensemble connected correlator\""
                <<",\"candidate_max_rule\":\"LEGACY: mixed grid+off-grid raw mz2 maximum over listed non-Q=0 points; for a uniform configuration off-grid raw mz2 can be entirely uniform-background leakage (Xbar^2 |W_L(Q)|^2); superseded by vbs.candidates.groups (grid_raw / offgrid_raw / offgrid_projected are reported separately and never mixed)\"";
            json<<",\"primitive\":[";
            for(int i=0;i<2;++i) { if(i)json<<',';const auto a=model.lattice().primitive(i);json<<'['<<a.x<<','<<a.y<<']'; }
            json<<"],\"reciprocal\":[";
            for(int i=0;i<2;++i) { if(i)json<<',';const auto b=model.lattice().reciprocal(i);json<<'['<<b.x<<','<<b.y<<']'; }
            json<<"]";
            if(gridBlock) {
                json<<",\"momenta\":[";
                for(std::size_t q=0;q<vbs.momenta().size();++q) {
                    if(q)json<<',';const auto &k=vbs.momenta()[q];
                    json<<"{\"index\":"<<q<<",\"n\":["<<k.n1<<','<<k.n2<<"],\"reciprocal_fraction\":["<<k.reciprocal_x<<','<<k.reciprocal_y
                        <<"],\"cartesian\":["<<k.cartesian.x<<','<<k.cartesian.y<<"]}";
                }
                json<<"]";
            }
            if(selBlock) {
                json<<",\"qpoints\":[";
                for(int j=0;j<vbs.nq();++j) { if(j)json<<','; qPointJson(json,vbs.qpoints()[j],j); }
                json<<"],\"candidate_set\":{\"excludes_Q0\":true,\"participant_indices\":[";
                const auto &ci=vbs.candidateInfo();
                for(std::size_t t=0;t<ci.participants.size();++t) {
                    if(t)json<<','; json<<ci.participants[t];
                }
                json<<"],\"participant_labels\":[";
                for(std::size_t t=0;t<ci.participants.size();++t) {
                    if(t)json<<','; json<<'"'<<vbs.qpoints()[ci.participants[t]].label<<'"';
                }
                json<<"],\"excluded_zero_indices\":[";
                for(std::size_t t=0;t<ci.excludedZero.size();++t) {
                    if(t)json<<','; json<<ci.excludedZero[t];
                }
                json<<"]}";
            }
            json<<",\"blocking\":[";
            auto merged=statistics.blocks(); int level=0;
            while(!merged.empty()) {
                if(level++)json<<',';
                const auto vbsBlockingBegin=std::chrono::steady_clock::now();
                const auto raw=summarizeBlocks(merged);
                const auto derived=summarizeDerivedBlocks(merged,[&](const auto &r){
                    const std::vector<Complex> original(r.begin()+vbsOffset,r.begin()+susceptibilityOffset);
                    return vbs.derive(original,0);});
                vbsBlockingSeconds+=std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-vbsBlockingBegin).count();
                json<<"{\"blocks\":"<<merged.size()<<",\"block_samples\":"<<merged.front().count
                    <<",\"block_sweeps\":"<<merged.front().count*schedule.stride
                    <<",\"jackknife_defined\":"<<raw.jackknife_defined<<",\"phase_resolved_5sigma\":"<<raw.phase_resolved
                    <<",\"average_phase\":";estimate(json,raw.phase);
                json<<",\"phase_diagnostics\":";phaseDiagnostics(json,raw);
                json<<",\"raw\":";estimates(json,names,raw);
                json<<",\"derived\":";estimatesDerived(json,vbs.derivedNames(),derived);json<<'}';
                if(merged.size()%2 || merged.size()<4)break;
                merged=mergeAdjacentBlocks(merged);
            }
            json<<"]";
            if(selBlock) {
                // Categorical peak-selection diagnostics at the finest block
                // level: same phase blocks, same complex reweighting as the
                // joint jackknife; no RNG, no layout change.
                const auto &finest=statistics.blocks();
                const auto width=statistics.names().size();
                Complex phaseTotal=0.; std::vector<Complex> sums(width,0.);
                std::uint64_t countTotal=0;
                for(const auto &b:finest) {
                    phaseTotal+=b.phase; countTotal+=b.count;
                    for(std::size_t jj=0;jj<width;++jj) sums[jj]+=b.numerator[jj];
                }
                std::vector<Complex> fullRatios(width);
                const bool ratioDefined=std::abs(phaseTotal)>1e-12*double(countTotal);
                if(ratioDefined) for(std::size_t jj=0;jj<width;++jj) fullRatios[jj]=sums[jj]/phaseTotal;
                const std::size_t selOffset=vbsOffset+std::size_t(
                    vbsMode==VBSMode::Both?VBSObservables::allStride(model.lattice().cells()):0);
                const auto candidateDiagnosticsBegin=std::chrono::steady_clock::now();
                const auto replicas=candidateReplicas(finest);
                const auto diag=runCandidateDiagnostics(vbs,fullRatios,selOffset,replicas,ratioDefined);
                candidateDiagnosticsSeconds+=std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-candidateDiagnosticsBegin).count();
                writeCandidatesJson(json,vbs,diag);
            }
            json<<"}";
        }
        json<<",\"fermions_enabled\":"<<bool(fermionsEnabled);
        if(fermionsEnabled) {json<<",\"fermions\":";fermions.writeMetadata(json);}
        json<<",\"measure_energy\":"<<bool(measureEnergy)
            <<",\"measurement_green_enabled\":"<<bool(fermionsEnabled || measureEnergy || bondDiagnosticsEnabled)
            <<",\"measure_total_energy\":"<<bool(measureTotalEnergy)
            <<",\"total_energy_estimator_available\":"<<bool(measureTotalEnergy)<<",\"vbs_conclusion\":\"not_assessed\"";
        if(measureTotalEnergy) {
            json<<",\"total_energy\":{\"schema\":\"phonon-link-total-energy-v1\",\"raw_offset\":"<<totalEnergyOffset
                <<",\"normalization\":\"extensive; no zero-point subtraction\",\"cuts\":["<<finalPath.intervals()-1<<','<<finalPath.intervals()+1<<']'
                <<",\"times\":["<<(finalPath.center()-.5)*finalPath.dt()<<','<<(finalPath.center()+.5)*finalPath.dt()<<']'
                <<",\"combination\":\"equal average of adjacent central phonon-link midpoints\""
                <<",\"legacy_majorana_energy\":\"unchanged central matter half-step cut; not a component of this link estimator\""
                <<",\"phonon_estimator\":\"exact forced Mehler bridge first/second moments; kinetic=-d_dt log K_force-potential-source\""
                <<",\"matter_estimator\":\"certified odd-cut Green; bridge mean for linear coupling\""
                <<",\"external_source\":\"-sum force_b*X_b, included in energy_total\""
                <<",\"finite_dt_location_differs_from_node_center\":true,\"variational_bound_claimed\":false"
                <<",\"rng_consumed\":false,\"timing_scope\":\"this process only\"}";
        }
        const double jsonPreTimingSeconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-jsonBuildBegin).count();
        json<<",\"timing\":{\"initialization_seconds\":"<<initializationSeconds
            <<",\"sampling_seconds\":"<<samplingWallSeconds
            <<",\"warmup_seconds\":"<<samplingTimings.warmup_seconds
            <<",\"measurement_seconds\":"<<samplingTimings.measurement_seconds
            <<",\"measurement_materialization_seconds\":"<<samplingTimings.measurement_materialization_seconds
            <<",\"measurement_callback_seconds\":"<<samplingTimings.measurement_callback_seconds
            <<",\"total_energy_seconds_included_in_measurement_callback\":"<<totalEnergySeconds
            <<",\"total_energy_measurements_this_process\":"<<totalEnergyMeasurements
            <<",\"statistics_add_seconds\":"<<samplingTimings.statistics_add_seconds
            <<",\"block_callback_seconds\":"<<samplingTimings.block_callback_seconds
            <<",\"statistics_finalize_seconds\":"<<statisticsFinalizeSeconds
            <<",\"vbs_blocking_seconds\":"<<vbsBlockingSeconds
            <<",\"candidate_diagnostics_seconds\":"<<candidateDiagnosticsSeconds
            <<",\"block_file_output_seconds\":"<<blockFileOutputSeconds
            <<",\"json_pre_timing_seconds\":"<<jsonPreTimingSeconds
            <<",\"elapsed_excludes_json_build_and_output\":true}"
            <<"}\n";
        const std::string jsonText=json.str();
        const double jsonBuildSeconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-jsonBuildBegin).count();
        const auto jsonOutputBegin=std::chrono::steady_clock::now();
        // Final summaries can take appreciable time. A request received after
        // the last sweep still saves that last complete state before output.
        const auto finalRequest=checkpointSignals.consume();
        if(finalRequest.first)saveCheckpoint(finalRequest.second);
        if(finalRequest.second)throw CheckpointStopped{};
        writeResult(prefix,jsonText);
        const double jsonOutputSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-jsonOutputBegin).count();
        const auto stdoutOutputBegin=std::chrono::steady_clock::now();
        std::cout<<jsonText; std::cout.flush();
        const double stdoutOutputSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-stdoutOutputBegin).count();
        std::cerr<<"TIMING_VBS_BLOCKING_SECONDS="<<vbsBlockingSeconds<<"\n"
                 <<"TIMING_CANDIDATE_DIAGNOSTICS_SECONDS="<<candidateDiagnosticsSeconds<<"\n"
                 <<"TIMING_JSON_BUILD_SECONDS="<<jsonBuildSeconds<<"\n"
                 <<"TIMING_BLOCK_FILE_OUTPUT_SECONDS="<<blockFileOutputSeconds<<"\n"
                 <<"TIMING_RESULTS_FILE_OUTPUT_SECONDS="<<jsonOutputSeconds<<"\n"
                 <<"TIMING_STDOUT_OUTPUT_SECONDS="<<stdoutOutputSeconds<<"\n"
                 <<"JSON output seconds (excluded from elapsed_seconds): "<<jsonOutputSeconds<<"\n";
        std::cerr<<"Statistics are finite-run block estimates; block independence/projection/size convergence not established.\n";
        if(!summary.phase_resolved || summary.too_few_blocks || !summary.jackknife_defined)
            std::cerr<<"STATISTICAL LIMITATION: unresolved phase denominator, insufficient blocks or undefined jackknife.\n";
        return 0;
    } catch(const CheckpointStopped &) {
        numerics.status(prefix,"checkpointed","saved at completed sweep boundary; resume with a fresh output prefix");
        std::cerr<<"Saved checkpoint and stopped cleanly (exit75); this is not a completed production result.\n";
        return 75;
    } catch(const std::exception &e) {
        if(numerics.owns_output) {
            try { numerics.status(prefix,"aborted",e.what()); }
            catch(const std::exception &io) {std::cerr<<"Could not write aborted status: "<<io.what()<<'\n';}
        }
        std::cerr<<"kitaev_ssh_mc FAILED (not a normal rejection): "<<e.what()
                 <<"\nCompleted raw blocks, if any: "<<prefix<<"_blocks.csv\n"; return 1;
    }
}
