#ifndef KITAEV_SSH_VBS_CANDIDATES_H
#define KITAEV_SSH_VBS_CANDIDATES_H

#include "vbs_observables.h"
#include "phase_statistics.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace kitaev_ssh {

// Categorical peak-selection diagnostics for the selected-mode candidate
// list. These are diagnostics, NOT continuous observables: they never enter
// the numeric derived vector, they define no "standard error of an index",
// they consume no RNG and they change no sampling quantity. The numeric
// derived layout of vbs_observables.h (including the legacy candidate_max_*
// entries) is untouched.
//
// Scoring rule (fixed, documented): the REAL part of the reweighted ensemble
//   mz2_raw(Q)       = Re [ Ssel_s{j}_zz ratio ] / Nc,
//   mz2_projected(Q) = Re [ Sprj_s{j}_zz ratio ] / Nc.
// Imaginary parts are kept per point as numerical diagnostics. Peak selection
// never takes |value| and never clips negative estimates to zero.
//
// Q=0 exclusion uses the unified normalization of vbs_momenta.h: fractional
// coordinates reduced mod the reciprocal lattice, on-grid test by
// makeQPoint, and (nx,ny)==(0,0) exactly. A momentum congruent to 0 mod the
// reciprocal lattice can never be off-grid, so the rule is single-valued.
//
// Winner-SET tie tolerance (documented constants, no random tie-break):
//   |a-b| <= candidateTieAbs + candidateTieRel*|group max score|.
// Replica validity mirrors summarizeBlocks: a delete-one-block ratio is
// invalid when |phase_total-block_phase| <= candidateReplicaDenominatorSlack
// * (remaining count) or when a ratio is non-finite. Invalid replicas are
// counted. They may be used for a conditional selection diagnostic, but they
// never become a smaller ordinary jackknife ensemble.

inline constexpr double candidateTieAbs=1e-12;
inline constexpr double candidateTieRel=1e-9;
inline constexpr double candidateReplicaDenominatorSlack=1e-12;

struct CandidatePointDiag {
    int index=-1;
    std::string label;
    std::string group;            // "grid_nonzero" / "offgrid" / "zero_excluded"
    double u=0.,v=0.;             // normalized fractional coordinates
    bool on_grid=false;
    double W_abs2=0.;             // |W_L(Q)|^2, uniform-background leakage probe
    double mz2_raw_real=std::numeric_limits<double>::quiet_NaN(),
           mz2_raw_imag=std::numeric_limits<double>::quiet_NaN();
    double mz2_proj_real=std::numeric_limits<double>::quiet_NaN(),
           mz2_proj_imag=std::numeric_limits<double>::quiet_NaN();
    double mz2_raw_error_real=std::numeric_limits<double>::quiet_NaN(),
           mz2_raw_error_imag=std::numeric_limits<double>::quiet_NaN();
    double mz2_raw_covariance=std::numeric_limits<double>::quiet_NaN();
    double mz2_proj_error_real=std::numeric_limits<double>::quiet_NaN(),
           mz2_proj_error_imag=std::numeric_limits<double>::quiet_NaN();
    double mz2_proj_covariance=std::numeric_limits<double>::quiet_NaN();
    bool raw_center_defined=false;
    std::string raw_center_reason="not_evaluated";
    bool projected_center_defined=false;
    std::string projected_center_reason="not_evaluated";
    bool raw_jackknife_defined=false;
    std::string raw_jackknife_reason="not_evaluated";
    bool projected_jackknife_defined=false;
    std::string projected_jackknife_reason="not_evaluated";
    int plus_minus_class=-1;
    bool class_self_conjugate=false;
};

struct CandidateGroupDiag {
    std::string name;             // "grid_raw" / "offgrid_raw" / "offgrid_projected"
    std::string estimator;        // "raw" / "projected"
    bool available=false;         // non-empty participant set in the q list
    std::vector<int> participants;
    std::vector<std::string> participant_labels;
    bool center_defined=false;    // full-sample reweighted ratios defined
    std::string center_reason="not_evaluated";
    double max_mz2_real=std::numeric_limits<double>::quiet_NaN();
    std::vector<int> winner_set;  // tied winners (empty when undefined)
    std::vector<std::string> winner_labels;
    int winner_display_index=-1;  // legacy input-order display number only
    std::uint64_t total_replicas=0,valid_replicas=0,invalid_replicas=0;
    bool jackknife_defined=false;
    std::string jackknife_reason="not_evaluated";
    bool stability_defined=false; // conditional selection diagnostic
    std::string stability_reason="not_evaluated";
    bool conditional_on_valid_replicas=false;
    std::vector<std::uint64_t> unique_win_counts,tie_membership_counts;
    std::vector<double> unique_win_frequency,tie_membership_frequency;
    // ±Q class view of the same replicas: a class is "tied-in" when any of
    // its members is in the replica winner set, and a "full-class win" when
    // the winner set equals exactly that class's member set. Because the
    // tie-membership indicator is an OR over members, class and member
    // frequencies do not sum to one in general.
    std::vector<std::uint64_t> class_tie_counts,class_full_set_counts;
    std::vector<double> class_tie_frequency,class_full_set_frequency;
};

struct CandidateDiagnostics {
    bool ratio_defined=false;
    bool center_defined=false;
    std::string center_reason="not_evaluated";
    bool jackknife_defined=false;
    std::string jackknife_reason="not_evaluated";
    bool stability_defined=false;
    std::string stability_reason="not_evaluated";
    std::uint64_t total_replicas=0,valid_replicas=0,invalid_replicas=0;
    std::vector<std::string> invalid_replica_reasons;
    std::vector<CandidatePointDiag> points;
    std::vector<CandidateGroupDiag> groups; // fixed order: grid_raw, offgrid_raw, offgrid_projected
    struct PlusMinusClass {
        std::vector<int> members;     // participant indices (never Q=0)
        bool self_conjugate=false;    // Q == -Q mod reciprocal lattice
        double u=0.,v=0.;             // canonical representative (min of Q,-Q)
    };
    std::vector<PlusMinusClass> classes;
};

namespace candidate_detail {
// All accessors take the absolute offset of the selected block inside the
// full-width ratio vector (vbsOffset in selected mode, vbsOffset+allStride in
// both mode) and index the selected-block layout of VBSObservables.
inline double rawScore(const std::vector<Complex> &r,std::size_t off,int j,int nc) {
    return r[off+VBSObservables::selS(j)+3*2+2].real()/double(nc);
}
inline double projScore(const std::vector<Complex> &r,std::size_t off,int j,int nq,int nc) {
    return r[off+VBSObservables::selPS(j,nq)+3*2+2].real()/double(nc);
}
inline Complex rawZz(const std::vector<Complex> &r,std::size_t off,int j) {
    return r[off+VBSObservables::selS(j)+3*2+2];
}
inline Complex projZz(const std::vector<Complex> &r,std::size_t off,int j,int nq) {
    return r[off+VBSObservables::selPS(j,nq)+3*2+2];
}
inline bool finiteScores(const std::vector<double> &s) {
    for(double v:s) if(!std::isfinite(v)) return false;
    return true;
}
// Winner set within the recorded tie tolerance; deterministic, no RNG.
inline std::vector<int> winnerSet(const std::vector<int> &participants,
                                  const std::vector<double> &score) {
    double best=std::numeric_limits<double>::quiet_NaN();
    if(score.size()!=participants.size())
        throw std::invalid_argument("candidate winner score dimension");
    for(std::size_t t=0;t<participants.size();++t)
        if(std::isnan(best)||score[t]>best) best=score[t];
    const double tol=candidateTieAbs+candidateTieRel*std::abs(best);
    std::vector<int> set;
    for(std::size_t t=0;t<participants.size();++t)
        if(score[t]>=best-tol) set.push_back(participants[t]);
    return set;
}
inline bool jackErrors(const std::vector<Complex> &v,double &errReal,double &errImag,
                       double &covariance) {
    errReal=errImag=covariance=std::numeric_limits<double>::quiet_NaN();
    const std::size_t r=v.size();
    if(r<2) return false;
    Complex mean=0.;
    for(auto z:v) {
        if(!finiteComplex(z)) return false;
        mean+=z/double(r);
    }
    double vr=0.,vi=0.,cov=0.;
    for(auto z:v) {
        z-=mean;
        vr+=z.real()*z.real(); vi+=z.imag()*z.imag(); cov+=z.real()*z.imag();
    }
    const double factor=(r-1.)/double(r);
    const double er=std::sqrt(factor*vr),ei=std::sqrt(factor*vi),co=factor*cov;
    if(!std::isfinite(er)||!std::isfinite(ei)||!std::isfinite(co)) return false;
    errReal=er; errImag=ei; covariance=co;
    return true;
}
} // namespace candidate_detail

// Categorical diagnostics from the same phase blocks the joint jackknife
// uses. `fullRatios` is the FULL-width ratio vector of all observables;
// `offset` is the absolute start of the selected block inside it (in both
// mode this is vbsOffset+allStride(Nc)). `replicas` hold the same full-width
// leave-one-block ratio vectors; entries flagged invalid are counted and
// excluded. Ratios are (sums - block.numerator)/(phase - block.phase),
// exactly as in summarizeDerivedBlocks.
struct CandidateReplica {
    std::vector<Complex> ratios;
    bool valid=false;
    std::string invalid_reason;
};
inline CandidateDiagnostics runCandidateDiagnostics(const VBSObservables &vbs,
                                                    const std::vector<Complex> &fullRatios,
                                                    std::size_t offset,
                                                    const std::vector<CandidateReplica> &replicas,
                                                    bool ratioDefined) {
    using namespace candidate_detail;
    CandidateDiagnostics out;
    out.ratio_defined=ratioDefined;
    out.total_replicas=replicas.size();
    const int nq=vbs.nq(),nc=vbs.lattice().cells();
    if(fullRatios.size()!=offset+std::size_t(vbs.selectedStride()))
        throw std::invalid_argument("candidate diagnostics ratio dimension");
    for(const auto &rep:replicas)
        if(rep.valid && rep.ratios.size()!=fullRatios.size())
            throw std::invalid_argument("candidate diagnostics replica dimension");
    const auto &qs=vbs.qpoints();
    // Classification and ±Q classes over the whole candidate list.
    out.points.resize(nq);
    for(int j=0;j<nq;++j) {
        auto &p=out.points[j];
        const auto &q=qs[j];
        p.index=j; p.label=q.label; p.u=q.u; p.v=q.v; p.on_grid=q.on_grid;
        p.W_abs2=std::norm(q.window);
        const bool zero=q.on_grid && q.n1==0 && q.n2==0;
        p.group=zero?"zero_excluded":(q.on_grid?"grid_nonzero":"offgrid");
        const VBSQPoint m=minusQPoint(q,vbs.lattice());
        const bool self=fractionalDelta(q.u,m.u)==0. && fractionalDelta(q.v,m.v)==0.;
        if(!zero) {
            double cu=q.u,cv=q.v;
            if(!self && (m.u<q.u-1e-12 || (std::abs(m.u-q.u)<=1e-12 && m.v<q.v-1e-12))) {cu=m.u;cv=m.v;}
            int cls=-1;
            for(std::size_t c=0;c<out.classes.size();++c) {
                const auto &k=out.classes[c];
                if(std::abs(fractionalDelta(k.u,cu))<=1e-9 && std::abs(fractionalDelta(k.v,cv))<=1e-9) {cls=int(c);break;}
            }
            if(cls<0) {
                CandidateDiagnostics::PlusMinusClass k;
                k.u=cu; k.v=cv; k.self_conjugate=self;
                out.classes.push_back(k);
                cls=int(out.classes.size())-1;
            } else if(self) out.classes[cls].self_conjugate=true;
            p.plus_minus_class=cls;
            out.classes[cls].members.push_back(j);
        }
        p.class_self_conjugate=p.plus_minus_class>=0 && out.classes[p.plus_minus_class].self_conjugate;
    }
    // Groups: fixed order and fixed semantics, never mixed:
    //   grid_raw         = on-grid non-zero candidates, raw mz2;
    //   offgrid_raw      = off-grid candidates ONLY, raw mz2 (finite-window
    //                      search diagnostic, may contain uniform-background
    //                      leakage Xbar^2|W_L(Q)|^2);
    //   offgrid_projected= off-grid candidates ONLY, projected mz2 (uniform
    //                      mode removed per configuration; window effects and
    //                      the incommensurability mismatch are NOT removed).
    const char *groupNames[3]={"grid_raw","offgrid_raw","offgrid_projected"};
    const char *groupEst[3]={"raw","raw","projected"};
    std::vector<std::vector<int>> groupParticipants(3);
    for(int j=0;j<nq;++j) {
        if(out.points[j].group=="grid_nonzero") groupParticipants[0].push_back(j);
        if(out.points[j].group=="offgrid")     {groupParticipants[1].push_back(j);groupParticipants[2].push_back(j);}
    }
    out.groups.resize(3);
    for(int g=0;g<3;++g) {
        auto &G=out.groups[g];
        G.name=groupNames[g]; G.estimator=groupEst[g];
        G.participants=groupParticipants[g];
        G.available=!G.participants.empty();
        for(int j:G.participants) G.participant_labels.push_back(qs[j].label);
        G.unique_win_counts.assign(G.participants.size(),0);
        G.tie_membership_counts.assign(G.participants.size(),0);
        G.class_tie_counts.assign(out.classes.size(),0);
        G.class_full_set_counts.assign(out.classes.size(),0);
    }
    // Center values are checked per estimator. A usable full denominator is
    // necessary but not sufficient: every complex ratio used by a point must
    // also be finite. This preserves a valid center estimate while exposing a
    // local failure instead of turning it into zero.
    for(int j=0;j<nq;++j) {
        auto &p=out.points[j];
        if(!out.ratio_defined) {
            p.raw_center_reason="full_sample_ratio_undefined";
            p.projected_center_reason="full_sample_ratio_undefined";
            continue;
        }
        const Complex raw=rawZz(fullRatios,offset,j)/double(nc);
        const Complex projected=projZz(fullRatios,offset,j,nq)/double(nc);
        if(finiteComplex(raw)) {
            p.raw_center_defined=true;
            p.raw_center_reason="defined";
            p.mz2_raw_real=raw.real(); p.mz2_raw_imag=raw.imag();
        } else {
            p.raw_center_reason="nonfinite_raw_estimate";
        }
        if(finiteComplex(projected)) {
            p.projected_center_defined=true;
            p.projected_center_reason="defined";
            p.mz2_proj_real=projected.real(); p.mz2_proj_imag=projected.imag();
        } else {
            p.projected_center_reason="nonfinite_projected_estimate";
        }
    }
    for(int g=0;g<3;++g) {
        auto &G=out.groups[g];
        if(!G.available) { G.center_reason="empty_group"; continue; }
        std::vector<double> score;
        bool defined=true;
        for(int j:G.participants) {
            const auto &p=out.points[j];
            const bool pointDefined=g==2?p.projected_center_defined:p.raw_center_defined;
            if(!pointDefined) {
                G.center_reason="point_"+std::to_string(j)+"_"+
                    (g==2?p.projected_center_reason:p.raw_center_reason);
                defined=false; break;
            }
            score.push_back(g==2?p.mz2_proj_real:p.mz2_raw_real);
        }
        if(!defined) continue;
        G.center_defined=true; G.center_reason="defined";
        G.max_mz2_real=score[0];
        for(double s:score) G.max_mz2_real=std::max(G.max_mz2_real,s);
        G.winner_set=winnerSet(G.participants,score);
        G.winner_display_index=G.winner_set.empty()?-1:G.winner_set.front(); // display only
        for(int j:G.winner_set) G.winner_labels.push_back(qs[j].label);
    }
    // Replica pass: re-selection stability and per-point jackknife errors.
    std::vector<std::vector<Complex>> rawRep(nq),projRep(nq);
    for(const auto &rep:replicas) {
        bool replicaValid=rep.valid;
        std::string invalidReason=rep.invalid_reason;
        if(!replicaValid && invalidReason.empty()) invalidReason="invalid_replica";
        if(replicaValid && rep.ratios.size()!=fullRatios.size())
            throw std::invalid_argument("candidate diagnostics replica dimension");
        if(replicaValid) {
            for(const auto z:rep.ratios) if(!finiteComplex(z)) {
                replicaValid=false; invalidReason="nonfinite_replica_ratio"; break;
            }
        }
        if(!replicaValid) {
            ++out.invalid_replicas;
            out.invalid_replica_reasons.push_back(invalidReason.empty()?"invalid_replica":invalidReason);
            continue;
        }
        bool finiteAll=true;
        std::vector<double> rawV(nq),projV(nq);
        for(int j=0;j<nq;++j) {
            rawV[j]=rawScore(rep.ratios,offset,j,nc);
            projV[j]=projScore(rep.ratios,offset,j,nq,nc);
            if(!std::isfinite(rawV[j])||!std::isfinite(projV[j])||
               !finiteComplex(rawZz(rep.ratios,offset,j))||
               !finiteComplex(projZz(rep.ratios,offset,j,nq))) finiteAll=false;
        }
        if(!finiteAll) {
            ++out.invalid_replicas;
            out.invalid_replica_reasons.push_back("nonfinite_candidate_estimate");
            continue;
        }
        ++out.valid_replicas;
        for(int j=0;j<nq;++j) {
            rawRep[j].push_back(Complex(rawV[j],rawZz(rep.ratios,offset,j).imag()/double(nc)));
            projRep[j].push_back(Complex(projV[j],projZz(rep.ratios,offset,j,nq).imag()/double(nc)));
        }
        for(int g=0;g<3;++g) {
            auto &G=out.groups[g];
            if(!G.available) continue;
            std::vector<double> score;
            for(int j:G.participants) score.push_back(g==2?projV[j]:rawV[j]);
            if(!finiteScores(score)) continue;
            const auto set=winnerSet(G.participants,score);
            for(std::size_t t=0;t<G.participants.size();++t) {
                const int j=G.participants[t];
                if(set.size()==1 && set.front()==j) ++G.unique_win_counts[t];
                for(int w:set) if(w==j) {++G.tie_membership_counts[t];break;}
            }
            for(std::size_t c=0;c<out.classes.size();++c) {
                const auto &members=out.classes[c].members;
                std::size_t hit=0;
                for(int m:members) for(int w:set) if(w==m) {++hit;break;}
                if(hit) ++G.class_tie_counts[c];
                if(hit==members.size() && set.size()==members.size()) ++G.class_full_set_counts[c];
            }
        }
    }
    const bool completeReplicas=out.total_replicas>=2 && out.invalid_replicas==0;
    for(int g=0;g<3;++g) {
        auto &G=out.groups[g];
        G.total_replicas=out.total_replicas;
        G.valid_replicas=out.valid_replicas; G.invalid_replicas=out.invalid_replicas;
        if(!G.available) G.stability_reason="empty_group";
        else if(!G.center_defined) G.stability_reason="center_undefined";
        else if(out.valid_replicas<2) G.stability_reason="fewer_than_two_valid_replicas";
        else {
            G.stability_defined=true;
            G.conditional_on_valid_replicas=out.invalid_replicas>0;
            G.stability_reason=G.conditional_on_valid_replicas?
                "defined_conditionally_on_valid_replicas":"defined";
        }
        const double Rv=double(out.valid_replicas);
        for(std::size_t t=0;t<G.participants.size();++t) {
            G.unique_win_frequency.push_back(G.stability_defined?G.unique_win_counts[t]/Rv:
                                             std::numeric_limits<double>::quiet_NaN());
            G.tie_membership_frequency.push_back(G.stability_defined?G.tie_membership_counts[t]/Rv:
                                                 std::numeric_limits<double>::quiet_NaN());
        }
        for(std::size_t c=0;c<out.classes.size();++c) {
            G.class_tie_frequency.push_back(G.stability_defined?G.class_tie_counts[c]/Rv:
                                            std::numeric_limits<double>::quiet_NaN());
            G.class_full_set_frequency.push_back(G.stability_defined?G.class_full_set_counts[c]/Rv:
                                                 std::numeric_limits<double>::quiet_NaN());
        }
    }
    for(int j=0;j<nq;++j) {
        auto &p=out.points[j];
        if(!p.raw_center_defined) p.raw_jackknife_reason="center_undefined";
        else if(out.total_replicas<2) p.raw_jackknife_reason="fewer_than_two_replicas";
        else if(out.invalid_replicas>0) p.raw_jackknife_reason="incomplete_delete_one_block_replicas";
        else if(jackErrors(rawRep[j],p.mz2_raw_error_real,p.mz2_raw_error_imag,p.mz2_raw_covariance)) {
            p.raw_jackknife_defined=true; p.raw_jackknife_reason="defined";
        } else p.raw_jackknife_reason="nonfinite_or_undefined_jackknife";
        if(!p.projected_center_defined) p.projected_jackknife_reason="center_undefined";
        else if(out.total_replicas<2) p.projected_jackknife_reason="fewer_than_two_replicas";
        else if(out.invalid_replicas>0) p.projected_jackknife_reason="incomplete_delete_one_block_replicas";
        else if(jackErrors(projRep[j],p.mz2_proj_error_real,p.mz2_proj_error_imag,p.mz2_proj_covariance)) {
            p.projected_jackknife_defined=true; p.projected_jackknife_reason="defined";
        } else p.projected_jackknife_reason="nonfinite_or_undefined_jackknife";
    }
    for(int g=0;g<3;++g) {
        auto &G=out.groups[g];
        if(!G.available) G.jackknife_reason="empty_group";
        else if(!G.center_defined) G.jackknife_reason="center_undefined";
        else if(out.total_replicas<2) G.jackknife_reason="fewer_than_two_replicas";
        else if(!completeReplicas) G.jackknife_reason="incomplete_delete_one_block_replicas";
        else {
            G.jackknife_defined=true;
            for(int j:G.participants) {
                const auto &p=out.points[j];
                const bool defined=g==2?p.projected_jackknife_defined:p.raw_jackknife_defined;
                if(!defined) { G.jackknife_defined=false; G.jackknife_reason="point_jackknife_undefined"; break; }
            }
            if(G.jackknife_defined) G.jackknife_reason="defined";
        }
    }
    bool anyGroup=false,allCenter=true,allJack=true,allStability=true;
    for(const auto &G:out.groups) if(G.available) {
        anyGroup=true; allCenter=allCenter&&G.center_defined;
        allJack=allJack&&G.jackknife_defined; allStability=allStability&&G.stability_defined;
    }
    out.center_defined=anyGroup&&allCenter;
    out.jackknife_defined=anyGroup&&allJack;
    out.stability_defined=anyGroup&&allStability;
    if(!anyGroup) out.center_reason="no_candidate_groups";
    else if(!out.center_defined) out.center_reason="one_or_more_groups_center_undefined";
    else out.center_reason="defined";
    out.jackknife_reason=out.jackknife_defined?"defined":
        (out.center_defined?"one_or_more_groups_jackknife_undefined":"center_undefined");
    out.stability_reason=out.stability_defined?"defined":
        (out.center_defined?"one_or_more_groups_stability_undefined":"center_undefined");
    return out;
}

// Build the replica list from phase blocks: ratios and validity per block.
inline std::vector<CandidateReplica> candidateReplicas(const std::vector<PhaseBlock> &blocks) {
    if(blocks.empty()) return {};
    const auto width=blocks.front().numerator.size();
    PhaseBlock total;
    total.numerator.assign(width,0.);
    bool allBlockDataFinite=true;
    for(const auto &b:blocks) {
        if(!b.count || b.count!=blocks.front().count || b.numerator.size()!=width)
            throw std::invalid_argument("candidate diagnostics require equal nonempty blocks");
        if(!finiteComplex(b.phase)) allBlockDataFinite=false;
        for(const auto z:b.numerator) if(!finiteComplex(z)) allBlockDataFinite=false;
        total.count+=b.count; total.phase+=b.phase;
        for(std::size_t j=0;j<width;++j) total.numerator[j]+=b.numerator[j];
    }
    std::vector<CandidateReplica> out;
    out.reserve(blocks.size());
    for(const auto &b:blocks) {
        CandidateReplica rep;
        const std::uint64_t remaining=total.count-b.count;
        const Complex denom=total.phase-b.phase;
        rep.valid=remaining>0;
        if(!rep.valid) rep.invalid_reason="fewer_than_two_blocks";
        else if(!allBlockDataFinite) {
            rep.valid=false; rep.invalid_reason="nonfinite_block_data";
        }
        else if(!finiteComplex(denom) ||
                !(std::abs(denom)>candidateReplicaDenominatorSlack*double(remaining))) {
            rep.valid=false; rep.invalid_reason="invalid_delete_one_denominator";
        }
        if(rep.valid) {
            rep.ratios.resize(width);
            for(std::size_t j=0;j<width;++j) {
                rep.ratios[j]=(total.numerator[j]-b.numerator[j])/denom;
                if(!finiteComplex(rep.ratios[j])) {
                    rep.valid=false; rep.invalid_reason="nonfinite_replica_ratio"; break;
                }
            }
        }
        out.push_back(std::move(rep));
    }
    return out;
}

} // namespace kitaev_ssh
#endif
