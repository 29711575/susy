#ifndef KITAEV_SSH_PHONON_MC_H
#define KITAEV_SSH_PHONON_MC_H

#include "phonon_projector.h"
#include "phase_statistics.h"
#include "phonon_q_mode.h"
#include "phonon_guided_proposal.h"
#include <chrono>
#include <numeric>
#include <functional>

namespace kitaev_ssh {
struct ConditionalGaussian {
    double mean,variance;
    double logDensity(double x) const { return -.5*(std::log(2*pi*variance)+(x-mean)*(x-mean)/variance); }
};
inline ConditionalGaussian freeConditional(const PhononPath &path,int node,int bond) {
    path.coordinate(node,bond);
    const double r=path.decay(),alpha=path.alpha();
    ConditionalGaussian q;
    if(node==0 || node==path.intervals())
        q={r*path.coordinate(node==0?1:node-1,bond)+path.noiseFactor()*path.trialShift(bond)/2,path.noiseFactor()/(2*alpha)};
    else q={r/(1+r*r)*(path.coordinate(node-1,bond)+path.coordinate(node+1,bond)),path.noiseFactor()/(2*alpha*(1+r*r))};
    if(path.force(bond)!=0.) q.mean+=q.variance*path.force(bond)*path.sourceNodeWeight(node);
    if(!std::isfinite(q.mean) || !std::isfinite(q.variance) || q.variance<=0)
        throw std::runtime_error("invalid free-phonon conditional");
    return q;
}
struct MoveResult { double log_ratio, general_log_ratio; bool accepted; };
struct MoveCounts {
    std::uint64_t proposed=0,accepted=0,untrusted_rejected=0;
    std::uint64_t first_failures=0,reference_recovered=0,reference_failed=0,metropolis_rejected=0;
    MoveCounts operator-(const MoveCounts &b) const {
        return {proposed-b.proposed,accepted-b.accepted,untrusted_rejected-b.untrusted_rejected,
            first_failures-b.first_failures,reference_recovered-b.reference_recovered,
            reference_failed-b.reference_failed,metropolis_rejected-b.metropolis_rejected};
    }
    double acceptance() const { return proposed?double(accepted)/proposed:0.; }
};
// An uncertifiable numerical weight is not a certified zero weight. Production
// transitions retry the SAME candidate once, then fail closed. No rejection
// compatibility switch is provided.
inline constexpr const char *numericalFailurePolicy="uncached_reference_once_then_abort";
struct NumericalEvent {
    std::string stage,operation,first_error,reference_error;
    bool recovered=false;
    std::uint64_t proposal_index=0;
    int node=-1,bond=-1;
    double u=0.;
    std::vector<double> values;
};
using NumericalObserver=std::function<void(const NumericalEvent&,const PhononPath&)>;
enum class PhononProposal { LocalConditional, BondRefresh, GuidedBond };
// Normalized free single-bond path density, evaluated independently via AR
// transitions of X-mu. Equals log Wph_b - log Zph_b, including shifted endpoints.
inline double freeBondLogDensity(const PhononPath &path,int bond) {
    const auto &p=path.parameters();const double variance=1/(2*path.alpha());
    const double r=path.decay(),noise=variance*path.noiseFactor();
    auto y=[&](int j){ return path.coordinate(j,bond)-path.freeMean(j,bond); };
    double result=ConditionalGaussian{0.,variance}.logDensity(y(0));
    for(int j=1;j<path.nodes();++j)
        result+=ConditionalGaussian{r*y(j-1),noise}.logDensity(y(j));
    if(!std::isfinite(result))throw std::runtime_error("nonfinite free bond proposal density");
    return result;
}
class PhononSampler {
public:
    PhononSampler(CoupledPhononPath &state,std::mt19937_64 rng,PhononProposal mode=PhononProposal::LocalConditional,
                  NumericalObserver observer={},QModeOptions qMode={},GuidedBondOptions guide={},
                  std::vector<QModeOptions> qPool={})
        : state_(state),rng_(std::move(rng)),mode_(mode),observer_(std::move(observer)),qMode_(std::move(qMode)),qPool_(std::move(qPool)) {
        if(mode!=PhononProposal::LocalConditional && mode!=PhononProposal::BondRefresh && mode!=PhononProposal::GuidedBond)
            throw std::invalid_argument("unknown phonon proposal mode");
        if(qMode_.every<0 || !std::isfinite(qMode_.beta) || qMode_.beta<=0 || qMode_.beta>1 ||
           (qMode_.every>0 && qMode_.basis.empty())) throw std::invalid_argument("invalid sampler q-mode options");
        const auto &path=state_.state().path;
        if(mode==PhononProposal::GuidedBond)guide_=std::make_unique<GuidedBondProposal>(path,std::move(guide));
        if(mode!=PhononProposal::GuidedBond && (!guide.forces.empty() || guide.beta!=1.))
            throw std::invalid_argument("guided options require GuidedBond proposal");
        // Orthonormality is required for the matter-only acceptance rule.
        for(std::size_t k=0;k<qMode_.basis.size();++k) {
            if(qMode_.basis[k].size()!=std::size_t(path.bonds()))throw std::invalid_argument("q-mode basis bond count");
            for(std::size_t l=0;l<=k;++l) {
                double dot=0.;for(int b=0;b<path.bonds();++b)dot+=qMode_.basis[k][b]*qMode_.basis[l][b];
                if(!std::isfinite(dot)||std::abs(dot-(k==l?1.:0.))>1e-10)throw std::invalid_argument("q-mode basis must be orthonormal");
            }
        }
        for(const auto &q:qPool_) {
            if(q.every!=qMode_.every || q.beta<=0 || q.beta>1 || !std::isfinite(q.beta) || q.basis.empty())
                throw std::invalid_argument("q pool needs fixed common cadence and nonempty modes");
            for(std::size_t k=0;k<q.basis.size();++k) {
                if(q.basis[k].size()!=std::size_t(path.bonds()))throw std::invalid_argument("q pool bond count");
                for(std::size_t l=0;l<=k;++l) {
                    double dot=0.;for(int b=0;b<path.bonds();++b)dot+=q.basis[k][b]*q.basis[l][b];
                    if(!std::isfinite(dot)||std::abs(dot-(k==l?1.:0.))>1e-10)throw std::invalid_argument("q pool basis not orthonormal");
                }
            }
        }
        qPoolCounts_.resize(qPool_.size());
        const long long size=static_cast<long long>(mode==PhononProposal::LocalConditional?path.nodes():1)*path.bonds();
        if(size>std::numeric_limits<int>::max()) throw std::invalid_argument("too many coordinates per sweep");
        order_.resize(size); std::iota(order_.begin(),order_.end(),0);
    }
    // Also permits deterministic decisions in the Fock/transaction test.
    MoveResult attempt(int node,int bond,double value,double u) {
        if(!(u>0 && u<1)) throw std::invalid_argument("Metropolis u must be strictly inside (0,1)");
        try {
            const auto &old=state_.state(); const auto q=freeConditional(old.path,node,bond);
            auto proposal=certified([&](bool reference){return state_.propose(node,bond,value,reference);},
                NumericalEvent{stage_,"local","","",false,counts_.proposed+1,node,bond,u,{value}});
            const double logRatio=proposal.candidate().weight.log_abs_matter-old.weight.log_abs_matter;
            const double general=proposal.logAbsRatio()+q.logDensity(old.path.coordinate(node,bond))-q.logDensity(value);
            if(!std::isfinite(logRatio) || !std::isfinite(general)) throw std::runtime_error("nonfinite acceptance ratio");
            const bool accept=std::log(u)<std::min(0.,logRatio);
            if(accept) state_.accept(proposal); else state_.reject(proposal);
            ++counts_.proposed; counts_.accepted+=accept; counts_.metropolis_rejected+=!accept;
            return {logRatio,general,accept};
        } catch(const std::exception &e) {
            throw std::runtime_error("MC numerical/proposal failure at node="+std::to_string(node)+
                " bond="+std::to_string(bond)+" completed_proposals="+std::to_string(counts_.proposed)+": "+e.what());
        }
    }
    MoveResult attemptBond(int bond,const std::vector<double> &values,double u) {
        if(!(u>0 && u<1))throw std::invalid_argument("Metropolis u must be strictly inside (0,1)");
        try {
            const auto &old=state_.state();
            auto proposal=profilePath("bond_refresh_candidate",[&] {
                return certified([&](bool reference){return state_.proposeBondPath(bond,values,reference);},
                    NumericalEvent{stage_,"bond-refresh","","",false,counts_.proposed+1,-1,bond,u,values});
            });
            const double logRatio=proposal.candidate().weight.log_abs_matter-old.weight.log_abs_matter;
            const double general=proposal.logAbsRatio()+freeBondLogDensity(old.path,bond)
                -freeBondLogDensity(proposal.candidate().path,bond);
            if(!std::isfinite(logRatio) || !std::isfinite(general))throw std::runtime_error("nonfinite bond acceptance ratio");
            const bool accept=std::log(u)<std::min(0.,logRatio);
            if(accept)state_.accept(proposal);else state_.reject(proposal);
            ++counts_.proposed;counts_.accepted+=accept; counts_.metropolis_rejected+=!accept;
            return {logRatio,general,accept};
        }catch(const std::exception &e) {
            throw std::runtime_error("MC bond-refresh failure at bond="+std::to_string(bond)+
                " completed_proposals="+std::to_string(counts_.proposed)+": "+e.what());
        }
    }
    MoveResult attemptGuidedBond(int bond,const std::vector<double> &values,double u) {
        if(!(u>0 && u<1) || !guide_)throw std::invalid_argument("invalid guided bond attempt");
        try {
            const auto &old=state_.state();
            auto proposal=profilePath("guided_bond_candidate",[&] {
                return certified([&](bool reference){return state_.proposeBondPath(bond,values,reference);},
                    NumericalEvent{stage_,"guided-bond","","",false,counts_.proposed+1,-1,bond,u,values});
            });
            const double ratio=proposal.candidate().weight.log_abs_matter-old.weight.log_abs_matter
                -guide_->referenceScoreDifference(old.path,bond,values);
            std::vector<double> previous(old.path.nodes());
            for(int j=0;j<old.path.nodes();++j)previous[j]=old.path.coordinate(j,bond);
            const double general=proposal.logAbsRatio()+guide_->logTransition(proposal.candidate().path,bond,previous)
                -guide_->logTransition(old.path,bond,values);
            if(!std::isfinite(ratio)||!std::isfinite(general))throw std::runtime_error("nonfinite guided MH ratio");
            const bool accept=std::log(u)<std::min(0.,ratio);
            if(accept)state_.accept(proposal);else state_.reject(proposal);
            ++counts_.proposed;counts_.accepted+=accept;counts_.metropolis_rejected+=!accept;
            return {ratio,general,accept};
        }catch(const std::exception &e) {
            throw std::runtime_error("guided bond failure at bond="+std::to_string(bond)+
                " completed_proposals="+std::to_string(counts_.proposed)+": "+e.what());
        }
    }
    MoveResult attemptQMode(const QModeAmplitudes &amplitudes,double u) {
        return attemptQMode(qMode_,amplitudes,u);
    }
    MoveResult attemptQMode(const QModeOptions &mode,const QModeAmplitudes &amplitudes,double u) {
        if(!(u>0 && u<1) || mode.basis.empty())throw std::invalid_argument("invalid q-mode attempt");
        const auto before=counts_;
        try {
            const auto &old=state_.state();
            const auto from=qModeAmplitudes(old.path,mode);
            const auto candidate=qModeCandidate(old.path,mode,amplitudes);
            std::vector<double> values;for(const auto &row:amplitudes)values.insert(values.end(),row.begin(),row.end());
            auto proposal=profilePath("q_mode_candidate",[&] {
                return certified([&](bool reference){return state_.proposeQModePath(candidate,reference);},
                    NumericalEvent{stage_,"q-mode","","",false,counts_.proposed+1,-1,-1,u,values});
            });
            const double ratio=proposal.candidate().weight.log_abs_matter-old.weight.log_abs_matter;
            const double general=proposal.logAbsRatio()+qModeTransitionLogDensity(old.path,amplitudes,from,mode.beta)
                -qModeTransitionLogDensity(old.path,from,amplitudes,mode.beta);
            if(!std::isfinite(ratio)||!std::isfinite(general))throw std::runtime_error("nonfinite q-mode acceptance ratio");
            const bool accept=std::log(u)<std::min(0.,ratio);
            if(accept)state_.accept(proposal);else state_.reject(proposal);
            ++counts_.proposed;counts_.accepted+=accept;counts_.metropolis_rejected+=!accept;
            addQCounts(counts_-before);
            return {ratio,general,accept};
        }catch(...) {addQCounts(counts_-before);throw;}
    }
    void sweep() {
        // Random permutation, covering each endpoint and the shared center once.
        std::shuffle(order_.begin(),order_.end(),rng_);
        const int bonds=state_.state().path.bonds();
        for(int index:order_) {
            if(mode_==PhononProposal::GuidedBond) {
                const auto values=guide_->draw(state_.state().path,index,rng_);
                attemptGuidedBond(index,values,(double(rng_()>>12)+.5)*0x1p-52);
                continue;
            }
            if(mode_==PhononProposal::BondRefresh) {
                const auto &p=state_.state().path;
                const auto drawn=PhononPath(1,p.intervals(),p.dt(),p.parameters(),{p.trialShift(index)},{p.force(index)}).freeSample(rng_);
                std::vector<double> values(p.nodes());for(int j=0;j<p.nodes();++j)values[j]=drawn.coordinate(j,0);
                attemptBond(index,values,(double(rng_()>>12)+.5)*0x1p-52);
                continue;
            }
            const int node=index/bonds,bond=index%bonds;
            const auto q=freeConditional(state_.state().path,node,bond);
            const double next=q.mean+std::sqrt(q.variance)*normal_(rng_);
            // 52 bits + half a bin gives representable values strictly in (0,1).
            const double u=(double(rng_()>>12)+.5)*0x1p-52;
            attempt(node,bond,next,u);
        }
        // Disabled mode consumes exactly zero additional RNG draws.
        if(qMode_.every>0 && ++qSweep_%qMode_.every==0) {
            // Fixed, state-independent uniform mixture: no adaptive peak feedback.
            std::size_t index=0;
            if(qPool_.size()>1)index=std::uniform_int_distribution<std::size_t>(0,qPool_.size()-1)(rng_);
            const auto &mode=qPool_.empty()?qMode_:qPool_[index];
            const auto &p=state_.state().path;
            auto a=qModeAmplitudes(p,mode);
            const double rho=std::sqrt((1-mode.beta)*(1+mode.beta));
            const double sigma=std::sqrt(1/(2*p.alpha())),noise=sigma*std::sqrt(p.noiseFactor());
            for(auto &row:a) {
                double xi=sigma*normal_(rng_);row[0]=rho*row[0]+mode.beta*xi;
                for(int j=1;j<p.nodes();++j){xi=p.decay()*xi+noise*normal_(rng_);row[j]=rho*row[j]+mode.beta*xi;}
            }
            const auto before=counts_;
            try {attemptQMode(mode,a,(double(rng_()>>12)+.5)*0x1p-52);}
            catch(...) {if(!qPool_.empty())addCounts(qPoolCounts_[index],counts_-before);throw;}
            if(!qPool_.empty())addCounts(qPoolCounts_[index],counts_-before);
        }
    }
    const std::vector<MoveCounts> &qPoolCounts()const{return qPoolCounts_;}
    const MoveCounts &counts() const { return counts_; }
    const MoveCounts &qCounts() const { return qCounts_; }
    const CoupledPhononPath::Snapshot &state() const { return state_.state(); }
    // Task B: no-op for eager engines; for deferred engines it fills the
    // signed phase and, when requested by an observable, center Green before
    // a measurement. Consumes no proposal RNG.
    void ensureMeasured(){
        certified([&](bool reference){state_.ensureMeasured(reference);return true;},
            NumericalEvent{stage_,"measurement","","",false,counts_.proposed,-1,-1,0.,{}});
    }
    void setStage(std::string stage){if(stage!=stage_)qSweep_=0;stage_=std::move(stage);}
private:
    template<class F> auto certified(F evaluate,NumericalEvent event)->decltype(evaluate(false)) {
        try { return evaluate(false); }
        catch(const std::runtime_error &e) { event.first_error=e.what(); }
        ++counts_.first_failures;
        // No random calls here. Neither failed attempt may commit a snapshot.
        // Keep observer I/O outside evaluation catches: a diagnostics write
        // failure must itself terminate, never trigger another computation.
        using Result=decltype(evaluate(false));
        std::unique_ptr<Result> recovered;
        try { recovered=std::make_unique<Result>(evaluate(true)); }
        catch(const std::exception &e) { event.reference_error=e.what(); }
        if(!recovered) {
            ++counts_.reference_failed;
            if(event.operation!="measurement") ++counts_.proposed;
            if(observer_) observer_(event,state_.state().path);
            throw std::runtime_error("uncertified "+event.operation+" in "+event.stage+
                "; first: "+event.first_error+"; uncached reference: "+event.reference_error);
        }
        ++counts_.reference_recovered; event.recovered=true;
        if(observer_) observer_(event,state_.state().path);
        return std::move(*recovered);
    }
    static void addCounts(MoveCounts &out,const MoveCounts &c) {
        out.proposed+=c.proposed;out.accepted+=c.accepted;out.untrusted_rejected+=c.untrusted_rejected;
        out.first_failures+=c.first_failures;out.reference_recovered+=c.reference_recovered;
        out.reference_failed+=c.reference_failed;out.metropolis_rejected+=c.metropolis_rejected;
    }
    void addQCounts(const MoveCounts &c) {addCounts(qCounts_,c);}
    NumericalObserver observer_;
    QModeOptions qMode_;
    std::unique_ptr<GuidedBondProposal> guide_;
    std::vector<QModeOptions> qPool_;
    std::vector<MoveCounts> qPoolCounts_;
    std::uint64_t qSweep_=0;
    MoveCounts qCounts_;
    std::string stage_="unspecified";
    CoupledPhononPath &state_;
    std::mt19937_64 rng_;
    PhononProposal mode_;
    std::normal_distribution<double> normal_;
    std::vector<int> order_;
    MoveCounts counts_;
};
struct SamplingSchedule {
    int warmup,measurement_sweeps,stride,block_samples;
    void validate() const {
        if(warmup<0 || measurement_sweeps<=0 || stride<=0 || block_samples<=0 ||
            measurement_sweeps%stride || (measurement_sweeps/stride)%block_samples)
            throw std::invalid_argument("measurement sweeps must contain whole strides and equal complete blocks");
    }
};
struct SamplingTimings {
    double warmup_seconds=0.;
    double measurement_seconds=0.;
    double measurement_materialization_seconds=0.;
    double measurement_callback_seconds=0.;
    double statistics_add_seconds=0.;
    double block_callback_seconds=0.;
};
template<class Measure,class OnBlock>
MoveCounts sample(PhononSampler &sampler,const SamplingSchedule &schedule,PhaseStatistics &statistics,
                  Measure measure,OnBlock onBlock,SamplingTimings *timings=nullptr,
                  const std::function<void(int,const CoupledPhononPath::Snapshot&)> &onWarmup={},
                  const std::function<void(const char*,int,const MoveCounts&)> &onSweep={}) {
    schedule.validate();
    if(statistics.blockSamples()!=schedule.block_samples)
        throw std::invalid_argument("schedule/statistics block sizes disagree");
    sampler.setStage("warmup");
    const auto warmupBegin=std::chrono::steady_clock::now();
    for(int i=0;i<schedule.warmup;++i) {
        sampler.sweep();
        if(onWarmup) {
            sampler.ensureMeasured();
            onWarmup(i+1,sampler.state());
        }
        if(onSweep) onSweep("warmup",i+1,sampler.counts());
    }
    if(timings) timings->warmup_seconds=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-warmupBegin).count();
    const auto warm=sampler.counts();
    sampler.setStage("measurement");
    const auto measurementBegin=std::chrono::steady_clock::now();
    for(int i=0;i<schedule.measurement_sweeps;++i) {
        sampler.sweep();
        if((i+1)%schedule.stride==0) {
            const auto materializationBegin=std::chrono::steady_clock::now();
            sampler.ensureMeasured();
            if(timings) timings->measurement_materialization_seconds+=
                std::chrono::duration<double>(std::chrono::steady_clock::now()-materializationBegin).count();
            const auto &s=sampler.state();
            const auto measureBegin=std::chrono::steady_clock::now();
            const auto values=measure(s);
            if(timings) timings->measurement_callback_seconds+=
                std::chrono::duration<double>(std::chrono::steady_clock::now()-measureBegin).count();
            const auto addBegin=std::chrono::steady_clock::now();
            const bool completed=statistics.add(s.weight.matter_phase,values);
            if(timings) timings->statistics_add_seconds+=
                std::chrono::duration<double>(std::chrono::steady_clock::now()-addBegin).count();
            if(completed) {
                const auto blockBegin=std::chrono::steady_clock::now();
                onBlock(statistics);
                if(timings) timings->block_callback_seconds+=
                    std::chrono::duration<double>(std::chrono::steady_clock::now()-blockBegin).count();
            }
        }
        if(onSweep) onSweep("measurement",i+1,sampler.counts()-warm);
    }
    if(timings) timings->measurement_seconds=
        std::chrono::duration<double>(std::chrono::steady_clock::now()-measurementBegin).count();
    return warm;
}
inline std::vector<std::string> basicObservableNames(bool measureEnergy=true) {
    std::vector<std::string> result={"X_x","X_y","X_z","X2_x","X2_y","X2_z"};
    if(measureEnergy) result.push_back("majorana_energy");
    return result;
}
inline std::vector<Complex> basicObservables(const KitaevSSHModel &model,
                                             const CoupledPhononPath::Snapshot &s,
                                             bool measureEnergy=true) {
    std::vector<Complex> result(measureEnergy?7:6,0.);
    const auto &x=s.path.node(s.path.center());
    for(const auto &b:model.lattice().bonds()) {
        const int mu=static_cast<int>(b.type); result[mu]+=x[b.index]/model.lattice().cells();
        result[3+mu]+=x[b.index]*x[b.index]/model.lattice().cells();
    }
    if(measureEnergy)
        result[6]=matterEnergyAtDisplacement(model,x,s.weight.center_green); // not total phonon+matter energy
    return result;
}
}
#endif
