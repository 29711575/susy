#ifndef KITAEV_SSH_PHASE_STATISTICS_H
#define KITAEV_SSH_PHASE_STATISTICS_H

#include <complex>
#include <cmath>
#include <cstdint>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "checkpoint_io.h"

namespace kitaev_ssh {
using Complex=std::complex<double>;
inline bool finiteComplex(Complex z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }
// Classification only: retain the original complex phase in all sums.
inline constexpr double phaseRealityTolerance=1e-10;
struct PhaseBlock {
    std::uint64_t count=0;
    Complex phase=0.;
    std::vector<Complex> numerator;
    std::uint64_t positive=0,negative=0,nonreal=0;
};
struct ComplexEstimate {
    Complex value=0.;
    double error_real=std::numeric_limits<double>::quiet_NaN();
    double error_imag=std::numeric_limits<double>::quiet_NaN();
    double covariance=std::numeric_limits<double>::quiet_NaN();
};
struct PhaseSummary {
    std::uint64_t samples=0;
    std::size_t blocks=0;
    ComplexEstimate phase;
    std::uint64_t positive=0,negative=0,nonreal=0;
    double phase_magnitude_error=std::numeric_limits<double>::quiet_NaN();
    std::vector<ComplexEstimate> observables;
    bool ratio_defined=false, jackknife_defined=false, phase_resolved=false;
    bool too_few_blocks=true;
};

inline void jackknifeErrors(ComplexEstimate &out,const std::vector<Complex> &leaveOne) {
    if(leaveOne.size()<2) return;
    Complex mean=0.; for(auto z:leaveOne) mean+=z/double(leaveOne.size());
    double vr=0.,vi=0.,cov=0.;
    for(auto z:leaveOne) { z-=mean; vr+=z.real()*z.real(); vi+=z.imag()*z.imag(); cov+=z.real()*z.imag(); }
    const double factor=(leaveOne.size()-1.)/leaveOne.size();
    if(!std::isfinite(vr) || !std::isfinite(vi) || !std::isfinite(cov))
        throw std::overflow_error("nonfinite jackknife variance");
    out.error_real=std::sqrt(factor*vr); out.error_imag=std::sqrt(factor*vi); out.covariance=factor*cov;
}

inline PhaseSummary summarizeBlocks(const std::vector<PhaseBlock> &blocks) {
    PhaseSummary result; result.blocks=blocks.size();
    if(blocks.empty()) return result;
    const auto width=blocks[0].numerator.size();
    PhaseBlock total; total.numerator.resize(width,0.);
    for(const auto &b:blocks) {
        if(!b.count || b.count!=blocks[0].count || b.numerator.size()!=width)
            throw std::invalid_argument("statistics require equal nonempty blocks");
        total.count+=b.count; total.phase+=b.phase;
        result.positive+=b.positive; result.negative+=b.negative; result.nonreal+=b.nonreal;
        for(std::size_t j=0;j<width;++j) total.numerator[j]+=b.numerator[j];
    }
    if(!finiteComplex(total.phase)) throw std::overflow_error("nonfinite phase sum");
    for(auto z:total.numerator) if(!finiteComplex(z)) throw std::overflow_error("nonfinite numerator sum");
    result.samples=total.count; result.phase.value=total.phase/double(total.count);
    result.too_few_blocks=blocks.size()<16;
    result.observables.resize(width);
    result.ratio_defined=std::abs(total.phase)>1e-12*total.count;
    if(result.ratio_defined) for(std::size_t j=0;j<width;++j) {
        result.observables[j].value=total.numerator[j]/total.phase;
        if(!finiteComplex(result.observables[j].value)) throw std::overflow_error("nonfinite reweighted estimate");
    }
    if(blocks.size()<2) return result;
    std::vector<Complex> phaseReplicas;
    for(const auto &b:blocks) phaseReplicas.push_back((total.phase-b.phase)/double(total.count-b.count));
    jackknifeErrors(result.phase,phaseReplicas);
    std::vector<Complex> magnitudeReplicas;
    for(auto z:phaseReplicas) magnitudeReplicas.emplace_back(std::abs(z),0.);
    ComplexEstimate magnitude; jackknifeErrors(magnitude,magnitudeReplicas);
    result.phase_magnitude_error=magnitude.error_real;
    result.phase_resolved=result.ratio_defined && std::abs(result.phase.value)>
        5*std::hypot(result.phase.error_real,result.phase.error_imag);
    if(!result.ratio_defined) return result;
    for(const auto &b:blocks) if(std::abs(total.phase-b.phase)<=1e-12*(total.count-b.count)) return result;
    std::vector<Complex> replicas(blocks.size());
    for(std::size_t j=0;j<width;++j) {
        for(std::size_t i=0;i<blocks.size();++i) {
            const auto &b=blocks[i];
            replicas[i]=(total.numerator[j]-b.numerator[j])/(total.phase-b.phase);
        }
        jackknifeErrors(result.observables[j],replicas);
    }
    result.jackknife_defined=true;
    return result;
}

// Transform ensemble ratios, NOT individual measurements. Every nonlinear
// operation is repeated with the SAME block removed from all s and sO sums.
// The point estimate stays the total-sum ratio (no jackknife bias correction).
template<class Transform>
PhaseSummary summarizeDerivedBlocks(const std::vector<PhaseBlock> &blocks,Transform transform) {
    auto result=summarizeBlocks(blocks);
    result.observables.clear();
    if(!result.ratio_defined) return result;
    Complex phase=0.; std::vector<Complex> sums(blocks.front().numerator.size(),0.);
    for(const auto &b:blocks) { phase+=b.phase; for(std::size_t j=0;j<sums.size();++j)sums[j]+=b.numerator[j]; }
    auto ratios=sums; for(auto &z:ratios)z/=phase;
    const auto values=transform(ratios); result.observables.resize(values.size());
    for(std::size_t j=0;j<values.size();++j) {
        if(!finiteComplex(values[j])) throw std::runtime_error("nonfinite derived observable");
        result.observables[j].value=values[j];
    }
    if(!result.jackknife_defined) return result;
    std::vector<std::vector<Complex>> replicas(values.size());
    for(const auto &b:blocks) {
        for(std::size_t j=0;j<sums.size();++j)ratios[j]=(sums[j]-b.numerator[j])/(phase-b.phase);
        const auto v=transform(ratios);
        if(v.size()!=values.size()) throw std::logic_error("derived observable dimension changed");
        for(std::size_t j=0;j<v.size();++j) replicas[j].push_back(v[j]);
    }
    for(std::size_t j=0;j<values.size();++j) jackknifeErrors(result.observables[j],replicas[j]);
    return result;
}

inline std::vector<PhaseBlock> mergeAdjacentBlocks(const std::vector<PhaseBlock> &blocks) {
    if(blocks.empty() || blocks.size()%2) return {}; // never discard a remainder
    std::vector<PhaseBlock> merged;
    for(std::size_t i=0;i<blocks.size();i+=2) {
        auto b=blocks[i]; const auto &c=blocks[i+1];
        if(b.count!=c.count || b.numerator.size()!=c.numerator.size())
            throw std::invalid_argument("unequal blocks");
        b.count+=c.count; b.phase+=c.phase;
        b.positive+=c.positive; b.negative+=c.negative; b.nonreal+=c.nonreal;
        for(std::size_t j=0;j<b.numerator.size();++j)b.numerator[j]+=c.numerator[j];
        merged.push_back(std::move(b));
    }
    return merged;
}

class PhaseStatistics {
public:
    PhaseStatistics(std::vector<std::string> names,int blockSamples)
        : names_(std::move(names)),block_samples_(blockSamples) {
        if(names_.empty() || blockSamples<=0) throw std::invalid_argument("invalid observable/block count");
        pending_.numerator.resize(names_.size(),0.);
    }
    // Called at fixed measurement times, whether the preceding moves accepted or not.
    bool add(Complex phase,const std::vector<Complex> &observables) {
        if(!finiteComplex(phase) || std::abs(std::abs(phase)-1.)>1e-10 || observables.size()!=names_.size())
            throw std::invalid_argument("invalid phase or observable dimension");
        for(auto z:observables) if(!finiteComplex(z)) throw std::invalid_argument("nonfinite observable");
        if(!finiteComplex(pending_.phase+phase)) throw std::overflow_error("phase accumulation overflow");
        for(std::size_t j=0;j<names_.size();++j)
            if(!finiteComplex(pending_.numerator[j]+phase*observables[j]))
                throw std::overflow_error("numerator accumulation overflow");
        // Mutually exclusive measured-configuration counts, not proposal counts.
        if(std::abs(phase.imag())>phaseRealityTolerance) ++pending_.nonreal;
        else if(phase.real()<0.) ++pending_.negative;
        else ++pending_.positive;
        pending_.phase+=phase; ++pending_.count;
        for(std::size_t j=0;j<names_.size();++j) pending_.numerator[j]+=phase*observables[j];
        if(pending_.count!=std::uint64_t(block_samples_)) return false;
        blocks_.push_back(pending_); pending_={0,0.,std::vector<Complex>(names_.size(),0.)}; return true;
    }
    const std::vector<PhaseBlock> &blocks() const { return blocks_; }
    int blockSamples() const { return block_samples_; }
    const std::vector<std::string> &names() const { return names_; }
    std::uint64_t samplesIncludingPending() const {
        return std::uint64_t(blocks_.size())*block_samples_+pending_.count;
    }
    void saveCheckpoint(std::ostream &out) const {
        out<<std::setprecision(17)<<block_samples_<<' '<<names_.size()<<' '<<blocks_.size()<<'\n';
        for(const auto &name:names_)out<<std::quoted(name)<<'\n';
        auto write=[&](const PhaseBlock &b) {
            out<<b.count<<' '<<b.phase.real()<<' '<<b.phase.imag()<<' '
               <<b.positive<<' '<<b.negative<<' '<<b.nonreal<<'\n';
            for(auto z:b.numerator)out<<z.real()<<' '<<z.imag()<<'\n';
        };
        for(const auto &b:blocks_)write(b);write(pending_);
        if(!out)throw std::runtime_error("statistics checkpoint write failed");
    }
    void loadCheckpoint(std::istream &in) {
        int bs;std::size_t width,count;
        checkpointRequire(bool(in>>bs>>width>>count) && bs==block_samples_ && width==names_.size(),"statistics schema mismatch");
        for(const auto &name:names_){std::string got;checkpointRequire(bool(in>>std::quoted(got)) && got==name,"observable names mismatch");}
        auto read=[&](bool pending) {
            PhaseBlock b;double re,im;
            checkpointRequire(bool(in>>b.count>>re>>im>>b.positive>>b.negative>>b.nonreal),"truncated phase block");
            b.phase={re,im};
            checkpointRequire(finiteComplex(b.phase) && (pending?b.count<std::uint64_t(bs):b.count==std::uint64_t(bs)) &&
                b.positive+b.negative+b.nonreal==b.count && std::abs(b.phase)<=double(b.count)*(1+1e-10),"invalid phase block");
            b.numerator.resize(width);
            for(auto &z:b.numerator){checkpointRequire(bool(in>>re>>im),"truncated observable sums");z={re,im};checkpointRequire(finiteComplex(z),"nonfinite observable sum");}
            return b;
        };
        // Bound allocation by the actual text remaining, including a minimal block header.
        const auto here=in.tellg();in.seekg(0,std::ios::end);const auto end=in.tellg();in.seekg(here);
        checkpointRequire(here>=0 && end>=here && count<=std::uint64_t(end-here)/12,"invalid block count");
        std::vector<PhaseBlock> blocks;for(std::size_t i=0;i<count;++i)blocks.push_back(read(false));
        auto pending=read(true);checkpointEnd(in);blocks_=std::move(blocks);pending_=std::move(pending);
    }
    void writeAllBlocks(std::ostream &out) const {
        for(std::size_t i=0;i<blocks_.size();++i)writeBlock(out,i);
    }
    PhaseSummary summary() const {
        if(pending_.count) throw std::logic_error("incomplete final block: no samples silently discarded");
        return summarizeBlocks(blocks_);
    }
    // Pair adjacent equal blocks. A comparison is a diagnostic, not proof of independence.
    PhaseSummary doubledBlockSummary() const {
        if(pending_.count || blocks_.size()%2) return {};
        return summarizeBlocks(mergeAdjacentBlocks(blocks_));
    }
    void writeHeader(std::ostream &out) const {
        out<<"block,count,sum_phase_real,sum_phase_imag,positive_count,negative_count,complex_count";
        for(const auto &name:names_) out<<",sum_phase_"<<name<<"_real,sum_phase_"<<name<<"_imag";
        out<<'\n';
    }
    void writeLastBlock(std::ostream &out) const {
        writeBlock(out,blocks_.size()-1);out.flush();
        if(!out) throw std::runtime_error("block statistics write failed");
    }
private:
    void writeBlock(std::ostream &out,std::size_t index) const {
        const auto &b=blocks_.at(index); out<<index<<','<<b.count<<','<<b.phase.real()<<','<<b.phase.imag();
        out<<','<<b.positive<<','<<b.negative<<','<<b.nonreal;
        for(auto z:b.numerator) out<<','<<z.real()<<','<<z.imag();
        out<<'\n';
        if(!out) throw std::runtime_error("block statistics write failed");
    }
    std::vector<std::string> names_;
    int block_samples_;
    PhaseBlock pending_;
    std::vector<PhaseBlock> blocks_;
};
}
#endif
