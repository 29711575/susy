#ifndef PURE_PROJECTOR_FAILURE_SNAPSHOT_H
#define PURE_PROJECTOR_FAILURE_SNAPSHOT_H
#include "pure_projector_fast.h"
#include <fstream>
#include <random>
#include <cstdio>

// Versioned native little-endian diagnostic snapshot (x86). No restart of
// sampler caches is implied. Exact double bits, including signed zeros, survive.
namespace pure_failure_snapshot {
inline void require(bool ok,const char *message){if(!ok)throw std::runtime_error(message);}
template<class T> void put(std::ostream &out,const T &value){
    out.write(reinterpret_cast<const char *>(&value),sizeof(value));}
template<class T> T take(std::istream &in){T value{};
    in.read(reinterpret_cast<char *>(&value),sizeof(value));require(bool(in),"truncated failure snapshot");return value;}
inline void matrix(std::ostream &out,const MatType &m){
    put(out,int(m.rows()));put(out,int(m.cols()));
    const DataType base=m.size()?m(0,m.cols()-1):DataType(0,0);put(out,base);
    int count=0;for(int k=0;k<m.size();++k)if(std::memcmp(&m.data()[k],&base,sizeof(base)))++count;
    put(out,count);for(int k=0;k<m.size();++k)if(std::memcmp(&m.data()[k],&base,sizeof(base))){put(out,k);put(out,m.data()[k]);}
}
inline MatType matrix(std::istream &in){
    const int r=take<int>(in),c=take<int>(in);require(r>0&&r<=4096&&c>0&&c<=4096,"invalid failure matrix shape");
    MatType m=MatType::Constant(r,c,take<DataType>(in));
    const int n=take<int>(in);require(n>=0&&n<=m.size(),"invalid failure matrix count");int prior=-1;
    for(int k=0;k<n;++k){const int at=take<int>(in);require(at>prior&&at<m.size(),"invalid failure matrix index");
        m.data()[at]=take<DataType>(in);prior=at;}
    return m;
}
struct Snapshot {
    MatType phi;
    PureFastConfiguration configuration;
    PureFastProposal proposal;
    std::uint64_t hash=0,proposal_hash=0;
    int sign=0,block=0,ratio_status=0;
    double log_abs=0;
    DataType observed_ratio=0;
    std::string rng;
};
inline void write(const std::string &path,const GaussianTrialState &trial,
                  const PureProjectorFastWalker &walker,const PureFastProposal &proposal,
                  const PureFastProposalResult &result,int block,const std::mt19937_64 &rng){
    require(!std::ifstream(path).good(),"refusing to overwrite failure snapshot");
    std::ofstream out(path+".tmp",std::ios::binary);require(bool(out),"cannot open failure snapshot");
    put(out,std::uint64_t(0x5046514641494c01ULL));put(out,walker.configurationHash());
    put(out,result.snapshot.configuration_hash);put(out,walker.z2Sign());put(out,block);
    put(out,walker.currentWeight().log_abs_weight);put(out,int(result.ratio.status));put(out,result.snapshot.fast_ratio);
    matrix(out,trial.Phi);const auto &cfg=walker.configuration();put(out,int(cfg.slices.size()));
    for(int k=0;k<int(cfg.slices.size());++k){
        put(out,cfg.hs_fields[k]);const auto &l=cfg.locations[k];put(out,int(l.branch));
        put(out,l.slice);put(out,l.factor);put(out,l.bond);put(out,l.aux);
        put(out,cfg.slices[k].eta);matrix(out,MatType(cfg.slices[k].matrix));
    }
    put(out,proposal.index);put(out,proposal.new_hs);put(out,proposal.new_eta);put(out,proposal.uniform);
    put(out,proposal.positive_ratio_log_bound);matrix(out,proposal.materialized().new_factor);
    std::ostringstream state;state<<rng;const auto value=state.str();put(out,int(value.size()));out.write(value.data(),value.size());
    out.close();require(bool(out),"failure snapshot write failed");
    require(std::rename((path+".tmp").c_str(),path.c_str())==0,"failure snapshot rename failed");
}
inline Snapshot read(const std::string &path){
    std::ifstream in(path,std::ios::binary);Snapshot s;
    require(take<std::uint64_t>(in)==0x5046514641494c01ULL,"failure snapshot version mismatch");
    s.hash=take<std::uint64_t>(in);s.proposal_hash=take<std::uint64_t>(in);s.sign=take<int>(in);s.block=take<int>(in);
    s.log_abs=take<double>(in);s.ratio_status=take<int>(in);s.observed_ratio=take<DataType>(in);s.phi=matrix(in);
    const int count=take<int>(in);require(count>0&&count<1000000,"invalid failure factor count");
    for(int k=0;k<count;++k){
        s.configuration.hs_fields.push_back(take<int>(in));PureSliceLocation l;l.branch=PureBranch(take<int>(in));
        l.slice=take<int>(in);l.factor=take<int>(in);l.bond=take<int>(in);l.aux=take<int>(in);s.configuration.locations.push_back(l);
        const auto eta=take<DataType>(in);const auto factor=matrix(in);
        require(factor.rows()==s.phi.rows()&&factor.cols()==s.phi.rows(),"failure factor shape mismatch");
        s.configuration.slices.emplace_back(factor,eta,"failure_snapshot");
    }
    s.proposal.index=take<int>(in);s.proposal.new_hs=take<int>(in);s.proposal.new_eta=take<DataType>(in);
    s.proposal.uniform=take<double>(in);s.proposal.positive_ratio_log_bound=take<double>(in);s.proposal.new_factor=matrix(in);
    const int n=take<int>(in);require(n>0&&n<100000,"invalid failure RNG length");s.rng.resize(n);in.read(&s.rng[0],n);
    require(bool(in)&&in.peek()==std::char_traits<char>::eof(),"invalid failure snapshot trailer");
    require(s.hash==pureFastConfigurationHash(s.configuration),"failure snapshot hash mismatch");
    require(s.proposal.index>=0&&s.proposal.index<count,"invalid failure proposal index");
    return s;
}
}
#endif
