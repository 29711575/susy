#ifndef KITAEV_SSH_MC_CHECKPOINT_H
#define KITAEV_SSH_MC_CHECKPOINT_H

#include "phonon_mc.h"
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <fstream>
#include <map>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace kitaev_ssh {
// Signal handlers do no I/O, allocation, sampling, or numerical work.
inline volatile std::sig_atomic_t checkpointSaveRequested=0,checkpointStopRequested=0;
inline void checkpointSignalHandler(int signal) {
    checkpointSaveRequested=1;
    if(signal==SIGUSR2 || signal==SIGTERM)checkpointStopRequested=1;
}
class CheckpointSignals {
public:
    CheckpointSignals() {
        sigemptyset(&mask_);for(int s:{SIGUSR1,SIGUSR2,SIGTERM})sigaddset(&mask_,s);
        struct sigaction action{};action.sa_handler=checkpointSignalHandler;action.sa_mask=mask_;action.sa_flags=SA_RESTART;
        checkpointSaveRequested=checkpointStopRequested=0;
        for(int i=0;i<3;++i) {
            if(sigaction(signals_[i],&action,&old_[i])!=0) {
                while(i--)sigaction(signals_[i],&old_[i],nullptr);
                throw std::runtime_error("install checkpoint signal handler failed");
            }
        }
    }
    ~CheckpointSignals(){for(int i=0;i<3;++i)sigaction(signals_[i],&old_[i],nullptr);}
    std::pair<bool,bool> consume() {
        sigset_t before; if(sigprocmask(SIG_BLOCK,&mask_,&before)!=0)throw std::runtime_error("block checkpoint signals failed");
        const bool save=checkpointSaveRequested,stop=checkpointStopRequested;
        checkpointSaveRequested=checkpointStopRequested=0;
        if(sigprocmask(SIG_SETMASK,&before,nullptr)!=0)throw std::runtime_error("restore signal mask failed");
        return {save,stop};
    }
private:
    int signals_[3]={SIGUSR1,SIGUSR2,SIGTERM};struct sigaction old_[3]{};sigset_t mask_{};
};
struct CheckpointStopped:std::exception {
    const char *what()const noexcept override{return "checkpoint saved; requested stop";}
};
inline std::string checkpointReadFile(const std::string &name) {
    std::ifstream in(name,std::ios::binary);if(!in)throw std::invalid_argument("cannot open checkpoint/input: "+name);
    std::ostringstream out;out<<in.rdbuf();if(in.bad())throw std::runtime_error("checkpoint/input read failed: "+name);return out.str();
}
inline std::string checkpointExecutableIdentity() {
    // Read executable bytes, not the path, so copying the identical binary is allowed.
    std::ifstream in("/proc/self/exe",std::ios::binary);if(!in)throw std::runtime_error("cannot identify checkpoint executable");
    std::uint64_t hash=14695981039346656037ull,size=0;char buf[65536];
    while(in){in.read(buf,sizeof(buf));const auto n=in.gcount();size+=n;for(std::streamsize i=0;i<n;++i){hash^=static_cast<unsigned char>(buf[i]);hash*=1099511628211ull;}}
    if(!in.eof())throw std::runtime_error("executable fingerprint read failed");
    std::ostringstream s;s<<size<<':'<<std::hex<<hash;return s.str();
}
inline void checkpointAtomicWrite(const std::string &name,const std::string &data) {
    const std::string temp=name+".tmp."+std::to_string(getpid());
    int fd=open(temp.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);
    if(fd<0)throw std::runtime_error("cannot create checkpoint temporary file: "+temp);
    bool renamed=false;
    try {
        std::size_t done=0;
        while(done<data.size()) {
            const auto n=write(fd,data.data()+done,data.size()-done);
            if(n<0 && errno==EINTR)continue;
            if(n<=0)throw std::runtime_error("checkpoint write failed");done+=n;
        }
        if(fsync(fd)!=0)throw std::runtime_error("checkpoint fsync failed");
        if(close(fd)!=0){fd=-1;throw std::runtime_error("checkpoint close failed");}fd=-1;
        if(rename(temp.c_str(),name.c_str())!=0)throw std::runtime_error("checkpoint rename failed");renamed=true;
        const auto slash=name.find_last_of('/');const auto dir=slash==std::string::npos?".":(slash==0?"/":name.substr(0,slash));
        const int d=open(dir.c_str(),O_RDONLY|O_DIRECTORY);
        if(d<0)throw std::runtime_error("checkpoint directory open failed");
        const int rc=fsync(d);const int error=errno;close(d);
        if(rc!=0 && error!=EINVAL && error!=EROFS)throw std::runtime_error("checkpoint directory fsync failed");
    }catch(...){if(fd>=0)close(fd);if(!renamed)unlink(temp.c_str());throw;}
}
inline void checkpointWriteMap(std::ostream &out,const std::map<std::string,std::string> &items) {
    out<<items.size()<<'\n';for(const auto &p:items){checkpointBlob(out,p.first);checkpointBlob(out,p.second);}
}
inline std::map<std::string,std::string> checkpointReadMap(std::istream &in) {
    std::size_t n;checkpointRequire(bool(in>>n) && n<10000,"invalid map length");
    checkpointRequire(in.get()=='\n',"invalid map separator");
    std::map<std::string,std::string> m;
    while(n--){auto k=checkpointBlob(in);auto v=checkpointBlob(in);checkpointRequire(m.emplace(k,v).second,"duplicate map key");}return m;
}
struct MCCheckpoint {
    std::string executable,origin,path,sampler,statistics,runner;
    std::map<std::string,std::string> args,inputs;
    SamplingProgress progress;
    std::uint64_t generation=0;
    std::string payload()const {
        std::ostringstream out;out.imbue(std::locale::classic());
        checkpointBlob(out,executable);checkpointBlob(out,origin);
        checkpointWriteMap(out,args);checkpointWriteMap(out,inputs);
        out<<generation<<' '<<progress.warmup_done<<' '<<progress.measurement_done<<'\n';
        PhononSampler::writeCounts(out,progress.warm_counts);
        for(const auto &s:{path,sampler,statistics,runner})checkpointBlob(out,s);
        return out.str();
    }
    void save(const std::string &file)const {
        const auto body=payload();std::ostringstream out;
        out<<"KITAEV_SSH_MC_CHECKPOINT 1\n"<<body.size()<<' '<<checkpointChecksum(body)<<'\n'<<body;
        checkpointAtomicWrite(file,out.str());
    }
    static MCCheckpoint load(const std::string &file) {
        const auto data=checkpointReadFile(file);std::istringstream in(data);in.imbue(std::locale::classic());
        std::string magic;int version;std::uint64_t size,hash;
        checkpointRequire(bool(in>>magic>>version>>size>>hash) && magic=="KITAEV_SSH_MC_CHECKPOINT" && version==1,"unknown format/version");
        checkpointRequire(in.get()=='\n',"invalid payload separator");
        const auto at=in.tellg();checkpointRequire(at>=0 && size==data.size()-std::size_t(at),"truncated/trailing payload");
        const auto body=data.substr(std::size_t(at));checkpointRequire(checkpointChecksum(body)==hash,"checksum mismatch");
        std::istringstream s(body);s.imbue(std::locale::classic());MCCheckpoint c;
        c.executable=checkpointBlob(s);c.origin=checkpointBlob(s);c.args=checkpointReadMap(s);c.inputs=checkpointReadMap(s);
        checkpointRequire(bool(s>>c.generation>>c.progress.warmup_done>>c.progress.measurement_done),"missing progress");
        c.progress.warm_counts=PhononSampler::readCounts(s);
        c.path=checkpointBlob(s);c.sampler=checkpointBlob(s);c.statistics=checkpointBlob(s);c.runner=checkpointBlob(s);checkpointEnd(s);return c;
    }
};
inline bool checkpointOperationalArgument(const std::string &k) {
    return k=="output"||k=="resume"||k=="checkpoint-file"||k=="checkpoint-every";
}
inline void checkpointMergeArguments(std::map<std::string,std::string> &args,const MCCheckpoint &c) {
    for(const auto &p:args)if(!checkpointOperationalArgument(p.first) && p.first!="sweeps") {
        const auto old=c.args.find(p.first);
        checkpointRequire(old!=c.args.end() && old->second==p.second,"changed physical, sampling, measurement or numerical option");
    }
    for(const auto &p:c.args)if(!checkpointOperationalArgument(p.first) && !args.count(p.first))args.emplace(p);
    checkpointRequire(args.count("output"),"resume requires a fresh --output prefix");
    checkpointRequire(args.at("output")!=c.origin,"resume output must differ from original prefix");
    if(args.count("sweeps"))checkpointRequire(std::stoll(args.at("sweeps"))>=c.progress.measurement_done,"sweeps below completed progress");
}
}
#endif
