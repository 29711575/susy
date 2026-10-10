#define KITAEV_SSH_TEST_FAILURE_INJECTION
#include "phonon_mc.h"
#include <iostream>
#include <sstream>
using namespace kitaev_ssh;
namespace {
void demand(bool ok,const char *s){if(!ok)throw std::runtime_error(s);}
void close(double a,double b,double tol,const char *s){demand(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<tol,s);}
void equal(const CoupledPhononPath::Snapshot &a,const CoupledPhononPath::Snapshot &b,bool exact=false) {
    for(int j=0;j<a.path.nodes();++j)for(int k=0;k<a.path.bonds();++k)
        demand(a.path.coordinate(j,k)==b.path.coordinate(j,k),"coordinates changed");
    close(a.weight.logAbsTotal(),b.weight.logAbsTotal(),1e-9,"weight mismatch");
    if(a.weight.phase_valid&&b.weight.phase_valid)close(std::abs(a.weight.matter_phase-b.weight.matter_phase),0.,1e-9,"phase mismatch");
    if(a.weight.green_valid&&b.weight.green_valid)close((a.weight.center_green-b.weight.center_green).norm(),0.,1e-9,"Green mismatch");
    demand(a.configuration_hash==b.configuration_hash,"configuration hash mismatch");
    if(exact) {
        demand(a.signed_prefix_cache==b.signed_prefix_cache&&a.lu_cayley_cache==b.lu_cayley_cache,"rejection cache pointer changed");
        for(std::size_t i=0;i<a.factors.size();++i)close((MatType(a.factors[i].matrix)-MatType(b.factors[i].matrix)).norm(),0.,1e-30,"rejection factors changed");
    }
}
void basisAndBalance() {
    std::mt19937_64 rng(123);
    for(const auto &lat:{HoneycombLattice(9,9),HoneycombLattice(4,4),HoneycombLattice(4,4,2)}) {
        const int n1=lat.Lx()==9?5:2,n2=lat.Lx()==9?4:(lat.shear()?1:0);
        const auto q=makeQModeOptions(lat,n1,n2,1,.7);
        demand(q.basis.size()==(lat.Lx()==9?2:1),"conjugate rank");
        for(std::size_t k=0;k<q.basis.size();++k)for(std::size_t l=0;l<q.basis.size();++l) {
            double dot=0.,sum=0.;for(const auto &b:lat.bonds()) {
                dot+=q.basis[k][b.index]*q.basis[l][b.index];sum+=q.basis[k][b.index];
                if(b.type!=BondType::Z)demand(q.basis[k][b.index]==0.,"non-z basis");
            }
            close(dot,k==l?1.:0.,1e-12,"basis orthonormality");close(sum,0.,1e-12,"zero mode contamination");
        }
        std::vector<double>d(lat.bonds().size()),f(d.size());
        for(std::size_t b=0;b<d.size();++b){d[b]=.2*std::sin(b);f[b]=.3*std::cos(b);}
        const auto p=PhononPath(d.size(),8,.13,{1.7,.8},d,f).freeSample(rng);
        const auto from=qModeAmplitudes(p,q);
        for(double beta:{.2,.7,1.}) {
            auto a=from;for(auto &row:a)for(auto &v:row)v+=.15*std::normal_distribution<double>{}(rng);
            const auto next=qModeCandidate(p,q,a);
            const double balance=next.logWeight()-p.logWeight()+qModeTransitionLogDensity(p,a,from,beta)-qModeTransitionLogDensity(p,from,a,beta);
            close(balance,0.,2e-10,"physical Mehler/proposal cancellation");
            for(int j=0;j<p.nodes();++j)for(const auto &b:lat.bonds())if(b.type!=BondType::Z)
                demand(p.coordinate(j,b.index)==next.coordinate(j,b.index),"non-z coordinate changed");
        }
    }
    // Non-self-conjugate sheared character as well.
    const auto q=makeQModeOptions(HoneycombLattice(5,4,2),2,1,1,1.);
    double dot=0.;for(std::size_t b=0;b<q.basis[0].size();++b)dot+=q.basis[0][b]*q.basis[1][b];
    close(dot,0.,1e-12,"sheared pair");
    bool failed=false;try{makeQModeOptions(HoneycombLattice(4,4),0,0,1,1.);}catch(const std::invalid_argument&){failed=true;}
    demand(failed,"enabled zero mode allowed");
}
// Independent two-state Fock oracle for a general linear bond coupling.
MajoranaGenerator generator(double g,const std::vector<double> &c) {
    return [g,c](const std::vector<double>&x){double t=.5;for(std::size_t b=0;b<c.size();++b)t+=g*c[b]*x[b];
        MatType k=MatType::Zero(2,2);k(0,1)=DataType(0,-2*t);k(1,0)=-k(0,1);return k;};
}
double fock(const PhononPath &p,double g,const std::vector<double>&c) {
    double log=0.;for(int j=1;j<p.nodes();++j){double t=1.;for(std::size_t b=0;b<c.size();++b)t+=g*c[b]*(p.coordinate(j,b)+p.coordinate(j-1,b));log+=p.dt()/2*t;}return log;
}
void productionOptions() {
    KitaevSSHModel model(HoneycombLattice(3,3,1),p2Couplings(.7));
    const auto &lat=model.lattice();const auto q=makeQModeOptions(lat,2,1,1,.65);
    demand(q.basis.size()==2,"production conjugate pair");
    std::vector<double> shifts(lat.bondCount()),forces(lat.bondCount());
    for(const auto &b:lat.bonds()) {
        shifts[b.index]=.13*std::cos(.7*b.index);
        forces[b.index]=.07*std::sin(.3*b.index)+.09*q.basis[0][b.index];
    }
    std::mt19937_64 rng(8291);
    const auto path=PhononPath(lat.bondCount(),6,.13,{1.7,.8},shifts,forces).freeSample(rng);
    const auto gen=honeycombGenerator(model);const auto trial=freeGroundState(model);
    const MatType k0=gen(path.node(0)),k1=gen(path.node(1));
    demand((k0*k1-k1*k0).norm()>1e-5,"production generators must not commute");
    PhononWeightOptions options;options.singular_values_only=true;options.merge_node_half_steps=true;
    options.defer_phase_measurement=true;options.reuse_signed_prefix_checkpoints=true;options.reuse_lu_cayley=true;
    CoupledPhononPath normal(path,trial,gen,4,options),recovered(path,trial,gen,4,options);
    std::vector<NumericalEvent> events;
    PhononSampler a(normal,rng,PhononProposal::BondRefresh,{},q);
    PhononSampler b(recovered,rng,PhononProposal::BondRefresh,[&](const auto&e,const auto&){events.push_back(e);},q);
    auto amplitude=qModeAmplitudes(path,q);
    for(std::size_t k=0;k<amplitude.size();++k)for(int j=0;j<path.nodes();++j)
        amplitude[k][j]+=.17*std::cos(.5*j+.8*k);
    const auto candidate=qModeCandidate(path,q,amplitude);
    auto fast=normal.proposeQModePath(candidate),reference=normal.proposeQModePath(candidate,true);
    equal(fast.candidate(),reference.candidate());
    const double ratio=fast.candidate().weight.log_abs_matter-normal.state().weight.log_abs_matter;
    normal.reject(fast);normal.reject(reference);
    const double acceptU=std::exp(std::max(-600.,std::min(0.,ratio)-2.));
    const auto accepted=a.attemptQMode(amplitude,acceptU);
    demand(accepted.accepted,"production deterministic acceptance");
    close(accepted.log_ratio,ratio,2e-9,"production matter ratio");
    close(accepted.general_log_ratio,ratio,2e-9,"production general MH");
    int tries=0;recovered.failureInjection=[&](const char *op,bool isReference){
        if(std::string(op)=="q-mode"&&!isReference&&tries++==0)throw std::runtime_error("production injected q failure");
    };
    const auto retried=b.attemptQMode(amplitude,acceptU);
    demand(retried.accepted&&events.size()==1&&events[0].recovered&&events[0].u==acceptU,"production retry decision");
    std::vector<double> flat;for(const auto &row:amplitude)flat.insert(flat.end(),row.begin(),row.end());
    demand(events[0].values==flat&&b.qCounts().reference_recovered==1,"production same candidate retry");
    equal(normal.state(),recovered.state());recovered.failureInjection={};
    a.ensureMeasured();b.ensureMeasured();
    demand(normal.state().weight.phase_valid&&normal.state().weight.green_valid,"production signed measurement materialization");
    equal(normal.state(),recovered.state());
    equal(normal.state(),normal.uncachedReference(normal.state().path,true,true));
    // Find a genuine unfavorable collective candidate without assuming a sign
    // of the interacting response, then certify atomic rejection in both paths.
    bool rejected=false;
    for(double delta:{-.2,.2,-.7,.7,-1.5,1.5}) {
        const auto saved=normal.state();auto next=qModeAmplitudes(saved.path,q);
        for(int j=0;j<saved.path.nodes();++j)next[0][j]+=delta;
        const auto nextPath=qModeCandidate(saved.path,q,next);
        auto probe=normal.proposeQModePath(nextPath);
        const double logRatio=probe.candidate().weight.log_abs_matter-saved.weight.log_abs_matter;
        normal.reject(probe);equal(saved,normal.state(),true);
        if(logRatio>=-1e-5)continue;
        const auto denied=a.attemptQMode(next,.999999);
        demand(!denied.accepted,"production rejection");close(denied.general_log_ratio,denied.log_ratio,2e-9,"production rejection general MH");
        equal(saved,normal.state(),true);
        const auto retrySaved=recovered.state();
        recovered.failureInjection=[](const char *op,bool isReference){if(std::string(op)=="q-mode"&&!isReference)throw std::runtime_error("production rejection retry");};
        const auto deniedRetry=b.attemptQMode(next,.999999);
        demand(!deniedRetry.accepted,"production recovered rejection");equal(retrySaved,recovered.state(),true);
        demand(b.qCounts().reference_recovered==2,"production recovered rejection counter");rejected=true;break;
    }
    demand(rejected,"production unfavorable candidate missing");
}
void samplerChecks() {
    const HoneycombLattice lat(2,2);auto q=makeQModeOptions(lat,1,0,2,.6);
    const auto c=q.basis[0];std::mt19937_64 rng(81);
    const auto p=PhononPath(c.size(),4,.12,{1.7,.8}).freeSample(rng);
    for(double g:{0.,.4}) {
        const auto gen=generator(g,c),trialGen=generator(0.,c);
        const auto trial=GaussianTrialState::fromMajoranaHamiltonian(trialGen(p.node(0)));
        CoupledPhononPath s(p,trial,gen),retry(p,trial,gen);
        std::vector<NumericalEvent> events;
        PhononSampler a(s,rng,PhononProposal::BondRefresh,{},q),b(retry,rng,PhononProposal::BondRefresh,[&](const auto&e,const auto&){events.push_back(e);},q);
        auto mode=qModeAmplitudes(p,q);for(auto &v:mode[0])v+=.3;
        const auto candidate=qModeCandidate(p,q,mode);
        auto proposal=s.proposeQModePath(candidate);equal(proposal.candidate(),s.uncachedReference(candidate));s.reject(proposal);
        equal(s.fullReference(candidate),s.uncachedReference(candidate));
        const auto result=a.attemptQMode(mode,.5);
        close(result.log_ratio,fock(candidate,g,c)-fock(p,g,c),1e-10,"Fock ratio");
        close(result.log_ratio,result.general_log_ratio,1e-10,"general MH");demand(result.accepted,"positive Fock proposal rejected");
        int first=0,ref=0;retry.failureInjection=[&](const char *op,bool reference){if(std::string(op)=="q-mode"){if(reference)++ref;else if(first++==0)throw std::runtime_error("injected");}};
        b.attemptQMode(mode,.5);equal(s.state(),retry.state());
        demand(ref==1&&events.size()==1&&events[0].recovered&&events[0].u==.5,"same candidate recovery");
        demand(events[0].values==mode[0]&&b.qCounts().reference_recovered==1,"retry values and counters");
        retry.failureInjection={};
        if(g>0.) {
            const auto saved=s.state();auto down=qModeAmplitudes(saved.path,q);for(auto &v:down[0])v-=5.;
            const auto rejected=a.attemptQMode(down,.999);demand(!rejected.accepted,"Fock rejection");equal(saved,s.state(),true);
        }
        const auto saved=retry.state();retry.failureInjection=[](const char*,bool){throw std::runtime_error("persistent");};
        bool failed=false;try{b.attemptQMode(mode,.5);}catch(const std::runtime_error&){failed=true;}
        demand(failed&&b.qCounts().reference_failed==1,"persistent failure must abort");equal(saved,retry.state(),true);
    }
    // Disabled options with a diagnostic basis preserve the default RNG stream.
    const auto gen=generator(0.,c);const auto trial=GaussianTrialState::fromMajoranaHamiltonian(gen(p.node(0)));
    for(auto mode:{PhononProposal::LocalConditional,PhononProposal::BondRefresh}) {
        CoupledPhononPath s(p,trial,gen),t(p,trial,gen);q.every=0;
        PhononSampler a(s,rng,mode),b(t,rng,mode,{},q);
        for(int i=0;i<3;++i){a.sweep();b.sweep();equal(s.state(),t.state());}
        demand(b.qCounts().proposed==0&&a.counts().proposed==b.counts().proposed,"disabled counters");
    }
    CoupledPhononPath s(p,trial,gen);q.every=2;PhononSampler a(s,rng,PhononProposal::BondRefresh,{},q);
    a.setStage("warmup");a.sweep();demand(a.qCounts().proposed==0,"interval start");a.setStage("warmup");a.sweep();
    demand(a.qCounts().proposed==1&&a.qCounts().accepted==1,"g0 collective acceptance");
    a.setStage("measurement");a.sweep();demand(a.qCounts().proposed==1,"stage interval reset");a.sweep();demand(a.qCounts().proposed==2,"stage interval cadence");
    // Small stationary g=0 ensemble exercises the actual sweep RNG/AR kernel.
    CoupledPhononPath freeState(p,trial,gen);q.every=1;q.beta=1.;
    PhononSampler freeSampler(freeState,std::mt19937_64(772),PhononProposal::BondRefresh,{},q);
    double sum=0.,squares=0.,ends=0.;constexpr int samples=256;
    for(int i=0;i<samples;++i) {
        freeSampler.sweep();const auto m=qModeAmplitudes(freeSampler.state().path,q)[0];
        sum+=m[p.center()];squares+=m[p.center()]*m[p.center()];ends+=m.front()*m.back();
    }
    const double variance=1/(2*p.alpha());
    close(sum/samples,0.,7*std::sqrt(variance/samples),"g0 mode mean");
    close(squares/samples,variance,7*variance*std::sqrt(2./samples),"g0 mode variance");
    close(ends/samples,variance*std::pow(p.decay(),p.intervals()),7*variance*std::sqrt(2./samples),"g0 temporal covariance");
    demand(freeSampler.qCounts().accepted==samples,"g0 Gaussian kernel acceptance");

}
}
int main(){try{basisAndBalance();samplerChecks();productionOptions();std::cout<<"q-mode validation passed\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
