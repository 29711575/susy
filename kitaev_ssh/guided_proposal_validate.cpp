#define KITAEV_SSH_TEST_FAILURE_INJECTION
// Fixed-reference pCN correctness tests, including deliberately wrong guides.
#include "phonon_mc.h"
#include <iostream>
#include <sstream>
using namespace kitaev_ssh;
namespace {
double worst=0.,worstStat=0.;
void check(double x,double y,double tol,const std::string &what) {
 const double err=std::abs(x-y);
 const bool statistical=what=="exact heatbath mean" || what=="exact heatbath variance";
 if(statistical)worstStat=std::max(worstStat,err);else worst=std::max(worst,err);
 if(!std::isfinite(err)||err>tol)throw std::runtime_error(what+" err="+std::to_string(err));
}
void demand(bool v,const std::string&w){if(!v)throw std::runtime_error(w);}
template<class F>void rejects(F f){bool threw=false;try{f();}catch(const std::exception&){threw=true;}demand(threw,"invalid input accepted");}
MajoranaGenerator gen(double g) {
 return [g](const std::vector<double>&x){MatType k=MatType::Zero(2,2);k(0,1)=DataType(0,-2*(.5+g*x.at(0)));k(1,0)=-k(0,1);return k;};
}
double fockLog(const PhononPath&p,double g) {
 double out=.5*p.totalTau();
 for(int j=0;j<p.nodes();++j)out+=g*p.dt()*(j==0||j==p.intervals()?.5:1.)*p.coordinate(j,0);
 return out;
}
void same(const CoupledPhononPath::Snapshot&a,const CoupledPhononPath::Snapshot&b) {
 for(int j=0;j<a.path.nodes();++j)for(int k=0;k<a.path.bonds();++k)check(a.path.coordinate(j,k),b.path.coordinate(j,k),0.,"rejected coordinate mutation");
 check(a.weight.log_abs_matter,b.weight.log_abs_matter,0.,"rejected weight mutation");
 demand(a.factors.size()==b.factors.size(),"factor count mutation");
 for(std::size_t i=0;i<a.factors.size();++i)check((MatType(a.factors[i].matrix)-MatType(b.factors[i].matrix)).norm(),0.,0.,"factor mutation");
}
void algebra() {
 std::mt19937_64 rng(581);
 for(double beta:{.17,.65,1.})for(double f:{-.9,0.,1.1}) {
  const auto p=PhononPath(1,8,.2,{2.3,.7},{.21},{-.13}).freeSample(rng);
  GuidedBondOptions o;o.beta=beta;o.forces={f};GuidedBondProposal q(p,o);
  const auto to=q.draw(p,0,rng);auto n=p.withBondCoordinates(0,to);
  std::vector<double> old(p.nodes());for(int j=0;j<p.nodes();++j)old[j]=p.coordinate(j,0);
  const double reverse=q.logTransition(n,0,old)-q.logTransition(p,0,to);
  check(n.logWeight()-p.logWeight()+reverse,-q.referenceScoreDifference(p,0,to),3e-12,"Gaussian reference cancellation");
  rejects([&]{q.logTransition(p,0,{});});
  rejects([&]{q.referenceScoreDifference(p,0,{});});
 }
 rejects([&]{GuidedBondOptions o;o.beta=0;GuidedBondProposal q(PhononPath(1,4,.1,{1.,1.}),o);});
 // Proposal mean is C times the reference source, including shifted trial and
 // physical force; compare to an independent dense covariance sum.
 const PhononPath p(2,10,.07,{1.7,.8},{.2,-.1},{.3,-.2});GuidedBondOptions o;o.forces={-.7,.9};GuidedBondProposal q(p,o);
 for(int j=0;j<p.nodes();++j)for(int b=0;b<p.bonds();++b) {
  double expected=p.freeMean(j,b);
  for(int l=0;l<p.nodes();++l)expected+=std::pow(p.decay(),std::abs(j-l))/(2*p.alpha())*p.sourceNodeWeight(l)*o.forces[b];
  check(q.mean(b)[j],expected,1e-13,"guided open endpoint covariance");
 }
}
void transactions() {
 std::mt19937_64 rng(3367);
 for(bool deferred:{false,true})for(double beta:{.3,1.}) {
  auto path=PhononPath(1,8,.2,{2.,.7},{.2},{.1}).freeSample(rng);
  const double g=.8;auto generator=gen(g);auto trial=GaussianTrialState::fromMajoranaHamiltonian(generator({0.}));
  PhononWeightOptions options;options.defer_phase_measurement=deferred;options.singular_values_only=true;options.merge_node_half_steps=true;
  CoupledPhononPath state(path,trial,generator,4,options),retry(path,trial,generator,4,options);
  GuidedBondOptions guide;guide.beta=beta;guide.forces={.3};
  GuidedBondProposal q(path,guide);std::vector<NumericalEvent> events;
  PhononSampler a(state,rng,PhononProposal::GuidedBond,{}, {},guide);
  PhononSampler b(retry,rng,PhononProposal::GuidedBond,[&](const auto&e,const auto&){events.push_back(e);}, {},guide);
  std::vector<double> values(path.nodes());for(int j=0;j<path.nodes();++j)values[j]=path.coordinate(j,0)+.3;
  auto candidate=path.withBondCoordinates(0,values);const double expected=fockLog(candidate,g)-fockLog(path,g)-q.referenceScoreDifference(path,0,values);
  auto move=a.attemptGuidedBond(0,values,.5);
  demand(move.accepted,"positive guided move refused");check(move.log_ratio,expected,2e-11,"Fock guided ratio");check(move.general_log_ratio,expected,2e-11,"full guided MH ratio");
  // Proposal density + target density equality, not only simplified ratio.
  std::vector<double> old(path.nodes());for(int j=0;j<path.nodes();++j)old[j]=path.coordinate(j,0);
  check(path.logWeight()+fockLog(path,g)+q.logTransition(path,0,values)+std::min(0.,expected),
        candidate.logWeight()+fockLog(candidate,g)+q.logTransition(candidate,0,old)+std::min(0.,-expected),2e-11,"full detailed balance");
  retry.failureInjection=[](const char*,bool reference){if(!reference)throw std::runtime_error("injected first failure");};
  auto recovered=b.attemptGuidedBond(0,values,.5);
  demand(recovered.accepted&&events.size()==1&&events[0].recovered&&events[0].values==values&&events[0].u==.5,"same candidate retry");
  retry.failureInjection={};a.ensureMeasured();b.ensureMeasured();check(state.state().weight.log_abs_matter,retry.state().weight.log_abs_matter,1e-12,"retry accepted weight");
  auto before=state.state();for(int j=0;j<path.nodes();++j)values[j]=before.path.coordinate(j,0)-2.;
  move=a.attemptGuidedBond(0,values,.999);demand(!move.accepted,"guided unfavorable proposal accepted");same(before,state.state());
  before=retry.state();retry.failureInjection=[](const char*,bool){throw std::runtime_error("persistent");};
  rejects([&]{b.attemptGuidedBond(0,values,.5);});same(before,retry.state());demand(b.counts().reference_failed==1,"persistent failure did not count");
 }
}
void noncommuting() {
 const KitaevSSHModel model(HoneycombLattice(3,3,1),p2Couplings(.7));
 const auto &lat=model.lattice();std::vector<double> shift(lat.bondCount()),force(lat.bondCount());
 GuidedBondOptions o;o.beta=.6;o.forces.resize(lat.bondCount());
 for(const auto &b:lat.bonds()) {shift[b.index]=.1*std::cos(b.index);force[b.index]=.08*std::sin(b.index);o.forces[b.index]=.4+.1*std::sin(.7*b.index);}
 std::mt19937_64 rng(923);const auto path=PhononPath(lat.bondCount(),6,.2,{1.7,.8},shift,force).freeSample(rng);
 auto generator=honeycombGenerator(model);auto trial=freeGroundState(model);
 demand((generator(path.node(0))*generator(path.node(1))-generator(path.node(1))*generator(path.node(0))).norm()>1e-6,"test accidentally commuting");
 PhononWeightOptions options;options.defer_phase_measurement=true;options.singular_values_only=true;options.merge_node_half_steps=true;
 CoupledPhononPath state(path,trial,generator,4,options);PhononSampler sampler(state,rng,PhononProposal::GuidedBond,{}, {},o);GuidedBondProposal q(path,o);
 for(int b:{0,5,10,20}) {
  const auto before=state.state();const auto values=q.draw(before.path,b,rng);const auto candidate=before.path.withBondCoordinates(b,values);
  const auto reference=state.uncachedReference(candidate,true,true);
  std::vector<double> old(before.path.nodes());for(int j=0;j<before.path.nodes();++j)old[j]=before.path.coordinate(j,b);
  const double full=reference.weight.logAbsTotal()-before.weight.logAbsTotal()+q.logTransition(candidate,b,old)-q.logTransition(before.path,b,values);
  const double u=std::exp(std::max(-600.,std::min(0.,full)-1.));
  auto result=sampler.attemptGuidedBond(b,values,u);
  demand(result.accepted,"noncommuting deterministic accepted proposal rejected");
  check(result.log_ratio,full,3e-9,"noncommuting full reference ratio");check(result.general_log_ratio,full,3e-9,"noncommuting full MH ratio");
  sampler.ensureMeasured();check(state.state().weight.log_abs_matter,reference.weight.log_abs_matter,3e-9,"noncommuting signed log weight");
  check(std::abs(state.state().weight.matter_phase-reference.weight.matter_phase),0.,3e-9,"noncommuting phase");
 }
}
void exactLinearReference() {
 // For a single conserved Majorana bilinear, the matter action is linear in X.
 // At this particular guide force, beta=1 is an exact independent heatbath.
 const PhononPath path(1,12,.2,{2.,.7});const double g=.8;
 const double f=g*(path.dt()/2)/path.sourceNodeWeight(0);
 GuidedBondOptions o;o.beta=1.;o.forces={f};GuidedBondProposal q(path,o);
 auto generator=gen(g);auto trial=GaussianTrialState::fromMajoranaHamiltonian(generator({0.}));
 PhononWeightOptions opts;opts.defer_phase_measurement=true;opts.compute_measurement_green=false;opts.singular_values_only=true;opts.merge_node_half_steps=true;
 CoupledPhononPath state(path,trial,generator,4,opts);PhononSampler sampler(state,std::mt19937_64(1181),PhononProposal::GuidedBond,{}, {},o);
 double sum=0.,square=0.;const int N=1500;const double mean=q.mean(0)[path.center()],var=1/(2*path.alpha());
 for(int i=0;i<N;++i){sampler.sweep();const double y=sampler.state().path.coordinate(path.center(),0)-mean;sum+=y;square+=y*y;}
 demand(sampler.counts().accepted==N,"exact linear-reference heatbath must accept all");
 check(sum/N,0.,7*std::sqrt(var/N),"exact heatbath mean");check(square/N,var,7*var*std::sqrt(2./N),"exact heatbath variance");
}
void pool() {
 for(const HoneycombLattice lat:{HoneycombLattice(2,2),HoneycombLattice(3,3,1),HoneycombLattice(9,9)}) {
  for(const std::string types:{"x","y","z","xyz","yx"}) {
   auto q=makeQModeOptions(lat,1,1,1,.5,types);
   for(std::size_t a=0;a<q.basis.size();++a)for(std::size_t b=0;b<q.basis.size();++b) {
    double dot=0;for(std::size_t i=0;i<q.basis[a].size();++i)dot+=q.basis[a][i]*q.basis[b][i];
    check(dot,a==b?1.:0.,2e-14,"multi-bond q orthonormality");
   }
   std::mt19937_64 rng(551);auto p=PhononPath(lat.bondCount(),4,.13,{2.,.7}).freeSample(rng);auto a=qModeAmplitudes(p,q),b=a;
   for(auto &row:b)for(auto&v:row)v+=.12;
   auto n=qModeCandidate(p,q,b);
   check(n.logWeight()-p.logWeight()+qModeTransitionLogDensity(p,b,a,q.beta)-qModeTransitionLogDensity(p,a,b,q.beta),0.,5e-10,"multi-bond q MH");
  }
 }
 const HoneycombLattice lat(9,9);
 std::istringstream in("# q pool\n5 4\n3 5\n4 6 # mirror\n");auto pool=readQModePool(in,lat,1,.4,"xyz");demand(pool.size()==3,"pool parser");
 for(const auto &text:{"5 4\n4 5\n","5 4 extra\n","0 0\n","9 3\n",""})rejects([&]{std::istringstream bad(text);readQModePool(bad,lat,1,.4);});
 rejects([&]{makeQModeOptions(lat,5,4,1,.4,"zz");});
 // Exercise state-independent pool selection and per-component counts at g=0.
 const HoneycombLattice small(2,2);std::istringstream modes("1 0\n0 1\n");auto opts=readQModePool(modes,small,1,.6,"xyz");
 auto path=PhononPath(small.bondCount(),4,.2,{2.,.7});KitaevSSHModel model(small,p2Couplings(0.));
 CoupledPhononPath state(path,freeGroundState(model),honeycombGenerator(model));
 PhononSampler sampler(state,std::mt19937_64(17),PhononProposal::BondRefresh,{},opts[0],{},opts);
 for(int n=0;n<40;++n)sampler.sweep();
 demand(sampler.qCounts().proposed==40&&sampler.qCounts().accepted==40,"g0 pool acceptance");
 demand(sampler.qPoolCounts()[0].proposed+sampler.qPoolCounts()[1].proposed==40,"pool accounting");
 demand(sampler.qPoolCounts()[0].proposed>0&&sampler.qPoolCounts()[1].proposed>0,"pool coverage");
}
}
int main(){try{algebra();transactions();noncommuting();pool();exactLinearReference();std::cout<<"guided proposal validation PASS; worst deterministic residual="<<std::setprecision(17)<<worst<<"; largest finite-sample difference="<<worstStat<<'\n';return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
