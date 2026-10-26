#define main original_q_mode_validation_main
#include "q_mode_validate.cpp"
#undef main
#include "physical_matter_trial.h"
namespace {
void gammaBalance() {
 std::mt19937_64 rng(20261020);std::normal_distribution<double> normal;
 for(auto lat:{HoneycombLattice(2,2),HoneycombLattice(6,6),HoneycombLattice(5,4,2)}) {
  std::vector<double>d(lat.bondCount()),f(d.size());for(std::size_t b=0;b<d.size();++b){d[b]=.3*std::sin(b);f[b]=.2*std::cos(b);}
  auto p=PhononPath(d.size(),8,.13,{1.7,.8},d,f).freeSample(rng);
  for(double beta:{.2,.5,1.})for(const auto&q:makeGammaOptions(lat,1,beta).modes) {
   close(std::inner_product(q.basis[0].begin(),q.basis[0].end(),q.basis[0].begin(),0.),1.,1e-12,"Gamma real norm");
   auto from=qModeAmplitudes(p,q),a=drawProjectedPCN(p,q,rng,normal);auto next=qModeCandidate(p,q,a);
   close(next.logWeight()-p.logWeight()+qModeTransitionLogDensity(p,a,from,beta)-qModeTransitionLogDensity(p,from,a,beta),0.,2e-9,"Gamma shifted/source full MH cancellation");
   for(int j=0;j<p.nodes();++j) {
    double delta=0.;bool first=true;
    for(const auto&b:lat.bonds()) {
     const double diff=next.coordinate(j,b.index)-p.coordinate(j,b.index);
     if(q.bond_types[0]!="xyz"[int(b.type)])close(diff,0.,1e-30,"Gamma preserves other bonds");
     else {if(first){delta=diff;first=false;}close(diff,delta,1e-14,"Gamma preserves nonzero modes");}
    }
   }
  }
 }
 std::cout<<"PASS Gamma real normalization; open endpoints, arbitrary shifts/sources, complementary modes, full MH beta .2/.5/1\n";
 // g=0 exact stationary open law; independent draws at beta1.
 HoneycombLattice lat(2,2);const auto gamma=makeGammaOptions(lat,1,1.);auto q=gamma.modes[0];
 auto p=PhononPath(12,8,.13,{1.7,.8}).freeSample(rng);auto gen=generator(0.,q.basis[0]);auto trial=GaussianTrialState::fromMajoranaHamiltonian(gen(p.node(0)));
 CoupledPhononPath state(p,trial,gen);PhononSampler sampler(state,rng,PhononProposal::BondRefresh,{}, {},{}, {},gamma);
 sampler.setStage("measurement");double mean=0,sq=0,cov=0;const int N=512;
 for(int i=0;i<N;++i){sampler.sweep();auto a=qModeAmplitudes(sampler.state().path,q)[0];mean+=a[p.center()];sq+=a[p.center()]*a[p.center()];cov+=a.front()*a.back();}
 const double var=1/(2*p.alpha());close(mean/N,0.,7*std::sqrt(var/N),"Gamma g0 mean");close(sq/N,var,7*var*std::sqrt(2./N),"Gamma g0 variance");close(cov/N,var*std::pow(p.decay(),p.intervals()),7*var*std::sqrt(2./N),"Gamma endpoint covariance");
 demand(sampler.gammaCounts().accepted==3*N,"Gamma g0 all accept");
 std::cout<<"PASS Gamma g0 analytic moments, endpoint covariance and acceptance; samples="<<N<<"\n";
}
std::vector<MatType> gammaClifford() {
 std::vector<MatType>c(8,MatType::Zero(16,16));for(int j=0;j<4;++j)for(int n=0;n<16;++n){double s=__builtin_popcount(unsigned(n&((1<<j)-1)))%2?-1.:1.;c[j](n^(1<<j),n)=s;c[4+j](n^(1<<j),n)=DataType(0,(n&(1<<j))?s:-s);}return c;
}
MatType gammaH(const std::vector<MatType>&c,const std::vector<double>&x,bool synthetic) {
 MatType h=MatType::Zero(16,16);
 for(int xx=0;xx<2;++xx)for(int yy=0;yy<2;++yy){int r=2*xx+yy;h+=DataType(0,1+.4*x[3*r])*c[r]*c[4+2*(1-xx)+yy];h+=DataType(0,1+.4*x[3*r+1])*c[r]*c[4+2*xx+1-yy];h+=DataType(0,.5+.4*x[3*r+2])*c[r]*c[4+r];}
 if(synthetic)h+=DataType(0,.7*x[0])*c[0]*c[1]+DataType(0,.3*x[1])*c[0]*c[2];return h;
}
MatType gammaExp(const MatType&h,double dt){Eigen::SelfAdjointEigenSolver<MatType>es(h);return es.eigenvectors()*(-dt*es.eigenvalues().array()).exp().matrix().asDiagonal()*es.eigenvectors().adjoint();}
struct GammaFock {
 DataType w;MatType green;
 GammaFock(const PhononPath&p,const std::vector<MatType>&c,const cVecType&psi,bool synthetic){
  std::vector<MatType>u;for(int j=0;j<p.intervals();++j){u.push_back(gammaExp(gammaH(c,p.node(j),synthetic),p.dt()/2));u.push_back(gammaExp(gammaH(c,p.node(j+1),synthetic),p.dt()/2));}
  cVecType ket=psi;for(int i=0;i<p.intervals();++i)ket=u[i]*ket;
  Eigen::RowVectorXcd bra=psi.adjoint();for(int i=int(u.size())-1;i>=p.intervals();--i)bra=bra*u[i];w=(bra*ket)(0);
  green=MatType::Identity(8,8);for(int i=0;i<8;++i)for(int j=0;j<8;++j)green(i,j)-=(bra*c[i]*c[j]*ket)(0)/w;
 }
};
void gammaFockTransactions() {
 const HoneycombLattice lat(2,2);KitaevSSHModel model(lat,{1,1,.5,.4});auto chosen=selectFreeMatterTrial(model);const auto c=gammaClifford();
 MatType parent=MatType::Zero(16,16);for(int i=0;i<8;++i)for(int j=i+1;j<8;++j)parent-=chosen.trial.G_T(i,j)*c[i]*c[j];Eigen::SelfAdjointEigenSolver<MatType>es(parent);cVecType psi=es.eigenvectors().col(0);
 bool phaseWitness=false;int n=0;
 for(bool synthetic:{false,true})for(bool production:{false,true})for(bool recover:{false,true})for(bool accept:{false,true})for(int mu=0;mu<3;++mu) {
  auto base=honeycombGenerator(model);MajoranaGenerator gen=[base,synthetic](const std::vector<double>&x){MatType k=base(x);if(synthetic){k(0,1)+=DataType(0,-1.4*x[0]);k(1,0)-=DataType(0,-1.4*x[0]);k(0,2)+=DataType(0,-.6*x[1]);k(2,0)-=DataType(0,-.6*x[1]);}return k;};
  auto gamma=makeGammaOptions(lat,1,.5);auto q=gamma.modes[mu];
  PhononPath p(12,6,.11,{1.7,.8},std::vector<double>(12,.13));
  for(int j=0;j<p.nodes();++j)for(int b=0;b<12;++b)p=p.withCoordinate(j,b,.4*std::sin(.73*b+.91*j));
  auto a=qModeAmplitudes(p,q);for(int j=0;j<p.nodes();++j)a[0][j]+=.4+.2*std::cos(j*.6);
  auto next=qModeCandidate(p,q,a);GammaFock oldRef(p,c,psi,synthetic),newRef(next,c,psi,synthetic);
  if(std::log(std::abs(newRef.w/oldRef.w))>0)std::swap(p,next);
  oldRef=GammaFock(p,c,psi,synthetic);newRef=GammaFock(next,c,psi,synthetic);
  const double expected=std::log(std::abs(newRef.w/oldRef.w));demand(expected<-1e-8,"Gamma nontrivial reject witness");
  PhononWeightOptions opts{production,production,production,production,production,true,production};
  CoupledPhononPath state(p,chosen.trial,gen,2,opts);std::vector<NumericalEvent>events;
  PhononSampler sampler(state,std::mt19937_64(7),PhononProposal::BondRefresh,[&](const auto&e,const auto&){events.push_back(e);},{},{},{},gamma);
  const auto saved=state.state();if(recover)state.failureInjection=[](const char*op,bool reference){if(std::string(op)=="q-mode"&&!reference)throw std::runtime_error("Gamma candidate recovery");};
  const double u=accept?std::exp(std::max(-600.,expected-1.)):(1+std::exp(expected))/2;
  const auto result=sampler.attemptGamma(mu,qModeAmplitudes(next,q),u);close(result.general_log_ratio,expected,3e-9,"Gamma independent 16-state Fock MH ratio");demand(result.accepted==accept,"Gamma Fock MH decision");
  if(recover)demand(events.size()==1&&events[0].recovered&&events[0].u==u,"Gamma same-candidate retry");
  if(!accept)equal(saved,state.state(),true);
  state.failureInjection={};
  if(production && recover && accept)state.failureInjection=[](const char*op,bool reference){if(std::string(op)=="measurement"&&!reference)throw std::runtime_error("Gamma deferred phase recovery");};
  sampler.ensureMeasured();state.failureInjection={};
  if(production && recover && accept)demand(events.size()==2&&events.back().operation=="measurement"&&events.back().recovered,"Gamma phase recovery observer");
  const auto&ref=accept?newRef:oldRef;
  close(std::abs(std::exp(state.state().weight.log_abs_matter)*state.state().weight.matter_phase/ref.w-DataType(1)),0.,3e-9,"Gamma signed Fock weight after decision");
  close((state.state().weight.center_green-ref.green).norm(),0.,3e-9,"Gamma Fock Green after decision");
  equal(state.state(),state.uncachedReference(state.state().path,true,true));
  phaseWitness=phaseWitness||std::abs(ref.w.imag()/std::abs(ref.w))>1e-6;
  // Persistent failure: no uncertified rejection, no changed path/cache.
  const auto before=state.state();state.failureInjection=[](const char*,bool){throw std::runtime_error("Gamma persistent failure");};bool failed=false;
  try{sampler.attemptGamma(mu,qModeAmplitudes(next,q),u);}catch(const std::runtime_error&){failed=true;}demand(failed,"Gamma must abort persistent failure");equal(before,state.state(),true);++n;
 }
 demand(phaseWitness,"Gamma nontrivial complex phase witness");std::cout<<"PASS Gamma independent physical-parity 16-state Fock, signed phase, accepted/rejected caches, same-candidate recovery and fail-closed; cases="<<n<<"\n";
}
void gammaOffAndPool() {
 std::mt19937_64 rng(123);HoneycombLattice lat(2,2);auto q=makeQModeOptions(lat,1,1,1,.5);auto p=PhononPath(12,4,.1,{1.7,.8}).freeSample(rng);auto gen=generator(0.,q.basis[0]);auto trial=GaussianTrialState::fromMajoranaHamiltonian(gen(p.node(0)));
 CoupledPhononPath a(p,trial,gen),b(p,trial,gen);PhononSampler x(a,rng,PhononProposal::BondRefresh,{},q,{}, {q});PhononSampler y(b,rng,PhononProposal::BondRefresh,{},q,{}, {q},makeGammaOptions(lat,0,.3));
 for(int i=0;i<4;++i){x.sweep();y.sweep();equal(a.state(),b.state(),true);}demand(!y.gammaCounts().proposed,"Gamma off counters");
 auto extra=makeQModeOptions(lat,1,0,1,.5);CoupledPhononPath s(p,trial,gen);PhononSampler z(s,rng,PhononProposal::BondRefresh,{},q,{}, {q},makeGammaOptions(lat,1,.5),{extra});
 z.setStage("warmup");z.sweep();z.setStage("measurement");for(int i=0;i<5;++i)z.sweep();
 demand(z.qMeasuredCounts()[0].proposed==5&&z.qMeasuredCounts()[1].proposed==5,"pool original cadence preserved");for(const auto&v:z.gammaMeasuredCounts())demand(v.proposed==5,"Gamma formal counts exclude warmup");
 std::cout<<"PASS Gamma off exact paths/caches/RNG; additional pool cadence and formal-stage counters\n";
}
}
int main(){try{basisAndBalance();samplerChecks();productionOptions();gammaBalance();gammaFockTransactions();gammaOffAndPool();return 0;}catch(const std::exception&e){std::cerr<<"FAIL Gamma "<<e.what()<<'\n';return 1;}}
