#include "phonon_projector.h"
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
using namespace kitaev_ssh;
using Complex = DataType;
namespace {
struct Metric { double error=0,tolerance=0; int count=0; };
std::map<std::string,Metric> metrics;
void check(const std::string &name,double error,double tolerance=2e-9) {
 auto &m=metrics[name];m.error=std::max(m.error,error);++m.count;m.tolerance=tolerance;
 if(!std::isfinite(error)||error>tolerance){std::ostringstream out;out<<std::scientific<<std::setprecision(17)<<name<<" error="<<error<<" tolerance="<<tolerance;throw std::runtime_error(out.str());}
}
std::vector<MatType> operators(int modes) {
 const int d=1<<modes;std::vector<MatType> c(2*modes,MatType::Zero(d,d));
 for(int j=0;j<modes;++j)for(int n=0;n<d;++n){double sign=__builtin_popcount(unsigned(n&((1<<j)-1)))%2?-1.:1.;c[j](n^(1<<j),n)=sign;c[modes+j](n^(1<<j),n)=Complex(0,(n&(1<<j))?sign:-sign);}
 for(int i=0;i<2*modes;++i)for(int j=0;j<2*modes;++j)check("Clifford_algebra",(c[i]*c[j]+c[j]*c[i]-(i==j?2.:0.)*MatType::Identity(d,d)).norm(),0.);
 return c;
}
// Independent physical bonds: do not convert production K back into H.
struct Edge {int a,b,index;double j;};
std::vector<Edge> edges(){std::vector<Edge> e;for(int x=0;x<2;++x)for(int y=0;y<2;++y){int r=2*x+y;e.push_back({r,4+2*(1-x)+y,3*r,1});e.push_back({r,4+2*x+1-y,3*r+1,1});e.push_back({r,4+r,3*r+2,.5});}return e;}
MatType hamiltonian(const std::vector<MatType>&c,const std::vector<double>&x,bool synthetic) {
 MatType h=MatType::Zero(c[0].rows(),c[0].cols());
 for(auto e:edges())h+=Complex(0,e.j+.4*x[e.index])*c[e.a]*c[e.b];
 if(synthetic)h+=Complex(0,.7*x[0])*c[0]*c[1]+Complex(0,.3*x[1])*c[0]*c[2];
 return h;
}
MatType exponential(const MatType&h,double dtau){Eigen::SelfAdjointEigenSolver<MatType> es(h);if(es.info()!=Eigen::Success)throw std::runtime_error("Fock eigensolver");return es.eigenvectors()*(-dtau*es.eigenvalues().array()).exp().matrix().asDiagonal()*es.eigenvectors().adjoint();}
struct Oracle {
 std::vector<MatType> u;std::vector<cVecType> right;std::vector<Eigen::RowVectorXcd> left;Complex weight;
 Oracle(const PhononPath&p,const std::vector<MatType>&c,const cVecType&psi,bool synthetic){
  // Explicit endpoint half-node rule, with interior nodes occurring twice.
  for(int link=0;link<p.intervals();++link){u.push_back(exponential(hamiltonian(c,p.node(link),synthetic),p.dt()/2));u.push_back(exponential(hamiltonian(c,p.node(link+1),synthetic),p.dt()/2));}
  right.resize(u.size()+1);left.resize(u.size()+1);right[0]=psi;for(int f=0;f<int(u.size());++f)right[f+1]=u[f]*right[f];left.back()=psi.adjoint();for(int f=int(u.size())-1;f>=0;--f)left[f]=left[f+1]*u[f];weight=(left[0]*right[0])(0);
 }
 Complex expectation(const MatType&o,int cut)const{return (left[cut]*o*right[cut])(0)/weight;}
 MatType covariance(const std::vector<MatType>&c,int cut)const{MatType out(c.size(),c.size());for(int i=0;i<int(c.size());++i)for(int j=0;j<int(c.size());++j)out(i,j)=expectation(c[i]*c[j],cut);return out;}
};
void run(bool synthetic){
 const int Ntau=6;const double dtau=.11,tau=Ntau*dtau,tau_mid=tau/2;
 auto c=operators(4);Couplings j{1,1,.5,.4};KitaevSSHModel model(HoneycombLattice(2,2),j);auto trial=freeGroundState(model);auto base=honeycombGenerator(model);
 MajoranaGenerator generator=synthetic?MajoranaGenerator([base](const std::vector<double>&x){MatType k=base(x);k(0,1)+=Complex(0,-1.4*x[0]);k(1,0)-=Complex(0,-1.4*x[0]);k(0,2)+=Complex(0,-.6*x[1]);k(2,0)-=Complex(0,-.6*x[1]);return k;}):base;
 Eigen::SelfAdjointEigenSolver<MatType> initial(hamiltonian(c,std::vector<double>(12,0.),false));if(initial.info()!=Eigen::Success)throw std::runtime_error("trial ED");cVecType psi=initial.eigenvectors().col(0);
 std::mt19937_64 rng(1997);std::uniform_real_distribution<double> random(-.8,.8);Complex previousReference=1,previousMeasured=1;bool phaseWitness=false;
 for(int kind=0;kind<3;++kind){PhononPath p(12,Ntau,dtau,{1.7,.8},std::vector<double>(12,.1378676470588235));
  if(kind==0){for(int node=0;node<=Ntau;++node)for(int b=0;b<12;++b)p=p.withCoordinate(node,b,.2);}
  else if(kind==1){for(int b=0;b<12;++b){double v=random(rng);for(int node=0;node<=Ntau;++node)p=p.withCoordinate(node,b,v);}}
  else p=p.freeSample(rng);
  Oracle ref(p,c,psi,synthetic);auto factors=pathFactors(p,8,generator);auto measured=evaluatePath(p,trial,factors,3);Complex weight=std::exp(measured.log_abs_matter)*measured.matter_phase;
  check("signed_complex_weight",std::abs(weight/ref.weight-Complex(1)));check("complex_weight_ratio",std::abs((weight/previousMeasured)/(ref.weight/previousReference)-Complex(1)));previousMeasured=weight;previousReference=ref.weight;
  MatType centerC=ref.covariance(c,Ntau);check("center_Green",(measured.center_green-(MatType::Identity(8,8)-centerC)).cwiseAbs().maxCoeff());
  if(!synthetic)check("production_matter_energy",std::abs(matterEnergyAtDisplacement(model,p.node(Ntau/2),measured.center_green)-ref.expectation(hamiltonian(c,p.node(Ntau/2),false),Ntau)));
  MatType covariance=MatType::Identity(8,8)-measured.center_green;
  for(auto e:edges()){double hopping=e.j+.4*p.coordinate(Ntau/2,e.index);Complex predicted=Complex(0,hopping)*covariance(e.a,e.b);check("bond_energy",std::abs(predicted-ref.expectation(Complex(0,hopping)*c[e.a]*c[e.b],Ntau)));
   for(auto f:edges()){double other=f.j+.4*p.coordinate(Ntau/2,f.index);Complex wick=covariance(e.a,e.b)*covariance(f.a,f.b)-covariance(e.a,f.a)*covariance(e.b,f.b)+covariance(e.a,f.b)*covariance(e.b,f.a);check("bond_energy_four_point",std::abs(-hopping*other*wick-ref.expectation(-hopping*other*c[e.a]*c[e.b]*c[f.a]*c[f.b],Ntau)));}}
  // No production unequal-time measurement API exists. Test its implied
  // one-body transport using production propagators and arbitrary-cut Green.
  for(auto cuts:{std::pair<int,int>{0,Ntau},std::pair<int,int>{2,Ntau+2},std::pair<int,int>{Ntau,2*Ntau}}){int early=cuts.first,late=cuts.second;auto endpoint=pureProjectorEndpointRebuild(trial,factors,late,3,PureProjectorOptions{});if(!endpoint.ok())throw std::runtime_error("cut Green rebuild");MatType lateC=MatType::Identity(8,8)-endpoint.green.adjoint();MatType b=MatType::Identity(8,8),v=MatType::Identity(16,16);for(int f=early;f<late;++f){b=(b*MatType(factors[f].matrix)).eval();v=(ref.u[f]*v).eval();}MatType predicted=lateC*b.transpose();for(int i=0;i<8;++i)for(int k=0;k<8;++k){Complex exact=(ref.left[late]*c[i]*v*c[k]*ref.right[early])(0)/ref.weight;check("unequal_time_transport",std::abs(predicted(i,k)-exact));}}
  if(kind==2){cVecType wrong=psi;
   // Reversing all Hermitian factors conjugates the amplitude; for real
   // paths exchange two interior factors to test a genuine order change.
   auto reordered=ref.u;std::swap(reordered[2],reordered[4]);wrong=psi;for(const auto &u:reordered)wrong=u*wrong;Complex wrongWeight=(psi.adjoint()*wrong)(0);check("slice_order_sensitivity",std::abs((wrongWeight-ref.weight)/ref.weight)>1e-6?0.:1.,0.);}
  if(kind==2)check("noncommuting_path_witness",(hamiltonian(c,p.node(1),synthetic)*hamiltonian(c,p.node(2),synthetic)-hamiltonian(c,p.node(2),synthetic)*hamiltonian(c,p.node(1),synthetic)).norm()>1e-3?0.:1.,0.);
  phaseWitness=phaseWitness||std::abs(ref.weight.imag()/std::abs(ref.weight))>1e-5;
  std::cout<<"REFERENCE "<<(synthetic?"synthetic_same_sublattice":"honeycomb")<<" path="<<kind<<" tau="<<tau<<" dtau="<<dtau<<" Ntau="<<Ntau<<" tau_mid="<<tau_mid<<" Fock_weight="<<ref.weight<<" PfQMC_weight="<<weight<<" relative_error="<<std::abs(weight/ref.weight-Complex(1))<<" PASS\n";
 }
 if(synthetic)check("nontrivial_complex_phase_witness",phaseWitness?0.:1.,0.);
}
void negativeRealWitness(){
 // Four Majoranas, real Hermitian Fock Hamiltonians. Products of positive
 // factors need not have positive trial overlaps. No scalar eta is injected.
 const auto c=operators(2);const int Ntau=6;const double dtau=.3;
 auto h=[&](const std::vector<double>&x){MatType out=MatType::Zero(4,4);for(int a=0;a<2;++a)for(int b=0;b<2;++b)out+=Complex(0,x[2*a+b])*c[a]*c[2+b];return out;};
 MajoranaGenerator generator=[](const std::vector<double>&x){MatType k=MatType::Zero(4,4);for(int a=0;a<2;++a)for(int b=0;b<2;++b){k(a,2+b)=Complex(0,-2*x[2*a+b]);k(2+b,a)=-k(a,2+b);}return k;};
 std::vector<double> initialX{1,0,0,1};Eigen::SelfAdjointEigenSolver<MatType> initial(h(initialX));cVecType psi=initial.eigenvectors().col(0);auto trial=GaussianTrialState::fromMajoranaHamiltonian(generator(initialX));
 std::mt19937_64 rng(7262026);std::normal_distribution<double> random(0,3.);
 for(int attempt=0;attempt<256;++attempt){PhononPath p(4,Ntau,dtau,{1,1});for(int node=0;node<=Ntau;++node)for(int b=0;b<4;++b)p=p.withCoordinate(node,b,random(rng));cVecType right=psi;for(int link=0;link<Ntau;++link){right=exponential(h(p.node(link)),dtau/2)*right;right=exponential(h(p.node(link+1)),dtau/2)*right;}Complex exact=(psi.adjoint()*right)(0);
  if(exact.real() < -1e-5 && std::abs(exact.imag())<1e-10*std::abs(exact)){auto measured=evaluatePath(p,trial,pathFactors(p,4,generator),2);Complex weight=std::exp(measured.log_abs_matter)*measured.matter_phase;check("genuine_negative_real_weight",std::abs(weight/exact-Complex(1)));check("negative_phase",std::abs(measured.matter_phase+Complex(1)));std::cout<<"REFERENCE negative_real attempt="<<attempt<<" seed=7262026 tau="<<Ntau*dtau<<" dtau="<<dtau<<" Ntau="<<Ntau<<" tau_mid="<<Ntau*dtau/2<<" Fock_weight="<<exact<<" PfQMC_weight="<<weight<<" PASS\n";return;}
 }
 throw std::runtime_error("deterministic negative real witness not found");
}

}
int main(){std::cout<<std::setprecision(17);try{run(false);run(true);negativeRealWitness();for(auto &m:metrics)std::cout<<"PASS "<<m.first<<" comparisons="<<m.second.count<<" max_error="<<m.second.error<<" tolerance="<<m.second.tolerance<<'\n';return 0;}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
