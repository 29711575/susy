#include "physical_matter_trial.h"
#include "phonon_projector.h"
#include <iostream>
#include <iomanip>
#include <numeric>
using namespace kitaev_ssh;
namespace {
void require(bool v,const char *s){if(!v)throw std::runtime_error(s);}
void near(double x,double tol,const char *s){require(std::isfinite(x)&&std::abs(x)<=tol,s);}
int cliffordParity(const KitaevSSHModel &model) {
 const int n=model.lattice().cells();std::vector<int>w;
 for(const auto &b:model.lattice().bonds()){int mu=int(b.type);w.push_back(4*b.a+mu);w.push_back(4*b.b+mu);}
 for(int r=0;r<n;++r){w.push_back(4*r+3);w.push_back(4*(n+r)+3);}
 int sign=1;for(std::size_t i=0;i<w.size();++i)for(std::size_t j=i+1;j<w.size();++j)if(w[i]>w[j])sign=-sign;
 if(n%2)sign=-sign;for(int u:model.gauge())sign*=u;return sign;
}
void geometry() {
 int count=0,blocked=0;
 for(auto shape:std::vector<std::array<int,3>>{{2,2,0},{4,4,0},{6,6,0},{9,9,0},{3,4,1}})
 for(auto sector:std::vector<std::array<int,2>>{{1,1},{-1,1},{1,-1},{-1,-1}}) {
  HoneycombLattice lat(shape[0],shape[1],shape[2]);KitaevSSHModel model(lat,{1,1,.5,.4},{},lat.zeroFluxGauge(sector[0],sector[1]));
  auto base=freeGroundState(model);auto chosen=selectFreeMatterTrial(model);auto legacy=selectFreeMatterTrial(model,"unrestricted");
  require(chosen.required_parity==cliffordParity(model),"direct Clifford physical condition");
  require(chosen.trial.fermionParity()==chosen.required_parity&&chosen.physicalParityCompatible(),"physical parity");
  near((base.Phi-legacy.trial.Phi).norm(),0.,"unrestricted Phi compatibility");
  near((base.G_T-legacy.trial.G_T).norm(),0.,"unrestricted Green compatibility");
  Eigen::SelfAdjointEigenSolver<MatType> es(model.frameworkGenerator());
  const double e0=-.25*es.eigenvalues().cwiseAbs().sum();
  const double expected=e0+(base.fermionParity()==chosen.required_parity?0.:es.eigenvalues()[lat.cells()]);
  near(matterEnergy(model,chosen.trial.G_T).real()-expected,2e-10,"lowest physical energy");
  if(chosen.blocked)++blocked;else near((base.Phi-chosen.trial.Phi).norm(),0.,"compatible trial unchanged");
  auto repeated=selectFreeMatterTrial(model);near((repeated.trial.G_T-chosen.trial.G_T).norm(),0.,"repeatable trial");
  std::vector<int>d(model.majoranas(),1);d[0]=-1;d[lat.cells()+1]=-1;d[2]=-1;
  auto gauged=model.gaugeTransformed(d);auto equivalent=selectFreeMatterTrial(gauged);
  MatType green=chosen.trial.G_T;
  for(int i=0;i<model.majoranas();++i)for(int j=0;j<model.majoranas();++j)green(i,j)*=d[i]*d[j];
  near((green-equivalent.trial.G_T).norm(),3e-9,"local gauge covariance");
  require(equivalent.required_parity==cliffordParity(gauged),"transformed physical condition");
  if(shape==std::array<int,3>{{6,6,0}}&&sector==std::array<int,2>{{1,1}}) {
   require(chosen.blocked&&chosen.actual_parity==1,"L6 physical default");
   require(chosen.momentum_n1==0&&chosen.momentum_n2==3,"L6 deterministic momentum");
   // This gauge is uniform. Every AA, BB and AB covariance is invariant
   // under simultaneous translations of the two unit-cell coordinates.
   for(int x=0;x<6;++x)for(int y=0;y<6;++y)for(int a=0;a<2;++a)for(int b=0;b<2;++b) {
    const int i=a*36+lat.cellIndex(x,y),j=b*36+lat.cellIndex((x+2)%6,(y+1)%6);
    near(std::abs(chosen.trial.G_T(i,j)-chosen.trial.G_T(a*36,b*36+lat.cellIndex(2,1))),2e-10,"L6 translation invariant");
   }
  }
  if(shape==std::array<int,3>{{9,9,0}}&&sector==std::array<int,2>{{1,1}})require(!chosen.blocked,"L9 default unmodified");
  std::cout<<"CASE geometry "<<shape[0]<<' '<<shape[1]<<' '<<shape[2]<<' '<<sector[0]<<' '<<sector[1]<<' ';writeMatterTrialSelection(std::cout,chosen);std::cout<<'\n';++count;
 }
 require(blocked>0,"no blocked validation cases");
 bool bad=false;try{selectFreeMatterTrial(KitaevSSHModel(HoneycombLattice(2,2)),"bad");}catch(const std::invalid_argument&){bad=true;}require(bad,"invalid policy rejected");
 bool zero=false;try{selectFreeMatterTrial(KitaevSSHModel(HoneycombLattice(6,6),{1,1,2,0}));}catch(const GaussianTrialStateError&e){zero=e.code()==GaussianTrialStateStatus::zero_mode_ambiguous;}require(zero,"zero mode gate preserved");
 std::cout<<"PASS geometry cases="<<count<<" blocked="<<blocked<<" plus local gauges\n";
}
std::vector<MatType> majoranas() {
 std::vector<MatType>c(8,MatType::Zero(16,16));
 for(int j=0;j<4;++j)for(int n=0;n<16;++n){double s=__builtin_popcount(unsigned(n&((1<<j)-1)))%2?-1.:1.;c[j](n^(1<<j),n)=s;c[4+j](n^(1<<j),n)=DataType(0,(n&(1<<j))?s:-s);}
 return c;
}
MatType fockH(const std::vector<MatType>&c,const std::vector<double>&x) {
 MatType h=MatType::Zero(16,16);
 for(int xx=0;xx<2;++xx)for(int yy=0;yy<2;++yy){int r=2*xx+yy;h+=DataType(0,1+.4*x[3*r])*c[r]*c[4+2*(1-xx)+yy];h+=DataType(0,1+.4*x[3*r+1])*c[r]*c[4+2*xx+1-yy];h+=DataType(0,.5+.4*x[3*r+2])*c[r]*c[4+r];}
 return h;
}
MatType exponential(const MatType&h,double dt){Eigen::SelfAdjointEigenSolver<MatType>e(h);return e.eigenvectors()*(-dt*e.eigenvalues().array()).exp().matrix().asDiagonal()*e.eigenvectors().adjoint();}
void fock() {
 auto c=majoranas();KitaevSSHModel model(HoneycombLattice(2,2),{1,1,.5,.4});auto chosen=selectFreeMatterTrial(model);require(chosen.blocked,"L2 blocked Fock case");
 auto h=fockH(c,std::vector<double>(12,0));MatType parity=MatType::Identity(16,16);
 for(int j=0;j<4;++j)parity=parity*DataType(0,1)*c[j]*c[4+j];
 std::vector<int>even;for(int j=0;j<16;++j)if(parity(j,j).real()>.5)even.push_back(j);
 MatType hp(even.size(),even.size());for(int i=0;i<int(even.size());++i)for(int j=0;j<int(even.size());++j)hp(i,j)=h(even[i],even[j]);
 Eigen::SelfAdjointEigenSolver<MatType>physical(hp);near(physical.eigenvalues()[0]+4.,1e-12,"independent physical Fock ground energy");
 // The selected Gaussian must lie inside this lowest physical eigenspace.
 // Its flattened covariance defines a nondegenerate parent only to select
 // the same vector inside the *independently* computed degenerate eigenspace.
 MatType parent=MatType::Zero(16,16);for(int i=0;i<8;++i)for(int j=i+1;j<8;++j)parent-=chosen.trial.G_T(i,j)*c[i]*c[j];
 Eigen::SelfAdjointEigenSolver<MatType>p(parent);cVecType psi=p.eigenvectors().col(0);
 near(std::abs((parity*psi-psi).norm()),2e-11,"independent Fock parity");
 near(((h*psi)+4.*psi).norm(),2e-11,"state lies in physical ground space");
 MatType target(8,8);for(int i=0;i<8;++i)for(int j=0;j<8;++j)target(i,j)=(i==j?1.:0.)-(psi.adjoint()*c[i]*c[j]*psi)(0);
 near((target-chosen.trial.G_T).norm(),2e-11,"independent Fock trial covariance");
 for(int kind=0;kind<3;++kind) {
  PhononPath path(12,4,.1,{1,1});for(int node=0;node<=4;++node)for(int b=0;b<12;++b)path=path.withCoordinate(node,b,kind==0?0.:.2*std::sin(.3*b+.2*node+kind));
  cVecType right=psi;std::vector<MatType>u;for(int link=0;link<4;++link){u.push_back(exponential(fockH(c,path.node(link)),.05));u.push_back(exponential(fockH(c,path.node(link+1)),.05));}
  for(auto &gate:u)right=gate*right;const DataType exact=(psi.adjoint()*right)(0);
  const auto gen=honeycombGenerator(model);const auto measured=evaluatePath(path,chosen.trial,pathFactors(path,8,gen),2);
  const DataType weight=std::exp(measured.log_abs_matter)*measured.matter_phase;
  near(std::abs(weight/exact-DataType(1)),3e-10,"blocked signed Pfaffian/Fock weight");
  cVecType ket=psi;for(int i=0;i<4;++i)ket=u[i]*ket;
  Eigen::RowVectorXcd bra=psi.adjoint();for(int i=7;i>=4;--i)bra=bra*u[i];
  MatType g(8,8);for(int i=0;i<8;++i)for(int j=0;j<8;++j)g(i,j)=(i==j?1.:0.)-(bra*c[i]*c[j]*ket)(0)/exact;
  near((g-measured.center_green).norm(),3e-10,"blocked Fock center Green");
  for(bool production:{false,true}) {
   PhononWeightOptions opts{production,production,production,production,production,true,production};
   CoupledPhononPath state(path,chosen.trial,gen,2,opts);state.ensureMeasured();
   near(std::abs(std::exp(state.state().weight.log_abs_matter)*state.state().weight.matter_phase/exact-DataType(1)),3e-10,"blocked production/reference weight");
   near((state.state().weight.center_green-g).norm(),3e-10,"blocked production/reference Green");
  }
 }
 std::cout<<"PASS blocked independent16-dimensional Fock parity/energy/weights/Green; reference and production\n";
}
}
int main(){try{std::cout<<std::setprecision(17);geometry();fock();return 0;}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
