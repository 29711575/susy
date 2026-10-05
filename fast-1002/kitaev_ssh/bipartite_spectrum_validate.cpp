#include "phonon_projector.h"
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace kitaev_ssh;
namespace {
void require(bool b,const char *s){if(!b)throw std::runtime_error(s);}
double error(const MatType&a,const MatType&b){return (a-b).cwiseAbs().maxCoeff();}
std::vector<MatType> majoranas(){
 std::vector<MatType> c(8,MatType::Zero(16,16));
 for(int a=0;a<4;++a)for(int n=0;n<16;++n){double s=__builtin_popcount(unsigned(n&((1<<a)-1)))%2?-1.:1.;c[a](n^(1<<a),n)=s;c[4+a](n^(1<<a),n)=DataType(0,(n&(1<<a))?s:-s);}return c;
}
MatType hamiltonian(const std::vector<MatType>&c,const std::vector<double>&x,double g){
 MatType h=MatType::Zero(16,16);
 for(int a=0;a<2;++a)for(int b=0;b<2;++b){int i=2*a+b;
  h+=DataType(0,1+g*x[3*i])*c[i]*c[4+2*(1-a)+b];
  h+=DataType(0,1+g*x[3*i+1])*c[i]*c[4+2*a+1-b];
  h+=DataType(0,.5+g*x[3*i+2])*c[i]*c[4+i];}return h;
}
MatType fockExp(const MatType&h,double t){
 Eigen::SelfAdjointEigenSolver<MatType> e(h);require(e.info()==Eigen::Success,"Fock solver");
 return e.eigenvectors()*(-t*e.eigenvalues().array()).exp().matrix().asDiagonal()*e.eigenvectors().adjoint();
}
}
int main(){try{
 std::mt19937_64 rng(10022026);BipartiteSpectrum reduced;
 double maxExp=0.,maxCayley=0.,maxLog=0.,maxPhase=0.,maxGreen=0.;int cases=0;
 for(int L:{2,6,9,12})for(int shear:{0,1}){
  HoneycombLattice lattice(L,L,shear);Couplings c;c.Jz=.5;c.g=1.2;
  auto gauge=lattice.zeroFluxGauge();
  for(const auto &b:lattice.bonds())if((b.a%3==0)!=(b.b%3==0))gauge[b.index]*=-1;
  std::vector<double>x(lattice.bondCount());std::normal_distribution<double> normal;
  for(auto &v:x)v=normal(rng);
  KitaevSSHModel model(lattice,c,x,gauge);MatType k=model.frameworkGenerator();
  require(reduced.compute(k),"eligible bipartite generator rejected");
  auto reference=selfAdjointSpectral(k);require(reference.ok,"complex reference spectrum");
  for(double a:{.05,.1,.2}){
   MatType b,cy;reduced.exponential(a,b);reduced.cayley(a,cy);
   cVecType f=(-a*reference.D.array()).exp().cast<DataType>();
   MatType expected=reference.V*f.asDiagonal()*reference.V.adjoint();
   maxExp=std::max(maxExp,error(b,expected));
   f=(.5*a*reference.D.array()).tanh().cast<DataType>();
   expected=reference.V*f.asDiagonal()*reference.V.adjoint();
   maxCayley=std::max(maxCayley,error(cy,expected));
   require((b-b.adjoint()).norm()<2e-12*std::max(1.,b.norm()),"half-step Hermiticity");
   ++cases;
  }
  BipartiteSpectrum workspace;PhononNodeSpectrum node;
  require(node.compute(k,true,workspace)&&node.reduced()==(k.rows()/2>=2),"production dimension dispatch");
  require(node.compute(k,false,workspace)&&!node.reduced(),"generic caller dispatch");
  MatType wrong=k;wrong(0,0)=1e-30;require(!reduced.compute(wrong),"same-sublattice coefficient discarded");
  wrong=k;wrong(0,L*L)+=DataType(1e-30,0.);require(!reduced.compute(wrong),"real cross coefficient discarded");
 }
 MatType zero=MatType::Zero(72,72),b;
 require(reduced.compute(zero),"zero hopping spectrum");reduced.exponential(.1,b);
 require(error(b,MatType::Identity(72,72))==0.,"zero generator exponential");
 // Independent many-body oracle with noncommuting, non-mirrored node fields.
 const auto c=majoranas();KitaevSSHModel model(HoneycombLattice(2,2),p2Couplings(.4));
 const auto trial=freeGroundState(model);
 Eigen::SelfAdjointEigenSolver<MatType> initial(hamiltonian(c,std::vector<double>(12,0.),.4));
 require(initial.info()==Eigen::Success,"Fock initial state");const cVecType psi=initial.eigenvectors().col(0);
 for(int repeat=0;repeat<3;++repeat){
  auto path=PhononPath(12,4,.13,{1.7,.8}).freeSample(rng);
  std::vector<PureProjectorSlice> factors(8);
  for(int j=0;j<path.nodes();++j){auto k=model.frameworkGeneratorWith(path.node(j));
   require(reduced.compute(k),"small Fock reduced spectrum");reduced.exponential(path.dt()/2,b);
   for(int f:nodeFactors(j,path.intervals()))factors[f]=PureProjectorSlice(b,1.,"oracle_node");}
  const auto w=evaluatePath(path,trial,factors,4);
  cVecType right=psi;Eigen::RowVectorXcd left=psi.adjoint();
  for(int l=1;l<=path.center();++l)right=fockExp(hamiltonian(c,path.node(l),.4),path.dt()/2)*fockExp(hamiltonian(c,path.node(l-1),.4),path.dt()/2)*right;
  for(int l=path.intervals();l>path.center();--l)left=left*fockExp(hamiltonian(c,path.node(l),.4),path.dt()/2)*fockExp(hamiltonian(c,path.node(l-1),.4),path.dt()/2);
  DataType amp=(left*right)(0);maxLog=std::max(maxLog,std::abs(w.log_abs_matter-std::log(std::abs(amp))));
  maxPhase=std::max(maxPhase,std::abs(w.matter_phase-amp/std::abs(amp)));
  MatType green=MatType::Zero(8,8);
  for(int i=0;i<8;++i)for(int j=0;j<8;++j)if(i!=j)green(i,j)=-(left*c[i]*c[j]*right)(0)/amp;
  maxGreen=std::max(maxGreen,error(w.center_green,green));
 }
 std::cout<<std::setprecision(17)<<"{\"status\":\"PASS\",\"matrix_cases\":"<<cases<<",\"fock_cases\":3,\"max_exponential_error\":"<<maxExp<<",\"max_cayley_error\":"<<maxCayley<<",\"max_fock_log_error\":"<<maxLog<<",\"max_fock_phase_error\":"<<maxPhase<<",\"max_fock_green_error\":"<<maxGreen<<"}\n";
 require(maxExp<2e-12&&maxCayley<2e-12,"spectral matrix reference tolerance");
 require(maxLog<2e-12&&maxPhase<2e-12&&maxGreen<2e-12,"independent Fock reference tolerance");
}catch(const std::exception&e){std::cerr<<"bipartite_spectrum_validate: "<<e.what()<<'\n';return 1;}}
