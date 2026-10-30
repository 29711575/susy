// Deterministic dense-Fock oracle for mode-count parity and scalar phase branches.
#include "pure_projector_fast.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>

namespace {
int contextN=0,contextBlocked=0,contextBackend=-1,contextCut=-1;double contextBranch=0.;
const char *contextOrder="setup";
void require(bool ok,const char *s){if(!ok){std::ostringstream out;out<<std::setprecision(17)<<s<<" n="<<contextN<<" blocked="<<contextBlocked<<" branch="<<contextBranch<<" backend="<<contextBackend<<" order="<<contextOrder<<" cut="<<contextCut;throw std::runtime_error(out.str());}}
void mismatch(double ew,double eg,DataType expectedWeight,DataType measuredWeight,
              const MatType &expectedGreen,const MatType &measuredGreen) {
 if(ew<2e-9&&eg<2e-9)return;
 Eigen::Index row=0,col=0;(expectedGreen-measuredGreen).cwiseAbs().maxCoeff(&row,&col);
 std::cerr<<std::setprecision(17)<<"MISMATCH order="<<contextOrder<<" cut="<<contextCut
     <<" ew="<<ew<<" eg="<<eg<<" expected_weight="<<expectedWeight
     <<" measured_weight="<<measuredWeight<<" green_row="<<row<<" green_col="<<col
     <<" expected_Green="<<expectedGreen(row,col)<<" measured_Green="<<measuredGreen(row,col)<<'\n';
}
MatType exponential(const MatType &h,double t){
 Eigen::SelfAdjointEigenSolver<MatType> es(h);
 require(es.info()==Eigen::Success,"independent Hermitian exponential");
 return es.eigenvectors()*(t*es.eigenvalues().array()).exp().matrix().asDiagonal()*es.eigenvectors().adjoint();
}
std::vector<MatType> clifford(int n){
 const int d=1<<n;std::vector<MatType> g(2*n,MatType::Zero(d,d));
 for(int j=0;j<n;++j)for(int v=0;v<d;++v){double s=__builtin_popcount(unsigned(v&((1<<j)-1)))%2?-1.:1.;
  g[j](v^(1<<j),v)=s;g[n+j](v^(1<<j),v)=DataType(0,(v&(1<<j))?s:-s);}
 return g;
}
MatType quadratic(const MatType &h,const std::vector<MatType> &g){
 MatType out=MatType::Zero(g[0].rows(),g[0].cols());
 // Production one-body K=-iA; physical H=i gamma^T A gamma/4.
 // Thus H=-sum_{a<b} K_ab gamma_a gamma_b/2.
 for(int a=0;a<h.rows();++a)for(int b=a+1;b<h.cols();++b)out-=.5*h(a,b)*g[a]*g[b];
 return out;
}
}
int main(){try{
 mkl_set_num_threads(1);int cases=0;double worstWeight=0,worstGreen=0,worstWick=0,worstEndpoint=0;
 const double pi=std::acos(-1.);
 for(int n:{1,2,3,5})for(int blocked:{0,1})for(double branch:{-pi-1e-7,-pi+1e-7,0.,pi-1e-7,pi+1e-7}){
  contextN=n;contextBlocked=blocked;contextBranch=branch;contextBackend=-1;contextOrder="setup";contextCut=-1;
  auto g=clifford(n);MatType ht=MatType::Zero(2*n,2*n);
  const int m=2*n,dimension=1<<n;
  std::vector<MatType> pairs(m*m);
  for(int a=0;a<m;++a)for(int b=0;b<m;++b)pairs[a*m+b]=g[a]*g[b];
  for(int a=0;a<m;++a)for(int b=0;b<m;++b)
   require((pairs[a*m+b]+pairs[b*m+a]-(a==b?2.:0.)*MatType::Identity(dimension,dimension)).norm()<1e-13,"CAR algebra");
  for(int j=0;j<n;++j){ht(j,n+j)=DataType(0,(j==0&&blocked?-1.:1.)*(1.+.17*j));ht(n+j,j)=-ht(j,n+j);}
  auto trial=GaussianTrialState::fromMajoranaHamiltonian(ht);
  Eigen::SelfAdjointEigenSolver<MatType> es(quadratic(ht,g));require(es.info()==Eigen::Success,"trial ED");
  require(es.eigenvalues()[1]-es.eigenvalues()[0]>1e-5,"unique trial required");
  cVecType psi=es.eigenvectors().col(0),right=psi;std::vector<PureProjectorSlice> slices;
  std::vector<MatType> fockFactors;
  for(int step=0;step<5;++step){MatType h=MatType::Zero(2*n,2*n),fockH=MatType::Zero(1<<n,1<<n);
   // Both representations start from physical scalar couplings. In particular
   // the oracle never reconstructs its Hamiltonian from a production propagator.
   for(int a=0;a<2*n;++a)for(int b=a+1;b<2*n;++b){
    const double coupling=.4*std::sin(.37*(a+1)*(b+2)+.83*step);
    h(a,b)=DataType(0,coupling);h(b,a)=-h(a,b);
    fockH+=DataType(0,-.5*coupling)*pairs[a*m+b];}
   const DataType eta=std::exp(DataType(0,step==2?branch:0.));
   slices.emplace_back(exponential(h,-.23),eta,"independent-Fock");
   fockFactors.push_back(eta*exponential(fockH,-.23));
   right=(fockFactors.back()*right).eval();
  }
  const DataType exact=(psi.adjoint()*right)(0);require(std::abs(exact)>1e-6,"resolved overlap required");
  Eigen::RowVectorXcd left=psi.adjoint();
  for(auto it=fockFactors.rbegin();it!=fockFactors.rend();++it)left=(left*(*it)).eval();
  require(std::abs((left*psi)(0)/exact-DataType(1))<2e-12,"left/right Fock product order differs");
  std::vector<Eigen::RowVectorXcd> leftPair0(m*m),leftPairN(m*m);
  std::vector<cVecType> rightPair0(m*m),rightPairN(m*m);
  MatType cov0(m,m),covN(m,m);for(int a=0;a<m;++a)for(int b=0;b<m;++b){
   const int index=a*m+b;
   leftPair0[index]=left*pairs[index];rightPair0[index]=pairs[index]*psi;
   leftPairN[index]=psi.adjoint()*pairs[index];rightPairN[index]=pairs[index]*right;
   cov0(a,b)=(leftPair0[index]*psi)(0)/exact;
   covN(a,b)=(leftPairN[index]*right)(0)/exact;}
  // The legacy weight core consumes reversed chronological Hermitian factors.
  // Scalar eta phases are reordered, never conjugated. Its final right frame
  // is B0...BN Phi, which is the physical left frame at canonical cut zero.
  // Therefore its Green already equals I-cov0; taking another adjoint is wrong.
  const std::vector<PureProjectorSlice> weightOrder(slices.rbegin(),slices.rend());
  const MatType expected0=MatType::Identity(m,m)-cov0;
  for(int backend:{0,1,2}){
   contextBackend=backend;contextOrder="reversed-core";contextCut=0;
   PureProjectorWeightResult result;
   if(backend==0)result=PureProjectorWeightEvaluator().evaluate(trial,weightOrder);
   else {PureFastOptions opts;opts.accelerated_linear_algebra=backend==2;result=pureProjectorStableReferenceWeight(trial,weightOrder,opts);}
   require(result.ok(),"weight backend rejected regular oracle");
   const DataType measured=std::exp(result.log_abs_weight)*result.complex_phase;
   double ew=std::abs(measured/exact-DataType(1));
   double eg=(result.green-expected0).cwiseAbs().maxCoeff();
   worstWeight=std::max(worstWeight,ew);worstGreen=std::max(worstGreen,eg);
   mismatch(ew,eg,exact,measured,expected0,result.green);
   require(ew<2e-9&&eg<2e-9,"Fock weight/Green mismatch");
  }
  // Endpoint reconstruction keeps the original physical action order. Its
  // core Green requires the physical adjoint, here at the opposite endpoint.
  contextBackend=-1;contextOrder="chronological-endpoint-adjoint";contextCut=int(slices.size());
  const auto endpoint=pureProjectorEndpointRebuild(trial,slices,contextCut,2);
  require(endpoint.ok(),"physical endpoint rejected regular oracle");
  const MatType expectedN=MatType::Identity(m,m)-covN,measuredN=endpoint.green.adjoint();
  const double endpointWeightError=std::abs(endpoint.log_abs_weight-std::log(std::abs(exact)));
  const double endpointGreenError=(measuredN-expectedN).cwiseAbs().maxCoeff();
  worstEndpoint=std::max(worstEndpoint,endpointGreenError);
  mismatch(endpointWeightError,endpointGreenError,DataType(std::abs(exact),0.),
      DataType(std::exp(endpoint.log_abs_weight),0.),expectedN,measuredN);
  require(endpointWeightError<2e-9&&endpointGreenError<2e-9,"physical endpoint Fock mismatch");
  contextOrder="chronological-Fock-four-point";
  for(int cut:{0,int(slices.size())}) {
   contextCut=cut;
   const MatType &cov=cut==0?cov0:covN;
   const auto &leftPair=cut==0?leftPair0:leftPairN;
   const auto &rightPair=cut==0?rightPair0:rightPairN;
   for(int a=0;a<m;++a)for(int b=0;b<m;++b)for(int c=0;c<m;++c)for(int d=0;d<m;++d){
    DataType direct=(leftPair[a*m+b]*rightPair[c*m+d])(0)/exact;
    double err=std::abs(direct-(cov(a,b)*cov(c,d)-cov(a,c)*cov(b,d)+cov(a,d)*cov(b,c)));
    worstWick=std::max(worstWick,err);require(err<2e-9,"direct Fock four-point mismatch");
   }
  }
  ++cases;
 }
 std::cout<<std::setprecision(17)<<"PASS cases="<<cases<<" modes=1,2,3,5 blocked=0,1 backends=3 core_order=reversed core_cut=0 endpoint_order=chronological endpoint_cut=N max_weight="<<worstWeight<<" max_Green="<<worstGreen<<" max_endpoint_Green="<<worstEndpoint<<" max_direct_four_point="<<worstWick<<'\n';return 0;
}catch(const std::exception &e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
