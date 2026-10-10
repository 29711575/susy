// Full odd-parity Hilbert-space imaginary-time evolution, no QMC or HS.
// For even L, PBC, t=Delta=1, mu=0, the occupation gauge
// d(n)=i^(sum_even n_j) gives K=-sum_j X_j X_{j+1}.
// Density observables are invariant under this diagonal gauge.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using Vec=std::vector<double>;
struct Model {
 int L,N; std::vector<int> masks; Vec interaction,spi,sq;
 Model(int l):L(l),N(1<<(l-1)),interaction(N),spi(N),sq(N) {
  for(int j=0;j<L;++j) masks.push_back(((1<<j)^(1<<((j+1)%L)))>>1);
  for(int i=0;i<N;++i) {
   unsigned n=(i<<1)|((__builtin_popcount(unsigned(i))&1)^1);
   double a=0; std::complex<double> b=0;
   for(int j=0;j<L;++j) {
    double d=((n>>j)&1)-.5,e=((n>>((j+1)%L))&1)-.5;
    interaction[i]+=d*e; a+=(j%2?-1:1)*d;
    b+=d*std::polar(1., (std::acos(-1.)-2*std::acos(-1.)/L)*j);
   }
   spi[i]=a*a/(L*L); sq[i]=std::norm(b)/(L*L);
  }
 }
 Vec initial() const {return Vec(N,1/std::sqrt(double(N)));}
 void H(const Vec& x,Vec& y,double v) const {
  for(int i=0;i<N;++i)y[i]=v*interaction[i]*x[i];
  for(int m:masks)for(int i=0;i<N;++i)y[i]-=x[i^m];
 }
 static void normalize(Vec& x) { double n=0;for(double a:x)n+=a*a;n=std::sqrt(n);if(!std::isfinite(n)||n==0)throw std::runtime_error("norm");for(double& a:x)a/=n; }
 void kinetic(Vec& x,double h) const {
  // All XX bonds commute; each exact two-state gate is cosh(h)+sinh(h)XX.
  double c=std::cosh(h),s=std::sinh(h);
  for(int m:masks)for(int i=0;i<N;++i)if(i<(i^m)){
   int j=i^m;double a=x[i],b=x[j];x[i]=c*a+s*b;x[j]=s*a+c*b;
  }
 }
 void trotter(Vec& x,double v,double h) const {
  kinetic(x,h/2);for(int i=0;i<N;++i)x[i]*=std::exp(-h*v*interaction[i]);
  kinetic(x,h/2);normalize(x);
 }
 void continuous(Vec& x,double v,double h) const {
  // Converged exponential action, Taylor degree selected by term norm.
  // Small steps avoid overflow; halving h is checked separately at L=18.
  Vec term=x,sum=x,next(N);bool converged=false;
  for(int k=1;k<=80;++k) {
   H(term,next,v);double a=0,b=0;
   for(int i=0;i<N;++i){term[i]=-h*next[i]/k;sum[i]+=term[i];a+=term[i]*term[i];b+=sum[i]*sum[i];}
   if(k>5 && a<1e-32*b){converged=true;break;}
  }
  if(!converged)throw std::runtime_error("exponential not converged");
  x.swap(sum);normalize(x);
 }
 Vec measure(const Vec& x,double v)const {
  Vec y(N);H(x,y,v);double e=0,a=0,b=0;
  for(int i=0;i<N;++i){e+=x[i]*y[i];a+=x[i]*x[i]*spi[i];b+=x[i]*x[i]*sq[i];}
  return {e/L,a,b,1-b/a};
 }
};
using C=std::complex<double>;
std::pair<unsigned,C> gamma(int L,int k,unsigned n){int j=k%L;double s=(__builtin_popcount(n&((1u<<j)-1))&1)?-1:1;return {n^(1u<<j),k<L?C(s,0):C(0,(n&(1u<<j))?s:-s)};}
void checkGauge(int L){
 // Independently assemble original Majorana generator convention, including PBC.
 std::vector<C>A(4*L*L);int D=2*L;
 auto add=[&](int i,int j,C v){A[i*D+j]+=v;A[j*D+i]-=v;};
 for(int j=0;j<L-1;++j){for(int k=0;k<2;++k)add(j+k*L,j+1+k*L,C(0,j%2?-1:1));add(j,j+1,C(0,1));add(j+L,j+1+L,C(0,-1));}
 add(2*L-1,L,C(0,-2));
 auto phase=[&](unsigned n){int a=0;for(int j=0;j<L;j+=2)a+=(n>>j)&1;return std::pow(C(0,1),a);};
 for(unsigned n=0;n<(1u<<L);++n)if(__builtin_popcount(n)&1){
  std::vector<C> col(1<<L);
  for(int i=0;i<D;++i)for(int j=i+1;j<D;++j)if(std::abs(A[i*D+j])>0){auto b=gamma(L,j,n);auto a=gamma(L,i,b.first);col[a.first]+=-.5*A[i*D+j]*a.second*b.second*std::conj(phase(a.first))*phase(n);}
  for(int j=0;j<L;++j)col[n^(1u<<j)^(1u<<((j+1)%L))]+=1.;
  for(C z:col)if(std::abs(z)>1e-12)throw std::runtime_error("original Majorana gauge mismatch");
 }
}
