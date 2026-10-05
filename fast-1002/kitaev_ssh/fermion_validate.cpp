#include "fermion_observables.h"
#include <iostream>
using namespace kitaev_ssh;
void check(Complex a,Complex b,double tol,const char *label){if(std::abs(a-b)>tol)throw std::runtime_error(label);}
int main(){try{
 const KitaevSSHModel model(HoneycombLattice(4,4),{1,1,.5,0});
 const auto &lat=model.lattice(); FermionObservables obs(model);
 const auto green=freeGroundState(model).G_T; const auto values=obs.measure(green);
 const double pi=std::acos(-1.); double maxerr=0;
 for(std::size_t n=0;n<values.size();++n){const auto &c=obs.channels()[n]; Complex expected=0.;
  if(c.alpha==c.beta)expected=(c.dx==0&&c.dy==0)?1.:0.;
  else for(int kx=0;kx<4;++kx)for(int ky=0;ky<4;++ky){double qx=2*pi*kx/4,qy=2*pi*ky/4;
   Complex f=.5+std::exp(Complex(0,-qx))+std::exp(Complex(0,-qy));
   expected+=std::exp(Complex(0,-qx*c.dx-qy*c.dy))*(c.alpha==0?Complex(0,1)*f:Complex(0,-1)*std::conj(f))/std::abs(f)/16.;
  }
  if(c.string)expected*=((c.dx+c.dy+(c.alpha==1&&c.beta==0))%2?-1.:1.);
  maxerr=std::max(maxerr,std::abs(values[n]-expected));check(values[n],expected,3e-12,"Bloch Fourier mismatch");
 }
 MatType arbitrary=MatType::Zero(32,32);for(int i=0;i<32;++i)for(int j=i+1;j<32;++j){arbitrary(i,j)=Complex(std::sin(i+3*j),std::cos(2*i-j));arbitrary(j,i)=-arbitrary(i,j);}
 auto measured=obs.measure(arbitrary);std::vector<int>d(32);MatType dg=arbitrary;
 for(int i=0;i<32;++i)d[i]=((i*i+3*i+1)%7<3?-1:1);
 for(int i=0;i<32;++i)for(int j=0;j<32;++j)dg(i,j)*=d[i]*d[j];
 auto transformed=model.gaugeTransformed(d);auto gauged=FermionObservables(transformed).measure(dg);double bareChange=0;
 for(std::size_t n=0;n<measured.size();++n){const auto &c=obs.channels()[n];Complex ref=0.;
  for(int x=0;x<4;++x)for(int y=0;y<4;++y){int i=lat.siteIndex(x,y,static_cast<Sublattice>(c.alpha)),j=lat.siteIndex((x+c.dx)%4,(y+c.dy)%4,static_cast<Sublattice>(c.beta));
   auto path=obs.pathFor(x,y,c);if(lat.pathEnd(path)!=j)throw std::runtime_error("path endpoint");
   ref+=(c.string?gaugeStringGreen(model,arbitrary,path):(i==j?Complex(1.):Complex(-arbitrary(i,j))))/16.;
  }check(measured[n],ref,3e-14,"translation average mismatch");
  if(c.string)check(measured[n],gauged[n],3e-14,"gauge invariance mismatch");else bareChange=std::max(bareChange,std::abs(measured[n]-gauged[n]));
 }
 if(bareChange<1e-3)throw std::runtime_error("gauge test vacuous");
 std::cout<<"PASS 128 complex channels: analytic Bloch, translations, PBC path endpoints, directed signs, gauge invariance; max Bloch error="<<maxerr<<'\n';
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
