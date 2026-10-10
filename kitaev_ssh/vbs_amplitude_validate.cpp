#include "vbs_amplitude.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
using namespace kitaev_ssh;
namespace {
void require(bool ok,const char *message) {if(!ok)throw std::runtime_error(message);}
Complex direct(const HoneycombLattice &l,const VBSObservables &v,const PhononPath &p,int j,int q,int mu) {
    Complex out=0.;const auto m=v.momenta()[q];
    for(const auto &b:l.bonds())if(static_cast<int>(b.type)==mu)
        out+=p.coordinate(j,b.index)*std::polar(1.,-2*pi*(m.reciprocal_x*b.cell.x+m.reciprocal_y*b.cell.y));
    return out/double(l.cells());
}
}
int main() {
    try {
        double maxDftError=0.,maxSourceError=0.;
        for(int shear:{0,1}) {
            HoneycombLattice lattice(3,2,shear);VBSObservables vbs(lattice);
            PhononPath path(lattice.bondCount(),6,.17,PhononParameters(4.,.5));
            for(int b=0;b<path.bonds();++b) {
                std::vector<double> x;
                for(int j=0;j<path.nodes();++j)x.push_back(.3*std::sin(.43*b+.71*j)+.2*std::cos(.7*b-.11*j));
                path=path.withBondCoordinates(b,x);
            }
            const VBSAmplitude estimator(vbs,true);const auto r=estimator.measure(path);
            require(r.size()==15*std::size_t(lattice.cells()) && r.size()==estimator.names().size(),"raw schema");
            const auto equal=VBSAmplitude(vbs).measure(path);
            for(int q=0;q<lattice.cells();++q)for(int mu=0;mu<3;++mu) {
                const double a=std::abs(direct(lattice,vbs,path,path.center(),q,mu));
                double integral=0.;
                // Sum over links independently, half of each link at either end.
                for(int j=0;j<path.intervals();++j)
                    integral+=.5*path.dt()*(std::abs(direct(lattice,vbs,path,j,q,mu))+
                                           std::abs(direct(lattice,vbs,path,j+1,q,mu)));
                const double ref[]={a,a*a,a*a*a*a,integral,a*integral};
                for(int k=0;k<5;++k)maxDftError=std::max(maxDftError,std::abs(r[15*q+3*k+mu]-ref[k]));
                for(int k=0;k<3;++k)require(equal[9*q+3*k+mu]==r[15*q+3*k+mu],"equal-only schema");
            }
            require(maxDftError<3e-13,"independent DFT/path quadrature");
        }
        // Rotate a non-self-conjugate D through a full circle in imaginary
        // time. |integral D| vanishes while integral |D| stays finite.
        {
            HoneycombLattice l(3,2);VBSObservables v(l);
            PhononPath p(l.bondCount(),6,.2,PhononParameters(4.,.5));
            for(const auto &b:l.bonds()) {
                std::vector<double> x;
                for(int j=0;j<p.nodes();++j)x.push_back(static_cast<int>(b.type)==2?
                    2*std::cos(2*pi*b.cell.x/3+2*pi*j/p.intervals()):0.);
                p=p.withBondCoordinates(b.index,x);
            }
            const int q=2;const auto r=VBSAmplitude(v,true).measure(p);
            Complex linear=0.;for(int j=0;j<p.nodes();++j)linear+=VBSAmplitude::nodeWeight(p,j)*direct(l,v,p,j,q,2);
            require(std::abs(linear)<1e-13 && std::abs(r[15*q+11]-p.totalTau())<1e-13,"modulus before integral");
        }
        // Exact finite-dt source reweighting of a finite complex-weight
        // ensemble, with fixed endpoints. Differentiate normalized <A>.
        {
            const double a[]={.2,.7,1.3,.4},i[]={.9,.4,1.7,.6};
            const Complex w[]={Complex(1.,.2),Complex(-.13,.07),Complex(.8,-.11),Complex(.3,.04)};
            Complex z=0.,sa=0.,si=0.,sai=0.;const double nc=6.;
            for(int k=0;k<4;++k){z+=w[k];sa+=w[k]*a[k];si+=w[k]*i[k];sai+=w[k]*a[k]*i[k];}
            const Complex covariance=nc*(sai/z-(sa/z)*(si/z));
            auto sourceMean=[&](double eps) {
                Complex num=0.,den=0.;
                for(int k=0;k<4;++k){const Complex weight=w[k]*std::exp(eps*nc*i[k]);den+=weight;num+=weight*a[k];}
                return num/den;
            };
            for(double eps:{1e-3,1e-4,1e-5}) {
                const double err=std::abs((sourceMean(eps)-sourceMean(-eps))/(2*eps)-covariance);
                if(eps==1e-5)maxSourceError=err;
                require(err<20*eps*eps+2e-10,"nonlinear source derivative");
            }
        }
        // Free stationary phonon vacuum: non-self q has a Rayleigh modulus
        // (B4=2), self-conjugate q a half-normal modulus (B4=3).
        {
            HoneycombLattice l(3,2);VBSObservables v(l);VBSAmplitude e(v);
            PhononPath p(l.bondCount(),2,.2,PhononParameters(4.,.5));std::mt19937_64 rng(831725);
            constexpr int samples=60000;double sums[2][3]{};
            for(int k=0;k<samples;++k) {
                const auto r=e.measure(p.freeSample(rng));
                for(int t=0;t<2;++t)for(int j=0;j<3;++j)sums[t][j]+=r[9*(t==0?2:1)+3*j+2].real();
            }
            const double variance=1./(2*4.*.5*l.cells());
            for(int t=0;t<2;++t) {
                const double mean=sums[t][0]/samples,m2=sums[t][1]/samples,m4=sums[t][2]/samples;
                const double refMean=std::sqrt(variance)*(t==0?std::sqrt(pi)/2:std::sqrt(2/pi));
                const double refBinder=t==0?2.:3.;
                require(std::abs(mean/refMean-1)<.025 && std::abs(m2/variance-1)<.035,"free amplitude moments");
                require(std::abs(m4/(m2*m2)-refBinder)<.09,"free Binder limit");
                const double fluct=l.cells()*(m2-mean*mean);
                const double refFluct=l.cells()*variance*(t==0?1-pi/4:1-2/pi);
                require(std::abs(fluct/refFluct-1)<.04,"free equal-time modulus variance");
                std::cout<<"free "<<(t==0?"complex":"real")<<" A="<<mean<<" VA="<<fluct<<" B4="<<m4/(m2*m2)<<'\n';
            }
        }
        std::cout<<std::setprecision(16)<<"PASS max_direct_DFT_error="<<maxDftError
                 <<" max_modulus_source_derivative_error="<<maxSourceError<<'\n';
        return 0;
    } catch(const std::exception &e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
