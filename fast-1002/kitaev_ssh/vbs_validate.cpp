#include "vbs_observables.h"
#include <iomanip>
#include <iostream>
#include <map>
#include <random>

using namespace kitaev_ssh;
namespace {
std::map<std::string,double> errors;
// O(1000) floating operations on O(1) fields: 5e-12 absolute bound,
// conservative roundoff allowance, not a statistical acceptance threshold.
void check(const std::string &name,double error,double tolerance=5e-12) {
    errors[name]=std::max(errors[name],error);
    if(!std::isfinite(error) || error>tolerance)throw std::runtime_error(name+" failed");
}
void fields() {
    const HoneycombLattice lattice(5,6); const VBSObservables vbs(lattice);
    const int nc=lattice.cells(); const double v[3]={1.,-.4,.7};
    std::vector<double> x(3*nc);
    for(int r=0;r<nc;++r)for(int mu=0;mu<3;++mu)x[3*r+mu]=v[mu];
    auto values=vbs.measure(x);
    for(int q=0;q<nc;++q)for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)
        check("uniform",std::abs(values[9*q+3*mu+nu]-(q==0?nc*v[mu]*v[nu]:0.)));
    for(int r=0;r<nc;++r) {
        const auto c=lattice.cellCoordinate(r);
        for(int mu=0;mu<3;++mu)x[3*r+mu]=v[mu]*std::cos(2*pi*(double(c.x)/5+2.*c.y/6)+.31);
    }
    values=vbs.measure(x);
    for(int q=0;q<nc;++q)for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)
        check("allowed_plane_wave",std::abs(values[9*q+3*mu+nu]-((q==8||q==28)?nc*v[mu]*v[nu]/4:0.)));
    std::mt19937_64 rng(7509);std::normal_distribution<double> normal;
    for(auto &a:x)a=normal(rng);
    const auto fft=vbs.fourier(x),direct=vbs.directFourier(x);
    for(int q=0;q<nc;++q)for(int mu=0;mu<3;++mu)check("FFT_direct",std::abs(fft[q][mu]-direct[q][mu]));
    values=vbs.measure(x); Complex reciprocal=0.;double real=0.;
    for(int q=0;q<nc;++q)for(int mu=0;mu<3;++mu)reciprocal+=values[9*q+4*mu];
    for(auto a:x)real+=a*a;
    check("Parseval",std::abs(reciprocal-real));
    std::vector<double> translated(x.size());
    for(int r=0;r<nc;++r) {
        const auto c=lattice.cellCoordinate(r);const int rr=lattice.cellIndex((c.x+1)%5,(c.y+2)%6);
        for(int mu=0;mu<3;++mu)translated[3*r+mu]=x[3*rr+mu];
    }
    const auto shifted=vbs.measure(translated);
    for(int j=0;j<9*nc;++j)check("translation_invariance",std::abs(values[j]-shifted[j]));
    const auto derived=vbs.derive(values);
    for(int d=0;d<nc;++d) {
        const auto dc=lattice.cellCoordinate(d);
        for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu) {
            double correlation=0.;
            for(int r=0;r<nc;++r) {
                const auto c=lattice.cellCoordinate(r);
                const int rr=lattice.cellIndexWrapped(
                    static_cast<std::int64_t>(c.x)+dc.x,
                    static_cast<std::int64_t>(c.y)+dc.y);
                correlation+=x[3*r+mu]*x[3*rr+nu]/nc;
            }
            check("inverse_Fourier_direct_correlation",std::abs(derived[4*nc+9*d+3*mu+nu]-correlation));
        }
    }
}
void freeCovariance() {
    for(int l:{2,4,6}) {
        const HoneycombLattice lattice(l,l);const VBSObservables vbs(lattice);
        const int nc=lattice.cells(),nb=lattice.bondCount();const double variance=1/(2*1.7*.8);
        // Exact second-moment quadrature of the independent Gaussian vacuum:
        // equal +/- basis impulses reproduce <Xi Xj>=variance*delta_ij.
        PhaseStatistics stats(vbs.names(),2);
        for(int b=0;b<nb;++b)for(int sign:{-1,1}) {
            std::vector<double>x(nb,0.);x[b]=sign*std::sqrt(nb*variance);
            stats.add(1.,vbs.measure(x));
        }
        const auto summary=stats.summary();
        for(int q=0;q<nc;++q)for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)
            check("g0_all_q_exact_Gaussian_covariance",std::abs(summary.observables[9*q+3*mu+nu].value-(mu==nu?variance:0.)));
        for(int i=9*nc;i<12*nc;++i)check("g0_zero_mean",std::abs(summary.observables[i].value));
    }
}
void jointJackknife() {
    const VBSObservables vbs(HoneycombLattice(2,3));const int nc=6;
    PhaseStatistics stats(vbs.names(),4);
    std::vector<std::vector<Complex>> measurements;std::vector<Complex> phases;
    std::vector<std::array<double,2>> direct; // real-space Czz(0), mean Xz
    for(int t=0;t<64;++t) {
        std::vector<double>x(3*nc);
        for(int j=0;j<3*nc;++j)x[j]=.4*std::sin(.71*j+.13*t)+.02*t;
        double xz=0.,zz=0.;for(int r=0;r<nc;++r) { xz+=x[3*r+2]/nc;zz+=x[3*r+2]*x[3*r+2]/nc; }
        direct.push_back({zz,xz});
        const auto phase=std::polar(1.,.8*std::sin(.29*t));const auto m=vbs.measure(x);
        phases.push_back(phase);measurements.push_back(m);stats.add(phase,m);
    }
    const auto result=summarizeDerivedBlocks(stats.blocks(),[&](const auto &r){return vbs.derive(r);});
    // Independent sample-level exclusion and direct real-space subtraction.
    auto reference=[&](int omitted) {
        Complex phase=0.,xx=0.,x=0.;
        for(int t=0;t<64;++t)if(t/4!=omitted) { phase+=phases[t];xx+=phases[t]*direct[t][0];x+=phases[t]*direct[t][1]; }
        return xx/phase-(x/phase)*(x/phase);
    };
    const auto &estimate=result.observables[4*nc+9*nc+8];
    check("joint_complex_connected_ratio",std::abs(estimate.value-reference(-1)));
    std::vector<Complex> replicas;Complex mean=0.;
    for(int b=0;b<16;++b){replicas.push_back(reference(b));mean+=replicas.back()/16.;}
    double vr=0.,vi=0.,cov=0.;
    for(auto z:replicas){z-=mean;vr+=z.real()*z.real()*15/16;vi+=z.imag()*z.imag()*15/16;cov+=z.real()*z.imag()*15/16;}
    check("joint_jackknife_real_error",std::abs(estimate.error_real-std::sqrt(vr)));
    check("joint_jackknife_imag_error",std::abs(estimate.error_imag-std::sqrt(vi)));
    check("joint_jackknife_covariance",std::abs(estimate.covariance-cov));
    const auto coarse=summarizeDerivedBlocks(mergeAdjacentBlocks(stats.blocks()),[&](const auto &r){return vbs.derive(r);});
    for(std::size_t j=0;j<result.observables.size();++j)
        check("blocking_same_total_ratio",std::abs(result.observables[j].value-coarse.observables[j].value));
    if(!result.phase_resolved || !result.jackknife_defined)throw std::runtime_error("test phase unexpectedly unresolved");
    PhaseStatistics singular(vbs.names(),1);
    singular.add(1.,measurements[0]);singular.add(-1.,measurements[1]);
    if(summarizeDerivedBlocks(singular.blocks(),[&](const auto &r){return vbs.derive(r);}).ratio_defined)
        throw std::runtime_error("singular phase ratio accepted");

    // Eigenvalue of the ensemble matrix is NOT mean(config eigenvalue).
    std::vector<double> x(3*nc,0.),y(3*nc,0.);
    for(int r=0;r<nc;++r){x[3*r]=1.;y[3*r+1]=1.;}
    auto mx=vbs.measure(x),my=vbs.measure(y);
    for(std::size_t j=0;j<mx.size();++j)mx[j]=(mx[j]+my[j])*.5;
    check("eigenvalue_after_ensemble_average",std::abs(vbs.derive(mx)[2]-double(nc)/2));
    std::vector<Complex> noisy(12*nc,0.);
    for(int q=0;q<nc;++q)for(int mu=0;mu<3;++mu)noisy[9*q+4*mu]=-1.-mu;
    check("negative_noise_not_clipped",std::abs(vbs.derive(noisy)[2]+1.));
}
}
int main() {
    try {
        fields();freeCovariance();jointJackknife();
        for(const auto &e:errors)std::cout<<e.first<<" max_error="<<std::setprecision(12)<<e.second<<" tolerance=5e-12\n";
        std::cout<<"P4 cheap observable tests PASS (synthetic fields are tests, not physical plots).\n";
        return 0;
    }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
