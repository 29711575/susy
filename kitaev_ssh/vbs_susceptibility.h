#ifndef KITAEV_SSH_VBS_SUSCEPTIBILITY_H
#define KITAEV_SSH_VBS_SUSCEPTIBILITY_H
#include "vbs_observables.h"
#include "phonon_path.h"

namespace kitaev_ssh {
// Finite-projection derivative for H_h=H-h sum_r f(r)X(r), FIXED trial
// endpoints. Differentiating the shifted exact Mehler kernel gives
// d_h log K_h(x,y)|0=tanh(omega*dt/2)/omega*(x+y).
// No explicit source derivative occurs in the matter factors.
class VBSSusceptibility {
public:
    static constexpr int stride=24;
    static double linkScoreCoefficient(double omega,double dt) {
        if(!std::isfinite(omega)||omega<=0||!std::isfinite(dt)||dt<=0)
            throw std::invalid_argument("invalid susceptibility omega/dt");
        return std::tanh(omega*dt/2)/omega;
    }
    static double nodeWeight(const PhononPath &p,int j) {
        if(j<0||j>=p.nodes())throw std::out_of_range("susceptibility node");
        return linkScoreCoefficient(p.parameters().omega,p.dt())*
            ((j==0||j==p.intervals())?1.:2.);
    }
    explicit VBSSusceptibility(const VBSObservables &vbs):vbs_(vbs) {
        if(vbs.mode()!=VBSMode::All)throw std::invalid_argument("--vbs-susceptibility requires --vbs-mode all");
    }
    std::vector<std::string> names() const {
        std::vector<std::string> out;
        for(int q=0;q<vbs_.lattice().cells();++q) {
            const std::string suffix="_q"+std::to_string(q)+"_";
            for(const auto *stem:{"ChiD","ChiI","ChiIconj"})
                for(const auto *mu:{"x","y","z"})out.push_back(stem+suffix+mu);
            for(const auto *mu:{"x","y","z"})for(const auto *nu:{"x","y","z"})
                out.push_back("ChiRaw"+suffix+mu+nu);
            for(const auto *stem:{"ChiDRe","ChiDIm","ChiIRe","ChiIIm"})out.push_back(stem+suffix+"z");
            out.push_back("ChiCosRaw"+suffix+"zz");out.push_back("ChiSinRaw"+suffix+"zz");
        }
        return out;
    }
    std::vector<Complex> measure(const PhononPath &p) const {
        if(p.bonds()!=vbs_.lattice().bondCount())throw std::invalid_argument("susceptibility path bonds");
        std::vector<double> integrated(p.bonds(),0.);
        for(int j=0;j<p.nodes();++j) {
            const double weight=nodeWeight(p,j);
            for(int b=0;b<p.bonds();++b)integrated[b]+=weight*p.coordinate(j,b);
        }
        // Fourier linearity: two transforms of three-channel fields, independent
        // of node count; the full real-space path is summed before the FFT.
        const auto d=vbs_.fourier(p.node(p.center())),i=vbs_.fourier(integrated);
        const double nc=vbs_.lattice().cells();
        std::vector<Complex> out;out.reserve(stride*d.size());
        for(std::size_t q=0;q<d.size();++q) {
            for(auto z:d[q])out.push_back(z/nc);
            for(auto z:i[q])out.push_back(z/nc);
            for(auto z:i[q])out.push_back(std::conj(z)/nc);
            for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)
                out.push_back(d[q][mu]*std::conj(i[q][nu])/nc);
            out.push_back(d[q][2].real()/nc);out.push_back(d[q][2].imag()/nc);
            out.push_back(i[q][2].real()/nc);out.push_back(i[q][2].imag()/nc);
            out.push_back(d[q][2].real()*i[q][2].real()/nc);
            out.push_back(d[q][2].imag()*i[q][2].imag()/nc);
        }
        return out;
    }
private:
    const VBSObservables &vbs_;
};
}
#endif
