#ifndef KITAEV_SSH_VBS_AMPLITUDE_H
#define KITAEV_SSH_VBS_AMPLITUDE_H
#include "vbs_observables.h"
#include "phonon_path.h"

namespace kitaev_ssh {
// Zero-source measurements. The finite-dt modulus source is DEFINED by
// W_epsilon=W_0 exp(epsilon*Nc*sum_j t_j |F_q(node_j)/Nc|), with
// trapezoidal t_j (dt/2 at fixed trial endpoints, dt in the interior).
// This is a symmetric diagonal insertion; its continuum limit corresponds
// to H_epsilon=H_0-epsilon*Nc*|D_q|. It is not the linear Mehler source score.
class VBSAmplitude {
public:
    explicit VBSAmplitude(const VBSObservables &vbs,bool response=false)
        :vbs_(vbs),response_(response) {
        if(vbs.mode()!=VBSMode::All)
            throw std::invalid_argument("--vbs-amplitude requires --vbs-mode all");
    }
    int stride() const {return response_?15:9;}
    static double nodeWeight(const PhononPath &p,int j) {
        if(j<0||j>=p.nodes())throw std::out_of_range("amplitude node");
        return p.dt()*((j==0||j==p.intervals())?.5:1.);
    }
    std::vector<std::string> names() const {
        std::vector<std::string> out;
        for(int q=0;q<vbs_.lattice().cells();++q) {
            const auto suffix="_q"+std::to_string(q)+"_";
            for(const auto *stem:{"AmpA","AmpA2","AmpA4"})
                for(const auto *mu:{"x","y","z"})out.push_back(stem+suffix+mu);
            if(response_)for(const auto *stem:{"AmpIA","AmpAIA"})
                for(const auto *mu:{"x","y","z"})out.push_back(stem+suffix+mu);
        }
        return out;
    }
    std::vector<Complex> measure(const PhononPath &p) const {
        if(p.bonds()!=vbs_.lattice().bondCount())
            throw std::invalid_argument("amplitude path bonds");
        const double nc=vbs_.lattice().cells();
        const auto d=vbs_.fourier(p.node(p.center()));
        std::vector<std::array<double,3>> integrals(d.size(),{0.,0.,0.});
        if(response_)for(int j=0;j<p.nodes();++j) {
            const auto f=(j==p.center())?d:vbs_.fourier(p.node(j));
            for(std::size_t q=0;q<f.size();++q)for(int mu=0;mu<3;++mu)
                integrals[q][mu]+=nodeWeight(p,j)*std::abs(f[q][mu])/nc;
        }
        std::vector<Complex> out;out.reserve(stride()*d.size());
        for(std::size_t q=0;q<d.size();++q) {
            std::array<double,3> a{},a2{};
            for(int mu=0;mu<3;++mu) {
                a2[mu]=std::norm(d[q][mu]/nc);a[mu]=std::sqrt(a2[mu]);
            }
            for(auto x:a)out.push_back(x);
            for(auto x:a2)out.push_back(x);
            for(auto x:a2)out.push_back(x*x); // Per-configuration fourth moment.
            if(response_) {
                for(auto x:integrals[q])out.push_back(x);
                for(int mu=0;mu<3;++mu)out.push_back(a[mu]*integrals[q][mu]);
            }
        }
        return out;
    }
private:
    const VBSObservables &vbs_;
    bool response_;
};
}
#endif
