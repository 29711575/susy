#ifndef KITAEV_SSH_TOTAL_ENERGY_H
#define KITAEV_SSH_TOTAL_ENERGY_H
#include "phonon_path.h"
#include "free_projector.h"
#include "phase_statistics.h"
#include <array>

namespace kitaev_ssh {
// Exact conditional insertions at the midpoint of ONE forced Mehler link.
// This is not the kinetic operator at the neighboring matter half-step cut.
struct PhononLinkEnergy {
    double mean=0., second=0., kinetic=0., potential=0., source=0.;
};
inline PhononLinkEnergy phononLinkEnergy(long double x,long double y,const PhononParameters &p,double dt,double force=0.) {
    if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(dt)||dt<=0||!std::isfinite(force))
        throw std::invalid_argument("invalid phonon energy link");
    using L=long double;
    const L m=p.mass,w=p.omega,k=m*w*w,z=w*L(dt),r=std::exp(-z),root=std::exp(-z/2),
        t=-std::expm1(-z)/(1+r),q=L(force)/k,
        shiftWeight=std::pow(-std::expm1(-z/2),2)/(1+r),
        mean=(L(x)+L(y))*root/(1+r)+q*shiftWeight,
        second=mean*mean+t/(2*m*w),potential=k*second/2,source=-L(force)*mean,
        difference=L(x)-L(y),sum=L(x)+L(y)-2*q,
        // -d_dt log K_force, using stable tanh/sech coefficients.
        hforced=w*(1+t*t)/(4*t)+k*r/(2*(1+r)*(1+r))*(sum*sum-(difference/t)*(difference/t))
            -L(force)*L(force)/(2*k),
        kinetic=hforced-potential-source;
    PhononLinkEnergy out{double(mean),double(second),double(kinetic),double(potential),double(source)};
    for(double value:{out.mean,out.second,out.kinetic,out.potential,out.source})
        if(!std::isfinite(value))throw std::runtime_error("nonfinite phonon link energy");
    return out;
}

inline std::vector<std::string> totalEnergyNames() {
    return {"energy_bare_matter","energy_matter_phonon","energy_phonon_kinetic",
            "energy_phonon_potential","energy_external_source","energy_total"};
}
// Extensively normalized, common-phase raw observables. Average the two link
// midpoints t=tau/2 +/- dt/2. All six components use these SAME cuts.
inline std::vector<Complex> totalEnergyObservables(const KitaevSSHModel &model,const PhononPath &path,
                                                  const std::array<MatType,2> &greens) {
    if(path.bonds()!=model.lattice().bondCount())throw std::invalid_argument("energy path/model dimension");
    std::vector<Complex> result(6,0.);
    for(int side=0;side<2;++side) {
        const int link=path.center()+side;
        const auto &green=greens[side];
        if(green.rows()!=model.majoranas()||green.cols()!=model.majoranas()||!green.allFinite())
            throw std::invalid_argument("invalid energy link Green");
        for(const auto &bond:model.lattice().bonds()) {
            const auto ph=phononLinkEnergy(path.coordinate(link-1,bond.index),path.coordinate(link,bond.index),
                                          path.parameters(),path.dt(),path.force(bond.index));
            const int mu=int(bond.type);
            const double j=mu==0?model.couplings().Jx:mu==1?model.couplings().Jy:model.couplings().Jz;
            const Complex bilinear=Complex(0,-1)*double(model.gauge()[bond.index])*green(bond.a,bond.b);
            result[0]+=.5*j*bilinear;
            result[1]+=.5*model.couplings().g*ph.mean*bilinear;
            result[2]+=.5*ph.kinetic;result[3]+=.5*ph.potential;result[4]+=.5*ph.source;
        }
    }
    for(int j=0;j<5;++j)result[5]+=result[j];
    for(auto value:result)if(!finiteComplex(value))throw std::runtime_error("nonfinite total energy observable");
    return result;
}
}
#endif
