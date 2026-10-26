#ifndef KITAEV_SSH_PHYSICAL_MATTER_TRIAL_H
#define KITAEV_SSH_PHYSICAL_MATTER_TRIAL_H

#include "free_projector.h"
#include <ostream>

namespace kitaev_ssh {

// Standard D_i=b_i^x b_i^y b_i^z c_i=+1 and positive spin Hamiltonian:
// u_code=-u_phys, P_code=prod_r(i*c_A,r*c_B,r). See ADR-023.
inline int requiredPhysicalMatterParity(const KitaevSSHModel &model) {
    const auto &lat=model.lattice();
    const long long theta=static_cast<long long>(lat.Lx())+lat.Ly()+
        static_cast<long long>(lat.shear())*(lat.Lx()-lat.shear());
    int result=(theta%2 ? -1 : 1)*(lat.cells()%2 ? -1 : 1);
    for(int u:model.gauge())result*=u;
    return result;
}

struct MatterTrialSelection {
    GaussianTrialState trial;
    std::string policy;
    int required_parity=0, unrestricted_parity=0, actual_parity=0;
    bool blocked=false;
    int momentum_n1=-1, momentum_n2=-1, minimum_mode_degeneracy=0;
    double momentum_u=0., momentum_v=0., blocking_energy_cost=0.;
    bool physicalParityCompatible() const { return actual_parity==required_parity; }
};

inline MatterTrialSelection selectFreeMatterTrial(
        const KitaevSSHModel &model,const std::string &policy="physical") {
    if(policy!="physical" && policy!="unrestricted")
        throw std::invalid_argument("--matter-trial must be physical or unrestricted");
    MatterTrialSelection result;
    result.policy=policy;
    result.required_parity=requiredPhysicalMatterParity(model);
    result.trial=freeGroundState(model);
    result.actual_parity=result.unrestricted_parity=result.trial.fermionParity();
    if(policy=="unrestricted" || result.physicalParityCompatible())return result;

    // Block a *fixed* lowest-cost quasiparticle once, before sampling. The
    // opposite parity free state need not be the interacting ground state.
    const auto &lat=model.lattice();const int n=lat.cells(),m=model.majoranas();
    const MatType k=model.frameworkGenerator();
    Eigen::SelfAdjointEigenSolver<MatType> es(k);
    if(es.info()!=Eigen::Success)throw std::runtime_error("physical trial eigensolver failed");
    const double gap=es.eigenvalues()(n);
    const double guard=512.*m*std::numeric_limits<double>::epsilon()*std::max(1.,k.norm());
    if(gap<=1e-10)throw std::runtime_error("physical trial has ambiguous zero mode");
    for(int j=0;j<n;++j)
        if(std::abs(es.eigenvalues()(j)+gap)<=guard)++result.minimum_mode_degeneracy;

    // Gauge-dressed Bloch modes preserve lattice translations for the uniform
    // reference Hamiltonian, even when the minimum is degenerate. Construct
    // the site gauge to the canonical twisted-boundary representative.
    const int wx=lat.wilsonLoop(0,model.gauge()),wy=lat.wilsonLoop(1,model.gauge());
    const auto canonical=lat.zeroFluxGauge(wx,wy);
    std::vector<int> d(m,0);d[0]=1;
    for(int pass=0;pass<m;++pass) {
        bool changed=false;
        for(const auto &b:lat.bonds()) {
            const int ratio=model.gauge()[b.index]*canonical[b.index];
            if(d[b.a] && !d[b.b]){d[b.b]=d[b.a]*ratio;changed=true;}
            else if(d[b.b] && !d[b.a]){d[b.a]=d[b.b]*ratio;changed=true;}
            if(d[b.a] && d[b.b] && d[b.a]*d[b.b]!=ratio)
                throw std::runtime_error("physical trial canonical gauge mismatch");
        }
        if(!changed)break;
    }
    for(int sign:d)if(!sign)throw std::runtime_error("disconnected physical trial gauge");
    cVecType mode;
    bool found=false;
    for(int n1=0;n1<lat.Lx() && !found;++n1)for(int n2=0;n2<lat.Ly() && !found;++n2) {
        const double u=(n1+(wx<0?.5:0.))/lat.Lx();
        const double v=(n2+(wy<0?.5:0.)-lat.shear()*u)/lat.Ly();
        cVecType seed=cVecType::Zero(m);
        for(int x=0;x<lat.Lx();++x)for(int y=0;y<lat.Ly();++y) {
            const int a=lat.cellIndex(x,y);
            seed[a]=double(d[a])*std::exp(DataType(0.,2*pi*(u*x+v*y)))/std::sqrt(double(n));
        }
        cVecType candidate=result.trial.Phi*(result.trial.Phi.adjoint()*seed);
        if(candidate.norm()<1e-10)continue;
        candidate.normalize();
        if((k*candidate+gap*candidate).norm()>guard)continue;
        mode=std::move(candidate);found=true;
        result.momentum_n1=n1;result.momentum_n2=n2;
        result.momentum_u=u;result.momentum_v=v;
    }
    if(!found)throw std::runtime_error(
        "physical trial needs a uniform reference with a lowest-energy Bloch mode; "
        "nonuniform trials require an explicit physical-state construction");

    // Rotate the occupied basis so column zero spans the selected mode, then
    // exchange it with its particle-hole conjugate. fromPhi rechecks purity.
    const cVecType coordinates=result.trial.Phi.adjoint()*mode;
    Eigen::HouseholderQR<MatType> qr(coordinates);
    const MatType rotation=qr.householderQ()*MatType::Identity(n,n);
    MatType phi=result.trial.Phi*rotation;
    phi.col(0)=phi.col(0).conjugate().eval();
    auto blocked=GaussianTrialState::fromPhi(phi);
    const MatType expected=result.trial.G_T+2.*mode*mode.adjoint()-
        2.*mode.conjugate()*mode.transpose();
    if((blocked.G_T-expected).norm()>guard)
        throw std::runtime_error("physical trial blocked covariance mismatch");
    const double cost=(matterEnergy(model,blocked.G_T)-matterEnergy(model,result.trial.G_T)).real();
    if(std::abs(cost-gap)>guard)
        throw std::runtime_error("physical trial is not the lowest-cost parity excitation");
    result.actual_parity=blocked.fermionParity();
    if(!result.physicalParityCompatible())
        throw std::runtime_error("physical trial parity correction failed");
    // A complex Bloch mode may lack a certified real representation; the
    // unchanged full-complex path is then used, never an approximate real one.
    certifyRealMagnitudeTrial(blocked);
    result.trial=std::move(blocked);result.blocked=true;result.blocking_energy_cost=gap;
    return result;
}

inline void writeMatterTrialSelection(std::ostream &out,const MatterTrialSelection &s) {
    out<<"{\"policy\":\""<<s.policy<<"\",\"required_code_parity\":"<<s.required_parity
       <<",\"unrestricted_code_parity\":"<<s.unrestricted_parity<<",\"actual_code_parity\":"<<s.actual_parity
       <<",\"physical_parity_compatible\":"<<(s.physicalParityCompatible()?"true":"false")
       <<",\"blocked\":"<<(s.blocked?"true":"false")<<",\"blocking_energy_cost\":"<<s.blocking_energy_cost
       <<",\"minimum_mode_degeneracy\":"<<s.minimum_mode_degeneracy
       <<",\"blocked_momentum_indices\":["<<s.momentum_n1<<','<<s.momentum_n2<<']'
       <<",\"blocked_momentum_fractional\":["<<s.momentum_u<<','<<s.momentum_v<<']'
       <<",\"spin_mapping\":\"positive_spin_H; standard_D_plus; u_code=-u_phys\""
       <<",\"full_physical_spin_projection\":false}";
}

} // namespace kitaev_ssh
#endif
