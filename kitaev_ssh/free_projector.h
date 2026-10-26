#ifndef KITAEV_SSH_FREE_PROJECTOR_H
#define KITAEV_SSH_FREE_PROJECTOR_H

#include "kitaev_ssh_model.h"
#include "gaussian_trial_state.h"
#include "pure_projector_endpoint.h"
#include "pure_projector_stack.h"

namespace kitaev_ssh {

inline void certifyRealMagnitudeTrial(GaussianTrialState &trial) {
    trial.real_magnitude_trial.reset();
    // D=diag(I,iI) makes this exact real bipartite Hamiltonian real symmetric.
    // A real basis is a unitary rotation within the occupied subspace; it
    // changes neither |overlap| nor the original trial used for phase/Green.
    const int m=trial.Phi.rows(),n=trial.Phi.cols();
    MatType transformed=trial.Phi;
    transformed.bottomRows(n)*=DataType(0.,-1.);
    Eigen::MatrixXd candidates(m,2*n);
    candidates.leftCols(n)=transformed.real();
    candidates.rightCols(n)=transformed.imag();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(candidates);
    Eigen::MatrixXd basis=qr.householderQ()*Eigen::MatrixXd::Identity(m,n);
    const MatType coordinates=basis.transpose().cast<DataType>()*transformed;
    const double subspace=(transformed-basis.cast<DataType>()*coordinates).norm();
    const double unitary=(coordinates.adjoint()*coordinates-MatType::Identity(n,n)).norm();
    const double orthogonal=(basis.transpose()*basis-Eigen::MatrixXd::Identity(n,n)).norm();
    const double guard=128.*m*std::numeric_limits<double>::epsilon();
    if(basis.allFinite()&&subspace<=guard&&unitary<=guard&&orthogonal<=guard) {
        auto representation=std::make_shared<GaussianRealMagnitudeTrial>();
        representation->source=trial.Phi;representation->basis=std::move(basis);
        trial.real_magnitude_trial=std::move(representation);
    }
}

// Unrestricted mathematical free vacuum, retained for independent numerical
// references. Production uses selectFreeMatterTrial's physical parity policy.
inline GaussianTrialState freeGroundState(const KitaevSSHModel &model) {
    // Exact/near-zero modes remain an error; never alter H or gauge to split them.
    auto trial=GaussianTrialState::fromMajoranaHamiltonian(model.frameworkGenerator());
    certifyRealMagnitudeTrial(trial);
    return trial;
}

struct StaticProjection {
    MatType green;
    double energy = 0., log_norm_squared = 0.;
    int matter_parity = 0;
    double green_residual = 0., solve_residual = 0., overlap_rcond = 0.;
};

inline DataType matterEnergy(const KitaevSSHModel &model, const MatType &green) {
    if (green.rows()!=model.majoranas() || green.cols()!=model.majoranas() || !green.allFinite())
        throw std::invalid_argument("invalid matter Green matrix");
    DataType energy = 0.;
    for (const auto &b : model.lattice().bonds())
        energy -= DataType(0,1)*double(model.gauge()[b.index])*model.bondHopping(b.index)*green(b.a,b.b);
    return energy;
}

// Same estimator for a sampled displacement field. Reuses the validated
// lattice, couplings, and gauge instead of constructing a temporary model.
inline DataType matterEnergyAtDisplacement(const KitaevSSHModel &model,
                                           const std::vector<double> &x,
                                           const MatType &green) {
    if (x.size()!=std::size_t(model.lattice().bondCount()) ||
        green.rows()!=model.majoranas() || green.cols()!=model.majoranas() || !green.allFinite())
        throw std::invalid_argument("invalid displacement or matter Green matrix");
    DataType energy=0.;
    for (const auto &b:model.lattice().bonds()) {
        const double j=b.type==BondType::X?model.couplings().Jx:
                       b.type==BondType::Y?model.couplings().Jy:model.couplings().Jz;
        energy-=DataType(0,1)*double(model.gauge()[b.index])*
            (j+model.couplings().g*x[b.index])*green(b.a,b.b);
    }
    return energy;
}

// An estimator for a specified path, not a production measurement accumulator.
inline DataType gaugeStringGreen(const KitaevSSHModel &model, const MatType &green, const BondPath &path) {
    if (green.rows()!=model.majoranas() || green.cols()!=model.majoranas() || !green.allFinite())
        throw std::invalid_argument("invalid matter Green matrix");
    const auto &lattice = model.lattice();
    const int end = lattice.pathEnd(path), string = lattice.gaugeProduct(path,model.gauge());
    return double(string)*(path.start==end ? DataType(1) : -green(path.start,end));
}

inline std::vector<PureProjectorSlice> staticContour(const KitaevSSHModel &model, double tauPerBranch, double dt) {
    if (!std::isfinite(tauPerBranch) || tauPerBranch<=0 || !std::isfinite(dt) || dt<=0)
        throw std::invalid_argument("tau per branch and dt must be positive and finite");
    const double steps = tauPerBranch/dt;
    if (!std::isfinite(steps) || steps>std::numeric_limits<int>::max()/2)
        throw std::invalid_argument("too many projection steps");
    const int count = int(std::llround(steps));
    if (count<1 || std::abs(steps-count)>1e-10*std::max(1.,steps))
        throw std::invalid_argument("tau-per-branch/dt must be a positive integer");
    const MatType b = model.propagator(dt);
    std::vector<PureProjectorSlice> contour;
    contour.reserve(2*count);
    // Static Hermitian factors: explicit reverse-adjoint bra and ket halves.
    // eta=1 is the canonical exponential lift, with no HS scalar prefactor.
    for (int n = count-1; n>=0; --n) contour.emplace_back(b.adjoint().eval(),1.,"static_bra:"+std::to_string(n));
    for (int n = 0; n<count; ++n) contour.emplace_back(b,1.,"static_ket:"+std::to_string(n));
    return contour;
}

inline StaticProjection projectStatic(const KitaevSSHModel &model, const GaussianTrialState &trial,
                                      double tauPerBranch, double dt, int block = 8) {
    if (trial.Phi.rows()!=model.majoranas() || block<=0) throw std::invalid_argument("trial dimension or stabilization block");
    const auto contour = staticContour(model,tauPerBranch,dt);
    const int center = int(contour.size()/2);
    PureProjectorOptions options;
    options.fail_on_green_structure = true;
    PureProjectorStackManager stack(trial,contour,block,options,center);
    if (!stack.ok()) throw std::runtime_error("static projector stack: "+stack.failureDetail());
    const auto green = stack.green();
    if (!green.ok()) throw std::runtime_error("static projector Green failed");
    const auto endpoint = pureProjectorEndpointRebuild(trial,contour,center,block,options);
    if (!endpoint.ok()) throw std::runtime_error("static projector norm: "+endpoint.message);
    const auto energy = matterEnergy(model,green.green);
    if (!std::isfinite(energy.real()) || std::abs(energy.imag())>1e-10*std::max(1.,std::abs(energy.real())))
        throw std::runtime_error("nonreal static matter energy");
    // This contour is exactly U^dagger U. Its norm is positive; no HS1
    // positivity contract is used. General MC signs are outside this adapter.
    return {green.green,energy.real(),endpoint.log_abs_weight,trial.fermionParity(),
            green.green_residual,green.solve_residual,green.overlap_rcond};
}

} // namespace kitaev_ssh
#endif
