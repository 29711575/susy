#include "gaussian_trial_state.h"
#include "pure_projector_green.h"
#include "pure_projector_mp.h"
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

namespace {
using MP = pure_projector_mp_detail::Number<50>;
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

MatType unitary(int n, std::mt19937_64 &rng) {
    std::normal_distribution<double> normal;
    MatType a(n,n);
    for (int j=0;j<n;++j) for (int i=0;i<n;++i)
        a(i,j)=DataType(normal(rng),normal(rng));
    Eigen::HouseholderQR<MatType> qr(a);
    return qr.householderQ()*MatType::Identity(n,n);
}

struct Bases { MatType right,left; };
Bases bases(int n, double smallest, std::uint64_t seed, bool dense=true) {
    std::mt19937_64 rng(seed);
    MatType phi=MatType::Zero(2*n,n);
    const double h=1/std::sqrt(2.0);
    for (int j=0;j<n;++j) { phi(j,j)=h; phi(n+j,j)=DataType(0,h); }
    MatType c=MatType::Zero(n,n),s=MatType::Zero(n,n);
    // Each pair has C=c I and skew S; C^T S+S^T C=0.
    for (int j=0;j<n;j+=2) {
        const double value=std::pow(smallest,double(j)/double(n-2));
        const double sine=std::sqrt(std::max(0.0,1-value*value));
        c(j,j)=c(j+1,j+1)=value;
        s(j,j+1)=sine; s(j+1,j)=-sine;
    }
    Bases result{phi*c+phi.conjugate()*s,phi};
    if (dense) { result.left*=unitary(n,rng); result.right*=unitary(n,rng); }
    return result;
}

void checkPure(const MatType &phi) {
    require((phi.adjoint()*phi-MatType::Identity(phi.cols(),phi.cols())).norm()<1e-12,
            "basis is not orthonormal");
    require((phi.transpose()*phi).norm()<1e-12,"basis is not pure/isotropic");
}

MatType reference(const Bases &b, bool roundedOverlap) {
    const auto left=MP::fromDouble(b.left),right=MP::fromDouble(b.right);
    const auto rhs=MP::adjoint(left);
    // Compare the rounded-overlap system and the actual supplied subspaces
    // independently; extended recovery recomputes the latter consistently.
    const auto overlap=roundedOverlap ? MP::fromDouble(MatType(b.left.adjoint()*b.right))
                                     : MP::multiply(rhs,right);
    const auto solved=MP::solve(overlap,rhs);
    require(solved.ok,"MP50 solve failed");
    require(solved.residual<MP::Real("1e-40"),"MP50 reference residual exceeded");
    return MP::toDouble(MP::sub(MP::Matrix::identity(b.right.rows()),
                              MP::scale(MP::multiply(right,solved.solution),MP::Real(2))));
}

double relativeError(const MatType &a,const MatType &b) {
    return (a-b).norm()/std::max(1.0,b.norm());
}
}

int main() {
    try {
        PureProjectorOptions options;
        require(options.solve_residual_tolerance==1e-10,"original solve gate changed");
        require(options.green_residual_tolerance==1e-10,"original Green gate changed");
        std::cout<<std::setprecision(17);
        int triggered=0,refined=0,accepted=0,greenRejected=0;
        double maximumError=0,maximumOverlapRoundingEffect=0;
        for (double smallest : {1e-5,3e-6,1e-6,3e-7}) {
            for (std::uint64_t seed=1;seed<=12;++seed) {
                const Bases b=bases(12,smallest,seed);
                checkPure(b.left); checkPure(b.right);
                const auto result=pureProjectorGreen(b.right,b.left,options);
                require(result.overlap_rank==12,"unexpected overlap rank failure");
                require(result.overlap_rcond>=options.minimum_overlap_rcond,"unexpected rcond failure");
                require(std::isfinite(result.solve_residual) && result.solve_residual<=1e-10,
                        "ill-conditioned solve failed original residual gate");
                require(result.solve_refinement_steps>=0 && result.solve_refinement_steps<=5,
                        "refinement step count outside limit");
                require(result.extended_precision_solve==(result.solve_refinement_steps==5),
                        "extended solve precision metadata inconsistent");
                const bool needs=result.initial_solve_residual>1e-10;
                triggered+=needs; refined+=result.solve_refinement_steps>0;
                require(!needs || result.solve_refinement_steps>0,"failed initial solve was not recovered");
                require(needs || result.solve_refinement_steps==0 || result.extended_precision_solve,
                        "solve changed without solve or Green certification failure");
                require(result.ok() || result.status==PureProjectorStatus::green_residual_exceeded,
                        "unexpected post-solve status");
                accepted+=result.ok(); greenRejected+=!result.ok();
                const MatType sameSystem=reference(b,true),physicalSystem=reference(b,false);
                const double error=relativeError(result.green,physicalSystem);
                require(error<=2e-10,"Green relative error versus original-subspace MP50 exceeded");
                maximumError=std::max(maximumError,error);
                maximumOverlapRoundingEffect=std::max(maximumOverlapRoundingEffect,
                    relativeError(sameSystem,physicalSystem));
                std::cout<<"case seed="<<seed<<" smallest="<<smallest
                         <<" initial="<<result.initial_solve_residual<<" final="<<result.solve_residual
                         <<" steps="<<result.solve_refinement_steps<<" mp_error="<<error
                         <<" status="<<pureProjectorStatusName(result.status)<<'\n';
            }
        }
        for (std::uint64_t seed=1;seed<=8;++seed) {
            const Bases b=bases(12,0.7,seed);
            checkPure(b.left); checkPure(b.right);
            const MatType overlap=b.left.adjoint()*b.right,rhs=b.left.adjoint();
            PureGreenSvd svd(overlap,options.accelerated_linear_algebra);
            svd.setThreshold(options.rank_tolerance);
            const MatType solved=svd.solve(rhs);
            const MatType oldGreen=MatType::Identity(24,24)-2.0*b.right*solved;
            const double oldResidual=(overlap*solved-rhs).norm()/std::max(1.0,rhs.norm());
            const auto result=pureProjectorGreen(b.right,b.left,options);
            require(result.ok(),"ordinary case rejected");
            require(result.solve_refinement_steps==0,"ordinary case refined");
            require(result.initial_solve_residual==oldResidual && result.solve_residual==oldResidual,
                    "ordinary residual changed");
            require((result.green.array()==oldGreen.array()).all(),"ordinary Green changed");
            require(relativeError(result.green,reference(b,true))<=2e-10,"ordinary MP50 error exceeded");
        }
        const Bases singular=bases(12,0,1,false);
        checkPure(singular.left); checkPure(singular.right);
        require(pureProjectorGreen(singular.right,singular.left,options).status==
                PureProjectorStatus::overlap_rank_deficient,"singular overlap accepted");
        const Bases ill=bases(12,1e-8,1,false);
        PureProjectorOptions rcondOptions=options;
        rcondOptions.minimum_overlap_rcond=1e-7;
        require(pureProjectorGreen(ill.right,ill.left,rcondOptions).status==
                PureProjectorStatus::overlap_ill_conditioned,"rcond gate bypassed");
        MatType deficient=ill.right; deficient.col(1)=deficient.col(0);
        require(pureProjectorGreenThinQr(deficient,ill.left,options).status==
                PureProjectorStatus::subspace_rank_deficient,"subspace rank gate bypassed");
        for (double bad : {std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
            for (bool changeLeft : {false,true}) {
                Bases b=ill;
                (changeLeft ? b.left : b.right)(0,0)=DataType(bad,0);
                require(pureProjectorGreen(b.right,b.left,options).status==
                        PureProjectorStatus::nonfinite_input,"nonfinite input accepted");
            }
        }
        std::cout<<"coverage initial_gate_exceeded="<<triggered<<" refined="<<refined
                 <<" accepted="<<accepted<<" original_green_gate_rejected="<<greenRejected
                 <<" max_original_subspace_mp_error="<<maximumError
                 <<" max_overlap_rounding_effect="<<maximumOverlapRoundingEffect<<'\n';
        // Coverage is reported rather than tied to one LAPACK/Eigen implementation.
        if (!triggered) std::cout<<"coverage_warning: refinement path not exercised on this platform\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr<<"green_refinement_validate: "<<error.what()<<'\n';
        return 1;
    }
}
