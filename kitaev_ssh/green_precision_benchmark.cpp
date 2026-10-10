#include "gaussian_trial_state.h"
#include "pure_projector_fast.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

MatType unitary(int n, std::mt19937_64 &rng) {
    std::normal_distribution<double> normal;
    MatType a(n,n);
    for (int j=0; j<n; ++j) for (int i=0; i<n; ++i)
        a(i,j)=DataType(normal(rng),normal(rng));
    Eigen::HouseholderQR<MatType> qr(a);
    return qr.householderQ()*MatType::Identity(n,n);
}

struct Bases { MatType right,left; };
Bases bases(int n, double smallest) {
    std::mt19937_64 rng(1);
    MatType phi=MatType::Zero(2*n,n);
    const double h=1/std::sqrt(2.0);
    for (int j=0; j<n; ++j) { phi(j,j)=h; phi(n+j,j)=DataType(0,h); }
    MatType c=MatType::Identity(n,n),s=MatType::Zero(n,n);
    const int paired=n-n%2;
    // An odd final column stays in phi; the other columns form pure pairs.
    for (int j=0; j<paired; j+=2) {
        const double value=std::pow(smallest,double(j)/double(paired-2));
        const double sine=std::sqrt(std::max(0.0,1-value*value));
        c(j,j)=c(j+1,j+1)=value;
        s(j,j+1)=sine; s(j+1,j)=-sine;
    }
    Bases result{phi*c+phi.conjugate()*s,phi};
    result.left*=unitary(n,rng);
    result.right*=unitary(n,rng);
    return result;
}

void checkPure(const MatType &phi) {
    require(phi.allFinite(),"nonfinite basis");
    require((phi.adjoint()*phi-MatType::Identity(phi.cols(),phi.cols())).norm()<1e-12,
            "basis is not orthonormal");
    require((phi.transpose()*phi).norm()<1e-12,"basis is not pure/isotropic");
}

PureProjectorGreenResult evaluate(const Bases &b, const PureProjectorOptions &options) {
    auto result=pureProjectorGreen(b.right,b.left,options);
    if (!result.ok()) result=pureProjectorCheckObliqueGreen(result,options);
    if (!result.ok()) throw std::runtime_error(pureProjectorStatusName(result.status));
    require(result.green.allFinite(),"nonfinite Green matrix");
    return result;
}
}

int main(int argc, char **argv) {
    try {
        require(argc<=2,"usage: kitaev_ssh_green_precision_benchmark [repeats]");
        int repeats=20;
        if (argc==2) {
            std::size_t consumed=0;
            repeats=std::stoi(argv[1],&consumed);
            require(consumed==std::string(argv[1]).size() && repeats>0,
                    "repeats must be a positive integer");
        }
        const auto options=pureFastGreenOptions(PureFastOptions{});
        std::cout<<std::setprecision(17);
        double totalSeconds=0;
        for (int n : {12,36,81}) for (double smallest : {0.7,1e-4,1e-6}) {
            const Bases b=bases(n,smallest);
            checkPure(b.left); checkPure(b.right);
            auto result=evaluate(b,options); // Untimed warm-up for each case.
            const auto start=std::chrono::steady_clock::now();
            for (int repeat=0; repeat<repeats; ++repeat) result=evaluate(b,options);
            const double seconds=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-start).count();
            totalSeconds+=seconds;
            // Successful strict checks do not populate this recovery-only field.
            const double dimensionScale=std::max(1.0,std::sqrt(double(result.green.rows())));
            const double relativeResidual=result.green_residual*
                (dimensionScale/std::max(dimensionScale,result.green.norm()));
            require(std::isfinite(relativeResidual),"nonfinite relative Green residual");
            std::cout<<"case n="<<n<<" dimension="<<2*n<<" smallest="<<smallest
                     <<" seed=1 repeats="<<repeats<<" rcond="<<result.overlap_rcond
                     <<" solve_residual="<<result.solve_residual
                     <<" green_relative_residual="<<relativeResidual
                     <<" extended_precision_solve="<<result.extended_precision_solve
                     <<" steps="<<result.solve_refinement_steps<<" seconds="<<seconds<<'\n';
        }
        std::cout<<"total_seconds="<<totalSeconds<<'\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr<<"green_precision_benchmark: "<<error.what()<<'\n';
        return 1;
    }
}
