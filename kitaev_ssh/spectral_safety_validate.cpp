// Deterministic condition-gate regression; no propagation or Monte Carlo.
#include "bipartite_spectrum.h"
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace kitaev_ssh;
namespace {
void require(bool ok,const char *message) {
    if(!ok)throw std::runtime_error(message);
}
MatType pairs(int n,double coupling) {
    MatType k=MatType::Zero(2*n,2*n);
    for(int j=0;j<n;++j) {
        k(j,n+j)=DataType(0.,-coupling);
        k(n+j,j)=DataType(0.,coupling);
    }
    return k;
}
}
int main() try {
    mkl_set_num_threads(1);
    int checks=0;
    for(int n:{1,2,3})for(double coupling:{0.,2.,800.}) {
        const MatType k=pairs(n,coupling);
        for(bool eligible:{false,true}) {
            BipartiteSpectrum workspace;
            PhononNodeSpectrum spectrum;
            require(spectrum.compute(k,eligible,workspace),"spectrum construction");
            require(spectrum.reduced()==(eligible&&n>=2),"spectral dispatch");
            for(double a:{0.,.001,1.,-.001,-1.}) {
                // For these exactly known paired eigenvalues, the condition
                // ratio of I+exp(-a K) is exp(-abs(a*coupling)).
                const double exponent=std::abs(a*coupling);
                const double exactRatio=std::exp(-exponent);
                const bool expected=std::isfinite(std::exp(exponent))&&
                    exactRatio>2.*double(2*n)*std::numeric_limits<double>::epsilon();
                require(spectrum.cayleySafe(a,2*n)==expected,
                        "paired-spectrum condition gate disagrees with analytic ratio");
                ++checks;
            }
            require(!spectrum.cayleySafe(std::numeric_limits<double>::quiet_NaN(),2*n),
                    "NaN scale must not certify a spectral Cayley");
            ++checks;
        }
    }
    std::cout<<"PASS spectral_safety checks="<<checks
             <<" analytic_paired_spectrum=1 overflow_rejected=1\n";
    return 0;
} catch(const std::exception &error) {
    std::cerr<<"FAIL spectral_safety: "<<error.what()<<'\n';return 1;
}
