#include "vbs_susceptibility.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace kitaev_ssh;
namespace {
void require(bool ok,const char *message) {
    if(!ok)throw std::runtime_error(message);
}
}
int main() {
    try {
        double maxScoreError=0.,maxFreeError=0.,maxFourierError=0.;
        // Differentiate the independently evaluated shifted Mehler kernel.
        for(double omega:{.05,.3,1.,4.})for(double dt:{.03,.1,.4}) {
            const double mass=1./(omega*omega);
            const PhononParameters parameters(mass,omega);
            const PhononPath kernel(1,4,dt,parameters);
            for(double x:{-.7,.2,1.3})for(double y:{-.4,.6})for(double f:{-.8,1.}) {
                // log K_h is exactly quadratic in h, so a symmetric finite
                // difference has zero truncation error even at this larger h.
                // Avoid subtracting nearly equal large small-omega actions.
                const double h=1e-2;
                auto forcedLog=[&](double field) {
                    const double shift=field*f/parameters.spring();
                    return -kernel.linkAction(x-shift,y-shift)+dt*field*field*f*f/(2*parameters.spring());
                };
                const double finiteDifference=(forcedLog(h)-forcedLog(-h))/(2*h);
                const double score=VBSSusceptibility::linkScoreCoefficient(omega,dt)*(x+y)*f;
                maxScoreError=std::max(maxScoreError,std::abs(finiteDifference-score));
                require(std::abs(finiteDifference-score)<2e-8,"forced-kernel score finite difference");
            }
            for(int intervals:{2,4,6,20,42,80,140}) {
                const PhononPath path(1,intervals,dt,parameters);
                double response=0.;
                for(int j=0;j<path.nodes();++j)
                    response+=VBSSusceptibility::nodeWeight(path,j)*
                        std::exp(-omega*dt*std::abs(j-path.center()))/(2*mass*omega);
                const double expected=-std::expm1(-omega*path.totalTau()/2)/parameters.spring();
                maxFreeError=std::max(maxFreeError,std::abs(response-expected));
                require(std::abs(response-expected)<2e-13,"finite-projection free-vacuum response");
            }
        }
        // Independent direct DFT and literal real sources on two geometries.
        for(int shear:{0,1}) {
            const HoneycombLattice lattice(3,2,shear);
            const VBSObservables vbs(lattice);
            const VBSSusceptibility susceptibility(vbs);
            PhononPath path(lattice.bondCount(),6,.1,PhononParameters(100./9,.3));
            for(int bond=0;bond<path.bonds();++bond) {
                std::vector<double> values;
                for(int j=0;j<path.nodes();++j)
                    values.push_back(.3*std::sin(.43*bond+.19*j)+.2*std::cos(.7*bond-.11*j));
                path=path.withBondCoordinates(bond,values);
            }
            const auto measured=susceptibility.measure(path);
            require(measured.size()==24*std::size_t(lattice.cells()),"raw moment size");
            require(susceptibility.names().size()==measured.size(),"raw names size");
            const double nc=lattice.cells();
            for(int q=0;q<lattice.cells();++q) {
                std::array<Complex,3> d{},integral{};
                const auto momentum=vbs.momenta()[q];
                for(const auto &bond:lattice.bonds()) {
                    const int mu=static_cast<int>(bond.type);
                    const auto phase=std::polar(1.,-2*pi*(momentum.reciprocal_x*bond.cell.x+momentum.reciprocal_y*bond.cell.y));
                    d[mu]+=phase*path.coordinate(path.center(),bond.index)/nc;
                    for(int j=0;j<path.nodes();++j) {
                        // Explicit link sum gives endpoints once/interior twice.
                        const double multiplicity=(j==0||j==path.intervals())?1.:2.;
                        const double w=multiplicity*std::tanh(.3*path.dt()/2)/.3;
                        integral[mu]+=w*phase*path.coordinate(j,bond.index)/nc;
                    }
                }
                const int offset=24*q;
                for(int mu=0;mu<3;++mu) {
                    maxFourierError=std::max(maxFourierError,std::abs(measured[offset+mu]-d[mu]));
                    maxFourierError=std::max(maxFourierError,std::abs(measured[offset+3+mu]-integral[mu]));
                    require(std::abs(measured[offset+6+mu]-std::conj(integral[mu]))<2e-14,"saved conjugate first moment");
                    for(int nu=0;nu<3;++nu)
                        require(std::abs(measured[offset+9+3*mu+nu]-nc*d[mu]*std::conj(integral[nu]))<2e-14,"cross response raw moment");
                }
                require(std::abs(measured[offset+22]-nc*d[2].real()*integral[2].real())<2e-14,"literal cosine source normalization");
                require(std::abs(measured[offset+23]-nc*d[2].imag()*integral[2].imag())<2e-14,"literal sine source normalization");
            }
        }
        require(maxFourierError<2e-14,"direct Fourier reference");
        std::cout<<std::setprecision(17)<<"{\"status\":\"passed\",\"max_score_finite_difference_error\":"<<maxScoreError
                 <<",\"max_free_response_error\":"<<maxFreeError
                 <<",\"max_direct_fourier_error\":"<<maxFourierError
                 <<",\"monte_carlo_used\":false}\n";
        return 0;
    } catch(const std::exception &error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
