#include "phonon_projector.h"
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>

using namespace kitaev_ssh;

namespace {
double maxEntryError(const MatType &a,const MatType &b) {
    if(a.rows()!=b.rows()||a.cols()!=b.cols()) return std::numeric_limits<double>::infinity();
    return (a-b).cwiseAbs().maxCoeff();
}
void require(bool condition,const char *message) {
    if(!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F &&function,const char *message) {
    bool rejected=false;
    try { function(); } catch(const std::exception &) { rejected=true; }
    require(rejected,message);
}
}

int main() {
    try {
        std::mt19937_64 rng(20260920);
        double maxGeneratorError=0.,maxIntoError=0.,maxHalfStepError=0.;
        long cases=0;
        for(const auto size:std::vector<std::pair<int,int>>{{2,3},{3,2},{4,4}}) {
            const double g=.37;
            const KitaevSSHModel model(HoneycombLattice(size.first,size.second),p2Couplings(g));
            const auto generator=honeycombGenerator(model);
            const auto path=PhononPath(model.lattice().bondCount(),6,.13,{1.7,.8}).freeSample(rng);
            MajoranaGenerationWorkspace workspace;
            MatType into;
            for(int j=0;j<path.nodes();++j) {
                const auto &x=path.node(j);
                const auto reference=KitaevSSHModel(model.lattice(),model.couplings(),x,model.gauge()).frameworkGenerator();
                const auto returned=model.frameworkGeneratorWith(x);
                model.frameworkGeneratorWith(x,into);
                maxGeneratorError=std::max(maxGeneratorError,maxEntryError(reference,returned));
                maxIntoError=std::max(maxIntoError,maxEntryError(reference,into));
                require(returned.rows()==model.majoranas() && into.rows()==model.majoranas(),
                        "generator dimension changed");

                MatType half;
                nodeHalfStepInto(path,j,model.majoranas(),generator,half,workspace);
                MatType referenceForExpm=reference;
                const auto referenceHalf=expm(referenceForExpm,-path.dt()/2);
                maxHalfStepError=std::max(maxHalfStepError,maxEntryError(referenceHalf,half));
                ++cases;
            }
            const auto factors=pathFactors(path,model.majoranas(),generator);
            require(int(factors.size())==2*path.intervals(),"factor count changed");
            for(int f=0;f<int(factors.size());++f) {
                const auto reference=KitaevSSHModel(model.lattice(),model.couplings(),
                    path.node(factorNode(f,path.intervals())),model.gauge()).frameworkGenerator();
                MatType referenceForExpm=reference;
                const auto expected=expm(referenceForExpm,-path.dt()/2);
                maxHalfStepError=std::max(maxHalfStepError,maxEntryError(expected,MatType(factors[f].matrix)));
            }
            auto badCount=path.node(0);
            badCount.pop_back();
            rejects([&] { (void)model.frameworkGeneratorWith(badCount); },"wrong X count accepted");
            auto badFinite=path.node(0);
            badFinite[0]=std::numeric_limits<double>::quiet_NaN();
            rejects([&] { (void)model.frameworkGeneratorWith(badFinite); },"nonfinite X accepted");
        }
        std::cout<<std::setprecision(17)
            <<"{\"cases\":"<<cases
            <<",\"max_generator_error\":"<<maxGeneratorError
            <<",\"max_into_error\":"<<maxIntoError
            <<",\"max_half_step_error\":"<<maxHalfStepError<<"}\n";
        require(maxGeneratorError<2e-13,"dynamic generator mismatch");
        require(maxIntoError<2e-13,"output generator mismatch");
        require(maxHalfStepError<2e-13,"half-step mismatch");
        std::cout<<"A3_GENERATOR_VALIDATE_OK\n";
        return 0;
    } catch(const std::exception &error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
