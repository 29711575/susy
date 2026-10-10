#include "phonon_mc.h"
#include "vbs_field.h"
#include <iostream>
#include <sstream>
using namespace kitaev_ssh;
void check(double a,double b,double tol,const char *name) {
    if(!std::isfinite(a) || !std::isfinite(b) || std::abs(a-b)>tol*(1+std::abs(b)))
        throw std::runtime_error(std::string(name)+": "+std::to_string(a)+" != "+std::to_string(b));
}
int main() {try {
    const PhononParameters p(100./9,.3);std::mt19937_64 rng(981734);
    for(int M:{4,80,120}) {
        auto zero=PhononPath(1,M,.1,p,{.17}).freeSample(rng);
        const double h=.05,dh=1e-5,k=p.spring();
        auto forced=[&](double f){auto q=PhononPath(1,M,.1,p,{.17},{f});for(int j=0;j<q.nodes();++j)q=q.withCoordinate(j,0,zero.coordinate(j,0));return q;};
        const auto path=forced(h);double exact=zero.logWeight();
        for(int j=1;j<path.nodes();++j) {
            const double x=path.coordinate(j,0),y=path.coordinate(j-1,0),s=h/k;
            exact+=zero.linkAction(x,y)-zero.linkAction(x-s,y-s)+path.dt()*h*h/(2*k);
        }
        check(path.logWeight(),exact,2e-12,"forced kernel weight");
        double score=h/k*M*(path.dt()-2*path.sourceNodeWeight(0));
        for(int j=0;j<path.nodes();++j)score+=path.sourceNodeWeight(j)*path.coordinate(j,0);
        check((forced(h+dh).logWeight()-forced(h-dh).logWeight())/(2*dh),score,2e-8,"finite field log derivative");
        for(int j:{0,M/2,M}) {
            const double next=path.coordinate(j,0)+.031;
            const auto changed=path.withCoordinate(j,0,next);const auto q=freeConditional(path,j,0);
            check(changed.logWeight()-path.logWeight(),path.localLogRatio(j,0,next),2e-12,"local ratio");
            check(path.localLogRatio(j,0,next)+q.logDensity(path.coordinate(j,0))-q.logDensity(next),0,2e-12,"forced conditional cancellation");
        }
        double cov=0;
        for(int j=0;j<path.nodes();++j)cov+=path.sourceNodeWeight(j)*std::exp(-p.omega*path.dt()*std::abs(j-path.center()))/(2*path.alpha());
        const double chi=-std::expm1(-p.omega*path.totalTau()/2)/k;
        check(cov,chi,2e-14,"free response covariance");
        check((forced(h+dh).freeMean(path.center(),0)-forced(h-dh).freeMean(path.center(),0))/(2*dh),chi,2e-12,"free finite difference");
        check(path.freeLogNorm()-zero.freeLogNorm(),h*.17*(-std::expm1(-p.omega*path.totalTau()))/p.omega+h*h/(2*k)*(path.totalTau()+std::expm1(-p.omega*path.totalTau())/p.omega),2e-14,"forced norm");
        std::stringstream ss;path.save(ss);const auto replay=PhononPath::load(ss);
        check(replay.logWeight(),path.logWeight(),0,"forced snapshot");
    }
    HoneycombLattice lat(5,5);VBSField field(lat,3,2,.4,.05);
    check(field.profileMeanSquare(),.5,2e-14,"cos normalization");
    auto fp=PhononPath(lat.bondCount(),80,.1,p,{},field.forces());
    for(int j=0;j<fp.nodes();++j)for(int b=0;b<fp.bonds();++b)fp=fp.withCoordinate(j,b,fp.freeMean(j,b));
    const auto moments=field.measure(fp);
    check(moments[0].real(),.05*.5*(-std::expm1(-p.omega*fp.totalTau()/2))/p.spring(),2e-14,"selected free mean");
    for(double g:{0.,.3}) {
        KitaevSSHModel model(HoneycombLattice(2,2),p2Couplings(g));VBSField f(model.lattice(),1,0,.3,.05);
        auto path=PhononPath(model.lattice().bondCount(),4,.1,p,{},f.forces()).freeSample(rng);
        CoupledPhononPath state(path,freeGroundState(model),honeycombGenerator(model),4);
        PhononSampler sampler(state,rng);
        const auto local=sampler.attempt(2,2,path.coordinate(2,2)+.02,.5);
        check(local.log_ratio,local.general_log_ratio,2e-10,"coupled local balance");
        const auto &current=sampler.state().path;
        auto drawn=PhononPath(1,4,.1,p,{current.trialShift(2)},{current.force(2)}).freeSample(rng);
        std::vector<double> values;for(int j=0;j<drawn.nodes();++j)values.push_back(drawn.coordinate(j,0));
        const auto bond=sampler.attemptBond(2,values,.5);
        check(bond.log_ratio,bond.general_log_ratio,2e-10,"coupled bond balance");
        if(g==0.) {check(local.log_ratio,0,2e-10,"free local acceptance");check(bond.log_ratio,0,2e-10,"free bond acceptance");}
    }
    std::cout<<"PASS forced kernel, finite-field score, free response/norm, local and bond detailed balance, source normalization and snapshot\n";
    return 0;
} catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
