#include "phonon_mc.h"
#include <iostream>
#include <map>
#include <sstream>
using namespace kitaev_ssh;
namespace {
std::map<std::string,std::pair<double,double>> errors;
void check(const std::string &name,double e,double tol) {
    auto &r=errors[name];r.first=std::max(r.first,e);r.second=tol;
    if(!std::isfinite(e)||e>tol)throw std::runtime_error(name+": "+std::to_string(e));
}
void require(bool b,const std::string &n){check(n,b?0.:1.,0.);}
template<class F> void rejects(F f) { bool yes=false;try{f();}catch(const std::exception&){yes=true;}require(yes,"invalid_input"); }
struct Gaussian { Eigen::MatrixXd q,c; Eigen::VectorXd s,mean; double logZ; };
// Extract quadratic coefficients from evaluations, independently of link formulas.
template<class F> Gaussian integrate(const PhononPath &p,F logWeight) {
    const int n=p.nodes(); const double l0=logWeight(p);
    Gaussian a; a.q=Eigen::MatrixXd::Zero(n,n);a.s.resize(n);
    Eigen::VectorXd plus(n);
    for(int j=0;j<n;++j) {
        plus(j)=logWeight(p.withCoordinate(j,0,1.));
        const double minus=logWeight(p.withCoordinate(j,0,-1.));
        a.s(j)=(plus(j)-minus)/2;a.q(j,j)=2*l0-plus(j)-minus;
    }
    for(int j=0;j<n;++j)for(int k=j+1;k<n;++k)
        a.q(j,k)=a.q(k,j)=-(logWeight(p.withCoordinate(j,0,1.).withCoordinate(k,0,1.))-plus(j)-plus(k)+l0);
    Eigen::LLT<Eigen::MatrixXd> llt(a.q);require(llt.info()==Eigen::Success,"precision_positive");
    a.c=llt.solve(Eigen::MatrixXd::Identity(n,n));a.mean=a.c*a.s;
    const Eigen::MatrixXd l=llt.matrixL();
    a.logZ=l0+n/2.*std::log(2*pi)-l.diagonal().array().log().sum()+a.s.dot(a.mean)/2;
    return a;
}
Eigen::MatrixXd covariance(const PhononPath &p) {
    Eigen::MatrixXd c(p.nodes(),p.nodes());
    for(int j=0;j<p.nodes();++j)for(int k=0;k<p.nodes();++k)
        c(j,k)=std::exp(-p.parameters().omega*p.dt()*std::abs(j-k))/(2*p.parameters().mass*p.parameters().omega);
    return c;
}
MajoranaGenerator generator(double g) {
    return [g](const std::vector<double>& x){MatType k=MatType::Zero(2,2);k(0,1)=DataType(0,-2*(.5+g*x[0]));k(1,0)=-k(0,1);return k;};
}
void freeChecks() {
    std::mt19937_64 r0(7211),r1(7211);
    PhononPath zero(2,8,.2,{1.7,.8}), shifted(2,8,.2,{1.7,.8},{.3,-.2});
    const auto a=zero.freeSample(r0),b=PhononPath(2,8,.2,{1.7,.8},{0,0}).freeSample(r1);
    check("d0_weight_compatibility",std::abs(a.logWeight()-b.logWeight()),0.);
    for(int j=0;j<a.nodes();++j)for(int k=0;k<2;++k)check("d0_rng_compatibility",std::abs(a.coordinate(j,k)-b.coordinate(j,k)),0.);
    const auto x=shifted.freeSample(r0);
    auto old=zero;for(int j=0;j<x.nodes();++j)for(int k=0;k<2;++k)old=old.withCoordinate(j,k,x.coordinate(j,k));
    double expected=0,density=0;
    for(int k=0;k<2;++k) {
        double d=x.trialShift(k);expected+=1.7*.8*(d*(x.coordinate(0,k)+x.coordinate(8,k))-d*d);
        density+=freeBondLogDensity(x,k);
    }
    check("endpoint_weight_shift",std::abs(x.logWeight()-old.logWeight()-expected),2e-13);
    check("normalized_AR_density",std::abs(density-(x.logWeight()-x.freeLogNorm())),2e-13);
    for(int j=0;j<x.nodes();++j)for(int k=0;k<2;++k) {
        const double v=x.coordinate(j,k)+.23;const auto y=x.withCoordinate(j,k,v);
        check("local_log_ratio",std::abs(x.localLogRatio(j,k,v)-(y.logWeight()-x.logWeight())),2e-13);
        const auto q=freeConditional(x,j,k);
        check("conditional_cancellation",std::abs(y.logWeight()-x.logWeight()+q.logDensity(x.coordinate(j,k))-q.logDensity(v)),2e-13);
    }
    std::stringstream ss; x.save(ss); const auto restored=PhononPath::load(ss);
    require(restored.trialShifts()==x.trialShifts(),"v2_shifts");
    check("v2_weight",std::abs(restored.logWeight()-x.logWeight()),0.);
    std::stringstream legacy;legacy<<std::setprecision(17)<<"KITAEV_SSH_PHONON_PATH 1\n2 8 .2 1.7 .8\n";
    for(int j=0;j<a.nodes();++j)for(double v:a.node(j))legacy<<v<<' ';
    const auto oldLoad=PhononPath::load(legacy);
    require(oldLoad.trialShifts()==std::vector<double>({0,0}),"v1_zero_shifts");
    check("v1_weight",std::abs(oldLoad.logWeight()-a.logWeight()),0.);
    rejects([]{PhononPath p(2,8,.2,{1.7,.8},{.1});});
    rejects([]{PhononPath p(1,8,.2,{1.7,.8},{std::numeric_limits<double>::infinity()});});
    for(double d:{0.,.3,-.7}) {
        const PhononPath p(1,8,.2,{1.7,.8},{d});const auto c=covariance(p);
        const auto integral=integrate(p,[](const PhononPath &z){return z.logWeight();});
        Eigen::VectorXd mu(p.nodes());for(int j=0;j<p.nodes();++j)mu(j)=p.freeMean(j,0);
        check("free_mean",(integral.mean-mu).cwiseAbs().maxCoeff(),2e-12);
        check("free_covariance",(integral.c-c).cwiseAbs().maxCoeff(),2e-12);
        check("free_raw_second_moment",(integral.c+integral.mean*integral.mean.transpose()-c-mu*mu.transpose()).cwiseAbs().maxCoeff(),2e-12);
        check("free_integrated_norm",std::abs(integral.logZ-p.freeLogNorm()),2e-12);
        for(int s=0;s<20;++s) {
            std::mt19937_64 ra(12+s),rb(12+s);
            auto px=p.freeSample(ra),zx=PhononPath(1,8,.2,{1.7,.8}).freeSample(rb);
            for(int j=0;j<p.nodes();++j)check("sampler_shifted_law",std::abs(px.coordinate(j,0)-zx.coordinate(j,0)-mu(j)),3e-16);
        }
    }
    // Reproduce the sweep RNG draws, including permutation; verifies each bond's d.
    auto gen=generator(0.);auto trial=GaussianTrialState::fromMajoranaHamiltonian(gen({0}));
    CoupledPhononPath state(shifted,trial,gen);PhononSampler sampler(state,std::mt19937_64(88),PhononProposal::BondRefresh);
    std::mt19937_64 reference(88);std::vector<int> order{0,1};std::shuffle(order.begin(),order.end(),reference);
    auto expectedPath=shifted;
    for(int bond:order) {
        auto p=PhononPath(1,8,.2,{1.7,.8},{shifted.trialShift(bond)}).freeSample(reference);
        std::vector<double> v(p.nodes());for(int j=0;j<p.nodes();++j)v[j]=p.coordinate(j,0);
        expectedPath=expectedPath.withBondCoordinates(bond,v);reference();
    }
    sampler.sweep();
    for(int j=0;j<shifted.nodes();++j)for(int k=0;k<2;++k)
        check("sweep_inherits_bond_shift",std::abs(sampler.state().path.coordinate(j,k)-expectedPath.coordinate(j,k)),0.);
}
void coupledChecks() {
    for(int eta:{-1,1}) {
        const double g=.6,J=.5;auto gen=generator(g);
        auto trial=GaussianTrialState::fromMajoranaHamiltonian(-double(eta)*gen({0}));
        const PhononPath p(1,8,.2,{1.7,.8},{.23});const auto c=covariance(p);
        auto weight=[&](const PhononPath &x){return evaluatePath(x,trial,pathFactors(x,2,gen),4);};
        auto integral=integrate(p,[&](const PhononPath &x){return weight(x).logAbsTotal();});
        Eigen::VectorXd mu(p.nodes()),w(p.nodes());
        for(int j=0;j<p.nodes();++j){mu(j)=p.freeMean(j,0);w(j)=p.dt()*((j==0||j==p.intervals())?.5:1.);}
        const double f=eta*g;
        const double logz=p.freeLogNorm()-eta*J*p.totalTau()-f*w.dot(mu)+f*f/2*w.dot(c*w);
        check("coupled_norm_eta_both",std::abs(integral.logZ-logz),5e-12);
        check("coupled_mean_eta_both",(integral.mean-(mu-f*c*w)).cwiseAbs().maxCoeff(),5e-12);
        check("coupled_covariance_eta_both",(integral.c-c).cwiseAbs().maxCoeff(),5e-12);
        CoupledPhononPath state(p,trial,gen);PhononSampler sampler(state,std::mt19937_64(9),PhononProposal::BondRefresh);
        for(bool accept:{false,true}) {
            const auto before=state.state();std::vector<double> values(p.nodes(),accept?-eta*.2:eta*2.);
            const double expectedMatter=-f*(w.sum()*values[0]-[&]{double s=0;for(int j=0;j<p.nodes();++j)s+=w(j)*before.path.coordinate(j,0);return s;}());
            const auto result=sampler.attemptBond(0,values,accept?.1:.99);
            require(result.accepted==accept,"accept_reject_decision");
            check("proposal_cancellation",std::abs(result.log_ratio-result.general_log_ratio),3e-12);
            check("single_Fock_ratio",std::abs(result.log_ratio-expectedMatter),3e-12);
            // Accepted eager bond proposals leave Green invalid until measured.
            sampler.ensureMeasured();
            require(state.state().weight.green_valid&&state.state().weight.phase_valid,
                    "committed_measurement_certified");
            const auto &now=state.state();require(now.path.trialShifts()==p.trialShifts(),"transaction_preserves_shift");
            for(int j=0;j<p.nodes();++j)check("commit_rollback_coordinates",std::abs(now.path.coordinate(j,0)-(accept?values[j]:before.path.coordinate(j,0))),0.);
            const auto full=state.fullReference(now.path);
            check("committed_weight_reference",std::abs(full.weight.logAbsTotal()-now.weight.logAbsTotal()),2e-12);
            check("committed_green_reference",(full.weight.center_green-now.weight.center_green).cwiseAbs().maxCoeff(),2e-12);
            check("physical_phase",std::abs(now.weight.matter_phase-DataType(1)),2e-12);
        }
    }
}
}
int main(){try{freeChecks();coupledChecks();for(auto &r:errors)std::cout<<std::setprecision(17)<<r.first<<" max_error="<<r.second.first<<" tolerance="<<r.second.second<<'\n';std::cout<<"SHIFTED_TRIAL_PASS\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
