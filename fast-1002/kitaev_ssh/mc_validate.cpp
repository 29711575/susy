#include "phonon_mc.h"
#include <iostream>
#include <map>

using namespace kitaev_ssh;
namespace {
struct Metric { double error=0.,tolerance=0.; };
std::map<std::string,Metric> metrics;
void check(const std::string &name,double error,double tolerance) {
    auto &m=metrics[name]; m.error=std::max(m.error,error); m.tolerance=tolerance;
    if(!std::isfinite(error) || error>tolerance) throw std::runtime_error(name+" error="+std::to_string(error));
}
void require(bool ok,const std::string &name) { check(name,ok?0.:1.,0.); }
template<class F> void rejects(F f) { bool failed=false; try{f();}catch(const std::exception&){failed=true;}require(failed,"expected_failure"); }
MajoranaGenerator singleGenerator(double g) {
    return [g](const std::vector<double> &x) {
        MatType k=MatType::Zero(2,2); k(0,1)=DataType(0,-2*(.5+g*x.at(0))); k(1,0)=-k(0,1); return k;
    };
}
// Explicit two-state Fock half-steps: cA=sigma_x, cB=-sigma_y,
// H=diag(J+gx,-J-gx). Occupied |1> has eta=-1.
double fockLog(const PhononPath &p,double g) {
    MatType u=MatType::Identity(2,2);
    auto half=[&](double x) { MatType f=MatType::Zero(2,2); const double t=p.dt()/2*(.5+g*x); f(0,0)=std::exp(-t);f(1,1)=std::exp(t);return f; };
    for(int l=1;l<=p.intervals();++l) u=half(p.coordinate(l,0))*half(p.coordinate(l-1,0))*u;
    return std::log(u(1,1).real());
}
void sameState(const CoupledPhononPath::Snapshot &a,const CoupledPhononPath::Snapshot &b) {
    for(int j=0;j<a.path.nodes();++j) check("state_coordinates",std::abs(a.path.coordinate(j,0)-b.path.coordinate(j,0)),0.);
    for(std::size_t f=0;f<a.factors.size();++f) check("state_factors",(MatType(a.factors[f].matrix)-MatType(b.factors[f].matrix)).norm(),0.);
    check("state_log_weight",std::abs(a.weight.logAbsTotal()-b.weight.logAbsTotal()),2e-12);
    check("state_phase",std::abs(a.weight.matter_phase-b.weight.matter_phase),2e-12);
    if(a.weight.green_valid && b.weight.green_valid)
        check("state_green",(a.weight.center_green-b.weight.center_green).norm(),2e-12);
    else require(!a.weight.green_valid && !b.weight.green_valid,"deferred_state_green");
}
void acceptance(PhononWeightOptions options={}) {
    const auto generator=singleGenerator(.6); const auto trial=GaussianTrialState::fromMajoranaHamiltonian(generator({0.}));
    std::mt19937_64 rng(33); const auto path=PhononPath(1,4,.4,{1.7,.8}).freeSample(rng);
    CoupledPhononPath state(path,trial,generator,4,options); PhononSampler sampler(state,rng);
    for(int j:{0,2,4}) for(bool accept:{false,true}) {
        const auto before=state.state(); const double next=before.path.coordinate(j,0)+(accept?.3:-2.);
        const auto candidate=before.path.withCoordinate(j,0,next); const auto q=freeConditional(before.path,j,0);
        const double reverseQ=q.logDensity(before.path.coordinate(j,0))-q.logDensity(next);
        check("conditional_phonon_cancellation",std::abs(before.path.localLogRatio(j,0,next)+reverseQ),2e-12);
        const auto move=sampler.attempt(j,0,next,accept?.5:.99);
        require(move.accepted==accept,"actual_accept_reject");
        check("general_vs_simplified_MH",std::abs(move.log_ratio-move.general_log_ratio),2e-12);
        check("Fock_acceptance_log_ratio",std::abs(move.log_ratio-(fockLog(candidate,.6)-fockLog(before.path,.6))),2e-12);
        const double oldRho=before.path.logWeight()+fockLog(before.path,.6),newRho=candidate.logWeight()+fockLog(candidate,.6);
        const double forward=oldRho+q.logDensity(next)+std::min(0.,move.log_ratio);
        const double backward=newRho+q.logDensity(before.path.coordinate(j,0))+std::min(0.,-move.log_ratio);
        check("MH_detailed_balance",std::abs(forward-backward),2e-12);
        sameState(state.state(),accept?state.fullReference(candidate):before);
    }
    require(sampler.counts().proposed==6 && sampler.counts().accepted==3,"proposal_counters");
    const auto before=state.state();
    rejects([&]{sampler.attempt(0,0,std::numeric_limits<double>::infinity(),.5);});
    sameState(state.state(),before);
    require(sampler.counts().proposed==6,"failure_not_rejection");
    rejects([&]{sampler.attempt(0,0,0.,1.);});
}
void refreshAcceptance(PhononWeightOptions options={}) {
    const auto generator=singleGenerator(.6); const auto trial=GaussianTrialState::fromMajoranaHamiltonian(generator({0.}));
    std::mt19937_64 rng(992);
    CoupledPhononPath state(PhononPath(1,8,.1,{1.7,.8}).freeSample(rng),trial,generator,4,options);
    PhononSampler sampler(state,rng,PhononProposal::BondRefresh);
    for(bool accept:{false,true}) {
        const auto before=state.state();std::vector<double> values(before.path.nodes());
        for(int j=0;j<before.path.nodes();++j)values[j]=before.path.coordinate(j,0)+(accept?.3:-2.);
        const auto candidate=before.path.withBondCoordinates(0,values);
        check("refresh_AR_endpoint_density",std::abs(freeBondLogDensity(before.path,0)-before.path.logWeight()
              -before.path.parameters().omega*before.path.totalTau()/2),2e-12);
        const auto move=sampler.attemptBond(0,values,accept?.5:.99);
        require(move.accepted==accept,"refresh_actual_accept_reject");
        check("refresh_general_MH",std::abs(move.log_ratio-move.general_log_ratio),2e-12);
        check("refresh_Fock_log_ratio",std::abs(move.log_ratio-(fockLog(candidate,.6)-fockLog(before.path,.6))),2e-12);
        const double forward=before.path.logWeight()+fockLog(before.path,.6)+freeBondLogDensity(candidate,0)+std::min(0.,move.log_ratio);
        const double backward=candidate.logWeight()+fockLog(candidate,.6)+freeBondLogDensity(before.path,0)+std::min(0.,-move.log_ratio);
        check("refresh_detailed_balance",std::abs(forward-backward),2e-12);
        sameState(state.state(),accept?state.fullReference(candidate,false,false):before);
        if(accept) {
            const auto &measured=state.ensureMeasured();
            const auto expected=state.fullReference(candidate);
            sameState(measured,expected);
        }
    }
    const auto before=state.state();
    rejects([&]{sampler.attemptBond(0,{1.,2.},.5);});sameState(state.state(),before);
    auto bad=std::vector<double>(before.path.nodes(),0.);bad.back()=std::numeric_limits<double>::infinity();
    rejects([&]{sampler.attemptBond(0,bad,.5);});sameState(state.state(),before);
    require(sampler.counts().proposed==2 && sampler.counts().accepted==1,"refresh_failures_not_rejections");
    const auto multi=PhononPath(3,4,.2,{1.7,.8}).freeSample(rng);
    const auto changed=multi.withBondCoordinates(1,std::vector<double>(multi.nodes(),.17));
    for(int j=0;j<multi.nodes();++j)for(int b=0;b<3;++b)
        check("refresh_other_bonds_untouched",std::abs(changed.coordinate(j,b)-(b==1?.17:multi.coordinate(j,b))),0.);
}
void sampledSingleBond(double g,PhononProposal mode=PhononProposal::LocalConditional) {
    const int m=4; const double h=.4,mass=1.7,omega=.8;
    const auto generator=singleGenerator(g); const auto trial=GaussianTrialState::fromMajoranaHamiltonian(generator({0.}));
    std::mt19937_64 rng(g==0?481516:234210);
    CoupledPhononPath state(PhononPath(1,m,h,{mass,omega}).freeSample(rng),trial,generator);
    PhononSampler sampler(state,rng,mode);
    Eigen::MatrixXd c(m+1,m+1); Eigen::VectorXd w(m+1);
    for(int j=0;j<=m;++j) { w(j)=h*((j==0||j==m)?.5:1.); for(int k=0;k<=m;++k)c(j,k)=std::exp(-omega*h*std::abs(j-k))/(2*mass*omega); }
    const Eigen::VectorXd mean=g*c*w; // eta=-1
    std::vector<std::pair<int,int>> pairs; for(int j=0;j<=m;++j)pairs.push_back({j,j});
    pairs.push_back({0,2});pairs.push_back({2,4});pairs.push_back({0,4});
    std::vector<std::string> names; for(int j=0;j<=m;++j)names.push_back("mean"+std::to_string(j));
    for(auto pair:pairs)names.push_back("cov"+std::to_string(pair.first)+std::to_string(pair.second));
    PhaseStatistics stats(names,100); const SamplingSchedule schedule{500,6400,2,100};
    sample(sampler,schedule,stats,[&](const auto &s) {
        std::vector<Complex> out;
        for(int j=0;j<=m;++j) out.push_back(s.path.coordinate(j,0));
        for(auto pair:pairs) out.push_back((s.path.coordinate(pair.first,0)-mean(pair.first))*(s.path.coordinate(pair.second,0)-mean(pair.second)));
        return out;
    },[](const auto&){});
    const auto result=stats.summary(),doubled=stats.doubledBlockSummary();
    require(result.samples==3200 && result.blocks==32 && result.jackknife_defined && result.phase_resolved,"measurement_schedule");
    require(sampler.counts().proposed==std::uint64_t(6900*(mode==PhononProposal::LocalConditional?5:1)),"sweep_covers_all_nodes");
    if(g==0) require(sampler.counts().accepted==sampler.counts().proposed,"g0_all_accepted");
    else require(sampler.counts().accepted<sampler.counts().proposed,"coupled_rejected_states_measured");
    for(std::size_t j=0;j<names.size();++j) {
        const double exact=j<=std::size_t(m)?mean(j):c(pairs[j-m-1].first,pairs[j-m-1].second);
        const auto &v=result.observables[j];
        require(v.error_real>0 && std::isfinite(v.error_real),"finite_sampling_error");
        check(g==0?"g0_sampling_zscore":"coupled_sampling_zscore",std::abs(v.value.real()-exact)/v.error_real,6.);
        check("sampled_imaginary_parts",std::abs(v.value.imag()),2e-12);
    }
    std::cout<<"SAMPLING g="<<g<<" acceptance="<<sampler.counts().acceptance()<<" center_X="<<result.observables[2].value.real()
        <<" stderr="<<result.observables[2].error_real<<" exact="<<mean(2)<<" center_cov="<<result.observables[7].value.real()
        <<" cov_stderr="<<result.observables[7].error_real<<" cov_exact="<<c(2,2)
        <<" doubled_block_center_stderr="<<doubled.observables[2].error_real<<'\n';
}
void phaseDiagnosticCounts() {
    PhaseStatistics stats({"identity"},2);
    const std::vector<Complex> phases{1.,-1.,std::polar(1.,.2),std::polar(1.,-2.),
        std::polar(1.,phaseRealityTolerance/2),std::polar(1.,std::acos(-1.)-phaseRealityTolerance/2),1.,1.};
    Complex total=0.;
    for(auto z:phases) {stats.add(z,{1.});total+=z;}
    const auto s=stats.summary(),merged=stats.doubledBlockSummary();
    require(s.samples==8 && s.positive==4 && s.negative==2 && s.nonreal==2,"phase_classification_counts");
    require(merged.positive==s.positive && merged.negative==s.negative && merged.nonreal==s.nonreal,"phase_counts_reblock");
    check("phase_not_rounded",std::abs(s.phase.value-total/8.),2e-15);
    check("phase_classification_identity_ratio",std::abs(s.observables[0].value-Complex(1.)),2e-15);
    std::vector<Complex> replicas;
    for(const auto &b:stats.blocks()) replicas.emplace_back(std::abs((total-b.phase)/6.),0.);
    ComplexEstimate expected;jackknifeErrors(expected,replicas);
    check("phase_modulus_joint_jackknife",std::abs(expected.error_real-s.phase_magnitude_error),2e-15);
}
void complexAccumulation() {
    const std::vector<double> angles{-.9,-.4,.2,.8,1.2};
    const std::vector<Complex> values{{-1.,.4},{.2,-.7},{.9,.3},{1.4,-.1},{-.5,.8}};
    PhaseStatistics stats({"complex_observable","identity"},4);
    std::vector<Complex> phases,weighted;
    Complex exactD=0.,exactN=0.;
    for(int j=0;j<5;++j){exactD+=std::polar(1.,angles[j])/5.;exactN+=std::polar(1.,angles[j])*values[j]/5.;}
    for(int r=0;r<40;++r) {
        int j=(r*3)%5; const auto s=std::polar(1.,angles[j]); phases.push_back(s);weighted.push_back(s*values[j]);stats.add(s,{values[j],1.});
    }
    const auto result=stats.summary();
    require(result.phase_resolved && result.jackknife_defined,"complex_denominator_resolved");
    check("complex_exact_phase",std::abs(result.phase.value-exactD),2e-14);
    check("complex_exact_ratio",std::abs(result.observables[0].value-exactN/exactD),2e-14);
    check("joint_jackknife_identity_error",std::hypot(result.observables[1].error_real,result.observables[1].error_imag),2e-14);
    // Independent delete-four enumeration from individual samples, not block API.
    std::vector<Complex> loo;
    for(int b=0;b<10;++b){Complex n=0.,d=0.;for(int r=0;r<40;++r)if(r/4!=b){n+=weighted[r];d+=phases[r];}loo.push_back(n/d);}
    Complex avg=0.;for(auto z:loo)avg+=z/10.;double vr=0.,vi=0.,cov=0.;
    for(auto z:loo){z-=avg;vr+=z.real()*z.real();vi+=z.imag()*z.imag();cov+=z.real()*z.imag();}
    check("complex_jackknife_error",std::abs(result.observables[0].error_real-std::sqrt(.9*vr)),2e-14);
    check("complex_jackknife_error",std::abs(result.observables[0].error_imag-std::sqrt(.9*vi)),2e-14);
    check("complex_jackknife_covariance",std::abs(result.observables[0].covariance-.9*cov),2e-14);
    Complex bad=0.;for(const auto &b:stats.blocks())bad+=b.numerator[0]/b.phase/10.;
    require(std::abs(bad-result.observables[0].value)>1e-4,"not_mean_of_block_ratios");
    PhaseStatistics zero({"x"},2);for(int n=0;n<8;++n)zero.add(n%2?Complex(-1):Complex(1),{1.});
    require(!zero.summary().ratio_defined && !zero.summary().phase_resolved,"zero_phase_denominator_flagged");
    PhaseStatistics few({"x"},2);few.add(1.,{2.});few.add(1.,{3.});
    require(few.summary().too_few_blocks && !few.summary().jackknife_defined,"insufficient_blocks_flagged");
    few.add(1.,{1.});rejects([&]{few.summary();});
    rejects([&]{few.add(Complex(2.),{1.});});
    PhaseStatistics unresolved({"x"},1);
    for(int n=0;n<8;++n) unresolved.add(Complex(.02,(n<4?1:-1)*std::sqrt(1-.02*.02)),{double(n)});
    require(unresolved.summary().ratio_defined && !unresolved.summary().phase_resolved,"unresolved_nonzero_phase_flagged");
    PhaseStatistics singularLeaveOne({"x"},2);
    for(double sign:{1.,-1.,1.}) for(int n=0;n<2;++n)singularLeaveOne.add(sign,{1.});
    require(singularLeaveOne.summary().ratio_defined && !singularLeaveOne.summary().jackknife_defined,"singular_leave_one_flagged");
    std::cout<<"COMPLEX exact_phase="<<exactD<<" exact_ratio="<<exactN/exactD<<'\n';
}
}
int main(){
    std::cout<<std::setprecision(17);
    try {
        for(const auto options:{PhononWeightOptions{false,false},PhononWeightOptions{true,false},
                                PhononWeightOptions{false,true},PhononWeightOptions{true,true}}) {
            acceptance(options);refreshAcceptance(options);
        }
        phaseDiagnosticCounts();complexAccumulation();sampledSingleBond(0.);sampledSingleBond(.6);
        sampledSingleBond(0.,PhononProposal::BondRefresh);sampledSingleBond(.6,PhononProposal::BondRefresh);
        for(const auto &m:metrics)std::cout<<"PASS "<<m.first<<" max_error="<<m.second.error<<" tolerance="<<m.second.tolerance<<'\n';
        return 0;
    }catch(const std::exception &e){std::cerr<<"FAIL P3: "<<e.what()<<'\n';return 1;}
}
