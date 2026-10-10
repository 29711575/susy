#include "phonon_mc.h"
#include <array>
#include <iomanip>
#include <iostream>

using namespace kitaev_ssh;

namespace {
struct ReplayMetrics {
    double maxLog=0.,maxPhase=0.,maxGreen=0.,maxRatio=0.,maxAcceptance=0.;
    int accepts=0,rejects=0,measurements=0;
};

void require(bool condition,const std::string &message) {
    if(!condition) throw std::runtime_error(message);
}

// Independent spin-lift oracle: rotate the first Majorana of each disjoint
// fermion pair by 2*theta. Each pair contributes cos(theta) to the vacuum
// amplitude; Q_trial^dag Q has two singular values cos(theta). Its determinant
// can be below 1e-14 while every singular value and the solve are well trusted.
void fullRankSmallDeterminantCheck() {
    const int n=36;
    MatType phi=MatType::Zero(2*n,n);
    phi.topRows(n)=MatType::Identity(n,n)/std::sqrt(2.);
    phi.bottomRows(n)=DataType(0.,1.)*MatType::Identity(n,n)/std::sqrt(2.);
    const auto trial=GaussianTrialState::fromPhi(phi);
    const double theta=std::acos(.3);
    MatType b=MatType::Identity(2*n,2*n);
    for(int j=0;j<n;j+=2) {
        b(j,j)=b(j+1,j+1)=std::cos(2*theta);
        b(j,j+1)=-std::sin(2*theta);b(j+1,j)=std::sin(2*theta);
    }
    require((b*b.transpose()-MatType::Identity(2*n,2*n)).norm()<1e-12,
            "small-determinant oracle must be a canonical Majorana rotation");
    const double expectedLog=(n/2)*std::log(std::cos(theta));
    for(bool accelerated:{false,true})for(bool valuesOnly:{false,true}) {
        PureFastOptions options;options.accelerated_linear_algebra=accelerated;
        options.reference_singular_values_only=valuesOnly;
        const std::vector<PureProjectorSlice> slices{{b,DataType(1.,0.),"pair-rotation"}};
        const auto magnitude=pureProjectorMagnitudeWeight(trial,slices,options);
        const auto signedWeight=pureProjectorStableReferenceWeight(trial,slices,options);
        require(magnitude.ok()&&signedWeight.ok(),"full-rank small determinant was rejected");
        require(std::abs(magnitude.log_abs_weight-expectedLog)<2e-10&&
                std::abs(signedWeight.log_abs_weight-expectedLog)<2e-10,
                "small-determinant magnitude disagrees with exact pair amplitude");
        require(std::abs(signedWeight.complex_phase-DataType(1.,0.))<2e-10,
                "small-determinant spin-lift phase");
        const std::vector<PureProjectorSlice> zero{{MatType::Zero(2*n,2*n),1.,"zero"}};
        require(!pureProjectorMagnitudeWeight(trial,zero,options).ok()&&
                !pureProjectorStableReferenceWeight(trial,zero,options).ok(),
                "singular propagation must still fail closed");
    }
}

void compareMagnitude(const PathWeight &a,const PathWeight &b,ReplayMetrics &metrics) {
    metrics.maxLog=std::max(metrics.maxLog,std::abs(a.log_abs_matter-b.log_abs_matter));
}

void compareMeasured(const PathWeight &a,const PathWeight &b,ReplayMetrics &metrics) {
    compareMagnitude(a,b,metrics);
    metrics.maxPhase=std::max(metrics.maxPhase,std::abs(a.matter_phase-b.matter_phase));
    if(a.center_green.size()!=0 && b.center_green.size()!=0)
        metrics.maxGreen=std::max(metrics.maxGreen,
            (a.center_green-b.center_green).cwiseAbs().maxCoeff());
}

void runReplay(std::uint64_t seed,double g,int intervals,int nprop,bool singularOnly,
               bool merge,ReplayMetrics &metrics) {
    const KitaevSSHModel model(HoneycombLattice(4,4),p2Couplings(g));
    std::mt19937_64 rng(seed);
    const auto path=PhononPath(48,intervals,.1,{1.7,.8},
        std::vector<double>(48,.1378676470588235)).freeSample(rng);
    const auto trial=freeGroundState(model);
    const auto generator=honeycombGenerator(model);
    const PhononWeightOptions eagerOptions{singularOnly,merge,false};
    const PhononWeightOptions deferredOptions{singularOnly,merge,true};
    CoupledPhononPath eager(path,trial,generator,4,eagerOptions);
    CoupledPhononPath deferred(path,trial,generator,4,deferredOptions);

    const auto &initial=deferred.state();
    require(eager.configurationHash()==deferred.configurationHash(),
            "eager/deferred initial configuration identity");
    require(initial.weight.magnitude_valid && !initial.weight.phase_valid &&
            !initial.weight.green_valid,"initial deferred validity flags");
    require(initial.weight.matter_phase==DataType(0) &&
            initial.weight.center_green.size()==0,"initial deferred phase is populated");
    compareMagnitude(eager.state().weight,initial.weight,metrics);
    deferred.ensureMeasured();
    compareMeasured(eager.state().weight,deferred.state().weight,metrics);
    ++metrics.measurements;

    // Keep the proposal sequence deterministic while forcing both transaction
    // outcomes. A tiny displacement is accepted with the smallest legal u;
    // the larger candidates are searched until one has a negative log ratio,
    // then u just below one forces a rejection.
    const std::array<double,8> amplitudes{{.001,-.001,.01,-.01,.1,-.1,.5,-.5}};
    for(int n=0;n<nprop;++n) {
        const bool bondProposal=(n%2)==1;
        const bool wantAccept=(n%2)==1;
        const int bond=n%48;
        const int node=(n%3==0)?intervals/2:(n%3==1?0:intervals);
        bool selected=false;
        for(double amplitude:amplitudes) {
            const auto &old=eager.state();
            std::vector<double> values;
            std::unique_ptr<CoupledPhononPath::Proposal> eagerProposal;
            std::unique_ptr<CoupledPhononPath::Proposal> deferredProposal;
            try {
                if(bondProposal) {
                    values.resize(old.path.nodes());
                    for(int j=0;j<old.path.nodes();++j)
                        values[j]=old.path.coordinate(j,bond)+amplitude;
                    eagerProposal=std::make_unique<CoupledPhononPath::Proposal>(
                        eager.proposeBondPath(bond,values));
                    deferredProposal=std::make_unique<CoupledPhononPath::Proposal>(
                        deferred.proposeBondPath(bond,values));
                } else {
                    const double value=old.path.coordinate(node,bond)+amplitude;
                    eagerProposal=std::make_unique<CoupledPhononPath::Proposal>(
                        eager.propose(node,bond,value));
                    deferredProposal=std::make_unique<CoupledPhononPath::Proposal>(
                        deferred.propose(node,bond,value));
                }
            } catch(...) {
                if(eagerProposal) eager.reject(*eagerProposal);
                throw;
            }
            const auto &a=eagerProposal->candidate().weight;
            const auto &b=deferredProposal->candidate().weight;
            require(eagerProposal->candidate().configuration_hash!=eager.configurationHash() &&
                    deferredProposal->candidate().configuration_hash!=deferred.configurationHash(),
                    "candidate configuration identity was not changed");
            require(eagerProposal->candidate().configuration_hash==
                    deferredProposal->candidate().configuration_hash,
                    "eager/deferred candidate configuration identity");
            require(b.magnitude_valid&&!b.phase_valid&&!b.green_valid,
                    "deferred candidate validity flags");
            compareMagnitude(a,b,metrics);
            const double ra=a.log_abs_matter-eager.state().weight.log_abs_matter;
            const double rb=b.log_abs_matter-deferred.state().weight.log_abs_matter;
            metrics.maxRatio=std::max(metrics.maxRatio,std::abs(ra-rb));
            metrics.maxAcceptance=std::max(metrics.maxAcceptance,
                std::abs(std::exp(std::min(0.,ra))-std::exp(std::min(0.,rb))));
            bool phaseReadFailed=false;
            try { (void)deferredProposal->phaseRatio(); }
            catch(const std::logic_error &) { phaseReadFailed=true; }
            require(phaseReadFailed,"deferred phaseRatio accepted invalid cache");

            const double u=wantAccept?.5*0x1p-52:1.-.5*0x1p-52;
            const bool acceptA=std::log(u)<std::min(0.,ra);
            const bool acceptB=std::log(u)<std::min(0.,rb);
            require(acceptA==acceptB,"eager/deferred decision mismatch");
            if(acceptA!=wantAccept) {
                eager.reject(*eagerProposal);
                deferred.reject(*deferredProposal);
                continue;
            }
            const bool wasMeasured=deferred.state().weight.phase_valid&&
                                   deferred.state().weight.green_valid;
            const DataType oldPhase=deferred.state().weight.matter_phase;
            const MatType oldGreen=deferred.state().weight.center_green;
            if(acceptA) {
                eager.accept(*eagerProposal);
                deferred.accept(*deferredProposal);
                ++metrics.accepts;
                require(!deferred.state().weight.phase_valid&&
                        !deferred.state().weight.green_valid,
                        "accept did not invalidate deferred caches");
                deferred.ensureMeasured();
                compareMeasured(eager.state().weight,deferred.state().weight,metrics);
                ++metrics.measurements;
                const DataType phase=deferred.state().weight.matter_phase;
                const MatType green=deferred.state().weight.center_green;
                deferred.ensureMeasured();
                require(deferred.state().weight.matter_phase==phase&&
                        (deferred.state().weight.center_green-green).norm()==0.,
                        "repeated measurement changed cache");
            } else {
                eager.reject(*eagerProposal);
                deferred.reject(*deferredProposal);
                ++metrics.rejects;
                require(eager.configurationHash()==deferred.configurationHash(),
                        "reject changed configuration identity");
                require(!wasMeasured || (deferred.state().weight.phase_valid&&
                        deferred.state().weight.green_valid),
                        "reject invalidated current measured cache");
                if(wasMeasured) {
                    require(deferred.state().weight.matter_phase==oldPhase&&
                        (deferred.state().weight.center_green-oldGreen).norm()==0.,
                        "reject changed current measured cache");
                }
            }
            compareMeasured(eager.state().weight,deferred.state().weight,metrics);
            selected=true;
            break;
        }
        require(selected,"could not force requested accept/reject outcome");
    }
    require(metrics.accepts>0&&metrics.rejects>0,"replay did not cover both transactions");
}

void runMeasurementCadence(std::uint64_t seed,int stride) {
    const KitaevSSHModel model(HoneycombLattice(2,2),p2Couplings(.3));
    std::mt19937_64 rng(seed);
    const auto path=PhononPath(12,6,.11,{1.7,.8}).freeSample(rng);
    const auto trial=freeGroundState(model);
    const auto generator=honeycombGenerator(model);
    const PhononWeightOptions deferredOptions{false,false,true};
    CoupledPhononPath eager(path,trial,generator,4);
    CoupledPhononPath deferred(path,trial,generator,4,deferredOptions);
    eager.ensureMeasured(); deferred.ensureMeasured();
    PhaseStatistics eagerStats({"identity"},1),deferredStats({"identity"},1);
    const auto initialHash=deferred.configurationHash();
    const int measurementSweeps=6;
    int measurements=0;
    ReplayMetrics cadenceMetrics;
    for(int sweep=0;sweep<measurementSweeps;++sweep) {
        // Explicitly reject every proposal. This keeps the accepted state
        // fixed while exercising the original fixed measurement schedule.
        const int node=sweep%path.nodes(),bond=(sweep+1)%path.bonds();
        const double value=deferred.state().path.coordinate(node,bond)+.021;
        auto eagerProposal=eager.propose(node,bond,value);
        auto deferredProposal=deferred.propose(node,bond,value);
        eager.reject(eagerProposal); deferred.reject(deferredProposal);
        require(deferred.configurationHash()==initialHash &&
                deferred.state().weight.phase_valid&&deferred.state().weight.green_valid,
                "consecutive reject changed measured deferred state");
        if((sweep+1)%stride==0) {
            eager.ensureMeasured(); deferred.ensureMeasured();
            require(eager.configurationHash()==deferred.configurationHash(),
                    "cadence configuration identity mismatch");
            compareMeasured(eager.state().weight,deferred.state().weight,cadenceMetrics);
            const std::vector<Complex> observable{DataType(1.)};
            require(eagerStats.add(eager.state().weight.matter_phase,observable)==
                    deferredStats.add(deferred.state().weight.matter_phase,observable),
                    "cadence block completion mismatch");
            ++measurements;
        }
    }
    require(measurements==measurementSweeps/stride,"measurement cadence count mismatch");
    const auto &a=eagerStats.blocks(); const auto &b=deferredStats.blocks();
    require(a.size()==b.size() && a.size()==std::size_t(measurements),
            "cadence block count mismatch");
    for(std::size_t i=0;i<a.size();++i) {
        require(a[i].count==b[i].count && a[i].positive==b[i].positive &&
                a[i].negative==b[i].negative && a[i].nonreal==b[i].nonreal,
                "cadence phase classification mismatch");
        require(std::abs(a[i].phase-b[i].phase)<2e-12 &&
                std::abs(a[i].numerator[0]-b[i].numerator[0])<2e-12,
                "cadence phase accumulation mismatch");
    }
    std::cout<<"CADENCE stride="<<stride<<" measurements="<<measurements
             <<" blocks="<<a.size()<<"\n";
}

void failureDoesNotCommit(double g) {
    const KitaevSSHModel model(HoneycombLattice(2,2),p2Couplings(g));
    std::mt19937_64 rng(6107);
    const auto path=PhononPath(12,6,.11,{1.7,.8}).freeSample(rng);
    const auto trial=freeGroundState(model);
    const auto fail=std::make_shared<bool>(false);
    MajoranaGenerator::ValueFunction value=[model,fail](const std::vector<double> &x) {
        if(*fail) throw std::runtime_error("injected deferred measurement failure");
        return model.frameworkGeneratorWith(x);
    };
    MajoranaGenerator::IntoFunction into=[model,fail](const std::vector<double> &x,MatType &out) {
        if(*fail) throw std::runtime_error("injected deferred measurement failure");
        model.frameworkGeneratorWith(x,out);
    };
    CoupledPhononPath state(path,trial,MajoranaGenerator(std::move(value),std::move(into)),4,
                            PhononWeightOptions{false,true,true});
    state.ensureMeasured();
    auto proposal=state.propose(path.center(),0,path.coordinate(path.center(),0)+.07);
    const auto candidateHash=proposal.candidate().configuration_hash;
    state.accept(proposal);
    const auto before=state.state();
    require(state.configurationHash()==candidateHash && !before.weight.phase_valid &&
            !before.weight.green_valid,"accepted deferred state validity before failure");
    *fail=true;
    bool failed=false;
    try { state.ensureMeasured(); } catch(const std::exception &) { failed=true; }
    require(failed,"injected deferred measurement failure was swallowed");
    require(state.configurationHash()==before.configuration_hash &&
            state.state().weight.log_abs_matter==before.weight.log_abs_matter &&
            !state.state().weight.phase_valid&&!state.state().weight.green_valid,
            "measurement failure committed a partial snapshot");
    *fail=false;
    state.ensureMeasured();
    require(state.configurationHash()==before.configuration_hash &&
            state.state().weight.phase_valid&&state.state().weight.green_valid,
            "deferred measurement did not recover after injected failure");
    std::cout<<"FAILURE_INJECTION_OK hash="<<candidateHash<<"\n";
}

void syntheticScalarPhaseCheck(double g) {
    const KitaevSSHModel model(HoneycombLattice(2,2),p2Couplings(g));
    std::mt19937_64 rng(4071);
    const auto path=PhononPath(12,6,.11,{1.7,.8}).freeSample(rng);
    const auto trial=freeGroundState(model);
    const auto generator=honeycombGenerator(model);
    const auto baseFactors=pathFactors(path,model.majoranas(),generator);
    const auto base=evaluatePath(path,trial,baseFactors,4);
    for(const DataType scalar:{DataType(-1),std::exp(DataType(0,.37))}) {
        auto factors=baseFactors;
        factors[0].eta=scalar;
        const auto signedWeight=evaluatePath(path,trial,factors,4);
        const auto magnitude=evaluatePathMagnitude(path,trial,factors,4);
        require(magnitude.magnitude_valid&&!magnitude.phase_valid&&
                !magnitude.green_valid,"synthetic magnitude flags");
        require(std::abs(magnitude.log_abs_matter-signedWeight.log_abs_matter)<2e-12,
                "synthetic scalar changed magnitude");
        require(std::abs(signedWeight.matter_phase-
                         base.matter_phase*pureProjectorUnitPhase(scalar))<2e-12,
                "synthetic scalar phase was erased");
        require(std::abs(signedWeight.matter_phase-DataType(1))>1e-6,
                "synthetic phase did not remain nontrivial");
    }
}
}

int main(int argc,char **argv) {
    try {
        fullRankSmallDeterminantCheck();
        const std::uint64_t seed=argc>1?std::stoull(argv[1]):7911;
        const double g=argc>2?std::stod(argv[2]):.3;
        const int intervals=argc>3?std::stoi(argv[3]):64;
        const int nprop=argc>4?std::stoi(argv[4]):12;
        ReplayMetrics merged,reference;
        runReplay(seed,g,intervals,nprop,false,true,merged);
        runReplay(seed,g,intervals,nprop,false,false,reference);
        syntheticScalarPhaseCheck(g);
        runMeasurementCadence(6121,1);
        runMeasurementCadence(6122,3);
        failureDoesNotCommit(g);
        std::cout<<std::setprecision(17)
            <<"{\"seed\":"<<seed<<",\"g\":"<<g<<",\"intervals\":"<<intervals
            <<",\"merged\":{\"max_log_abs_error\":"<<merged.maxLog
            <<",\"max_phase_error\":"<<merged.maxPhase
            <<",\"max_green_error\":"<<merged.maxGreen
            <<",\"max_log_ratio_error\":"<<merged.maxRatio
            <<",\"max_acceptance_error\":"<<merged.maxAcceptance
            <<",\"accepts\":"<<merged.accepts<<",\"rejects\":"<<merged.rejects
            <<",\"measurements\":"<<merged.measurements<<"},"
            <<"\"reference\":{\"max_log_abs_error\":"<<reference.maxLog
            <<",\"max_phase_error\":"<<reference.maxPhase
            <<",\"max_green_error\":"<<reference.maxGreen
            <<",\"max_log_ratio_error\":"<<reference.maxRatio
            <<",\"max_acceptance_error\":"<<reference.maxAcceptance
            <<",\"accepts\":"<<reference.accepts<<",\"rejects\":"<<reference.rejects
            <<",\"measurements\":"<<reference.measurements<<"},"
            <<"\"physical_negative_or_nonreal_samples\":0"
            <<",\"synthetic_nontrivial_phase_samples\":2"
            <<",\"cadence_measure_every\":[1,3]"
            <<",\"failure_injection\":true}\n";
        require(merged.maxLog<1e-8&&reference.maxLog<1e-8&&
                merged.maxPhase<2e-10&&reference.maxPhase<2e-10&&
                merged.maxGreen<2e-12&&reference.maxGreen<2e-12&&
                merged.maxRatio<2e-9&&reference.maxRatio<2e-9&&
                merged.maxAcceptance<2e-9&&reference.maxAcceptance<2e-9,
                "deferred-phase replay tolerance exceeded");
        std::cout<<"DEFERRED_PHASE_VALIDATE_OK\n";
        return 0;
    } catch(const std::exception &error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
