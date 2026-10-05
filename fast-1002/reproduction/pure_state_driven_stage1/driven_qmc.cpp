#include "stage1_common.h"
#include "pure_projector_failure_snapshot.h"

#include <fstream>
#include <iostream>

using namespace driven_stage1;

namespace {

std::map<std::string, std::string> arguments(int argc, char **argv) {
    std::map<std::string, std::string> result;
    for (int index = 1; index < argc; ++index) {
        const std::string key = argv[index];
        if (key.rfind("--", 0) != 0 || index + 1 >= argc)
            throw std::invalid_argument("arguments must be --name value pairs");
        result[key.substr(2)] = argv[++index];
    }
    return result;
}

std::string get(const std::map<std::string, std::string> &args,
                const std::string &name) {
    const auto found = args.find(name);
    if (found == args.end()) throw std::invalid_argument("missing --" + name);
    return found->second;
}

bool getBoolean(const std::map<std::string,std::string> &args,
                const std::string &name,bool fallback) {
    const auto found=args.find(name);
    if(found==args.end())return fallback;
    require(found->second=="true"||found->second=="false",
            "--"+name+" must be true or false");
    return found->second=="true";
}

PureImaginaryTimeProtocol protocolFromArguments(
        const std::map<std::string, std::string> &args, double dt) {
    const std::string kind = get(args, "protocol");
    const double initial = std::stod(get(args, "Vi"));
    const double tau = std::stod(get(args, "tau"));
    if (kind == "quench")
        return PureImaginaryTimeProtocol::suddenQuench(
            initial, std::stod(get(args, "Vf")), tau, dt);
    if (kind == "ramp")
        return PureImaginaryTimeProtocol::linearRamp(
            initial, std::stod(get(args, "R")), tau, dt);
    if (kind == "constant")
        return PureImaginaryTimeProtocol::constant(initial, tau, dt);
    throw std::invalid_argument("unknown protocol");
}

}  // namespace

int main(int argc, char **argv) {
    try {
        const auto args = arguments(argc, argv);
        ModelParameters p;
        p.L = std::stoi(get(args, "L"));
        p.boundary = get(args, "boundary") == "pbc" ? 0 : 1;
        p.hs = std::stoi(get(args, "hs"));
        p.dt = std::stod(get(args, "dt"));
        const int explicitParity = std::stoi(get(args, "parity"));
        const int burn = std::stoi(get(args, "burn"));
        const int measurements = std::stoi(get(args, "measurements"));
        const int block = args.count("block") ? std::stoi(get(args, "block")) : 16;
        const std::uint64_t seed = std::stoull(get(args, "seed"));
        const std::string outputPath = get(args, "output");
        const std::string protocolKind = get(args, "protocol");
        const std::string updatePolicy=args.count("update-policy")?
            get(args,"update-policy"):"stabilized-local";
        require(updatePolicy=="strict"||updatePolicy=="stabilized-local",
                "--update-policy must be strict or stabilized-local");
        const bool localUpdates=updatePolicy=="stabilized-local";
        p.analytic_hs_factors=getBoolean(args,"analytic-hs",localUpdates);
        const bool fullTraceHashes=getBoolean(args,"trace-hashes",!localUpdates);
        const auto initializationArgument = args.find("initialization-policy");
        const std::string initializationPolicyName =
            initializationArgument != args.end() ? initializationArgument->second :
            (protocolKind == "ramp" ? "independent" : "mirrored-theorem");
        require(p.L >= 2 && (p.boundary != 0 || (p.L % 2) == 0) &&
                    p.dt > 0.0 && burn >= 0 &&
                    measurements > 0 && block > 0 &&
                    (explicitParity == 1 || explicitParity == -1) &&
                    (initializationPolicyName == "independent" ||
                     initializationPolicyName == "mirrored-theorem"),
                "invalid driven QMC arguments");
        const PureImaginaryTimeProtocol protocol =
            protocolFromArguments(args, p.dt);
        const std::string initialState = args.count("initial-state") ?
            get(args,"initial-state") : "free";
        require(initialState == "free" || initialState == "cdw-even" ||
                initialState == "cdw-odd", "unknown --initial-state");
        const bool cdw = initialState != "free";
        const GaussianTrialState initial = cdw ?
            makeCdwInitialState(p, initialState == "cdw-odd", explicitParity) :
            makeInitialState(p, 0.0, 0.0, explicitParity);

        std::mt19937_64 rng(seed);
        const bool mirroredInitialization =
            initializationPolicyName == "mirrored-theorem";
        PureDrivenContourResult contour = mirroredInitialization ?
            makeMirroredDrivenContour(p, protocol, rng) :
            makeDrivenContour(p, protocol, rng);
        const bool braKetHsIndependent = contour.bra_ket_hs_independent;
        const PureFastInitializationPolicy walkerInitializationPolicy =
            mirroredInitialization ?
                PureFastInitializationPolicy::MirroredTheoremZ2Plus :
                PureFastInitializationPolicy::SequentialAudit;
        std::vector<int> indices = proposalIndices(contour.configuration);
        require(!indices.empty(), "driven contour has no HS fields");
        int ketCount = 0, braCount = 0;
        for (int index : indices) {
            if (contour.configuration.locations[index].branch == PureBranch::Ket)
                ++ketCount;
            else ++braCount;
        }
        require(ketCount == braCount && ketCount > 0,
                "bra/ket HS field counts are not balanced");

        PureFastOptions options;
        options.weight_mode = PureProjectorWeightMode::RealZ2;
        options.stabilized_local_updates=localUpdates;
        options.pair_block_factors=getBoolean(args,"pair-block-factors",localUpdates);
        options.direct_green_columns=getBoolean(args,"direct-green-columns",localUpdates);
        // Local stabilized updates reuse their reconstruction workspaces by
        // default; pass --reuse-rebuild-workspace false for an ablation.
        options.reuse_rebuild_workspace=localUpdates&&getBoolean(args,"reuse-rebuild-workspace",true);
        const bool compactDensityProposals=localUpdates&&getBoolean(args,"compact-density-proposals",true);
        PureDrivenProposalFactory proposalFactory(p,protocol,getBoolean(args,"cache-density-factors",localUpdates));
        options.positive_hs1_ratio_contract=positiveHs1RatioContract(p,initial);
        PureProjectorFastWalker walker(
            initial, std::move(contour.configuration), block,
            PureFastRunMode::FastStrict, options,
            walkerInitializationPolicy);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        long long accepted = 0, attempted = 0, sameProposalMismatch = 0;
        std::uint64_t trajectoryHash = 1469598103934665603ULL;
        auto mix = [&](std::uint64_t value) {
            trajectoryHash ^= value;
            trajectoryHash *= 1099511628211ULL;
        };
        int direction = 1;
        int completedSweeps=0; long long flips=0, proposalSignSum=0;
        long long measuredSweepSignSum=0;
        std::ofstream trace, sweepTrace;
        if(args.count("sign-trace")) {
            trace.open(get(args,"sign-trace")); require(bool(trace),"cannot open sign trace");
            trace<<std::setprecision(17)<<"event,proposal,sweep,phase,direction,index,branch,slice,bond,aux,uniform,ratio_real,ratio_imag,accepted,sign_before,sign_after,terminated,ratio_status\n";
            trace<<"initial,0,0,initial,0,-1,none,-1,-1,-1,0,1,0,0,"<<walker.z2Sign()<<','<<walker.z2Sign()<<",0,success\n";
        }
        if(args.count("sweep-sign-trace")) {
            sweepTrace.open(get(args,"sweep-sign-trace")); require(bool(sweepTrace),"cannot open sweep trace");
            sweepTrace<<std::setprecision(17)<<"sweep,phase,physical_sign,sign_flips,proposal_sign_mean,cumulative_measurement_sign\n";
        }
        auto sweep = [&]() {
            indices = proposalIndices(walker.configuration());
            flips=0; proposalSignSum=0;
            auto step = [&](int index) {
                const int signBefore=walker.z2Sign();
                const auto traceLocation=walker.configuration().locations[index];
                const double u = uniform(rng);
                const PureFastProposal proposal=proposalFactory.proposal(walker.configuration(),index,u,!compactDensityProposals);
                const PureFastProposalResult result = walker.propose(proposal);
                ++attempted;
                accepted += result.accepted;
                flips+=(signBefore!=walker.z2Sign()); proposalSignSum+=walker.z2Sign();
                if(trace.is_open()) {
                    trace<<"proposal,"<<attempted<<','<<completedSweeps+1<<','<<(completedSweeps<burn?"burn":"measurement")
                        <<','<<direction<<','<<index<<','<<(traceLocation.branch==PureBranch::Ket?"ket":"bra")
                        <<','<<traceLocation.slice<<','<<traceLocation.bond<<','<<traceLocation.aux<<','<<u
                        <<','<<result.ratio.ratio.real()<<','<<result.ratio.ratio.imag()<<','<<result.accepted
                        <<','<<signBefore<<','<<walker.z2Sign()<<','<<result.terminated<<','<<pureFastRatioStatusName(result.ratio.status)<<'\n';
                    require(bool(trace),"sign trace write failed");
                }
                if (result.ratio.used_reference)
                    sameProposalMismatch += result.z2_reference_mismatch;
                if (result.terminated || !result.ratio.ok()) {
                    const std::string failurePath=outputPath+".failure."+std::to_string(
                        std::chrono::system_clock::now().time_since_epoch().count())+".bin";
                    try {pure_failure_snapshot::write(failurePath,initial,walker,proposal.materialized(),result,block,rng);
                        std::cerr<<"failure_snapshot="<<failurePath<<'\n';
                    }catch(const std::exception &e){std::cerr<<"failure_snapshot_error="<<e.what()<<'\n';}
                    const PureFastDiagnostics &diagnostics = walker.diagnostics();
                    std::ostringstream message;
                    message << std::setprecision(17)
                            << "driven proposal reference failed closed"
                            << " attempted=" << attempted
                            << " index=" << index
                            << " stage="
                            << (result.failure_stage.empty() ?
                                "fast_ratio_untrusted" : result.failure_stage)
                            << " ratio_status="
                            << pureFastRatioStatusName(result.ratio.status)
                            << " branch="
                            << (result.snapshot.branch == PureBranch::Ket ?
                                "ket" : "bra")
                            << " slice=" << result.snapshot.slice
                            << " bond=" << result.snapshot.bond
                            << " aux=" << result.snapshot.aux
                            << " old_hs=" << result.snapshot.old_hs
                            << " new_hs=" << result.snapshot.new_hs
                            << " uniform=" << result.snapshot.uniform
                            << " configuration_hash="
                            << result.snapshot.configuration_hash
                            << " fast_ratio=" << result.snapshot.fast_ratio
                            << " overlap_rcond=" << result.ratio.overlap_rcond
                            << " solve_residual=" << result.ratio.solve_residual
                            << " trust_alarm=" << result.ratio.trust_alarm
                            << " used_reference=" << result.ratio.used_reference
                            << " mp_status="
                            << pureMpProposalStatusName(
                                   result.ratio.mp_reference.status)
                            << " mp_precision="
                            << result.ratio.mp_reference.precision_digits
                            << " mp_pre_rcond="
                            << result.ratio.mp_reference.pre_endpoint_rcond
                            << " mp_post_rcond="
                            << result.ratio.mp_reference.post_endpoint_rcond
                            << " mp_pre_residual="
                            << result.ratio.mp_reference.pre_endpoint_residual
                            << " mp_post_residual="
                            << result.ratio.mp_reference.post_endpoint_residual
                            << " mp_message=\""
                            << result.ratio.mp_reference.message << "\""
                            << " proposals=" << diagnostics.proposals
                            << " rebuilds=" << diagnostics.rebuild_count
                            << " trust_alarms="
                            << diagnostics.trust_alarm_count
                            << " mp_fallbacks="
                            << diagnostics.mp_fallback_count
                            << " mp_failures="
                            << diagnostics.mp_fallback_failure_count
                            << " mp_endpoint_fallbacks="
                            << diagnostics.mp_endpoint_fallback_count
                            << " mp_endpoint_failures="
                            << diagnostics.mp_endpoint_fallback_failure_count
                            << " mp_endpoint_rejection_recoveries="
                            << diagnostics.mp_endpoint_rejection_recovery_count
                            << " accepted_factor_center_recoveries="
                            << diagnostics.accepted_factor_center_recovery_count
                            << " post_decision_center_recoveries="
                            << diagnostics.post_decision_center_recovery_count
                            << " first_failure_proposal="
                            << diagnostics.first_failure_proposal
                            << " stack_detail=\"" << walker.stackFailureDetail() << "\""
                            << " local_scaled_green_rebuilds=" << walker.scaledGreenRebuilds();
                    throw std::runtime_error(message.str());
                }
                mix(index);
                mix(result.accepted);
                mix(std::uint64_t(result.ratio.status));
                if(fullTraceHashes)mix(walker.configurationHash());
                else mix(std::uint64_t(std::int64_t(result.snapshot.new_hs)));
                mix(std::uint64_t(std::int64_t(walker.z2Sign())));
            };
            if (direction > 0)
                for (int index : indices) step(index);
            else
                for (auto iterator = indices.rbegin(); iterator != indices.rend();
                     ++iterator) step(*iterator);
            if(completedSweeps>=burn) measuredSweepSignSum+=walker.z2Sign();
            if(sweepTrace.is_open()) {
                sweepTrace<<completedSweeps+1<<','<<(completedSweeps<burn?"burn":"measurement")<<','<<walker.z2Sign()
                    <<','<<flips<<','<<double(proposalSignSum)/indices.size()<<',';
                if(completedSweeps>=burn) sweepTrace<<double(measuredSweepSignSum)/(completedSweeps-burn+1);
                sweepTrace<<'\n'; sweepTrace.flush(); require(bool(sweepTrace),"sweep trace write failed");
            }
            if(trace.is_open()) { trace.flush(); require(bool(trace),"sign trace flush failed"); }
            ++completedSweeps;
            direction = -direction;
        };
        for (int index = 0; index < burn; ++index) sweep();

        std::ofstream csv(outputPath);
        require(bool(csv), "cannot open driven QMC measurement CSV");
        csv << "seed,measurement,physical_sign,energy,S_pi,S_pi_dq,R_CDW,"
               "fermion_parity,G_0_1_real,G_0_1_imag,G_0_L_real,G_0_L_imag,"
               "configuration_hash,acceptance";
        if (cdw) csv << ",m_pi_real,m_pi_imag";
        csv << '\n';
        double signDenominator = 0.0;
        int positive = 0, negative = 0;
        double endpointResidual = 0.0, phaseImaginary = 0.0;
        for (int measurement = 0; measurement < measurements; ++measurement) {
            sweep();
            const PureProjectorGreenResult green = walker.measurementGreen();
            require(green.ok(), "driven measurement Green failed");
            const int sign = walker.z2Sign();
            const Observables observable = observe(
                p, protocol.finalValue(), green.green,localUpdates?explicitParity:0);
            require(std::isfinite(observable.energy)&&std::isfinite(observable.S_pi)&&
                    std::isfinite(observable.S_pi_dq)&&std::isfinite(observable.R_CDW),
                    "nonfinite driven observable");
            signDenominator += sign;
            positive += sign > 0;
            negative += sign < 0;
            endpointResidual = std::max(endpointResidual, green.green_residual);
            phaseImaginary = std::max(
                phaseImaginary,
                std::abs(walker.currentWeight().complex_phase.imag()));
            csv << seed << ',' << measurement << ',' << sign << ','
                << std::setprecision(17) << observable.energy << ','
                << observable.S_pi << ',' << observable.S_pi_dq << ','
                << observable.R_CDW << ',' << observable.parity << ','
                << green.green(0, 1).real() << ',' << green.green(0, 1).imag()
                << ',' << green.green(0, p.L).real() << ','
                << green.green(0, p.L).imag() << ','
                << walker.configurationHash() << ','
                << double(accepted) / attempted;
            if (cdw) {
                const auto m = staggeredDensity(green.green);
                require(std::isfinite(m.real()) && std::isfinite(m.imag()),
                        "nonfinite staggered density");
                csv << ',' << m.real() << ',' << m.imag();
            }
            csv << '\n';
        }
        csv.flush();
        require(bool(csv), "driven QMC measurement CSV write failed");
        const PureFastDiagnostics &diagnostics = walker.diagnostics();
        std::cout << std::setprecision(17)
                  << "{\"status\":\"complete\",\"L\":" << p.L
                  << ",\"green_rebuild_check_policy\":\""
                  << (localUpdates?"absolute_then_relative_to_green":"absolute_strict") << "\""
                  << ",\"green_rebuild_tolerance\":"
                  << (localUpdates?options.local_green_rebuild_tolerance:options.green_update_tolerance)
                  << ",\"local_scaled_green_rebuilds\":" << walker.scaledGreenRebuilds()
                  << ",\"boundary\":\"" << (p.boundary ? "obc" : "pbc")
                  << "\",\"protocol\":\"" << protocol.name() << "\""
                  << ",\"V_i\":" << protocol.initialValue()
                  << ",\"V_f\":" << protocol.finalValue()
                  << ",\"R\":" << protocol.rate()
                  << ",\"tau_f\":" << protocol.tauFinal()
                  << ",\"Delta_tau\":" << p.dt
                  << ",\"update_policy\":\"" << updatePolicy << "\""
                  << ",\"analytic_hs_factors\":" << (p.analytic_hs_factors?"true":"false")
                  << ",\"pair_block_factors\":" << (options.pair_block_factors?"true":"false")
                  << ",\"reuse_rebuild_workspace\":" << (options.reuse_rebuild_workspace?"true":"false")
                  << ",\"direct_green_columns\":" << (options.direct_green_columns?"true":"false")
                  << ",\"compact_density_proposals\":" << (compactDensityProposals&&proposalFactory.enabled()?"true":"false")
                  << ",\"cached_density_factors\":" << (proposalFactory.enabled()?"true":"false")
                  << ",\"fermion_parity_source\":\""
                  << (localUpdates?"conserved_initial_sector":"green_pfaffian") << "\""
                  << ",\"stabilization_block\":" << block
                  << ",\"stabilization_unit\":\"factor\""
                  << ",\"trajectory_hash_schema\":\""
                  << (fullTraceHashes?"configuration_v1":"proposal_events_v2") << "\""
                  << ",\"local_double_retries\":" << diagnostics.local_double_retries
                  << ",\"local_ratio_imaginary_relative_max\":"
                  << diagnostics.maximum_local_ratio_imaginary
                  << ",\"Phi_0_parity\":" << explicitParity
                  << ",\"seed\":" << seed
                  << ",\"burn_sweeps\":" << burn
                  << ",\"measurements\":" << measurements
                  << ",\"initialization_policy\":\""
                  << walker.initializationPolicy() << "\""
                  << ",\"bra_ket_hs_independent\":"
                  << (braKetHsIndependent ? "true" : "false")
                  << ",\"contour_order\":\"" << contour.contour_order << "\""
                  << ",\"ket_hs_fields\":" << ketCount
                  << ",\"bra_hs_fields\":" << braCount
                  << ",\"positive_weights\":" << positive
                  << ",\"negative_weights\":" << negative
                  << ",\"sign_denominator\":" << signDenominator
                  << ",\"average_sign\":" << signDenominator / measurements
                  << ",\"acceptance\":" << double(accepted) / attempted
                  << ",\"trajectory_hash\":" << trajectoryHash
                  << ",\"final_configuration_hash\":"
                  << walker.configurationHash()
                  << ",\"rng_hash\":" << hashRng(rng)
                  << ",\"same_proposal_z2_mismatch\":"
                  << sameProposalMismatch
                  << ",\"mp_fallback_count\":"
                  << diagnostics.mp_fallback_count
                  << ",\"mp_fallback_failure_count\":"
                  << diagnostics.mp_fallback_failure_count
                  << ",\"mp_endpoint_fallback_count\":"
                  << diagnostics.mp_endpoint_fallback_count
                  << ",\"mp_endpoint_fallback_failure_count\":"
                  << diagnostics.mp_endpoint_fallback_failure_count
                  << ",\"mp_endpoint_rejection_recovery_count\":"
                  << diagnostics.mp_endpoint_rejection_recovery_count
                  << ",\"mp_subspace_legacy_fallback_count\":"
                  << diagnostics.mp_profile.subspace_legacy_fallback_count
                  << ",\"endpoint_green_residual_max\":"
                  << endpointResidual
                  << ",\"complex_phase_imaginary_max\":"
                  << phaseImaginary;
        if (cdw) std::cout << ",\"initial_state\":\"" << initialState << "\"";
        std::cout << "}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "driven_qmc: " << error.what() << '\n';
        return 1;
    }
}
