#include "phonon_mc.h"
#include <chrono>
#include <iomanip>
#include <iostream>
using namespace kitaev_ssh;

namespace {
double worst=0.;
void require(bool condition,const char *message) {
    if(!condition)throw std::runtime_error(message);
}
void compare(double error,double scale,const char *message) {
    const double scaled=error/std::max(1.,scale);
    worst=std::max(worst,scaled);
    require(std::isfinite(scaled)&&scaled<2e-9,message);
}
void kernelBoundaries() {
    const auto fast=phononCenterGreenOptions(true),ref=phononCenterGreenOptions(false);
    require(fast.rank_tolerance==ref.rank_tolerance&&
            fast.minimum_overlap_rcond==ref.minimum_overlap_rcond&&
            fast.solve_residual_tolerance==ref.solve_residual_tolerance&&
            fast.green_residual_tolerance==ref.green_residual_tolerance&&
            fast.well_conditioned_solve_tolerance==0.&&!fast.allow_oblique_green,
            "center optimization changed a strict numerical gate");
    for(int n:{4,32,36,81}) {
        MatType matrix=MatType::Zero(2*n,n);
        for(int j=0;j<n;++j)for(int i=0;i<2*n;++i)
            matrix(i,j)=DataType(.08*std::sin(.2*i+.3*j),.06*std::cos(.3*i-.17*j));
        matrix.topRows(n)+=MatType::Identity(n,n);
        const auto a=pure_projector_endpoint_detail::scaledThinQr(matrix,ref);
        const auto b=pure_projector_endpoint_detail::scaledThinQr(matrix,fast);
        require(a.ok()&&b.ok(),"generic complex endpoint QR failed");
        compare((a.q*a.q.adjoint()-b.q*b.q.adjoint()).norm(),a.q.norm(),"QR subspace mismatch");
        compare(std::abs(a.log_abs_determinant-b.log_abs_determinant),std::abs(a.log_abs_determinant),"QR determinant scale mismatch");
        matrix=MatType::Zero(2*n,n);
        matrix.topRows(n)=MatType::Identity(n,n);
        matrix(n-1,n-1)=.5e-12;
        require(!pure_projector_endpoint_detail::scaledThinQr(matrix,ref).ok()&&
                !pure_projector_endpoint_detail::scaledThinQr(matrix,fast).ok(),
                "rank boundary was relaxed");
        matrix(n-1,n-1)=2e-12;
        require(pure_projector_endpoint_detail::scaledThinQr(matrix,ref).ok()&&
                pure_projector_endpoint_detail::scaledThinQr(matrix,fast).ok(),
                "trusted diagonal QR rejected");
    }
    MatType left=MatType::Zero(64,32),right=left;
    left.topRows(32)=right.topRows(32)=MatType::Identity(32,32);
    left(31,31)=.5e-12;
    require(!pureProjectorGreen(right,left,ref).ok()&&!pureProjectorGreen(right,left,fast).ok(),
            "Green overlap boundary was relaxed");
}
double seconds(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}
}

int main(int argc,char **argv) {try {
    kernelBoundaries();
    const int L=argc>1?std::stoi(argv[1]):6,M=argc>2?std::stoi(argv[2]):12;
    const double omega=argc>3?std::stod(argv[3]):.3,lambda=argc>4?std::stod(argv[4]):2.4;
    const KitaevSSHModel model(HoneycombLattice(L,L,0),p2Couplings(std::sqrt(lambda)));
    PhononPath path(model.lattice().bondCount(),M,.1,{1/(omega*omega),omega});
    std::vector<std::vector<double>> x(path.nodes(),std::vector<double>(path.bonds()));
    for(int j=0;j<path.nodes();++j)for(int b=0;b<path.bonds();++b)
        x[j][b]=.17*std::sin(.7*b+.13*j)+.08*std::cos(.3*b-.21*j);
    path=path.withCoordinates(x);
    const auto generator=honeycombGenerator(model);const auto trial=freeGroundState(model);
    const auto factors=pathFactors(path,trial.Phi.rows(),generator);
    double referenceSeconds=0.,fastSeconds=0.;
    for(int repeat=0;repeat<4;++repeat) {
        PureEndpointRebuildResult a,b;
        auto reference=[&] {const auto start=std::chrono::steady_clock::now();
            a=pureProjectorEndpointRebuild(trial,factors,M,4,phononCenterGreenOptions(false));
            referenceSeconds+=seconds(start);};
        auto accelerated=[&] {const auto start=std::chrono::steady_clock::now();
            b=pureProjectorEndpointRebuild(trial,factors,M,4,phononCenterGreenOptions(true));
            fastSeconds+=seconds(start);};
        if(repeat%2){accelerated();reference();}else{reference();accelerated();}
        require(a.ok()&&b.ok(),"production-size endpoint comparison failed");
        compare((a.green-b.green).norm(),a.green.norm(),"complete endpoint Green mismatch");
        compare(std::abs(a.log_abs_weight-b.log_abs_weight),std::abs(a.log_abs_weight),"endpoint weight mismatch");
    }
    for(bool deferred:{false,true}) {
        PhononWeightOptions options{true,true,deferred,true,true,true,true};
        CoupledPhononPath state(path,trial,generator,4,options);
        const auto current=state.ensureMeasured();
        const auto reference=state.uncachedReference(path,true,true);
        compare(std::abs(current.weight.log_abs_matter-reference.weight.log_abs_matter),std::abs(reference.weight.log_abs_matter),"signed weight mismatch");
        compare(std::abs(current.weight.matter_phase-reference.weight.matter_phase),1.,"phase mismatch");
        compare((current.weight.center_green-reference.weight.center_green).norm(),reference.weight.center_green.norm(),"measurement Green mismatch");
        const auto a=basicObservables(model,current,true),b=basicObservables(model,reference,true);
        for(std::size_t k=0;k<a.size();++k)compare(std::abs(a[k]-b[k]),std::abs(b[k]),"matter energy or basic observable mismatch");
        auto accepted=state.propose(M/2,5,path.coordinate(M/2,5)+.03);
        state.accept(accepted);
        const auto updated=state.ensureMeasured();
        const auto full=state.uncachedReference(updated.path,true,true);
        compare((updated.weight.center_green-full.weight.center_green).norm(),full.weight.center_green.norm(),"accepted measurement mismatch");
        const auto before=state.state();
        auto rejected=state.propose(M/2,6,path.coordinate(M/2,6)-.02);state.reject(rejected);
        require(state.state().configuration_hash==before.configuration_hash&&
                (state.state().weight.center_green-before.weight.center_green).norm()==0.,"rejection changed accepted Green");
    }
    // A rank-two change in K does not produce a rank-two change in exp(-dt K/2).
    const int node=M/2,bond=5;
    const auto changed=path.withCoordinate(node,bond,path.coordinate(node,bond)+.15);
    const MatType generatorDelta=generator(changed.node(node))-generator(path.node(node));
    const MatType oldHalf=nodeHalfStep(path,node,trial.Phi.rows(),generator);
    const MatType newHalf=nodeHalfStep(changed,node,trial.Phi.rows(),generator);
    Eigen::JacobiSVD<MatType> ksvd(generatorDelta),bsvd(newHalf-oldHalf);
    ksvd.setThreshold(1e-12);bsvd.setThreshold(1e-12);
    const int generatorRank=ksvd.rank(),halfDifferenceRank=bsvd.rank();
    require(generatorRank==2&&halfDifferenceRank>2,"low-rank feasibility probe did not distinguish generator and propagator");
    std::cout<<std::setprecision(17)<<"PASS L="<<L<<" M="<<M<<" omega="<<omega<<" lambda="<<lambda
        <<" worst_scaled_error="<<worst<<" reference_endpoint_seconds="<<referenceSeconds
        <<" fast_endpoint_seconds="<<fastSeconds<<" fast_qr_calls="<<pure_projector_endpoint_detail::endpointAcceleratedQrCalls
        <<" reference_qr_calls="<<pure_projector_endpoint_detail::endpointReferenceQrCalls
        <<" generator_change_rank="<<generatorRank<<" half_step_difference_rank="<<halfDifferenceRank<<'\n';
    return 0;
}catch(const std::exception &error){std::cerr<<"center_green_validate FAIL: "<<error.what()<<'\n';return 1;}}
