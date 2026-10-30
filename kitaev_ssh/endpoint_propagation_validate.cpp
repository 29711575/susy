// Deterministic endpoint propagation checks. No sampler, MC, or proposal RNG.
#include "phonon_projector.h"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace kitaev_ssh;
namespace {
int checks=0,endpointCases=0;
double maximumGreenError=0.,maximumLogError=0.;
void require(bool ok,const std::string &what) {
    ++checks;if(!ok)throw std::runtime_error(what);
}
bool same(const MatType &a,const MatType &b) {
    return a.rows()==b.rows()&&a.cols()==b.cols()&&
        (a.size()==0||std::memcmp(a.data(),b.data(),sizeof(DataType)*a.size())==0);
}
std::uint64_t hash(const PureFactorMatrix &factor) {
    std::uint64_t result=1469598103934665603ULL;factor.mixHash(result);return result;
}
MatType input(int rows,int cols) {
    MatType result(rows,cols);
    for(int col=0;col<cols;++col)for(int row=0;row<rows;++row)
        result(row,col)=DataType(.04*std::sin(.37*row+.21*col),.03*std::cos(.19*row-.33*col));
    return result;
}
void kernelChecks() {
    PureEndpointRebuildWorkspace work;
    for(int dimension:{6,64,162,8,64})for(int kind:{0,1,2})for(bool adjoint:{false,true}) {
        MatType b=input(dimension,dimension);
        b.diagonal().array()+=DataType(1.,0.);
        // Preserve nontrivial complex entries and all four signed-zero pairs.
        b(0,0)=DataType(-0.,-0.);b(1,0)=DataType(0.,-0.);
        b(0,1)=DataType(-0.,0.);b(1,1)=DataType(0.,0.);
        if(kind==2) {
            b=MatType::Identity(dimension,dimension);
            b(0,1)=DataType(0.,-.07);b(1,0)=DataType(0.,.07);
            b(2,0)=DataType(-0.,-0.);b(0,2)=DataType(0.,-0.);
        }
        const PureFactorMatrix factor=kind==1?PureFactorMatrix::fromValidatedDense(b):PureFactorMatrix(b);
        require(bool(factor.denseStorage())==(kind!=2),"kernel storage did not exercise requested dispatch");
        const MatType savedFactor=factor;const auto savedHash=hash(factor);
        MatType frame=input(dimension,dimension/2),reference=frame;
        const MatType savedInput=frame;
        if(adjoint)reference=factor.adjoint()*reference;
        else reference=factor*reference;
        const auto oldRight=pure_projector_endpoint_detail::endpointDirectDenseRightCalls;
        const auto oldLeft=pure_projector_endpoint_detail::endpointDirectDenseLeftCalls;
        work.propagation=MatType::Constant(dimension,dimension/2,
            DataType(std::numeric_limits<double>::quiet_NaN(),0.));
        const auto *oldFrame=frame.data();
        pure_projector_endpoint_detail::propagateFactor(factor,frame,adjoint,&work);
        require(same(frame,reference),"kernel complex coefficients or signed zeros changed");
        require(same(MatType(factor),savedFactor)&&hash(factor)==savedHash,"kernel mutated immutable factor");
        require(frame.data()!=work.propagation.data(),"kernel output aliases propagation scratch");
        if(kind!=2) {
            require(work.propagation.data()==oldFrame,"kernel did not recycle the previous frame allocation");
            require(frame.data()!=factor.denseStorage()->data(),"kernel output aliases immutable factor");
            require(same(work.propagation,savedInput),"kernel modified its input during the dense product");
        }
        require(pure_projector_endpoint_detail::endpointDirectDenseRightCalls-oldRight==
                    std::uint64_t(kind!=2&&!adjoint),"right dense dispatch count differs");
        require(pure_projector_endpoint_detail::endpointDirectDenseLeftCalls-oldLeft==
                    std::uint64_t(kind!=2&&adjoint),"left dense dispatch count differs");
        MatType nullFrame=savedInput;
        pure_projector_endpoint_detail::propagateFactor(factor,nullFrame,adjoint,nullptr);
        require(same(nullFrame,reference),"null workspace reference changed");
        const MatType returned=frame;
        MatType next=input(dimension,dimension/2);
        pure_projector_endpoint_detail::propagateFactor(factor,next,!adjoint,&work);
        require(same(frame,returned),"later propagation mutated a returned frame");
        PureEndpointRebuildWorkspace copied=work,assigned;assigned=work;
        require(copied.propagation.data()!=work.propagation.data()&&
                assigned.propagation.data()!=work.propagation.data(),"workspace copies share writable propagation storage");
        const MatType savedScratch=work.propagation;
        next=input(dimension,dimension/2);
        pure_projector_endpoint_detail::propagateFactor(factor,next,adjoint,&copied);
        require(same(work.propagation,savedScratch),"copied workspace mutated original scratch");
    }
}
void compare(const PureEndpointRebuildResult &reference,const PureEndpointRebuildResult &candidate,
             const std::string &where,bool exact) {
    ++endpointCases;
    require(reference.status==candidate.status&&reference.cut==candidate.cut,where+": status or cut differs");
    if(!reference.ok())return;
    require(reference.overlap_rank==candidate.overlap_rank,where+": overlap rank differs");
    const double logError=std::abs(reference.log_abs_weight-candidate.log_abs_weight);
    const double greenError=(reference.green-candidate.green).cwiseAbs().maxCoeff();
    maximumGreenError=std::max(maximumGreenError,greenError);
    maximumLogError=std::max(maximumLogError,logError);
    require(std::isfinite(logError)&&logError<2e-10,where+": log magnitude differs");
    require(std::isfinite(greenError)&&greenError<2e-10,where+": full complex Green differs");
    require(reference.scalar_prefactor_phase==candidate.scalar_prefactor_phase&&
            reference.log_abs_scalar_prefactor==candidate.log_abs_scalar_prefactor,
            where+": complex scalar prefactor differs");
    if(exact) {
        require(same(reference.q_right,candidate.q_right)&&same(reference.q_left,candidate.q_left)&&
                same(reference.green,candidate.green),where+": same-backend coefficients or signed zeros differ");
        require(reference.log_abs_weight==candidate.log_abs_weight&&
                reference.overlap_rcond==candidate.overlap_rcond&&
                reference.solve_residual==candidate.solve_residual&&
                reference.green_residual==candidate.green_residual&&
                reference.green_skew_residual==candidate.green_skew_residual&&
                reference.green_diagonal_residual==candidate.green_diagonal_residual&&
                reference.long_double_green_residual==candidate.long_double_green_residual&&
                reference.canonical_frame_maintenance==candidate.canonical_frame_maintenance,
                where+": numerical certificate changed");
    }
}
void endpointChecks() {
    PureEndpointRebuildWorkspace work;
    for(int L:{2,6,9,2}) {
        const int intervals=6,block=3;
        const KitaevSSHModel model(HoneycombLattice(L,L,0),Couplings{1.,1.,.5,.4});
        PhononPath path(model.lattice().bondCount(),intervals,.04,{7.,.6});
        std::vector<std::vector<double>> coordinates(path.nodes(),std::vector<double>(path.bonds()));
        for(int j=0;j<path.nodes();++j)for(int b=0;b<path.bonds();++b)
            coordinates[j][b]=.07*std::sin(.21*j+.13*b)+.01*std::cos(.37*j-.09*b);
        path=path.withCoordinates(coordinates);
        const auto trial=freeGroundState(model);
        const auto factors=pathFactors(path,trial.Phi.rows(),honeycombGenerator(model));
        const auto fast=phononCenterGreenOptions(true),strict=phononCenterGreenOptions(false);
        for(int cut:{0,1,intervals,2*intervals-1,2*intervals}) {
            const auto original=pureProjectorEndpointRebuild(trial,factors,cut,block,fast);
            const auto candidate=pureProjectorEndpointRebuild(trial,factors,cut,block,fast,&work);
            require(original.ok()&&candidate.ok(),"valid lattice endpoint failed");
            compare(original,candidate,"L"+std::to_string(L)+" accelerated null/workspace",true);
            const auto independent=pureProjectorEndpointRebuild(trial,factors,cut,block,strict);
            require(independent.ok(),"independent Eigen endpoint failed");
            compare(independent,candidate,"L"+std::to_string(L)+" Eigen/workspace",false);
            const MatType retained=candidate.green;
            auto copied=work;
            const auto repeated=pureProjectorEndpointRebuild(trial,factors,cut,block,fast,&copied);
            compare(candidate,repeated,"copied endpoint workspace",true);
            require(same(retained,candidate.green),"subsequent endpoint changed returned Green");
        }
        auto mixed=factors;
        mixed[0]=PureProjectorSlice(MatType::Identity(trial.Phi.rows(),trial.Phi.rows()),
            std::polar(1.03,.17),"sparse identity with complex scalar");
        const auto mixedNull=pureProjectorEndpointRebuild(trial,mixed,intervals,block,fast);
        const auto mixedWork=pureProjectorEndpointRebuild(trial,mixed,intervals,block,fast,&work);
        require(mixedNull.ok()&&mixedWork.ok(),"sparse fallback endpoint failed");
        compare(mixedNull,mixedWork,"sparse fallback/complex eta",true);
        auto bad=factors;
        bad[0].eta=DataType(0,0);
        const auto zero=pureProjectorEndpointRebuild(trial,bad,intervals,block,fast,&work);
        require(zero.status==PureEndpointRebuildStatus::zero_scalar_prefactor,"zero eta gate changed");
        bad=factors;
        bad[0].matrix=MatType::Zero(trial.Phi.rows(),trial.Phi.rows());
        const auto rank=pureProjectorEndpointRebuild(trial,bad,intervals,block,fast,&work);
        const auto rankNull=pureProjectorEndpointRebuild(trial,bad,intervals,block,fast);
        require(!rank.ok(),"zero factor rank gate changed");
        compare(rankNull,rank,"rank failure",true);
        bad=factors;
        MatType nonfinite=bad[0].matrix;
        nonfinite(0,0)=DataType(std::numeric_limits<double>::infinity(),0.);
        bad[0].matrix=nonfinite;
        const auto invalid=pureProjectorEndpointRebuild(trial,bad,intervals,block,fast,&work);
        require(invalid.status==PureEndpointRebuildStatus::invalid_configuration,"nonfinite factor gate changed");
        const auto recovered=pureProjectorEndpointRebuild(trial,factors,intervals,block,fast,&work);
        const auto recoveredNull=pureProjectorEndpointRebuild(trial,factors,intervals,block,fast);
        require(recovered.ok()&&recoveredNull.ok(),"failed endpoint poisoned workspace");
        compare(recoveredNull,recovered,"failure recovery",true);
        const auto beforeRight=pure_projector_endpoint_detail::endpointDirectDenseRightCalls;
        const auto beforeLeft=pure_projector_endpoint_detail::endpointDirectDenseLeftCalls;
        const auto disabled=pureProjectorEndpointRebuild(trial,factors,intervals,block,strict,&work);
        require(disabled.ok(),"nonaccelerated workspace endpoint failed");
        require(beforeRight==pure_projector_endpoint_detail::endpointDirectDenseRightCalls&&
                beforeLeft==pure_projector_endpoint_detail::endpointDirectDenseLeftCalls,
                "nonaccelerated reference used new dense propagation");
    }
}
}
int main() {try {
    mkl_set_num_threads(1);
    kernelChecks();endpointChecks();
    require(pure_projector_endpoint_detail::endpointDirectDenseRightCalls>0&&
            pure_projector_endpoint_detail::endpointDirectDenseLeftCalls>0,"dense endpoint dispatch was not exercised");
    std::cout<<std::setprecision(17)<<"{\"status\":\"PASS\",\"checks\":"<<checks
        <<",\"endpoint_comparisons\":"<<endpointCases<<",\"maximum_green_error\":"<<maximumGreenError
        <<",\"maximum_log_error\":"<<maximumLogError
        <<",\"direct_dense_right_calls\":"<<pure_projector_endpoint_detail::endpointDirectDenseRightCalls
        <<",\"direct_dense_left_calls\":"<<pure_projector_endpoint_detail::endpointDirectDenseLeftCalls<<"}\n";
    return 0;
}catch(const std::exception &error){std::cerr<<"endpoint_propagation_validate: "<<error.what()<<'\n';return 1;}}
