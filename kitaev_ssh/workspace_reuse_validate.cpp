#include "../inc/gaussian_trial_state.h"
#include "../inc/pure_projector_green.h"
#include "../inc/pure_projector_real.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <limits>

namespace {
int checks=0;
double worst=0.;
void require(bool value,const char *message) {
    ++checks;if(!value)throw std::runtime_error(message);
}
void close(double error,double scale,const char *message,double tolerance=2e-12) {
    const double relative=error/std::max(1.,scale);
    worst=std::max(worst,relative);
    require(std::isfinite(relative)&&relative<tolerance,message);
}
MatType gaussian(int n) {
    MatType q=MatType::Zero(2*n,n);
    q.topRows(n)=MatType::Identity(n,n)/std::sqrt(2.);
    q.bottomRows(n)=DataType(0.,1.)*MatType::Identity(n,n)/std::sqrt(2.);
    return q;
}
MatType input(int n) {
    MatType a=MatType::Zero(2*n,n);
    for(int j=0;j<n;++j)for(int i=0;i<2*n;++i)
        a(i,j)=DataType(.015*std::sin(.37*i+.21*j),.012*std::cos(.19*i-.33*j));
    a.topRows(n)+=MatType::Identity(n,n);return a;
}
void qrChecks() {
    PureComplexQRWorkspace work;
    MatType empty;
    require(!work.formQ(empty),"QR formQ before factorize succeeded");
    require(!work.factorize(empty),"QR empty input accepted");
    for(int n:{32,81,36,32}) {
        const MatType a=input(n);MatType ref=a,reused=a;cVecType tau(n);
        require(LAPACKE_zgeqrf(LAPACK_COL_MAJOR,2*n,n,ref.data(),2*n,tau.data())==0,"reference QR failed");
        require(work.factorize(reused),"workspace QR failed");
        DataType detRef=1.,detReused=1.;
        for(int j=0;j<n;++j){detRef*=ref(j,j);detReused*=reused(j,j);}
        close((ref-reused).norm(),ref.norm(),"packed QR mismatch");
        close(std::abs(detRef-detReused),std::abs(detRef),"QR determinant mismatch");
        Eigen::ColPivHouseholderQR<MatType> independent(a);independent.setThreshold(1e-12);
        require(independent.rank()==n,"independent QR rank mismatch");
        require(LAPACKE_zungqr(LAPACK_COL_MAJOR,2*n,n,n,ref.data(),2*n,tau.data())==0,"reference formQ failed");
        require(work.formQ(reused),"workspace formQ failed");
        close((ref-reused).norm(),ref.norm(),"full thin Q mismatch");
        close((reused.adjoint()*reused-MatType::Identity(n,n)).norm(),std::sqrt(double(n)),"QR orthogonality mismatch");
        const MatType saved=reused;
        require(!work.formQ(reused),"QR formQ consumed twice");
        require((reused-saved).norm()==0.,"failed formQ changed result");
    }
    MatType bad=MatType::Zero(2,3);
    require(!work.factorize(bad),"wide QR accepted");
    bad=input(32);bad(0,0)=DataType(std::numeric_limits<double>::quiet_NaN(),0.);
    require(!work.factorize(bad)&&!work.formQ(bad),"nonfinite QR accepted or left stale state");
    bad=input(32);require(work.factorize(bad)&&work.formQ(bad),"QR failed to recover");
    bad=input(32);require(work.factorize(bad),"QR shape check setup failed");
    MatType wrongShape=MatType::Identity(8,4);
    require(!work.formQ(wrongShape),"QR formQ accepted wrong shape");
    require(!work.formQ(bad),"QR shape mismatch failed to consume pending factorization");
    bad=input(32);require(work.factorize(bad),"QR buffer check setup failed");
    MatType wrongBuffer=bad;
    const MatType wrongBufferSaved=wrongBuffer;
    require(!work.formQ(wrongBuffer),"QR formQ accepted unrelated same-size storage");
    require((wrongBuffer-wrongBufferSaved).norm()==0.,"QR wrong-buffer failure changed input");
    require(!work.formQ(bad),"QR wrong-buffer attempt failed to consume pending factorization");
    bad=input(32);require(work.factorize(bad)&&work.formQ(bad),"QR binding failure prevented fresh recovery");
    // A zero column is not rejected by non-pivoted LAPACK QR. Independently
    // verify its rank rather than inventing a new workspace rank gate.
    bad=input(32);bad.col(31).setZero();
    require(work.factorize(bad),"rank deficient QR status changed");
    require(std::abs(bad(31,31))<1e-14,"rank deficient QR lost zero diagonal");
    require(work.formQ(bad),"rank deficient formQ status changed");
}
void canonicalChecks() {
    PureCanonicalWorkspace work;
    for(int n:{16,32,81,36,32}) {
        MatType original=gaussian(n);original(0,0)*=1.+2.*std::numeric_limits<double>::epsilon();
        MatType ref=original,reused=original;DataType detRef,detReused;
        const bool a=pureProjectorCanonicalMicrocorrect(ref,detRef,1e-10,true);
        const bool b=pureProjectorCanonicalMicrocorrect(reused,detReused,1e-10,true,&work);
        require(a&&b,"canonical valid roundoff frame failed");
        close((ref-reused).norm(),ref.norm(),"canonical frame mismatch");
        close(std::abs(detRef-detReused),std::abs(detRef),"canonical determinant mismatch");
        require(pureProjectorCanonicalFrameTrusted(reused,1e-10),"canonical independent certificate failed");
        const MatType returned=reused;
        MatType another=gaussian(n);DataType detAnother;
        require(pureProjectorCanonicalMicrocorrect(another,detAnother,1e-10,true,&work),"canonical reuse failed");
        require((returned-reused).norm()==0.,"canonical returned frame aliases workspace");
        for(int mode:{0,1}) {
            MatType failure=gaussian(n);
            if(mode==0)failure*=1.01; // orthogonality failure
            else {const MatType conjugate=failure.conjugate();failure=(failure+1e-6*conjugate)/std::sqrt(1.+1e-12);} // bounded orth, excessive isotropy correction
            ref=failure;reused=failure;
            const bool oldOK=pureProjectorCanonicalMicrocorrect(ref,detRef,1e-10,true);
            const bool newOK=pureProjectorCanonicalMicrocorrect(reused,detReused,1e-10,true,&work);
            require(!oldOK&&!newOK,"canonical failure gate changed");
            require((ref-failure).norm()==0.&&(reused-failure).norm()==0.,"failed canonical correction committed q");
        }
        reused=original;
        require(pureProjectorCanonicalMicrocorrect(reused,detReused,1e-10,true,&work),"canonical failure poisoned workspace");
    }
    MatType bad=MatType::Zero(3,2);DataType det;
    require(!pureProjectorCanonicalMicrocorrect(bad,det,1e-10,true,&work),"canonical invalid dimensions accepted");
    bad=gaussian(32);bad(0,0)=DataType(std::numeric_limits<double>::infinity(),0.);
    require(!pureProjectorCanonicalMicrocorrect(bad,det,1e-10,true,&work),"canonical nonfinite input accepted");
    auto extended=gaussian(4).cast<PureCanonicalExtendedComplex>().eval();
    double relative=0.,isotropy=0.;
    require(pureProjectorCanonicalMicrocorrectExtended(extended,relative,isotropy,1e-10),"extended canonical helper failed");
    require(relative<=pureProjectorCanonicalCorrectionBound(8)&&std::isfinite(isotropy),"extended canonical certificate failed");
}
void realChecks() {
    PureRealMagnitudeWorkspace work;
    const auto initialStats=pureCanonicalCorrectionStats;
    int successes=0;
    for(int n:{16,32,81,36,32}) {
        Eigen::MatrixXd basis=Eigen::MatrixXd::Zero(2*n,n);
        basis.topRows(n)=basis.bottomRows(n)=Eigen::MatrixXd::Identity(n,n)/std::sqrt(2.);
        Eigen::MatrixXd triangular=Eigen::MatrixXd::Identity(n,n);
        double expectedDet=1.;
        for(int j=0;j<n;++j) {
            triangular(j,j)=1.+.001*j;expectedDet*=triangular(j,j);
            for(int i=0;i<j;++i)triangular(i,j)=.004*std::sin(.31*i+.27*j);
        }
        const Eigen::MatrixXd original=basis*triangular;
        Eigen::MatrixXd retained;
        for(int repeat=0;repeat<2;++repeat) {
            work.next_q=original;double det=0.;
            require(work.orthogonalize(det,1e-10),"real canonical QR failed");++successes;
            close(std::abs(std::abs(det)-expectedDet),expectedDet,"real QR determinant mismatch");
            close((work.next_q.transpose()*work.next_q-Eigen::MatrixXd::Identity(n,n)).norm(),std::sqrt(double(n)),"real orthogonality mismatch");
            Eigen::MatrixXd metric=work.next_q;metric.bottomRows(n)*=-1.;
            close((work.next_q.transpose()*metric).norm(),std::sqrt(double(n)),"real isotropy mismatch");
            close((work.next_q*work.next_q.transpose()-basis*basis.transpose()).norm(),basis.norm(),"real propagated span mismatch");
            if(repeat==0)retained=work.next_q;
            else close((retained-work.next_q).norm(),retained.norm(),"real repeated input mismatch");
        }
        work.next_q=original;work.next_q(0,0)=std::numeric_limits<double>::infinity();double det;
        require(!work.orthogonalize(det,1e-10),"real nonfinite input accepted");
        // Full-rank orthonormal but non-isotropic span: the mandatory
        // correction is too large and must fail after the initial QR.
        work.next_q=Eigen::MatrixXd::Zero(2*n,n);work.next_q.topRows(n).setIdentity();
        require(!work.orthogonalize(det,1e-10),"real excessive correction accepted");
        work.next_q=original;require(work.orthogonalize(det,1e-10),"real recovery failed");++successes;
    }
    require(pureCanonicalCorrectionStats.real_calls>=initialStats.real_calls+successes,"real call counters lost");
    require(pureCanonicalCorrectionStats.failures>=initialStats.failures+5,"real failure counters lost");
}
void svdChecks() {
    PureGreenSvd work;
    for(int n:{4,32,81,36,16,32})for(bool accelerated:{false,true})for(bool vectors:{true,false}) {
        MatType matrix=MatType::Identity(n,n);
        for(int i=0;i<n;++i)matrix(i,i)=1.+.01*i;
        work.compute(matrix,accelerated,vectors);work.setThreshold(1e-12);
        require(work.rank()==n,"reused SVD rank mismatch");
        close(std::abs(work.singularValues()(0)-(1.+.01*(n-1))),1.+.01*(n-1),"reused SVD spectrum mismatch");
        if(vectors) {
            const MatType solution=work.solve(MatType::Identity(n,n));
            close((matrix*solution-MatType::Identity(n,n)).norm(),std::sqrt(double(n)),"reused SVD solve mismatch");
        } else {
            bool threw=false;try {work.solve(MatType::Identity(n,n));}catch(const std::logic_error &){threw=true;}
            require(threw,"values-only SVD unexpectedly solved");
        }
        matrix(n-1,n-1)=.5e-12;work.compute(matrix,accelerated,vectors);work.setThreshold(1e-12);
        require(work.rank()==n-1,"SVD rank boundary relaxed");
    }
}
void greenChecks() {
    PureProjectorGreenWorkspace work;
    for(int n:{4,32,81,36,16,32})for(bool accelerated:{false,true})for(bool lu:{false,true}) {
        PureProjectorOptions options;options.accelerated_linear_algebra=accelerated;
        options.prefer_lu_green_solve=lu;
        const MatType right=gaussian(n),left=right;
        const auto ref=pureProjectorGreen(right,left,options);
        const auto reused=pureProjectorGreen(right,left,options,&work);
        require(ref.ok()&&reused.ok(),"Green valid frames failed");
        require(ref.overlap_rank==reused.overlap_rank,"Green rank changed");
        close(std::abs(ref.overlap_rcond-reused.overlap_rcond),1.,"Green rcond changed");
        close((ref.green-reused.green).norm(),ref.green.norm(),"Green workspace mismatch");
        const MatType expected=MatType::Identity(2*n,2*n)-2.*right*right.adjoint();
        close((reused.green-expected).norm(),expected.norm(),"independent Green expression mismatch");
        const MatType saved=reused.green;
        MatType wrongRows=MatType::Zero(2*n+1,n);
        require(pureProjectorGreen(right,wrongRows,options,&work).status==PureProjectorStatus::invalid_dimension,"Green invalid dimensions accepted");
        MatType changed=left;changed.col(n-1)*=.5e-12;
        const auto rankRef=pureProjectorGreen(right,changed,options);
        const auto rankReused=pureProjectorGreen(right,changed,options,&work);
        require(rankRef.status==PureProjectorStatus::overlap_rank_deficient&&rankReused.status==rankRef.status,"Green rank boundary changed");
        options.rank_tolerance=1e-14;options.minimum_overlap_rcond=1e-8;
        for(double scale:{.5e-8,2.e-8}) {
            changed=left;changed.col(n-1)*=scale;
            const auto a=pureProjectorGreen(right,changed,options);
            const auto b=pureProjectorGreen(right,changed,options,&work);
            require(a.status==b.status,"Green rcond boundary status changed");
            require(scale<1e-8?a.status==PureProjectorStatus::overlap_ill_conditioned:a.ok(),"Green rcond boundary unexpected");
            if(a.ok())close((a.green-b.green).norm(),a.green.norm(),"Green rcond trusted output mismatch");
        }
        require((reused.green-saved).norm()==0.,"returned Green aliases scratch");
        changed=left;changed(0,0)=DataType(std::numeric_limits<double>::quiet_NaN(),0.);
        require(pureProjectorGreen(right,changed,options,&work).status==PureProjectorStatus::nonfinite_input,"Green nonfinite accepted");
        require(pureProjectorGreen(right,left,options,&work).ok(),"Green failure poisoned workspace");
    }
}
void greenCopyChecks() {
    PureProjectorOptions options;options.accelerated_linear_algebra=true;
    PureProjectorGreenWorkspace source;
    const MatType frame=gaussian(32);
    const auto original=pureProjectorGreen(frame,frame,options,&source);
    require(original.ok(),"Green copy setup failed");
    MatType maintained=frame;DataType det;
    require(pureProjectorCanonicalMicrocorrect(maintained,det,1e-10,true,&source.canonical),"Green canonical copy setup failed");
    const MatType overlap=source.overlap,rhs=source.right_hand_side,solved=source.solved;
    const MatType canonicalGram=source.canonical.gram;
    const dVecType singular=source.svd.singularValues();
    PureProjectorGreenWorkspace constructed(source),assigned;
    const MatType small=gaussian(4);
    require(pureProjectorGreen(small,small,options,&assigned).ok(),"Green assignment destination setup failed");
    assigned=source;
    for(auto *copy:{&constructed,&assigned}) {
        require(copy->overlap.size()==0&&copy->right_hand_side.size()==0&&copy->solved.size()==0&&
                copy->canonical.gram.size()==0,"copied Green workspace retained source scratch");
        const MatType other=gaussian(36);
        const auto result=pureProjectorGreen(other,other,options,copy);
        require(result.ok(),"copied Green workspace unusable");
        require(copy->overlap.data()!=source.overlap.data()&&
                copy->right_hand_side.data()!=source.right_hand_side.data()&&
                copy->solved.data()!=source.solved.data(),"copied Green scratch aliases source");
        MatType next=other;
        require(pureProjectorCanonicalMicrocorrect(next,det,1e-10,true,&copy->canonical),"copied canonical workspace unusable");
        require(copy->canonical.gram.data()!=source.canonical.gram.data(),"copied canonical scratch aliases source");
        require((source.overlap-overlap).norm()==0.&&(source.right_hand_side-rhs).norm()==0.&&
                (source.solved-solved).norm()==0.&&(source.canonical.gram-canonicalGram).norm()==0.&&
                (source.svd.singularValues()-singular).norm()==0.,"copy use changed source decompositions");
        require((original.green-(MatType::Identity(64,64)-2.*frame*frame.adjoint())).norm()<2e-12,
                "copy use changed owned Green result");
    }
    source=source;
    require((source.overlap-overlap).norm()==0.&&(source.solved-solved).norm()==0.,"Green workspace self-assignment cleared live scratch");
    require(pureProjectorGreen(frame,frame,options,&source).ok(),"source Green workspace unusable after copies");
}
}
int main() {try {
    qrChecks();canonicalChecks();realChecks();svdChecks();greenChecks();greenCopyChecks();
    std::cout<<std::setprecision(17)<<"PASS workspace_reuse checks="<<checks<<" worst_scaled_error="<<worst
             <<" real_calls="<<pureCanonicalCorrectionStats.real_calls
             <<" complex_calls="<<pureCanonicalCorrectionStats.complex_calls
             <<" failures="<<pureCanonicalCorrectionStats.failures<<'\n';
    return 0;
}catch(const std::exception &error){std::cerr<<"workspace_reuse_validate FAIL: "<<error.what()<<'\n';return 1;}}
