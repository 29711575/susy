// Independent validation of the A2 spectral node identities used by the
// optimized merged-node path:
//   B(a)=U exp(-a eps) U^dag
//   C(a)=2(I+B)^{-1}-I = U tanh(a eps/2) U^dag
// against (half*half) construction and an independent FullPivLU solve, over
// random physical generators K=-i A (A real skew-symmetric, hence K Hermitian),
// including zero/degenerate/large spectra and identity (K=0) factors.
//
// The Cayley matrix C is Hermitian with purely imaginary entries, i.e.
// transpose-skew (C^T=-C); it is NOT an involution (the projector Green
// function is). In the extreme ill-conditioned tail where FullPivLU(I+B)
// declares non-invertibility, production code routes through FullPivLU
// (rcond sentinel); here that tail is recorded separately, not failed.
#include "skewMatUtils.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

using Complex=std::complex<double>;
using Mat=Eigen::MatrixXcd;
using RVec=Eigen::VectorXd;

static bool referenceCayley(const Mat &B,Mat &out) {
    const int d=B.rows();
    Mat I=Mat::Identity(d,d);
    Eigen::FullPivLU<Mat> lu(I+B);
    if(!lu.isInvertible()) return false;
    out=lu.solve(2.0*I)-I;
    return true;
}

int main() {
    std::mt19937_64 rng(20260920);
    std::normal_distribution<double> nd(0.,1.);
    double maxFullVsProduct=0,maxCayleyVsLu=0,maxResidual=0,maxHalfCayleyVsLu=0,maxSkew=0;
    long cases=0,failures=0,tail=0;
    auto check=[&](const Mat &K,double dt,const char *tag,bool expectTail=false){
        HermitianSpectral s=selfAdjointSpectral(K);
        if(!s.ok){ std::cerr<<"spectral failed "<<tag<<"\n"; ++failures; return; }
        const int d=K.rows();
        Mat I=Mat::Identity(d,d);
        RVec eh=(s.D.array()*(-dt/2)).exp();
        Mat h; h.noalias()=s.V*eh.cast<Complex>().asDiagonal()*s.V.adjoint();
        RVec ef=(s.D.array()*(-dt)).exp();
        Mat Bs; Bs.noalias()=s.V*ef.cast<Complex>().asDiagonal()*s.V.adjoint();
        Mat Bp=(h*h).eval();
        double relB=(Bs-Bp).norm()/std::max(1.0,Bp.norm());
        maxFullVsProduct=std::max(maxFullVsProduct,relB);
        RVec tf=(s.D.array()*(dt/2)).tanh();
        Mat Cf; Cf.noalias()=s.V*tf.cast<Complex>().asDiagonal()*s.V.adjoint();
        Mat Clu; bool invertible=referenceCayley(Bp,Clu);
        double relC=0,relCh=0;
        if(invertible){
            relC=(Cf-Clu).norm()/std::max(1.0,Clu.norm());
            maxCayleyVsLu=std::max(maxCayleyVsLu,relC);
        } else {
            ++tail;
            if(!expectTail){ std::cerr<<"unexpected singular tail "<<tag<<"\n"; ++failures; }
        }
        Mat R=(I+Bp)*Cf-(I-Bp);
        double res=R.norm()/std::max(1.0,2.0*std::sqrt(double(d)));
        // C Hermitian with purely imaginary entries <=> transpose-skew.
        double herm=(Cf-Cf.adjoint()).norm()/std::max(1.0,std::sqrt(double(d)));
        double skew=(Cf+Cf.transpose()).norm()/std::max(1.0,std::sqrt(double(d)));
        maxSkew=std::max(maxSkew,skew);
        RVec th=(s.D.array()*(dt/4)).tanh();
        Mat Ch; Ch.noalias()=s.V*th.cast<Complex>().asDiagonal()*s.V.adjoint();
        Mat Chlu;
        bool halfInvertible=referenceCayley(h,Chlu);
        if(halfInvertible){
            relCh=(Ch-Chlu).norm()/std::max(1.0,Chlu.norm());
            maxHalfCayleyVsLu=std::max(maxHalfCayleyVsLu,relCh);
        } else { ++tail; if(!expectTail){ std::cerr<<"unexpected half singular "<<tag<<"\n"; ++failures; } }
        ++cases;
        const bool tailCase=(!invertible)||(!halfInvertible);
        bool ok=(relB<1e-10&&herm<1e-10&&skew<1e-9&&Cf.allFinite());
        if(!tailCase){
            // Residual certificate and LU agreement are meaningful only in the
            // nonsingular regime; production uses the spectral C there too.
            maxResidual=std::max(maxResidual,res);
            ok=ok&&(res<1e-10&&relC<1e-9&&relCh<1e-9);
        }
        if(!ok){
            ++failures;
            std::cerr<<"FAIL "<<tag<<" relB="<<relB<<" relC="<<relC<<" res="<<res
                     <<" herm="<<herm<<" skew="<<skew<<" relCh="<<relCh<<"\n";
        }
    };
    for(int d:{4,8,16,32,64}){
        for(int rep=0;rep<24;++rep){
            Eigen::MatrixXd A(d,d);
            for(int a=0;a<d;++a)for(int b=a+1;b<d;++b){ A(a,b)=nd(rng); A(b,a)=-A(a,b); }
            A.diagonal().setZero();
            Mat K=Eigen::MatrixXcd::Zero(d,d);
            for(int a=0;a<d;++a)for(int b=0;b<d;++b) K(a,b)=Complex(0,-A(a,b));
            double scale=std::pow(4.0,std::uniform_real_distribution<double>(-1.5,1.5)(rng));
            check(K*scale,0.1,"random");
        }
    }
    {
        const int d=16;
        check(Mat::Zero(d,d),0.1,"zero",false);
        check(Mat::Zero(d,d),5.0,"zero-large-dt",false);
        check(Mat::Zero(d,d),0.1,"degenerate-zero",false);
        Eigen::MatrixXd A(d,d); A.setZero();
        for(int a=0;a+1<d;a+=2){ double v=(a<4?40.0:(a>10?1e-12:1.0)); A(a,a+1)=v; A(a+1,a)=-v; }
        Mat K2=Eigen::MatrixXcd::Zero(d,d);
        for(int a=0;a<d;++a)for(int b=0;b<d;++b) K2(a,b)=Complex(0,-A(a,b));
        check(K2,0.1,"mixed-spectrum",false);
        check(K2,2.0,"mixed-large-dt",true); // expected FullPivLU tail
    }
    std::cout<<"{\"cases\":"<<cases
             <<",\"failures\":"<<failures
             <<",\"singular_tail_cases\":"<<tail
             <<",\"max_full_vs_half_product\":"<<maxFullVsProduct
             <<",\"max_cayley_vs_lu\":"<<maxCayleyVsLu
             <<",\"max_residual\":"<<maxResidual
             <<",\"max_transpose_skew\":"<<maxSkew
             <<",\"max_half_cayley_vs_lu\":"<<maxHalfCayleyVsLu<<"}\n";
    return failures?1:0;
}
