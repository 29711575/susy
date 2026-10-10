#include "pure_projector_svd.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

namespace {
void require(bool ok,const char *message){if(!ok)throw std::runtime_error(message);}
MatType denseMatrix(int m,int n,std::mt19937_64 &rng){
    std::normal_distribution<double> normal;
    MatType a(m,n);
    for(int j=0;j<n;++j)for(int i=0;i<m;++i)a(i,j)=DataType(normal(rng),normal(rng));
    for(int i=0;i<std::min(m,n);++i)a(i,i)+=double(n);
    return a/double(n);
}
}
int main(){try{
    std::mt19937_64 rng(10022026);
    PureOverlapSingularValues accelerated,reference;
    double maxError=0.,maxLogError=0.;int cases=0;
    for(int n:{31,32,36,81,144,32})for(int extra:{0,3}){
        const MatType a=denseMatrix(n+extra,n,rng);
        accelerated.compute(a,true,1e-12,1e-12,1e-14);
        reference.compute(a,false,1e-12,1e-12,1e-14);
        require(accelerated.rank(1e-12)==reference.rank(1e-12),"dense overlap rank");
        maxError=std::max(maxError,(accelerated.singularValues()-reference.singularValues()).cwiseAbs().maxCoeff());
        maxLogError=std::max(maxLogError,std::abs(accelerated.singularValues().array().log().sum()-reference.singularValues().array().log().sum()));
        require(accelerated.accelerated()==(n>=32),"dimension dispatch");++cases;
    }
    // Exact threshold, deficient, ill-conditioned and determinant-zero inputs
    // must return to the original solver instead of relaxing a trust gate.
    for(double smallest:{0.,1e-14,1e-12,1.01e-12,1e-10}){
        MatType a=MatType::Identity(36,36);a(35,35)=smallest;
        accelerated.compute(a,true,1e-12,1e-12,1e-14);
        reference.compute(a,false,1e-12,1e-12,1e-14);
        require(accelerated.rank(1e-12)==reference.rank(1e-12),"boundary rank");
        require((accelerated.singularValues()-reference.singularValues()).norm()==0.,"diagonal singular values");
        if(smallest<=1.01e-12)require(!accelerated.accelerated(),"boundary reference fallback");++cases;
    }
    MatType a=.1*MatType::Identity(36,36);
    accelerated.compute(a,true,1e-12,1e-12,1e-14);
    require(!accelerated.accelerated(),"determinant-zero reference fallback");++cases;
    require(maxError<1e-12&&maxLogError<1e-10,"values/logdet reference tolerance");
    std::cout<<std::setprecision(17)<<"{\"status\":\"PASS\",\"cases\":"<<cases
        <<",\"max_singular_value_error\":"<<maxError<<",\"max_logdet_error\":"<<maxLogError
        <<",\"accelerated_calls\":"<<accelerated.acceleratedCalls()
        <<",\"reference_calls\":"<<accelerated.referenceCalls()<<"}\n";
}catch(const std::exception&e){std::cerr<<"overlap_svd_validate: "<<e.what()<<'\n';return 1;}}
