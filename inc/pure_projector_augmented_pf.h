#ifndef PURE_PROJECTOR_AUGMENTED_PF_H
#define PURE_PROJECTOR_AUGMENTED_PF_H
#include "types.h"
#include <Eigen/Dense>
#include <complex>
#include <limits>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>

namespace pure_projector_augmented_pf {

template<class T> using AugC=std::complex<T>;
template<class T> using AugM=Eigen::Matrix<AugC<T>,Eigen::Dynamic,Eigen::Dynamic>;
template<class T> struct AugLog {
 T logabs=0; AugC<T> phase=1; T minimum_pivot=std::numeric_limits<T>::infinity(),growth=0,factor_backward_bound=0; bool ok=true;
};
template<class T> AugLog<T> augLogPf(AugM<T> a,bool certificate=false) {
 AugLog<T> out;const int d=a.rows();T initial=0;
 AugM<T> original,l;std::vector<AugC<T>> pivots;
 if(certificate){original=a;l=AugM<T>::Identity(d,d);pivots.reserve(d/2);}
 for(int j=0;j<d;++j)for(int i=0;i<d;++i)initial=std::max(initial,std::abs(a(i,j)));
 T largest=initial;
 for(int k=0;k<d;k+=2) {
  int p=k,q=k+1;T best=0;
  for(int j=k+1;j<d;++j)for(int i=k;i<j;++i)if(std::abs(a(i,j))>best){best=std::abs(a(i,j));p=i;q=j;}
  if(!(best>0)||!std::isfinite(best)){out.ok=false;return out;}
  if(p!=k){a.row(p).swap(a.row(k));a.col(p).swap(a.col(k));out.phase=-out.phase;
   if(certificate){original.row(p).swap(original.row(k));original.col(p).swap(original.col(k));for(int j=0;j<k;++j)std::swap(l(p,j),l(k,j));}}
  if(q!=k+1){a.row(q).swap(a.row(k+1));a.col(q).swap(a.col(k+1));out.phase=-out.phase;
   if(certificate){original.row(q).swap(original.row(k+1));original.col(q).swap(original.col(k+1));for(int j=0;j<k;++j)std::swap(l(q,j),l(k+1,j));}}
  const AugC<T> pivot=a(k,k+1);const T ab=std::abs(pivot);
  out.minimum_pivot=std::min(out.minimum_pivot,ab);out.logabs+=std::log(ab);out.phase*=pivot/ab;
  if(certificate){pivots.push_back(pivot);for(int i=k+2;i<d;++i){l(i,k)=-a(k+1,i)/pivot;l(i,k+1)=a(k,i)/pivot;}}
  for(int j=k+2;j<d;++j)for(int i=k+2;i<j;++i) {
   a(i,j)+=(a(k+1,i)*a(k,j)-a(k,i)*a(k+1,j))/pivot;
   a(j,i)=-a(i,j);largest=std::max(largest,std::abs(a(i,j)));
  }
 }
 if(certificate) {
  // A posteriori certificate of P A P^T=L T L^T. Treat the stored L,T
  // as exact and bound reconstruction rounding conservatively, entrywise.
  // Each complex summand uses <=64*d real rounded operations. The modulus
  // sum gives a standard gamma_m absolute summation/product bound.
  const T gamma=(T(64*d)*std::numeric_limits<T>::epsilon())/(T(1)-T(64*d)*std::numeric_limits<T>::epsilon());
  T boundSquared=0;for(int i=0;i<d;++i)boundSquared+=std::norm(original(i,i));
  for(int j=0;j<d;++j)for(int i=0;i<j;++i){AugC<T> value=0;T magnitudeSum=0;
   for(int k=0;k<d;k+=2){value+=pivots[k/2]*(l(i,k)*l(j,k+1)-l(i,k+1)*l(j,k));magnitudeSum+=std::abs(pivots[k/2])*(std::abs(l(i,k))*std::abs(l(j,k+1))+std::abs(l(i,k+1))*std::abs(l(j,k)));}
   const T upper=std::abs(value-original(i,j))+gamma*magnitudeSum;
   const T lower=std::abs(-value-original(j,i))+gamma*magnitudeSum;
   boundSquared+=upper*upper+lower*lower;
  }
  out.factor_backward_bound=std::sqrt(boundSquared)*(T(1)+T(64*d)*std::numeric_limits<T>::epsilon());
 }
 out.phase/=std::abs(out.phase);out.growth=largest/initial;return out;
}
// Independent row-pivoted determinant elimination (no Pfaffian work reused).
template<class T> AugLog<T> augLogDet(AugM<T> a,bool certificate=false) {
 AugLog<T> out;const int d=a.rows();
 AugM<T> original,l;if(certificate){original=a;l=AugM<T>::Identity(d,d);}
 for(int k=0;k<d;++k) {
  int p=k;for(int i=k+1;i<d;++i)if(std::abs(a(i,k))>std::abs(a(p,k)))p=i;
  if(p!=k){a.row(p).swap(a.row(k));out.phase=-out.phase;
   if(certificate){original.row(p).swap(original.row(k));for(int j=0;j<k;++j)std::swap(l(p,j),l(k,j));}}
  const auto pivot=a(k,k);const T ab=std::abs(pivot);
  if(!(ab>0)||!std::isfinite(ab)){out.ok=false;return out;}
  out.minimum_pivot=std::min(out.minimum_pivot,ab);out.logabs+=std::log(ab);out.phase*=pivot/ab;
  for(int i=k+1;i<d;++i){const auto f=a(i,k)/pivot;if(certificate)l(i,k)=f;for(int j=k+1;j<d;++j)a(i,j)-=f*a(k,j);}
 }
 if(certificate){
  const T gamma=T(64*d)*std::numeric_limits<T>::epsilon()/(T(1)-T(64*d)*std::numeric_limits<T>::epsilon());T boundSquared=0;
  for(int j=0;j<d;++j)for(int i=0;i<d;++i){AugC<T> value=0;T magnitudeSum=0;for(int k=0;k<=std::min(i,j);++k){value+=l(i,k)*a(k,j);magnitudeSum+=std::abs(l(i,k))*std::abs(a(k,j));}const T bound=std::abs(value-original(i,j))+gamma*magnitudeSum;boundSquared+=bound*bound;}
  out.factor_backward_bound=std::sqrt(boundSquared)*(T(1)+T(64*d)*std::numeric_limits<T>::epsilon());
 }
 out.phase/=std::abs(out.phase);return out;
}
template<class T> T augCanonicalCorrect(AugM<T>&q) {
 const auto original=q.eval();
 // Only the frame is maintained. No Green matrix is formed or projected.
 for(int iteration=0;iteration<3;++iteration) {
  const AugM<T> defect=q.transpose()*q;
  q=(q-T(.5)*q.conjugate()*defect).eval();
 }
 return (q-original).norm()/original.norm();
}
struct AugCertificate {
 bool ok=false,cayley_ok=false,phase_ok=false; DataType phase=1.;double logabs=0.,log_identity_error=0.,phase_identity_error=0.;
 double double_log_disagreement=0.,double_phase_disagreement=0.,ld_isotropy=0.,relative_ld_correction=0.,skew=0.,max_entry=0.,growth=0.,min_pivot=0.,seconds=0.,predicted_log_ratio=0.;
 double cayley_factor_difference=0.,cayley_factor_residual=0.,cayley_log_identity_error=0.,cayley_phase_identity_error=0.,cayley_predicted_log_ratio=0.;
 bool condition_recorded=false;
 double k_norm=0.,k_inverse_norm=0.,k_condition=0.,k_inverse_residual=0.,k_inverse_backward_error=0.,m_condition=0.,m_inverse_residual=0.,conservative_phase_bound=0.;
 double factor_backward_bound=0.,relative_nonsingularity_bound=0.,posterior_phase_bound=0.,inverse_residual_upper=0.;
 double m_factor_backward_bound=0.,m_inverse_residual_upper=0.,m_relative_nonsingularity_bound=0.,m_posterior_phase_bound=0.,ratio_posterior_phase_bound=0.;
 double phase_branch_margin=0.,input_cayley_residual=0.,input_overlap_rcond=0.;int input_overlap_rank=0;
};
inline AugCertificate augmentedCertificate(const MatType&p0,const MatType&q0,const MatType&b0,const MatType&c0,double tolerance=1e-10,bool forceCondition=false) {
 const auto start=std::chrono::steady_clock::now();AugCertificate out;
 using T=long double;using C=AugC<T>;using M=AugM<T>;
 const int d=p0.rows(),n=p0.cols();M p=p0.cast<C>(),q=q0.cast<C>(),b=b0.cast<C>(),c=c0.cast<C>();
 out.relative_ld_correction=double(std::max(augCanonicalCorrect(p),augCanonicalCorrect(q)));
 out.ld_isotropy=double(std::max((p.transpose()*p).norm(),(q.transpose()*q).norm())/std::sqrt(T(n)));
 const M overlap=p.adjoint()*q;
 const M next=p.adjoint()*b*q;
 const auto dm=augLogDet<T>(overlap),dn=augLogDet<T>(next),db=augLogDet<T>(M::Identity(d,d)+b);
 M u(d,d);u.leftCols(n)=q;u.rightCols(n)=p.conjugate();
 M e=M::Zero(d,d);e.topRightCorner(n,n)=overlap.transpose();e.bottomLeftCorner(n,n)=-overlap;
 M k=M::Zero(3*d,3*d);k.topLeftCorner(d,d)=c;k.block(0,d,d,d)=-M::Identity(d,d);k.block(d,0,d,d)=M::Identity(d,d);
 k.block(d,2*d,d,d)=u;k.block(2*d,d,d,d)=-u.transpose();k.bottomRightCorner(d,d)=e;
 out.skew=double((k+k.transpose()).norm()/std::max(T(1),k.norm()));
 for(int j=0;j<3*d;++j)for(int i=0;i<3*d;++i)out.max_entry=std::max(out.max_entry,double(std::abs(k(i,j))));
 const auto pk=augLogPf<T>(k);const auto pd=augLogPf<double>(k.cast<DataType>());
 const T ordering=((n*(n-1)/2)%2)?T(-1):T(1);
 const C ratio=pk.phase/(ordering*dm.phase);const T logratio=pk.logabs-dm.logabs;
 const C expected=dn.phase/(db.phase*dm.phase);
 out.logabs=double(logratio);out.phase=DataType(double(ratio.real()),double(ratio.imag()));
 out.predicted_log_ratio=double(logratio+T(.5)*db.logabs-T(d/2)*std::log(T(2)));
 out.log_identity_error=double(std::abs(2*logratio-(T(d)*std::log(T(2))+dn.logabs-db.logabs-dm.logabs)));
 out.phase_identity_error=double(std::abs(ratio*ratio-expected));
 out.double_log_disagreement=double(std::abs(T(pd.logabs)-pk.logabs));
 out.double_phase_disagreement=double(std::abs(C(pd.phase.real(),pd.phase.imag())-pk.phase));
 out.growth=double(pk.growth);out.min_pivot=double(pk.minimum_pivot);
 out.ok=dm.ok&&dn.ok&&db.ok&&pk.ok&&pd.ok&&out.relative_ld_correction<=128.*d*std::numeric_limits<double>::epsilon()&&
  out.ld_isotropy<=128.*d*std::numeric_limits<T>::epsilon()&&out.skew<=tolerance&&
  out.log_identity_error<=tolerance&&out.phase_identity_error<=tolerance;
 // Diagnostic only: recover the unique factor represented by the actual C
 // in long double. Never replace the propagated B or alter its entries.
 const M eye=M::Identity(d,d),plus=eye+c,minus=eye-c;
 Eigen::PartialPivLU<M> lu(plus);const M bc=lu.solve(minus);
 out.cayley_factor_difference=double((bc-b).norm()/std::max(T(1),b.norm()));
 out.cayley_factor_residual=double((plus*bc-minus).norm()/std::max(T(1),minus.norm()));
 const auto dnc=augLogDet<T>((p.adjoint()*bc*q).eval());const auto dbc=augLogDet<T>((eye+bc).eval());
 out.cayley_log_identity_error=double(std::abs(2*logratio-(T(d)*std::log(T(2))+dnc.logabs-dbc.logabs-dm.logabs)));
 out.cayley_phase_identity_error=double(std::abs(ratio*ratio-dnc.phase/(dbc.phase*dm.phase)));
 out.cayley_predicted_log_ratio=double(logratio+T(.5)*dbc.logabs-T(d/2)*std::log(T(2)));
 out.cayley_ok=dm.ok&&dnc.ok&&dbc.ok&&pk.ok&&pd.ok&&out.relative_ld_correction<=128.*d*std::numeric_limits<double>::epsilon()&&
  out.ld_isotropy<=128.*d*std::numeric_limits<T>::epsilon()&&out.skew<=tolerance&&
  out.cayley_factor_difference<=128.*d*std::numeric_limits<double>::epsilon()&&out.cayley_factor_residual<=tolerance&&
  out.cayley_log_identity_error<=tolerance&&out.cayley_phase_identity_error<=tolerance;
 if(forceCondition||!out.ok||!out.cayley_ok||out.double_phase_disagreement>tolerance||out.double_log_disagreement>tolerance) {
  const M ki=Eigen::PartialPivLU<M>(k).solve(M::Identity(3*d,3*d));
  const M mi=Eigen::PartialPivLU<M>(overlap).solve(M::Identity(n,n));
  const T kn=k.norm(),kin=ki.norm(),mn=overlap.norm(),min=mi.norm();
  const T kres=(k*ki-M::Identity(3*d,3*d)).norm(),mres=(overlap*mi-M::Identity(n,n)).norm();
  out.condition_recorded=true;out.k_norm=double(kn);out.k_inverse_norm=double(kin);out.k_condition=double(kn*kin);
  out.k_inverse_residual=double(kres/std::sqrt(T(3*d)));out.k_inverse_backward_error=double(kres/(kn*kin+std::sqrt(T(3*d))));
  out.m_condition=double(mn*min);out.m_inverse_residual=double(mres/std::sqrt(T(n)));
  // Deliberately conservative screening bound, not a theorem-specific
  // validated Pfaffian error bound: complete-pivot growth and operation
  // count bound an O(eps) matrix perturbation, then kappa converts it.
  out.conservative_phase_bound=double(T(3*d/2)*(kn*kin)*T(64*(3*d)*(3*d))*std::numeric_limits<T>::epsilon()*pk.growth);
  const auto certifiedPk=augLogPf<T>(k,true);
  out.factor_backward_bound=double(certifiedPk.factor_backward_bound);
  // ||K^-1|| <= ||X||/(1-||I-KX||). This uses the unnormalized residual.
  // A norm below 1 excludes zero along the whole K -> LTL^T segment.
  const T eps=std::numeric_limits<T>::epsilon();
  const T inverseGamma=T(16*3*d)*eps/(T(1)-T(16*3*d)*eps);
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> absoluteProduct=k.cwiseAbs()*ki.cwiseAbs();
  const T residualUpper=(kres+inverseGamma*(absoluteProduct.norm()+std::sqrt(T(3*d))))*(T(1)+T(64*3*d)*eps);
  out.inverse_residual_upper=double(residualUpper);
  const T inverseUpper=residualUpper<T(1)?kin/(T(1)-residualUpper):std::numeric_limits<T>::infinity();
  const T perturbation=inverseUpper*certifiedPk.factor_backward_bound;
  out.relative_nonsingularity_bound=double(perturbation);
  out.posterior_phase_bound=perturbation<T(1)?double(-T(3*d/2)*std::log1p(-perturbation)):std::numeric_limits<double>::infinity();
  const auto certifiedM=augLogDet<T>(overlap,true);
  const T mGamma=T(16*n)*eps/(T(1)-T(16*n)*eps);
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> mAbsoluteProduct=overlap.cwiseAbs()*mi.cwiseAbs();
  const T mResidualUpper=(mres+mGamma*(mAbsoluteProduct.norm()+std::sqrt(T(n))))*(T(1)+T(64*n)*eps);
  const T mInverseUpper=mResidualUpper<T(1)?min/(T(1)-mResidualUpper):std::numeric_limits<T>::infinity();
  const T mPerturbation=mInverseUpper*certifiedM.factor_backward_bound;
  out.m_factor_backward_bound=double(certifiedM.factor_backward_bound);out.m_inverse_residual_upper=double(mResidualUpper);out.m_relative_nonsingularity_bound=double(mPerturbation);
  out.m_posterior_phase_bound=mPerturbation<T(1)?double(-T(n)*std::log1p(-mPerturbation)):std::numeric_limits<double>::infinity();
  out.ratio_posterior_phase_bound=out.posterior_phase_bound+out.m_posterior_phase_bound+double(T(64*3*d)*eps);
 }
 out.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();return out;
}
// Phase-only fallback candidate. It does not supply a Green matrix or alter
// any rank/rcond gate. `root` is the conventional determinant square root
// before its +/- ambiguity is resolved.
inline AugCertificate augmentedPhaseFallbackImpl(const MatType&p,const MatType&q,const MatType&b,const MatType&c,
 const DataType&root,double tolerance,double rankTolerance,double minimumRcond,double orthTolerance) {
 AugCertificate out;const int d=p.rows(),n=p.cols();
 if(d<=0||d!=2*n||q.rows()!=d||q.cols()!=n||b.rows()!=d||b.cols()!=d||c.rows()!=d||c.cols()!=d||
  !p.allFinite()||!q.allFinite()||!b.allFinite()||!c.allFinite()||!std::isfinite(root.real())||!std::isfinite(root.imag())||!(std::abs(root)>0))return out;
 const double normalization=std::sqrt(double(n)),bound=128.*d*std::numeric_limits<double>::epsilon();
 for(const MatType*frame:{&p,&q}) {
  const double orth=(frame->adjoint()*(*frame)-MatType::Identity(n,n)).norm()/normalization;
  const double isotropy=(frame->transpose()*(*frame)).norm()/normalization;
  if(!std::isfinite(orth)||!std::isfinite(isotropy)||orth>orthTolerance||isotropy>bound)return out;
 }
 const MatType overlap=p.adjoint()*q;Eigen::JacobiSVD<MatType> svd(overlap);svd.setThreshold(rankTolerance);
 if(!svd.singularValues().allFinite()||!(svd.singularValues()(0)>0.))return out;
 out.input_overlap_rank=svd.rank();out.input_overlap_rcond=svd.singularValues()(n-1)/svd.singularValues()(0);
 if(out.input_overlap_rank!=n||!std::isfinite(out.input_overlap_rcond)||out.input_overlap_rcond<minimumRcond)return out;
 const MatType eye=MatType::Identity(d,d);
 const double cayleyResidual=((eye+b)*c-(eye-b)).norm()/std::max(1.,2.*std::sqrt(double(d)));
 if(!std::isfinite(cayleyResidual)||cayleyResidual>tolerance)return out;
 out=augmentedCertificate(p,q,b,c,tolerance,true);out.input_overlap_rank=svd.rank();
 out.input_overlap_rcond=svd.singularValues()(n-1)/svd.singularValues()(0);out.input_cayley_residual=cayleyResidual;
 const DataType relative=std::conj(out.phase)*root/std::abs(root);
 const double acuteAngle=std::atan2(std::abs(relative.imag()),std::abs(relative.real()));
 out.phase_branch_margin=std::acos(-1.)/2.-acuteAngle;
 out.phase_ok=out.cayley_ok&&out.condition_recorded&&std::isfinite(out.phase.real())&&std::isfinite(out.phase.imag())&&
  std::isfinite(relative.real())&&std::isfinite(relative.imag())&&std::isfinite(out.phase_branch_margin)&&
  std::isfinite(out.ratio_posterior_phase_bound)&&
  out.relative_nonsingularity_bound<1.&&out.m_relative_nonsingularity_bound<1.&&
  out.phase_branch_margin>0.&&out.ratio_posterior_phase_bound<.5*out.phase_branch_margin;
 return out;
}
struct Statistics {
 std::uint64_t attempts=0,accepted=0,failed=0,raw_b_identity_failures=0,unbounded_phase=0;
 double max_phase_bound=0.,max_branch_fraction=0.,max_raw_b_log_error=0.,max_raw_b_phase_error=0.;
 double max_same_c_log_error=0.,max_same_c_phase_error=0.,max_relative_frame_correction=0.;
 double max_k_condition=0.,max_m_condition=0.,max_cayley_factor_difference=0.,total_seconds=0.;
};
inline thread_local Statistics statistics;
inline AugCertificate augmentedPhaseFallback(const MatType&p,const MatType&q,const MatType&b,const MatType&c,
 const DataType&root,double tolerance,double rankTolerance,double minimumRcond,double orthTolerance) {
 auto &stats=statistics;++stats.attempts;
 const auto out=augmentedPhaseFallbackImpl(p,q,b,c,root,tolerance,rankTolerance,minimumRcond,orthTolerance);
 if(out.phase_ok)++stats.accepted;else ++stats.failed;
 if(!out.ok)++stats.raw_b_identity_failures;
 if(std::isfinite(out.ratio_posterior_phase_bound))stats.max_phase_bound=std::max(stats.max_phase_bound,out.ratio_posterior_phase_bound);else ++stats.unbounded_phase;
 if(out.phase_branch_margin>0.&&std::isfinite(out.phase_branch_margin)&&std::isfinite(out.ratio_posterior_phase_bound))
  stats.max_branch_fraction=std::max(stats.max_branch_fraction,out.ratio_posterior_phase_bound/out.phase_branch_margin);
 stats.max_raw_b_log_error=std::max(stats.max_raw_b_log_error,out.log_identity_error);
 stats.max_raw_b_phase_error=std::max(stats.max_raw_b_phase_error,out.phase_identity_error);
 stats.max_same_c_log_error=std::max(stats.max_same_c_log_error,out.cayley_log_identity_error);
 stats.max_same_c_phase_error=std::max(stats.max_same_c_phase_error,out.cayley_phase_identity_error);
 stats.max_relative_frame_correction=std::max(stats.max_relative_frame_correction,out.relative_ld_correction);
 stats.max_k_condition=std::max(stats.max_k_condition,out.k_condition);stats.max_m_condition=std::max(stats.max_m_condition,out.m_condition);
 stats.max_cayley_factor_difference=std::max(stats.max_cayley_factor_difference,out.cayley_factor_difference);stats.total_seconds+=out.seconds;
 return out;
}
} // namespace pure_projector_augmented_pf
#endif
