#ifndef SkewMatUtils_H
#define SkewMatUtils_H

#include "types.h"
#include<iostream>
#include <limits>
#include "qmc_core_profile.h"
#include "../inc/pfapack/c_interface/pfapack.h"
enum class PfaffianStatus {
    success,
    untrusted_condition,
    untrusted_phase,
    invalid_dimension,
    lapack_failure,
    zero_pivot,
    nonfinite_pivot
};

struct PfaffianResult {
    PfaffianStatus status = PfaffianStatus::invalid_dimension;
    DataType value = DataType(0.0, 0.0);
    int lapack_info = 0;
    double min_pivot = std::numeric_limits<double>::infinity();
    double condition_proxy = 1.0;
    double phase_reality_error = 0.0;

    bool ok() const { return status == PfaffianStatus::success; }
    bool untrusted() const {
        return status == PfaffianStatus::untrusted_condition ||
               status == PfaffianStatus::untrusted_phase;
    }
};

const char *pfaffianStatusName(PfaffianStatus status);

// x^\dagger . x
void complexNorm2(const DataType* x, MKL_INT len, DataType* res);

/*
 * Householder is significantly slower than the Parlett-Reid approach
 * provided by the PFAPACK.
 * This routine is only for test purpose.
*/ 
void SkewMatHouseholder_PureMKL(const int N, DataType* A, DataType* temp, DataType* kVec);

DataType matDet(uint L, DataType* mat, lapack_int* temp);

DataType pfaf(const int N, MatType& A);
PfaffianResult pfafWithStatus(const int N, MatType& A);

inline MatType expm(MatType &H, double lambda)
{
    int N = H.rows();
    MatType V,D,expK(N,N);
    profileCore("expm_eigensolver",[&]{
        Eigen::SelfAdjointEigenSolver<MatType> es(H);
        V = es.eigenvectors();
        D = es.eigenvalues();
    });
    profileCore("expm_assemble",[&]{
        D = D.array() * lambda;
        D = D.array().exp();
        expK.noalias() = V * D.asDiagonal() * V.adjoint();
    });
    return expK;
}
// A2: reusable Hermitian spectral data. One SelfAdjointEigenSolver per node
// simultaneously produces the half-step exponential, the merged full-step
// exponential and the Cayley Green matrix 2(I+B)^{-1}-I = U tanh(a eps/2) U^dag.
struct HermitianSpectral {
    MatType V;
    dVecType D;
    bool ok = false;
};
inline HermitianSpectral selfAdjointSpectral(const MatType &H) {
    HermitianSpectral s;
    Eigen::SelfAdjointEigenSolver<MatType> es(H);
    if(es.info()==Eigen::Success) {
        s.V=es.eigenvectors();
        s.D=es.eigenvalues();
        s.ok=true;
    }
    return s;
}
// U diag(f) U^dag for a real spectral function evaluated on the real spectrum.
inline MatType spectralRecombination(const HermitianSpectral &s,const cVecType &f) {
    MatType r(s.V.rows(),s.V.rows());
    r.noalias()=s.V*f.asDiagonal()*s.V.adjoint();
    return r;
}

inline MatType sinhHQuarterSqrt2(MatType& H) {    Eigen::SelfAdjointEigenSolver<MatType> es(H);
    MatType V = es.eigenvectors();
    MatType D = es.eigenvalues();
    int N = H.rows();
    MatType sinhK(N, N);
    D = D.array() * 0.25;
    MatType Dinv = - D.array();
    D = (D.array().exp() - Dinv.array().exp()) / sqrt(2);
    sinhK.noalias() = V * D.asDiagonal() * V.adjoint();
    return sinhK;
}

// calculate eta directly using
// eta = (-2)^N Pf [...]
// should only be used for testing
void generateMatForEta(const MatType& H, MatType& A);

// note that this calculation happens in place
DataType signOfPfaf(MatType& A);
PfaffianResult signOfPfafWithStatus(MatType& A);
DataType pfaffianForEta(const MatType &H);
DataType pfaffianForSignOfEta(const MatType &H);
DataType pfaffianForSignOfProduct(const MatType &G1, const MatType &G2 /*, bool diagno=false*/);
PfaffianResult pfaffianForSignOfProductWithStatus(const MatType &G1,
                                                  const MatType &G2);

DataType signOfHamiltonian(const MatType &H);
// DataType normalizeToPlusMinus1(DataType x);

#endif
