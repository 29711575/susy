#ifndef PURE_PROJECTOR_QR_WORKSPACE_H
#define PURE_PROJECTOR_QR_WORKSPACE_H

#include "types.h"
#include <cmath>
#include <limits>
#include <vector>

// Explicitly owned by one caller/chain. No shared scratch: a nested QR must
// use a separate workspace. The packed matrix remains owned by the caller.
class PureComplexQRWorkspace {
    std::vector<DataType> tau_,qr_work_,q_work_;
    lapack_int rows_=0,cols_=0;
    bool factorized_=false;
    DataType *factorized_data_=nullptr;
    static bool resizeWork(const DataType &query,std::vector<DataType> &buffer) {
        if(!std::isfinite(query.real())||query.real()<1.||
           query.real()>double(std::numeric_limits<lapack_int>::max()))return false;
        buffer.resize(static_cast<lapack_int>(query.real()));return true;
    }
    bool prepare(MatType &packed) {
        const lapack_int m=packed.rows(),n=packed.cols();
        if(m==rows_&&n==cols_&&!qr_work_.empty()&&!q_work_.empty())return true;
        // A failed size-changing query must not leave a valid old size with
        // a newly resized tau buffer.
        rows_=cols_=0;
        tau_.resize(n);DataType query;
        if(LAPACKE_zgeqrf_work(LAPACK_COL_MAJOR,m,n,packed.data(),m,tau_.data(),
                              &query,-1)!=0||!resizeWork(query,qr_work_))return false;
        if(LAPACKE_zungqr_work(LAPACK_COL_MAJOR,m,n,n,packed.data(),m,tau_.data(),
                              &query,-1)!=0||!resizeWork(query,q_work_))return false;
        rows_=m;cols_=n;return true;
    }
public:
    bool factorize(MatType &packed) {
        factorized_=false;factorized_data_=nullptr;
        if(packed.cols()<=0||packed.rows()<packed.cols()||!packed.allFinite()||
           !prepare(packed))return false;
        factorized_=LAPACKE_zgeqrf_work(LAPACK_COL_MAJOR,rows_,cols_,packed.data(),rows_,
            tau_.data(),qr_work_.data(),static_cast<lapack_int>(qr_work_.size()))==0;
        if(factorized_)factorized_data_=packed.data();
        return factorized_;
    }
    bool formQ(MatType &packed) {
        const bool ready=factorized_&&packed.rows()==rows_&&packed.cols()==cols_&&
            packed.data()==factorized_data_;
        factorized_=false;factorized_data_=nullptr;
        if(!ready)return false;
        return LAPACKE_zungqr_work(LAPACK_COL_MAJOR,rows_,cols_,cols_,packed.data(),rows_,
            tau_.data(),q_work_.data(),static_cast<lapack_int>(q_work_.size()))==0;
    }
};

#endif
