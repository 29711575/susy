#ifndef PURE_PROJECTOR_STACK_H
#define PURE_PROJECTOR_STACK_H

#include "gaussian_trial_state.h"
#include "pure_projector_green.h"
#include "pure_projector_weight.h"
#include "pure_projector_factor_solve.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <sstream>
#include <vector>

enum class PureStackStatus {
    success,
    invalid_configuration,
    factor_solve_failure,
    subspace_failure,
    overlap_untrusted
};

inline const char *pureStackStatusName(PureStackStatus status) {
    switch (status) {
    case PureStackStatus::success: return "success";
    case PureStackStatus::invalid_configuration: return "invalid_configuration";
    case PureStackStatus::factor_solve_failure: return "factor_solve_failure";
    case PureStackStatus::subspace_failure: return "subspace_failure";
    case PureStackStatus::overlap_untrusted: return "overlap_untrusted";
    }
    return "unknown";
}

struct PureStackCheckpoint {
    int cut = 0;
    MatType q_right;
    MatType q_left;
};

class PureProjectorStackManager {
public:
    PureProjectorStackManager(const GaussianTrialState &trial,
                              std::vector<PureProjectorSlice> slices,
                              int blockSize,
                              PureProjectorOptions options = PureProjectorOptions(),
                              int initialCut = 0,bool localUpdates = false,
                              bool retainRebuildCache = true,bool reuseCheckedQr = false,
                              bool cachedFactorSupport = false,
                              bool borrowGreen = false,bool reuseBuffers = false,
                              bool reuseRebuildResult = false)
        : trial_(trial), slices_(std::move(slices)), block_size_(blockSize),
          options_(options),local_updates_(localUpdates),
          retain_rebuild_cache_(retainRebuildCache),reuse_checked_qr_(reuseCheckedQr) {
        cached_factor_support_=localUpdates&&cachedFactorSupport;
        borrow_green_=localUpdates&&borrowGreen;reuse_buffers_=localUpdates&&reuseBuffers;
        reuse_rebuild_result_=localUpdates&&retainRebuildCache&&reuseCheckedQr&&reuseRebuildResult;
        if (block_size_ <= 0 || initialCut < 0 || initialCut > int(slices_.size())) {
            status_ = PureStackStatus::invalid_configuration;
            return;
        }
        for (const auto &slice : slices_) {
            if (slice.matrix.rows() != trial_.Phi.rows() ||
                slice.matrix.cols() != trial_.Phi.rows() ||
                !pureProjectorMatrixFinite(slice.matrix)) {
                status_ = PureStackStatus::invalid_configuration;
                return;
            }
        }
        if(!local_updates_)buildCheckpoints();
        if (status_ != PureStackStatus::success) return;
        if(!rebuildAtCut(initialCut)&&local_updates_)checkObliqueGreenResidual();
    }

    bool ok() const { return status_ == PureStackStatus::success; }
    PureStackStatus status() const { return status_; }
    int cut() const { return cut_; }
    int length() const { return int(slices_.size()); }
    int blockSize() const { return block_size_; }
    int checkpointCount() const { return checkpoint_count_; }
    std::uint64_t rebuildReuseCount() const { return rebuild_reuse_count_; }
    std::uint64_t scaledGreenRebuilds() const { return scaled_green_rebuilds_; }
    std::string failureDetail() const {
        std::ostringstream out;out.precision(17);
        out<<"cut="<<cut_<<" stack_status="<<pureStackStatusName(status_);
        if(green_valid_)out<<" green_status="<<pureProjectorStatusName(cached_green_.status)
            <<" overlap_rank="<<cached_green_.overlap_rank
            <<" overlap_rcond="<<cached_green_.overlap_rcond
            <<" solve_residual="<<cached_green_.solve_residual
            <<" green_residual="<<cached_green_.green_residual
            <<" green_relative_residual="<<cached_green_.green_relative_residual;
        return out.str();
    }
    const MatType &phiRight() const { return phi_right_; }
    const MatType &phiLeft() const { return phi_left_; }
    const std::vector<PureProjectorSlice> &slices() const { return slices_; }
    const GaussianTrialState &trial() const { return trial_; }
    std::size_t rebuildCacheBytes() const {
        std::size_t bytes=0;
        for(const auto& entry:right_rebuild_cache_)bytes+=entry.second.size()*sizeof(DataType);
        for(const auto& entry:left_rebuild_cache_)bytes+=entry.second.size()*sizeof(DataType);
        return bytes;
    }

    // Recover the pre-factor right frame in the live cut's gauge without
    // transporting the left frame or Green. At cut index+1, R = B_index R0
    // (up to an invertible QR gauge), so replacing the factor uses B_new R0.
    // This read-only solve also works when the live cut is a QR checkpoint.
    bool rightBeforeFactorAtCut(int index,MatType &rightBefore) const {
        if(!ok()||index<0||index>=length()||cut_!=index+1)return false;
        PureFactorSolve solve(factorView(slices_[index].matrix),DataType(1,0));
        if(!solve.isInvertible())return false;
        rightBefore=solve.solve(phi_right_);
        return pureProjectorMatrixFinite(rightBefore);
    }

    bool acceptFactorAtCut(int index,const PureProjectorSlice &replacement,
                           const MatType &rightBefore,const MatType *updatedGreen=nullptr) {
        if(!ok()||index<0||index>=length()||cut_!=index+1||
           replacement.matrix.rows()!=trial_.Phi.rows()||
           replacement.matrix.cols()!=trial_.Phi.rows()||
           rightBefore.rows()!=trial_.Phi.rows()||rightBefore.cols()!=trial_.Phi.cols()||
           !pureProjectorMatrixFinite(replacement.matrix))return false;
        green_valid_=false;
        rebuild_result_cut_=-1;
        readonly_result_cut_=-1;
        // A changed factor invalidates only prefixes/suffixes containing it.
        right_rebuild_cache_.erase(right_rebuild_cache_.upper_bound(index),right_rebuild_cache_.end());
        left_rebuild_cache_.erase(left_rebuild_cache_.begin(),left_rebuild_cache_.upper_bound(index));
        slices_[index]=replacement;
        phi_right_=local_updates_?pureFactorApply(factorView(replacement.matrix),rightBefore):
            replacement.matrix*rightBefore;
        if(local_updates_&&updatedGreen){
            if(updatedGreen->rows()!=trial_.Phi.rows()||
               updatedGreen->cols()!=trial_.Phi.rows()||
               !pureProjectorMatrixFinite(*updatedGreen))return false;
            cached_green_.green=*updatedGreen;
            cached_green_.green_residual=std::numeric_limits<double>::infinity();
            cached_green_.green_relative_residual=std::numeric_limits<double>::infinity();
            green_valid_=true;
            return true;
        }
        if(!stabilizeCurrent())return false;
        if(!checkedGreenOk()){status_=PureStackStatus::overlap_untrusted;return false;}
        return true;
    }

    PureProjectorGreenResult green() const {
        return greenView();
    }

    // Valid only until the next manager mutation. Owning snapshots continue
    // to use green(); callers never retain this reference across cut motion.
    const PureProjectorGreenResult &greenView() const {
        if(!green_valid_){
            cached_green_=pureProjectorGreenThinQr(phi_right_, phi_left_, options_);
            green_valid_=true;
        }
        return cached_green_;
    }

    // A temporary measurement manager can transfer its owning result instead
    // of copying a full matrix just before destroying the manager.
    PureProjectorGreenResult takeGreen() && {
        greenView();
        PureProjectorGreenResult result=std::move(cached_green_);
        green_valid_=false;
        return result;
    }

    // Failure deliberately returns an untrusted result: the walker then runs
    // its original fresh measurement path, including its residual recovery.
    PureProjectorGreenResult readOnlyGreenAtCut(int target) {
        if(!ok()||!local_updates_||!retain_rebuild_cache_||!reuse_checked_qr_||
           target<0||target>length())return PureProjectorGreenResult();
        if(target==readonly_result_cut_)return readonly_result_;
        PureProjectorGreenResult result;
        if(target==rebuild_result_cut_){
            result=rebuild_result_;
        }else{
            MatType right,left;bool rightChecked=false,leftChecked=false;
            if(!reconstructFrames(target,right,left,rightChecked,leftChecked,false))return result;
            result=greenFromFrames(right,left,rightChecked,leftChecked);
        }
        if(result.ok()){
            readonly_result_=result;
            readonly_result_cut_=target;
        }
        return result;
    }

    bool moveToCut(int target) {
        if (target < 0 || target > length()) return false;
        while (cut_ < target) if (!forwardOne()) return false;
        while (cut_ > target) if (!backwardOne()) return false;
        return true;
    }

    bool forwardOne() {
        if (!ok() || cut_ >= length()) return false;
        if(local_updates_&&((cut_+1)%block_size_==0||cut_+1==length()))
            return rebuildAtCut(cut_+1);
        const auto &stored = slices_[cut_].matrix;
        MatType dense;
        if(!cached_factor_support_)dense=stored;
        const MatType &factor=dense;
        MatType temporaryRight,temporaryLeft,temporaryGreen;
        MatType &nextRight=reuse_buffers_?next_right_:temporaryRight;
        MatType &nextLeft=reuse_buffers_?next_left_:temporaryLeft;
        MatType &nextGreen=reuse_buffers_?next_green_:temporaryGreen;
        if(reuse_buffers_)pureFactorApplyInto(factorView(stored),phi_right_,nextRight,frame_apply_work_);
        else nextRight=local_updates_?pureFactorApply(factorView(stored),phi_right_):factor*phi_right_;
        PureFactorSolve solve(factorView(stored,PureFactorOperation::adjoint),DataType(1,0));
        if (!solve.isInvertible()) {
            status_ = PureStackStatus::factor_solve_failure;
            return false;
        }
        if(reuse_buffers_)solve.solveInto(phi_left_,nextLeft,frame_solve_work_);
        else nextLeft=solve.solve(phi_left_);
        if (!pureProjectorMatrixFinite(nextRight) || !pureProjectorMatrixFinite(nextLeft)) {
            status_ = PureStackStatus::factor_solve_failure;
            return false;
        }
        const bool propagatedGreen=local_updates_&&green_valid_;
        if(propagatedGreen){
            PureFactorSolve inverseTranspose(factorView(stored,PureFactorOperation::transpose),DataType(1,0));
            if(!inverseTranspose.isInvertible())return false;
            if(options_.direct_green_columns){
                if(reuse_buffers_)pureFactorApplyInto(factorView(stored),cached_green_.green,green_stage_,green_apply_work_);
                else green_stage_=pureFactorApply(factorView(stored),cached_green_.green);
                if(inverseTranspose.solveTransposedInPlace(green_stage_,green_solve_work_))nextGreen.swap(green_stage_);
                else nextGreen=inverseTranspose.solve(MatType(green_stage_.transpose())).transpose();
            }else if(reuse_buffers_){
                pureFactorApplyInto(factorView(stored),cached_green_.green,green_stage_,green_apply_work_);
                green_transposed_=green_stage_.transpose();
                inverseTranspose.solveInto(green_transposed_,green_stage_,green_solve_work_);
                nextGreen=green_stage_.transpose();
            }else nextGreen=inverseTranspose.solve(
                pureFactorApply(factorView(stored),cached_green_.green).transpose()).transpose();
        }
        green_valid_=false;
        ++cut_;
        phi_right_.swap(nextRight);
        phi_left_.swap(nextLeft);
        if ((cut_ % block_size_) == 0 || cut_ == length()) {
            if (!stabilizeCurrent()) return false;
        }else if(propagatedGreen&&nextGreen.size()){
            if(!pureProjectorMatrixFinite(nextGreen))return false;
            cached_green_.green.swap(nextGreen);
            cached_green_.green_residual=std::numeric_limits<double>::infinity();
            cached_green_.green_relative_residual=std::numeric_limits<double>::infinity();
            green_valid_=true;
        }
        if (!checkedGreenOk()) {
            status_ = PureStackStatus::overlap_untrusted;
            return false;
        }
        return true;
    }

    bool backwardOne() {
        if (!ok() || cut_ <= 0) return false;
        if(local_updates_&&(cut_-1)%block_size_==0)return rebuildAtCut(cut_-1);
        const auto &stored = slices_[cut_ - 1].matrix;
        MatType dense;
        if(!cached_factor_support_)dense=stored;
        const MatType &factor=dense;
        PureFactorSolve solve(factorView(stored),DataType(1,0));
        if (!solve.isInvertible()) {
            status_ = PureStackStatus::factor_solve_failure;
            return false;
        }
        MatType temporaryRight,temporaryLeft,temporaryGreen;
        MatType &nextRight=reuse_buffers_?next_right_:temporaryRight;
        MatType &nextLeft=reuse_buffers_?next_left_:temporaryLeft;
        MatType &nextGreen=reuse_buffers_?next_green_:temporaryGreen;
        if(reuse_buffers_){
            solve.solveInto(phi_right_,nextRight,frame_solve_work_);
            pureFactorApplyInto(factorView(stored,PureFactorOperation::adjoint),phi_left_,nextLeft,frame_apply_work_);
        }else{
            nextRight=solve.solve(phi_right_);
            nextLeft=local_updates_?pureFactorApply(factorView(stored,PureFactorOperation::adjoint),phi_left_):
                factor.adjoint()*phi_left_;
        }
        if (!pureProjectorMatrixFinite(nextRight) || !pureProjectorMatrixFinite(nextLeft)) {
            status_ = PureStackStatus::factor_solve_failure;
            return false;
        }
        const bool propagatedGreen=local_updates_&&green_valid_;
        if(propagatedGreen){
            if(options_.direct_green_columns){
                if(reuse_buffers_)solve.solveInto(cached_green_.green,green_stage_,green_solve_work_);
                else green_stage_=solve.solve(cached_green_.green);
                if(pureFactorApplyTransposedInPlace(factorView(stored,PureFactorOperation::transpose),green_stage_,green_apply_work_))
                    nextGreen.swap(green_stage_);
                else nextGreen=pureFactorApply(factorView(stored,PureFactorOperation::transpose),MatType(green_stage_.transpose())).transpose();
            }else if(reuse_buffers_){
                solve.solveInto(cached_green_.green,green_stage_,green_solve_work_);
                green_transposed_=green_stage_.transpose();
                pureFactorApplyInto(factorView(stored,PureFactorOperation::transpose),green_transposed_,green_stage_,green_apply_work_);
                nextGreen=green_stage_.transpose();
            }else nextGreen=pureFactorApply(factorView(stored,PureFactorOperation::transpose),
                solve.solve(cached_green_.green).transpose()).transpose();
        }
        green_valid_=false;
        --cut_;
        phi_right_.swap(nextRight);
        phi_left_.swap(nextLeft);
        if (cut_ == 0) {
            // Re-anchor the endpoint gauge deterministically so the Phase 3A
            // round-trip invariant remains exact.
            return rebuildAtCut(0);
        }
        if ((cut_ % block_size_) == 0) {
            if (!stabilizeCurrent()) return false;
        }else if(propagatedGreen&&nextGreen.size()){
            if(!pureProjectorMatrixFinite(nextGreen))return false;
            cached_green_.green.swap(nextGreen);
            cached_green_.green_residual=std::numeric_limits<double>::infinity();
            cached_green_.green_relative_residual=std::numeric_limits<double>::infinity();
            green_valid_=true;
        }
        if (!checkedGreenOk()) {
            status_ = PureStackStatus::overlap_untrusted;
            return false;
        }
        return true;
    }

    bool forwardToEnd() { return moveToCut(length()); }
    bool backwardToStart() { return moveToCut(0); }
    bool rebuildCurrent() { return rebuildAtCut(cut_); }
    bool rebuildToCut(int target) {
        if(target<0||target>length())return false;
        return rebuildAtCut(target);
    }
    bool rebuildToCutFresh(int target) {
        if(target<0||target>length())return false;
        // A failed cached rebuild can leave us retrying with exactly the same
        // prefix/suffix gauges.  Drop both caches so recovery really
        // re-anchors from the trial state without changing the configuration.
        readonly_result_cut_=-1;
        right_rebuild_cache_.clear();
        left_rebuild_cache_.clear();
        rebuild_result_cut_=-1;
        green_valid_=false;
        status_=PureStackStatus::success;
        if(rebuildAtCut(target))return true;
        return local_updates_&&checkObliqueGreenResidual();
    }

    std::uint64_t configurationHash() const {
        std::uint64_t hash = 1469598103934665603ULL;
        auto mix = [&hash](std::uint64_t value) {
            hash ^= value; hash *= 1099511628211ULL;
        };
        for (const auto &slice : slices_) {
            for (int col=0; col<slice.matrix.cols(); ++col)
                for (int row=0; row<slice.matrix.rows(); ++row) {
                    std::uint64_t re=0,im=0;
                    const double r=slice.matrix(row,col).real();
                    const double i=slice.matrix(row,col).imag();
                    std::memcpy(&re,&r,sizeof(double));
                    std::memcpy(&im,&i,sizeof(double));
                    mix(re); mix(im);
                }
        }
        return hash;
    }

private:
    GaussianTrialState trial_;
    std::vector<PureProjectorSlice> slices_;
    int block_size_ = 1;
    PureProjectorOptions options_;
    bool local_updates_=false;
    bool retain_rebuild_cache_=true,reuse_checked_qr_=false,cached_factor_support_=false;
    PureStackStatus status_ = PureStackStatus::success;
    int cut_ = 0;
    MatType phi_right_;
    MatType phi_left_;
    int checkpoint_count_=0;
    // These are exact snapshots from rebuildAtCut, with its original QR
    // cadence and gauge. Live-cut propagation never populates these caches.
    std::map<int,MatType> right_rebuild_cache_,left_rebuild_cache_;
    mutable bool green_valid_=false;
    mutable PureProjectorGreenResult cached_green_;
    bool borrow_green_=false,reuse_buffers_=false;
    MatType next_right_,next_left_,next_green_,green_stage_,green_transposed_;
    PureFactorApplyWorkspace frame_apply_work_,green_apply_work_;
    PureFactorSolveWorkspace frame_solve_work_,green_solve_work_;
    bool reuse_rebuild_result_=false;
    int rebuild_result_cut_=-1;
    PureProjectorGreenResult rebuild_result_;
    int readonly_result_cut_=-1;
    PureProjectorGreenResult readonly_result_;
    std::uint64_t rebuild_reuse_count_=0;
    std::uint64_t scaled_green_rebuilds_=0;

    bool checkObliqueGreenResidual() {
        if(status_!=PureStackStatus::overlap_untrusted||!green_valid_)return false;
        const auto checked=pureProjectorCheckObliqueGreen(cached_green_,options_);
        if(!checked.ok())return false;
        cached_green_=checked;
        status_=PureStackStatus::success;
        ++scaled_green_rebuilds_;
        return true;
    }

    void restoreRebuildResult() {
        // Both exact QR checkpoints still exist: no factor has changed since
        // this successful solve. Never substitute inverse-propagated frames.
        phi_right_=right_rebuild_cache_.at(rebuild_result_cut_);
        phi_left_=left_rebuild_cache_.at(rebuild_result_cut_);
        cut_=rebuild_result_cut_;
        cached_green_=rebuild_result_;
        green_valid_=true;
        status_=PureStackStatus::success;
        ++rebuild_reuse_count_;
    }

    bool checkedGreenOk()const{
        if(borrow_green_)return greenView().ok();
        const PureProjectorGreenResult checked=green();
        return checked.ok();
    }

    PureFactorView factorView(const PureFactorMatrix &factor,
                              PureFactorOperation operation=PureFactorOperation::normal)const{
        return PureFactorView(factor,operation,cached_factor_support_,options_.pair_block_factors);
    }

    bool stabilizeCurrent() {
        green_valid_=false;
        const ThinQrResult right = thinQrSubspace(phi_right_, options_);
        const ThinQrResult left = thinQrSubspace(phi_left_, options_);
        if (!right.ok() || !left.ok()) {
            status_ = PureStackStatus::subspace_failure;
            return false;
        }
        phi_right_ = right.q;
        phi_left_ = left.q;
        return true;
    }

    // Both callers use the same factor order and QR boundaries. Read-only
    // reconstruction consumes certified checkpoints without populating the
    // live caches, changing the sweep cut, or using inverse-propagated frames.
    PureFactorApplyWorkspace rebuild_apply_work_;
    MatType rebuild_frame_work_;
    void applyRebuildFactor(const PureFactorView &factor,MatType &frame){
        if(pureFactorApplyInPlace(factor,frame,rebuild_apply_work_))return;
        pureFactorApplyInto(factor,frame,rebuild_frame_work_,rebuild_apply_work_);
        frame.swap(rebuild_frame_work_);
    }
    bool reconstructFrames(int target,MatType &right,MatType &left,
                           bool &rightChecked,bool &leftChecked,bool retain) {
        const bool reuse=local_updates_&&options_.reuse_rebuild_workspace;
        if(!reuse)right = trial_.Phi;
        rightChecked=false;leftChecked=false;
        int rightStart=0;
        auto prefix=right_rebuild_cache_.upper_bound(target);
        if(prefix!=right_rebuild_cache_.begin()){
            --prefix;rightStart=prefix->first;right=prefix->second;rightChecked=true;
        }else if(reuse)right=trial_.Phi;
        for (int index=rightStart; index<target; ++index) {
            rightChecked=false;
            if(reuse)applyRebuildFactor(factorView(slices_[index].matrix),right);
            else right = local_updates_?pureFactorApply(factorView(slices_[index].matrix),right):
                slices_[index].matrix*right;
            if (((index+1)%block_size_)==0) {
                ThinQrResult qr=thinQrSubspace(right,options_);
                if(!qr.ok()){return false;}
                if(reuse)right.swap(qr.q);else right=qr.q;
                rightChecked=true;
                if(retain)right_rebuild_cache_[index+1]=right;
            }
        }
        if(!reuse)left = trial_.Phi;
        int leftStart=length();
        auto suffix=left_rebuild_cache_.lower_bound(target);
        if(suffix!=left_rebuild_cache_.end()){
            leftStart=suffix->first;left=suffix->second;leftChecked=true;
        }else if(reuse)left=trial_.Phi;
        int propagated=length()-leftStart;
        for (int index=leftStart-1; index>=target; --index) {
            leftChecked=false;
            if(reuse)applyRebuildFactor(factorView(slices_[index].matrix,PureFactorOperation::adjoint),left);
            else left = local_updates_?pureFactorApply(factorView(slices_[index].matrix,PureFactorOperation::adjoint),left):
                slices_[index].matrix.adjoint()*left;
            ++propagated;
            if ((propagated%block_size_)==0) {
                ThinQrResult qr=thinQrSubspace(left,options_);
                if(!qr.ok()){return false;}
                if(reuse)left.swap(qr.q);else left=qr.q;
                leftChecked=true;
                if(retain)left_rebuild_cache_[index]=left;
            }
        }
        return true;
    }

    PureProjectorGreenResult greenFromFrames(const MatType &right,const MatType &left,
                                             bool rightChecked,bool leftChecked) const {
        if(!local_updates_||!reuse_checked_qr_)
            return pureProjectorGreenThinQr(right,left,options_);
        MatType rightQ,leftQ;
        const MatType *r=&right,*l=&left;
        if(!rightChecked){
            auto qr=thinQrSubspace(right,options_);
            if(!qr.ok()){PureProjectorGreenResult result;result.status=qr.status;return result;}
            rightQ=std::move(qr.q);r=&rightQ;
        }
        if(!leftChecked){
            auto qr=thinQrSubspace(left,options_);
            if(!qr.ok()){PureProjectorGreenResult result;result.status=qr.status;return result;}
            leftQ=std::move(qr.q);l=&leftQ;
        }
        return pureProjectorGreen(*r,*l,options_);
    }

    bool rebuildAtCut(int target) {
        if(reuse_rebuild_result_&&ok()&&target==rebuild_result_cut_){
            restoreRebuildResult();
            return true;
        }
        rebuild_result_cut_=-1;
        green_valid_=false;
        MatType right,left;
        bool rightChecked=false,leftChecked=false;
        if(!reconstructFrames(target,right,left,rightChecked,leftChecked,retain_rebuild_cache_)){
            status_=PureStackStatus::subspace_failure;return false;
        }
        phi_right_.swap(right);
        phi_left_.swap(left);
        cut_=target;
        status_=PureStackStatus::success;
        if(local_updates_&&reuse_checked_qr_){
            cached_green_=greenFromFrames(phi_right_,phi_left_,rightChecked,leftChecked);
            green_valid_=true;
        }
        if(!checkedGreenOk()){status_=PureStackStatus::overlap_untrusted;return false;}
        // Cache only a solve whose input gauges are exactly the two retained,
        // certified checkpoints at this cut. Other cuts use the original path.
        if(reuse_rebuild_result_&&rightChecked&&leftChecked&&
           right_rebuild_cache_.count(target)&&left_rebuild_cache_.count(target)){
            rebuild_result_=cached_green_;
            rebuild_result_cut_=target;
        }
        return true;
    }

    void buildCheckpoints() {
        checkpoint_count_=0;
        std::vector<int> cuts;
        for (int cut=0; cut<=length(); cut+=block_size_) cuts.push_back(cut);
        if (cuts.empty() || cuts.back()!=length()) cuts.push_back(length());
        // Preserve all construction-time QR/rank/residual checks. The old
        // checkpoint matrices were never read; do not retain duplicate copies.

        MatType right=trial_.Phi;
        std::size_t next=0;
        for(int cut=0;cut<=length();++cut){
            if(next<cuts.size()&&cut==cuts[next]){
                ThinQrResult qr=thinQrSubspace(right,options_);
                if(!qr.ok()){status_=PureStackStatus::subspace_failure;return;}
                right=qr.q;++next;
            }
            if(cut<length())right=slices_[cut].matrix*right;
        }

        MatType left=trial_.Phi;
        int previous=int(cuts.size())-1;
        for(int cut=length();cut>=0;--cut){
            if(previous>=0&&cut==cuts[previous]){
                ThinQrResult ql=thinQrSubspace(left,options_);
                if(!ql.ok()){status_=PureStackStatus::subspace_failure;return;}
                left=ql.q;--previous;
            }
            if(cut>0)left=slices_[cut-1].matrix.adjoint()*left;
        }
        checkpoint_count_=int(cuts.size());
    }
};

#endif
