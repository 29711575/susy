#ifndef KITAEV_SSH_VBS_OBSERVABLES_H
#define KITAEV_SSH_VBS_OBSERVABLES_H

#include "honeycomb_lattice.h"
#include "vbs_momenta.h"
#include "phase_statistics.h"
#include <Eigen/Eigenvalues>
#include <mkl_dfti.h>

namespace kitaev_ssh {
struct Momentum { int n1,n2; double reciprocal_x,reciprocal_y; Position cartesian; };

// Cell Fourier convention: no bond-midpoint or sublattice form factor.
// Storage is cell=(x*Ly+y), mu=x,y,z, exactly as HoneycombLattice::bondIndex.
//
// Three measurement modes (only the set of evaluated momenta differs; the
// Hamiltonian, phonon PBC, QMC weight, proposals and RNG stream are untouched):
//   All      : the original full PBC momentum grid via MKL FFT (Nc momenta),
//              plus the full inverse-transform real-space correlations.
//   Selected : only the user-supplied q points (Nq momenta) via precomputed
//              phases, O(Nq*Nc) per configuration. No full FFT is run and no
//              full real-space correlation is produced. Every selected result
//              is a finite-window Fourier estimator; off-grid points do NOT
//              change the phonon boundary or the allowed momentum grid.
//   Both     : full grid block first (identical names/order/layout as All),
//              followed by the selected block, for cross-validation.
class VBSObservables {
public:
    struct FFTDescriptor {
        DFTI_DESCRIPTOR_HANDLE handle=nullptr;
        FFTDescriptor()=default;
        FFTDescriptor(const FFTDescriptor &)=delete;
        FFTDescriptor &operator=(const FFTDescriptor &)=delete;
        ~FFTDescriptor(){if(handle)DftiFreeDescriptor(&handle);}
    };
    explicit VBSObservables(HoneycombLattice lattice)
        : VBSObservables(std::move(lattice), VBSMode::All, {}) {}
    VBSObservables(HoneycombLattice lattice, VBSMode mode, std::vector<VBSQPoint> qpoints)
        : lattice_(std::move(lattice)), mode_(mode), qpoints_(std::move(qpoints)) {
        for (const auto &grid : allowedVBSMomenta(lattice_))
            momenta_.push_back({grid.n1,grid.n2,grid.u,grid.v,grid.cartesian});
        // Geometry-only inverse phases are reused by every jackknife replica.
        // Bound the cache to 16 MiB; larger lattices keep the direct path.
        const std::size_t nc=lattice_.cells();
        if(mode_!=VBSMode::Selected && nc<=1024) {
            inverse_phases_.reserve(nc*nc);
            for(int d=0;d<int(nc);++d) {
                const auto r=lattice_.cellCoordinate(d);
                for(const auto &k:momenta_)
                    inverse_phases_.push_back(std::polar(1./nc,
                        -2*pi*(k.reciprocal_x*r.x+k.reciprocal_y*r.y)));
            }
        }
        initializeFourierDescriptors();
        if(mode_==VBSMode::Selected || mode_==VBSMode::Both) {
            if(qpoints_.empty()) throw std::invalid_argument("selected/both VBS modes need --vbs-q-file points");
            phases_.reserve(qpoints_.size());
            for(const auto &q:qpoints_) phases_.push_back(qPointPhases(lattice_,q));
            for(std::size_t j=0;j<qpoints_.size();++j) {
                const bool zero=qpoints_[j].on_grid && qpoints_[j].n1==0 && qpoints_[j].n2==0;
                if(zero) candidate_.excludedZero.push_back(int(j));
                else candidate_.participants.push_back(int(j));
            }
        }
    }
    const auto &momenta() const { return momenta_; }
    const auto &lattice() const { return lattice_; }
    const auto &qpoints() const { return qpoints_; }
    const auto &selectedPhases() const { return phases_; }
    VBSMode mode() const { return mode_; }
    int nq() const { return int(qpoints_.size()); }
    struct CandidateInfo { std::vector<int> participants, excludedZero; };
    const CandidateInfo &candidateInfo() const { return candidate_; }

    using Fourier=std::vector<std::array<Complex,3>>;
    Fourier fourier(const std::vector<double> &x) const {
        validate(x);
        if (lattice_.shear() != 0) return shearedFourier(x);
        return rectangularFourier(x);
    }
    Fourier directFourier(const std::vector<double> &x) const {
        validate(x); Fourier out(lattice_.cells());
        for(int q=0;q<lattice_.cells();++q)for(int r=0;r<lattice_.cells();++r) {
            const auto c=lattice_.cellCoordinate(r);
            const auto a1=lattice_.primitive(0),a2=lattice_.primitive(1);
            const Position R{c.x*a1.x+c.y*a2.x,c.x*a1.y+c.y*a2.y};
            const auto k=momenta_[q].cartesian;
            const auto e=std::polar(1.,-k.x*R.x-k.y*R.y);
            for(int mu=0;mu<3;++mu)out[q][mu]+=e*x[3*r+mu];
        }
        return out;
    }

    // Direct quotient-translation check for the all-grid inverse transform:
    // C_mu,nu(d)=(1/Nc) sum_R X_mu(R) X_nu(R+d), with d reduced by the same
    // sheared quotient as the lattice bonds.
    std::vector<Eigen::Matrix3cd> directRealSpaceCorrelation(const std::vector<double> &x) const {
        validate(x);
        const double nc=double(lattice_.cells());
        std::vector<Eigen::Matrix3cd> result(lattice_.cells(),Eigen::Matrix3cd::Zero());
        for (int d=0; d<lattice_.cells(); ++d) {
            const auto delta=lattice_.cellCoordinate(d);
            for (int r=0; r<lattice_.cells(); ++r) {
                const auto source=lattice_.cellCoordinate(r);
                const int target=lattice_.cellIndexWrapped(
                    static_cast<std::int64_t>(source.x)+delta.x,
                    static_cast<std::int64_t>(source.y)+delta.y);
                for (int mu=0; mu<3; ++mu) for (int nu=0; nu<3; ++nu)
                    result[d](mu,nu)+=x[3*r+mu]*x[3*target+nu]/nc;
            }
        }
        return result;
    }

    // Inverse of the all-grid S(q) convention for a single configuration.
    std::vector<Eigen::Matrix3cd> inverseStructureFactor(const Fourier &f) const {
        if (f.size()!=std::size_t(lattice_.cells())) throw std::invalid_argument("VBS Fourier dimension");
        const double nc=double(lattice_.cells());
        std::vector<Eigen::Matrix3cd> result(lattice_.cells(),Eigen::Matrix3cd::Zero());
        for (int d=0; d<lattice_.cells(); ++d) {
            const auto delta=lattice_.cellCoordinate(d);
            for (int q=0; q<lattice_.cells(); ++q) {
                const auto &k=momenta_[q];
                const Complex phase=std::polar(1.,-2.*pi*(k.reciprocal_x*delta.x+
                                                           k.reciprocal_y*delta.y));
                for (int mu=0; mu<3; ++mu) for (int nu=0; nu<3; ++nu)
                    result[d](mu,nu)+=phase*f[q][mu]*std::conj(f[q][nu])/(nc*nc);
            }
        }
        return result;
    }

private:
    Fourier rectangularFourier(const std::vector<double> &x) const {
        const int nx=lattice_.Lx(),ny=lattice_.Ly();
        Fourier out(lattice_.cells());
        auto checked=[](MKL_LONG code) {
            if(code!=DFTI_NO_ERROR)throw std::runtime_error(std::string("VBS FFT: ")+DftiErrorMessage(code));
        };
        std::vector<Complex> line(lattice_.cells());
        for(int mu=0;mu<3;++mu) {
            for(int r=0;r<lattice_.cells();++r)line[r]=x[3*r+mu];
            checked(DftiComputeForward(rectangular_fft_.handle,line.data()));
            for(int q=0;q<lattice_.cells();++q)out[q][mu]=line[q];
        }
        return out;
    }

    Fourier shearedFourier(const std::vector<double> &x) const {
        const int nx=lattice_.Lx(),ny=lattice_.Ly();
        Fourier out(lattice_.cells());
        auto checked=[](MKL_LONG code,const char *stage) {
            if(code!=DFTI_NO_ERROR)throw std::runtime_error(std::string("VBS FFT ")+stage+": "+DftiErrorMessage(code));
        };
        // oneMKL 2026 rejects a tiny C2C rank-1 descriptor in this runtime;
        // rank-2 vectors with a singleton second dimension have the same
        // contiguous layout and commit reliably.
        std::vector<Complex> xLine(nx),yLine(ny),staged(lattice_.cells());
        for (int mu=0; mu<3; ++mu) {
            // G(n1,y): x FFT in storage order x*Ly+y.
            for (int y=0; y<ny; ++y) {
                for (int xCell=0; xCell<nx; ++xCell)
                    xLine[xCell]=x[3*(xCell*ny+y)+mu];
                checked(DftiComputeForward(sheared_x_fft_.handle,xLine.data())," forward x");
                for (int n1=0; n1<nx; ++n1) {
                    const double angle=2.*pi*double(lattice_.shear())*n1*y/
                        (double(nx)*ny);
                    staged[n1*ny+y]=xLine[n1]*std::polar(1.,angle);
                }
            }
            // F(n1,n2): y FFT of the phase-corrected columns.
            for (int n1=0; n1<nx; ++n1) {
                for (int y=0; y<ny; ++y) yLine[y]=staged[n1*ny+y];
                checked(DftiComputeForward(sheared_y_fft_.handle,yLine.data())," forward y");
                for (int n2=0; n2<ny; ++n2) out[n1*ny+n2][mu]=yLine[n2];
            }
        }
        return out;
    }

public:
    // Direct finite-window sums at the user-selected momenta, O(Nq*Nc).
    Fourier selectedFourier(const std::vector<double> &x) const {
        validate(x);
        Fourier f(nq());
        for(int j=0;j<nq();++j)for(int r=0;r<lattice_.cells();++r) {
            const Complex e=phases_[j][r];
            for(int mu=0;mu<3;++mu)f[j][mu]+=e*x[3*r+mu];
        }
        return f;
    }

    // ---- All-mode layout (unchanged): 9*Nc S entries, then 3*Nc Delta. ----
    static constexpr int allStride(int nc) { return 12*nc; }
    std::vector<Complex> measureAll(const std::vector<double> &x) const {
        const auto f=fourier(x); const double nc=lattice_.cells();
        std::vector<Complex> out;
        out.reserve(allStride(lattice_.cells()));
        for(const auto &v:f)for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)out.push_back(v[mu]*std::conj(v[nu])/nc);
        for(const auto &v:f)for(auto z:v)out.push_back(z/nc);
        return out;
    }
    std::vector<std::string> namesAll() const {
        std::vector<std::string> out;
        for(int q=0;q<lattice_.cells();++q)for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)
            out.push_back(matrixName("S",q,mu,nu));
        for(int q=0;q<lattice_.cells();++q)for(int mu=0;mu<3;++mu)
            out.push_back("Delta_q"+std::to_string(q)+"_"+axes[mu]);
        return out;
    }

    // ---- Selected layout: raw S, raw Delta, projected S, projected Delta,
    //      then the uniform background Delta_mu(0) used by the projection. ----
    int selectedStride() const { return 24*nq()+3; }
    static constexpr int selS(int j){return 9*j;}
    static constexpr int selD(int j,int nq){return 9*nq+3*j;}
    static constexpr int selPS(int j,int nq){return 12*nq+9*j;}
    static constexpr int selPD(int j,int nq){return 21*nq+3*j;}
    static constexpr int selBg(int nq){return 24*nq;}
    std::vector<Complex> measureSelected(const std::vector<double> &x) const {
        std::vector<Complex> out; out.reserve(selectedStride());
        appendSelected(x,lattice_.cells(),out);
        return out;
    }
private:
    static void appendGrid(const Fourier &f,double nc,std::vector<Complex> &out) {
        for(const auto &v:f)for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)
            out.push_back(v[mu]*std::conj(v[nu])/nc);
        for(const auto &v:f)for(auto z:v)out.push_back(z/nc);
    }
    void appendSelected(const std::vector<double> &x,double nc,std::vector<Complex> &out) const {
        const auto f=selectedFourier(x);
        std::array<Complex,3> f0{};
        for(int r=0;r<lattice_.cells();++r)for(int mu=0;mu<3;++mu)f0[mu]+=x[3*r+mu];
        auto block=[&](const Fourier &fv,bool delta) {
            if(!delta) for(const auto &v:fv)for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)
                out.push_back(v[mu]*std::conj(v[nu])/nc);
            else for(const auto &v:fv)for(auto z:v)out.push_back(z/nc);
        };
        Fourier proj(nq());
        for(int j=0;j<nq();++j)for(int mu=0;mu<3;++mu)
            proj[j][mu]=f[j][mu]-qpoints_[j].window*f0[mu];
        block(f,false); block(f,true);
        block(proj,false); block(proj,true);
        for(int mu=0;mu<3;++mu)out.push_back(f0[mu]/nc);
    }
public:
    std::vector<std::string> namesSelected() const {
        std::vector<std::string> out; out.reserve(selectedStride());
        for(int j=0;j<nq();++j)for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)
            out.push_back("Ssel_s"+std::to_string(j)+"_"+axes[mu]+axes[nu]);
        for(int j=0;j<nq();++j)for(int mu=0;mu<3;++mu)
            out.push_back("Dsel_s"+std::to_string(j)+"_"+axes[mu]);
        for(int j=0;j<nq();++j)for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)
            out.push_back("Sprj_s"+std::to_string(j)+"_"+axes[mu]+axes[nu]);
        for(int j=0;j<nq();++j)for(int mu=0;mu<3;++mu)
            out.push_back("Dprj_s"+std::to_string(j)+"_"+axes[mu]);
        for(int mu=0;mu<3;++mu)out.push_back(std::string("Dsel0_")+axes[mu]);
        return out;
    }

    // Dispatcher used by the MC at every measured configuration.
    std::vector<Complex> measure(const std::vector<double> &x) const {
        switch(mode_) {
            case VBSMode::All: return measureAll(x);
            case VBSMode::Selected: return measureSelected(x);
            case VBSMode::Both: {
                const auto f=fourier(x);
                std::vector<Complex> out;
                out.reserve(allStride(lattice_.cells())+selectedStride());
                appendGrid(f,lattice_.cells(),out);
                appendSelected(x,lattice_.cells(),out);
                return out;
            }
        }
        throw std::logic_error("unknown VBS mode");
    }
    // Configuration-level diagnostic used during warmup. This is the raw
    // per-configuration max of mz2(q), not an ensemble/reweighted estimator;
    // it must never be mixed into PhaseStatistics.
    std::pair<double,int> rawAllQMax(const std::vector<double> &x) const {
        if(mode_!=VBSMode::All) throw std::invalid_argument("warmup all-q-max needs VBS all mode");
        const auto raw=measureAll(x); const int nc=lattice_.cells();
        double best=-std::numeric_limits<double>::infinity(); int winner=-1;
        for(int q=0;q<nc;++q) {
            const double value=raw[9*q+8].real()/double(nc);
            if(value>best) { best=value; winner=q; }
        }
        return {best,winner};
    }
    std::vector<std::string> names() const {
        switch(mode_) {
            case VBSMode::All: return namesAll();
            case VBSMode::Selected: return namesSelected();
            case VBSMode::Both: {
                auto a=namesAll(); const auto b=namesSelected();
                a.insert(a.end(),b.begin(),b.end()); return a;
            }
        }
        throw std::logic_error("unknown VBS mode");
    }

    // Called only on ensemble ratios and joint jackknife replicas, never a
    // configuration. Keep raw complex S; Hermitian part used ONLY for eigensolve.
    // Selected mode has no grid-derived names: C/C_connected require the full
    // FFT block and must not be represented by placeholder entries.
    std::vector<std::string> derivedNames() const {
        if(mode_==VBSMode::Selected) return derivedNamesSelected();
        std::vector<std::string> out=derivedNamesAll();
        if(mode_==VBSMode::Both) {
            const auto s=derivedNamesSelected();
            out.insert(out.end(),s.begin(),s.end());
        }
        return out;
    }
    std::vector<std::string> derivedNamesAll() const {
        std::vector<std::string> out;
        for(int q=0;q<lattice_.cells();++q)for(const auto *name:{"TrS","mz2","leading_eigenvalue","hermitian_residual"})
            out.push_back(std::string(name)+"_q"+std::to_string(q));
        for(const auto *prefix:{"C","C_connected"})for(int d=0;d<lattice_.cells();++d)
            for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)out.push_back(matrixName(prefix,d,mu,nu));
        return out;
    }
    std::vector<std::string> derivedNamesSelected() const {
        std::vector<std::string> out;
        for(int j=0;j<nq();++j)for(const auto *name:{"TrS","mz2","leading_eigenvalue","hermitian_residual"})
            out.push_back(std::string(name)+"_s"+std::to_string(j));
        for(int j=0;j<nq();++j)for(const auto *name:{"TrS_proj","mz2_proj","leading_eigenvalue_proj","hermitian_residual_proj"})
            out.push_back(std::string(name)+"_s"+std::to_string(j));
        out.push_back("candidate_max_mz2");
        out.push_back("candidate_max_winner_index");
        out.push_back("candidate_participant_count");
        return out;
    }

    std::vector<Complex> derive(const std::vector<Complex> &ratios,std::size_t offset=0) const {
        std::vector<Complex> out;
        if(mode_==VBSMode::All) {
            deriveAll(ratios,offset,out);
        } else if(mode_==VBSMode::Selected) {
            deriveSelected(ratios,offset,out);
        } else {
            // Both mode: the ratio vector is the all block followed by the
            // selected block. Each helper keeps its exact-dimension contract,
            // so slice off exactly the grid block before deriveAll.
            const std::size_t grid=allStride(lattice_.cells());
            if(ratios.size()<offset+grid)
                throw std::invalid_argument("VBS ratio dimension");
            std::vector<Complex> gridRatios(ratios.begin()+offset,ratios.begin()+offset+grid);
            deriveAll(gridRatios,0,out);
            deriveSelected(ratios,offset+grid,out);
        }
        return out;
    }

    static std::string matrixName(const std::string &prefix,int q,int mu,int nu) {
        return prefix+"_q"+std::to_string(q)+"_"+axes[mu]+axes[nu];
    }
private:
    inline static const char *axes[3]={"x","y","z"};
    HoneycombLattice lattice_;
    VBSMode mode_;
    std::vector<Momentum> momenta_;
    std::vector<VBSQPoint> qpoints_;
    std::vector<std::vector<Complex>> phases_;
    CandidateInfo candidate_;
    std::vector<Complex> inverse_phases_;
    FFTDescriptor rectangular_fft_,sheared_x_fft_,sheared_y_fft_;

    void initializeFourierDescriptors() {
        auto create=[](FFTDescriptor &descriptor,int nx,int ny) {
            MKL_LONG shape[2]={nx,ny};
            const auto createStatus=DftiCreateDescriptor(&descriptor.handle,DFTI_DOUBLE,
                DFTI_COMPLEX,2,shape);
            if(createStatus!=DFTI_NO_ERROR)
                throw std::runtime_error(std::string("VBS FFT create: ")+DftiErrorMessage(createStatus));
            const auto commitStatus=DftiCommitDescriptor(descriptor.handle);
            if(commitStatus!=DFTI_NO_ERROR)
                throw std::runtime_error(std::string("VBS FFT commit: ")+DftiErrorMessage(commitStatus));
        };
        if(lattice_.shear()==0) create(rectangular_fft_,lattice_.Lx(),lattice_.Ly());
        else {
            // Singleton rank-2 descriptors avoid a oneMKL rank-1 commit issue.
            create(sheared_x_fft_,lattice_.Lx(),1);
            create(sheared_y_fft_,lattice_.Ly(),1);
        }
    }

    void validate(const std::vector<double> &x) const {
        if(x.size()!=std::size_t(lattice_.bondCount()))throw std::invalid_argument("VBS bond dimension");
        for(double v:x)if(!std::isfinite(v))throw std::invalid_argument("nonfinite VBS coordinate");
    }
    static void eigen3(const Eigen::Matrix3cd &s,Complex &tr,Complex &mz2,Complex &lambda,Complex &residual,double nc) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd> eig((s+s.adjoint())*.5,Eigen::EigenvaluesOnly);
        if(eig.info()!=Eigen::Success)throw std::runtime_error("VBS eigensolve failed");
        tr=s.trace(); mz2=s(2,2)/nc; lambda=eig.eigenvalues()(2);
        residual=Complex((s-s.adjoint()).norm(),0.);
    }
    void deriveAll(const std::vector<Complex> &ratios,std::size_t offset,std::vector<Complex> &out) const {
        const int nc=lattice_.cells();
        if(ratios.size()!=offset+12*nc)throw std::invalid_argument("VBS ratio dimension");
        std::vector<Eigen::Matrix3cd> s(nc),c(nc,Eigen::Matrix3cd::Zero());
        for(int q=0;q<nc;++q) {
            for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)s[q](mu,nu)=ratios[offset+9*q+3*mu+nu];
            Complex tr,mz2,lambda,residual;
            eigen3(s[q],tr,mz2,lambda,residual,double(nc));
            out.push_back(tr);out.push_back(mz2);out.push_back(lambda);out.push_back(residual);
        }
        for(int d=0;d<nc;++d) {
            const auto r=lattice_.cellCoordinate(d);
            for(int q=0;q<nc;++q) {
                const auto &k=momenta_[q];
                const Complex phase=inverse_phases_.empty()
                    ? std::polar(1./nc,-2*pi*(k.reciprocal_x*r.x+k.reciprocal_y*r.y))
                    : inverse_phases_[std::size_t(d)*nc+q];
                c[d]+=phase*s[q];
            }
        }
        for(bool connected:{false,true})for(const auto &matrix:c)for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu) {
            const Complex background=ratios[offset+9*nc+mu]*ratios[offset+9*nc+nu];
            out.push_back(matrix(mu,nu)-(connected?background:Complex(0.)));
        }
    }
    void deriveSelected(const std::vector<Complex> &ratios,std::size_t base,std::vector<Complex> &out) const {
        const int n=nq(), nc=lattice_.cells();
        if(ratios.size()!=base+std::size_t(selectedStride()))
            throw std::invalid_argument("selected VBS ratio dimension");
        auto matrix=[&](int block,int j,int mu,int nu){
            return ratios[base+block+9*j+3*mu+nu];
        };
        // raw ensemble matrices
        for(int j=0;j<n;++j) {
            Eigen::Matrix3cd s=Eigen::Matrix3cd::Zero();
            for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)s(mu,nu)=matrix(0,j,mu,nu);
            Complex tr,mz2,lambda,residual; eigen3(s,tr,mz2,lambda,residual,double(nc));
            out.push_back(tr);out.push_back(mz2);out.push_back(lambda);out.push_back(residual);
        }
        // per-configuration uniform-projected matrices (separate estimator)
        for(int j=0;j<n;++j) {
            Eigen::Matrix3cd s=Eigen::Matrix3cd::Zero();
            for(int mu=0;mu<3;++mu)for(int nu=0;nu<3;++nu)s(mu,nu)=matrix(12*n,j,mu,nu);
            Complex tr,mz2,lambda,residual; eigen3(s,tr,mz2,lambda,residual,double(nc));
            out.push_back(tr);out.push_back(mz2);out.push_back(lambda);out.push_back(residual);
        }
        // Candidate maximum over raw mz2; Q=0 excluded by construction. The
        // argmax is re-evaluated in every jackknife replica because this
        // transform is called per replica; the winner-index error reports
        // candidate switching. Selection bias is not removed by this error.
        if(candidate_.participants.empty()) {
            out.push_back(Complex(0.)); out.push_back(Complex(-1.)); out.push_back(Complex(0.));
        } else {
            int winner=candidate_.participants.front();
            double best=std::numeric_limits<double>::lowest();
            for(int j:candidate_.participants) {
                const double value=matrix(0,j,2,2).real()/double(nc);
                if(value>best){best=value;winner=j;}
            }
            out.push_back(Complex(best,0.));
            out.push_back(Complex(winner,0.));
            out.push_back(Complex(double(candidate_.participants.size()),0.));
        }
    }
};
}
#endif
