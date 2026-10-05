#ifndef KITAEV_SSH_MODEL_H
#define KITAEV_SSH_MODEL_H

#include "honeycomb_lattice.h"
#include "skewMatUtils.h"
#include <utility>

namespace kitaev_ssh {

struct Couplings { double Jx = 1., Jy = 1., Jz = .6, g = 0.; };

class KitaevSSHModel {
public:
    KitaevSSHModel(HoneycombLattice lattice, Couplings couplings = {},
                   std::vector<double> displacement = {}, std::vector<int> gauge = {})
        : lattice_(std::move(lattice)), couplings_(couplings), x_(std::move(displacement)), u_(std::move(gauge)) {
        for (double v : {couplings.Jx,couplings.Jy,couplings.Jz,couplings.g})
            if (!std::isfinite(v)) throw std::invalid_argument("nonfinite model coupling");
        if (x_.empty()) x_.assign(lattice_.bondCount(),0.);
        if (x_.size() != std::size_t(lattice_.bondCount())) throw std::invalid_argument("static X bond count");
        for (double v : x_) if (!std::isfinite(v)) throw std::invalid_argument("nonfinite static X");
        if (u_.empty()) u_ = lattice_.zeroFluxGauge();
        lattice_.validateGauge(u_);
        for (int p = 0; p < lattice_.cells(); ++p)
            if (lattice_.plaquetteFlux(p,u_) != 1) throw std::invalid_argument("P1 model requires zero plaquette flux");
        bond_terms_.reserve(lattice_.bonds().size());
        for (const auto &b : lattice_.bonds()) {
            const double bare = bareCoupling(b.type);
            const double gauge = 2.0*double(u_[b.index]);
            if (!std::isfinite(gauge*(bare+couplings_.g*x_[b.index])))
                throw std::invalid_argument("nonfinite effective hopping");
            bond_terms_.push_back({b.index,b.a,b.b,gauge,bare});
        }
    }
    const HoneycombLattice &lattice() const { return lattice_; }
    const Couplings &couplings() const { return couplings_; }
    const std::vector<double> &displacements() const { return x_; }
    const std::vector<int> &gauge() const { return u_; }
    int majoranas() const { return lattice_.sites(); }
    double bondHopping(int id) const {
        const auto type = lattice_.bond(id).type;
        const double j = type == BondType::X ? couplings_.Jx : type == BondType::Y ? couplings_.Jy : couplings_.Jz;
        return j+couplings_.g*x_.at(id);
    }
    // H = (i/4) c^T A c; each oriented physical bond occurs once.
    Eigen::MatrixXd realAntisymmetricA() const {
        Eigen::MatrixXd a = Eigen::MatrixXd::Zero(majoranas(),majoranas());
        for (const auto &b : lattice_.bonds()) {
            const double value = 2*u_[b.index]*bondHopping(b.index);
            a(b.a,b.b) += value; a(b.b,b.a) -= value;
        }
        return a;
    }
    // Existing framework: H = -(1/4) c^T K_code c, hence K_code = -i A.
    MatType frameworkGenerator() const {
        MatType k = DataType(0,-1)*realAntisymmetricA().cast<DataType>();
        validateGenerator(k,majoranas());
        return k;
    }
    // A3: generator for a dynamic phonon configuration x without rebuilding
    // the model. Lattice connectivity, couplings and gauge are immutable for
    // the engine lifetime and their checks run once at construction; every
    // dynamic gate (x size/finiteness, hopping finiteness, dimension,
    // transpose-skew and Hermiticity of K) is retained per call, and the
    // assembly arithmetic is identical to construction + frameworkGenerator.
    MatType frameworkGeneratorWith(const std::vector<double> &x) const {
        MatType k;
        frameworkGeneratorWith(x,couplings_.g,k);
        return k;
    }
    // Protocol-aware assembly: the instantaneous SSH coupling is supplied
    // explicitly while lattice connectivity, bare couplings and gauge signs
    // remain those of this immutable model.
    MatType frameworkGeneratorWith(const std::vector<double> &x,double g) const {
        MatType k;
        frameworkGeneratorWith(x,g,k);
        return k;
    }
    // A3 output-buffer form. The lattice connectivity, gauge signs and bare
    // bond couplings are immutable in this model, so only the X-dependent
    // hopping values are recomputed for each proposal.
    void frameworkGeneratorWith(const std::vector<double> &x,MatType &k) const {
        frameworkGeneratorWith(x,couplings_.g,k);
    }
    void frameworkGeneratorWith(const std::vector<double> &x,double g,MatType &k) const {
        if (x.size() != std::size_t(lattice_.bondCount()))
            throw std::invalid_argument("dynamic X bond count");
        if (!std::isfinite(g)) throw std::invalid_argument("nonfinite dynamic g");
        for (double v : x) if (!std::isfinite(v)) throw std::invalid_argument("nonfinite dynamic X");
        k.setZero(majoranas(),majoranas());
        for (const auto &term : bond_terms_) {
            const double value = term.gauge*(term.bare+g*x[term.index]);
            if (!std::isfinite(value)) throw std::invalid_argument("nonfinite effective hopping");
            k(term.a,term.b) += DataType(0.,-value);
            k(term.b,term.a) += DataType(0.,value);
        }
        validateGenerator(k,majoranas());
    }
    static void validateGenerator(const MatType &k, int dimension) {
        if (dimension <= 0 || dimension%2 || k.rows()!=dimension || k.cols()!=dimension || !k.allFinite())
            throw std::invalid_argument("invalid Majorana generator dimension or entries");
        const double scale = std::max(1.,k.norm());
        if ((k+k.transpose()).norm()>1e-12*scale || (k-k.adjoint()).norm()>1e-12*scale)
            throw std::invalid_argument("Majorana generator must be transpose-skew and Hermitian");
    }
    MatType propagator(double dt) const {
        if (!std::isfinite(dt) || dt<=0) throw std::invalid_argument("dt must be positive and finite");
        MatType k = frameworkGenerator(); // expm requires Hermitian input.
        MatType b = expm(k,-dt);
        if (!b.allFinite()) throw std::runtime_error("nonfinite static propagator");
        return b;
    }
    KitaevSSHModel gaugeTransformed(const std::vector<int> &d) const {
        if (d.size()!=std::size_t(majoranas())) throw std::invalid_argument("site gauge count");
        for (int v : d) if (v!=1 && v!=-1) throw std::invalid_argument("site gauge must be +/-1");
        auto u = u_;
        for (const auto &b : lattice_.bonds()) u[b.index] *= d[b.a]*d[b.b];
        return KitaevSSHModel(lattice_,couplings_,x_,u);
    }

private:
    struct BondTerm {
        int index, a, b;
        double gauge, bare;
    };
    double bareCoupling(BondType type) const {
        return type == BondType::X ? couplings_.Jx :
               type == BondType::Y ? couplings_.Jy : couplings_.Jz;
    }
    HoneycombLattice lattice_;
    Couplings couplings_;
    std::vector<double> x_;
    std::vector<int> u_;
    std::vector<BondTerm> bond_terms_;
};

} // namespace kitaev_ssh
#endif
