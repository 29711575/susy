#ifndef KITAEV_SSH_PHONON_PATH_H
#define KITAEV_SSH_PHONON_PATH_H

#include <cmath>
#include <iomanip>
#include <istream>
#include <limits>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

namespace kitaev_ssh {

// Unscaled displacement. g, mass and omega are independent; k/lambda derived.
struct PhononParameters {
    double mass, omega;
    PhononParameters(double m, double w) : mass(m), omega(w) {
        if (!std::isfinite(m) || !std::isfinite(w) || m<=0 || w<=0 ||
            !std::isfinite(m*w) || m*w<=0 || !std::isfinite(m*w*w) || m*w*w<=0)
            throw std::invalid_argument("phonon mass, omega and derived scales must be finite and positive");
    }
    double spring() const { return mass*omega*omega; }
    double lambda(double g) const {
        const double result=g*g/spring();
        if (!std::isfinite(g) || !std::isfinite(result)) throw std::invalid_argument("nonfinite g or lambda");
        return result;
    }
};

class PhononPath {
public:
    // M even intervals, M+1 independent integration coordinates (one center).
    PhononPath(int bonds, int intervals, double h, PhononParameters p, std::vector<double> shifts={})
        : p_(p.mass,p.omega), h_(h), x_(checkedNodes(bonds,intervals,h,p),std::vector<double>(bonds,0.)), d_(std::move(shifts)),
          z_(p.omega*h), decay_(std::exp(-z_)), noise_factor_(-std::expm1(-2*z_)), alpha_(p.mass*p.omega), tanh_half_(std::tanh(z_/2)) {
        if(d_.empty()) d_.assign(bonds,0.);
        if(d_.size()!=std::size_t(bonds)) throw std::invalid_argument("trial displacement bond count");
        for(double d:d_) if(!std::isfinite(d) || !std::isfinite(p_.mass*p_.omega*d*d))
            throw std::invalid_argument("invalid trial displacement");
    }
    int bonds() const { return int(x_[0].size()); }
    int intervals() const { return int(x_.size())-1; }
    int nodes() const { return int(x_.size()); }
    int center() const { return intervals()/2; }
    double dt() const { return h_; }
    double totalTau() const { return intervals()*h_; }
    double decay() const { return decay_; }
    double noiseFactor() const { return noise_factor_; }
    double alpha() const { return alpha_; }
    double tanhHalf() const { return tanh_half_; }
    const PhononParameters &parameters() const { return p_; }
    double trialShift(int b) const { return d_.at(b); }
    const std::vector<double> &trialShifts() const { return d_; }
    double freeMean(int j,int b) const {
        x_.at(j).at(b);
        return d_.at(b)/2*(std::exp(-p_.omega*h_*j)+std::exp(-p_.omega*h_*(intervals()-j)));
    }
    // Integrated normalized-endpoint free path norm, not kernel prefactors.
    double freeLogNorm() const {
        double value=-bonds()*p_.omega*totalTau()/2;
        for(double d:d_) value-=p_.mass*p_.omega*d*d*(-std::expm1(-p_.omega*totalTau()))/2;
        return value;
    }
    const std::vector<double> &node(int j) const { return x_.at(j); }
    double coordinate(int j, int b) const { return x_.at(j).at(b); }
    PhononPath withCoordinate(int j, int b, double value) const {
        coordinate(j,b);
        if (!std::isfinite(value)) throw std::invalid_argument("nonfinite phonon coordinate");
        PhononPath changed=*this; changed.x_.at(j).at(b)=value; return changed;
    }
    PhononPath withBondCoordinates(int bond,const std::vector<double> &values) const {
        coordinate(0,bond);
        if(values.size()!=x_.size())throw std::invalid_argument("bond path node count");
        for(double value:values)if(!std::isfinite(value))throw std::invalid_argument("nonfinite bond path coordinate");
        PhononPath changed=*this;
        for(int j=0;j<nodes();++j)changed.x_[j][bond]=values[j];
        return changed;
    }

    double linkAction(double x, double y) const {
        const double t=tanh_half_, alpha=alpha_;
        return alpha/4*((x-y)*(x-y)/t+t*(x+y)*(x+y));
    }
    double action() const {
        double result=0.;
        for (int b=0;b<bonds();++b) {
            result+=p_.mass*p_.omega/2*((x_[0][b]-d_[b])*(x_[0][b]-d_[b])+(x_.back()[b]-d_[b])*(x_.back()[b]-d_[b]));
            for (int j=1;j<nodes();++j) result+=linkAction(x_[j][b],x_[j-1][b]);
        }
        if (!std::isfinite(result)) throw std::overflow_error("nonfinite phonon action");
        return result;
    }
    double logNormalization() const {
        const double z=p_.omega*h_;
        // log(sinh z), no overflow for large z or cancellation for small z.
        const double logSinh=z<1 ? std::log(std::sinh(z)) : z+std::log1p(-std::exp(-2*z))-std::log(2.);
        const double logAlpha=std::log(p_.mass)+std::log(p_.omega), logPi=std::log(std::acos(-1.));
        return bonds()/2.*(logAlpha-logPi)+double(intervals())*bonds()/2.*(logAlpha-std::log(2.)-logPi-logSinh);
    }
    double logWeight() const {
        const double result=logNormalization()-action();
        if (!std::isfinite(result)) throw std::overflow_error("nonfinite phonon log weight");
        return result;
    }
    // log Wph(new) - log Wph(old); only neighbors and an endpoint if present.
    double localLogRatio(int j, int b, double value) const {
        const double old=coordinate(j,b), delta=value-old;
        if (!std::isfinite(value) || !std::isfinite(delta)) throw std::invalid_argument("nonfinite proposal");
        const double alpha=alpha_, t=tanh_half_;
        auto difference=[&](double y) {
            return alpha/4*(delta*(2*(old-y)+delta)/t+t*delta*(2*(old+y)+delta));
        };
        double ds=0;
        if (j>0) ds+=difference(x_[j-1][b]);
        if (j<intervals()) ds+=difference(x_[j+1][b]);
        if (j==0 || j==intervals()) ds+=alpha/2*delta*(2*(old-d_[b])+delta);
        if (!std::isfinite(ds)) throw std::overflow_error("nonfinite local phonon action change");
        return -ds;
    }

    // Exact normalized FREE open-path law, NOT a coupled equilibrium sampler.
    template<class Rng> PhononPath freeSample(Rng &rng) const {
        PhononPath result=*this;
        const double z=p_.omega*h_, r=std::exp(-z), sigma=std::sqrt(.5/(p_.mass*p_.omega));
        const double noise=sigma*std::sqrt(-std::expm1(-2*z));
        std::normal_distribution<double> normal;
        for (int b=0;b<bonds();++b) {
            result.x_[0][b]=sigma*normal(rng);
            for (int j=1;j<nodes();++j) result.x_[j][b]=r*result.x_[j-1][b]+noise*normal(rng);
            if(d_[b]!=0.) for(int j=0;j<nodes();++j) result.x_[j][b]+=freeMean(j,b);
        }
        return result;
    }
    // Text snapshot of coordinates and phonon parameters, not an RNG/MC checkpoint.
    void save(std::ostream &out) const {
        // Equivalent to defaultfloat, including clearing both floatfield bits;
        // supported by the older libstdc++ on the HPC Intel toolchain.
        out.unsetf(std::ios_base::floatfield);
        out<<std::setprecision(std::numeric_limits<double>::max_digits10)
           <<"KITAEV_SSH_PHONON_PATH 2\n"<<bonds()<<' '<<intervals()<<' '<<h_<<' '<<p_.mass<<' '<<p_.omega<<'\n';
        for(double d:d_) out<<d<<' '; out<<'\n';
        for (const auto &row:x_) { for (double value:row) out<<value<<' '; out<<'\n'; }
        if (!out) throw std::runtime_error("phonon path save failed");
    }
    static PhononPath load(std::istream &in) {
        std::string magic; int version,bonds,intervals; double h,m,w;
        if (!(in>>magic>>version>>bonds>>intervals>>h>>m>>w) || magic!="KITAEV_SSH_PHONON_PATH" || (version!=1 && version!=2))
            throw std::invalid_argument("invalid phonon path header");
        PhononPath result(bonds,intervals,h,{m,w});
        if(version==2) {
            std::vector<double> shifts(bonds);
            for(auto &d:shifts) if(!(in>>d)) throw std::invalid_argument("truncated trial shifts");
            result=PhononPath(bonds,intervals,h,{m,w},shifts);
        }
        for (auto &row:result.x_) for (double &value:row)
            if (!(in>>value) || !std::isfinite(value)) throw std::invalid_argument("invalid/truncated phonon coordinates");
        return result;
    }
private:
    PhononParameters p_;
    double h_;
    double z_, decay_, noise_factor_, alpha_, tanh_half_;
    std::vector<std::vector<double>> x_;
    std::vector<double> d_;
    static int checkedNodes(int bonds,int intervals,double h,PhononParameters p) {
        if (bonds<=0 || intervals<2 || intervals%2 || intervals>std::numeric_limits<int>::max()/2 ||
            !std::isfinite(h) || h<=0 || !std::isfinite(h*p.omega) || h*p.omega<=0 ||
            !std::isfinite(h*intervals) || std::tanh(h*p.omega/2)==0)
            throw std::invalid_argument("open contour needs bonds>0, even M>=2 and finite positive time scales");
        return intervals+1;
    }
};

} // namespace kitaev_ssh
#endif
