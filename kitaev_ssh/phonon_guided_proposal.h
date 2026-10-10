#ifndef KITAEV_SSH_PHONON_GUIDED_PROPOSAL_H
#define KITAEV_SSH_PHONON_GUIDED_PROPOSAL_H

#include "phonon_path.h"

namespace kitaev_ssh {
// A FIXED proposal reference, never a term in the sampled Hamiltonian or trial.
// The reference Gaussian is proportional to W_ph(X) exp(sum_b f_b I_b).
// Freeze forces and beta before production. A poor reference costs efficiency,
// not correctness, provided the reference-source correction below is retained.
struct GuidedBondOptions {
    double beta=1.;
    std::vector<double> forces;
};

class GuidedBondProposal {
public:
    GuidedBondProposal(const PhononPath &path,GuidedBondOptions options)
        : options_(std::move(options)),rho_(std::sqrt((1-options_.beta)*(1+options_.beta))),
          decay_(path.decay()),noise_factor_(path.noiseFactor()),variance_(1/(2*path.alpha())) {
        if(!std::isfinite(options_.beta) || options_.beta<=0 || options_.beta>1)
            throw std::invalid_argument("guided bond beta must be in (0,1]");
        if(options_.forces.empty())options_.forces.assign(path.bonds(),0.);
        if(options_.forces.size()!=std::size_t(path.bonds()))
            throw std::invalid_argument("guided bond reference force count");
        mean_.assign(path.bonds(),std::vector<double>(path.nodes()));
        weights_.resize(path.nodes());
        for(int j=0;j<path.nodes();++j) {
            weights_[j]=path.sourceNodeWeight(j);
            const double response=(1.-.5*(std::exp(-path.parameters().omega*path.dt()*j)
                +std::exp(-path.parameters().omega*path.dt()*(path.intervals()-j))))/path.parameters().spring();
            for(int b=0;b<path.bonds();++b) {
                if(!std::isfinite(options_.forces[b]))throw std::invalid_argument("nonfinite guided bond force");
                mean_[b][j]=path.freeMean(j,b)+options_.forces[b]*response;
                if(!std::isfinite(mean_[b][j]))throw std::invalid_argument("nonfinite guided bond reference mean");
            }
        }
    }
    const GuidedBondOptions &options()const{return options_;}
    const std::vector<double> &mean(int b)const{return mean_.at(b);}
    template<class Rng> std::vector<double> draw(const PhononPath &path,int b,Rng &rng)const {
        check(path,b);
        std::normal_distribution<double> normal;
        const double sigma=std::sqrt(variance_),noise=sigma*std::sqrt(noise_factor_);
        std::vector<double> out(path.nodes());
        double xi=sigma*normal(rng);
        for(int j=0;j<path.nodes();++j) {
            if(j)xi=decay_*xi+noise*normal(rng);
            out[j]=mean_[b][j]+rho_*(path.coordinate(j,b)-mean_[b][j])+options_.beta*xi;
        }
        return out;
    }
    // log r = log |W_f(new)/W_f(old)| - sum_j f_b w_j (X'_jb-X_jb).
    double referenceScoreDifference(const PhononPath &path,int b,const std::vector<double> &to)const {
        check(path,b,to.size());double result=0.;
        for(int j=0;j<path.nodes();++j)result+=weights_[j]*(to[j]-path.coordinate(j,b));
        result*=options_.forces[b];
        if(!std::isfinite(result))throw std::runtime_error("nonfinite guided reference score");
        return result;
    }
    // Independent normalized AR Gaussian density for the full MH cross-check.
    double logTransition(const PhononPath &from,int b,const std::vector<double> &to)const {
        check(from,b,to.size());
        const double var=options_.beta*options_.beta*variance_,noise=var*noise_factor_;
        auto ln=[](double x,double v){return -.5*(std::log(2*std::acos(-1.)*v)+x*x/v);};
        double previous=to[0]-mean_[b][0]-rho_*(from.coordinate(0,b)-mean_[b][0]);
        double result=ln(previous,var);
        for(int j=1;j<from.nodes();++j) {
            const double y=to[j]-mean_[b][j]-rho_*(from.coordinate(j,b)-mean_[b][j]);
            result+=ln(y-decay_*previous,noise);previous=y;
        }
        if(!std::isfinite(result))throw std::runtime_error("nonfinite guided proposal density");
        return result;
    }
private:
    void check(const PhononPath &path,int b,std::size_t n=std::numeric_limits<std::size_t>::max())const {
        if(b<0 || b>=int(mean_.size()) || path.bonds()!=int(mean_.size()) ||
           path.nodes()!=int(weights_.size()) || (n!=std::numeric_limits<std::size_t>::max() && n!=weights_.size()) ||
           path.decay()!=decay_ || path.noiseFactor()!=noise_factor_ || 1/(2*path.alpha())!=variance_)
            throw std::invalid_argument("guided proposal geometry changed");
    }
    GuidedBondOptions options_;
    double rho_,decay_,noise_factor_,variance_;
    std::vector<std::vector<double>> mean_;
    std::vector<double> weights_;
};
} // namespace kitaev_ssh
#endif
