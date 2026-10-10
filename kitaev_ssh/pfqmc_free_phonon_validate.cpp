// Independent finite-contour Gaussian reference and statistical production
// estimator checks. The zero-shift trial is the free oscillator ground state.
#include "phonon_mc.h"
#include "vbs_susceptibility.h"
#include <Eigen/Cholesky>
#include <iostream>
#include <iomanip>

using namespace kitaev_ssh;
namespace {
int failures=0;
void check(const std::string &name,double reference,double measured,double error,double tolerance) {
    const bool pass=std::isfinite(measured)&&std::isfinite(error)&&std::abs(measured-reference)<=tolerance;
    std::cout<<(pass?"PASS ":"FAIL ")<<name<<" reference="<<reference<<" measured="<<measured
             <<" abs_error="<<std::abs(measured-reference)<<" stderr="<<error<<" tolerance="<<tolerance<<'\n';
    failures+=!pass;
}
// S=1/2 X^T Q X, obtained from the exact oscillator propagator
// m*w/(2*sinh(w*dtau))*[(x*x+y*y)*cosh(w*dtau)-2*x*y]
// and the two normalized ground-state endpoint wavefunctions.
Eigen::MatrixXd covariance(int Ntau,double dtau,double mass,double omega_ph) {
    Eigen::MatrixXd Q=Eigen::MatrixXd::Zero(Ntau+1,Ntau+1);
    const double diagonal=mass*omega_ph/std::tanh(omega_ph*dtau);
    const double off=-mass*omega_ph/std::sinh(omega_ph*dtau);
    for(int j=0;j<Ntau;++j) {
        Q(j,j)+=diagonal; Q(j+1,j+1)+=diagonal;
        Q(j,j+1)+=off; Q(j+1,j)+=off;
    }
    Q(0,0)+=mass*omega_ph; Q(Ntau,Ntau)+=mass*omega_ph;
    Eigen::LDLT<Eigen::MatrixXd> factor(Q);
    if(factor.info()!=Eigen::Success || !factor.isPositive())throw std::runtime_error("Gaussian precision is not positive");
    return factor.solve(Eigen::MatrixXd::Identity(Ntau+1,Ntau+1));
}
struct Measurements {
    VBSObservables vbs;
    VBSSusceptibility susceptibility;
    double max_fourier_error=0.,max_response_error=0.;
    explicit Measurements(const HoneycombLattice &lattice):vbs(lattice),susceptibility(vbs){}
    std::vector<Complex> raw(const PhononPath &path,Complex phase) {
        const int nc=vbs.lattice().cells(),Ntau=path.intervals();
        const double dtau=path.dt(),omega_ph=path.parameters().omega;
        const int mid=Ntau/2;
        const auto values=vbs.measure(path.node(mid)),response=susceptibility.measure(path);
        std::vector<Complex> out(4+4*nc,0.);
        for(int b=0;b<path.bonds();++b) {
            const double x=path.coordinate(mid,b);
            out[0]+=x/path.bonds();out[1]+=x*x/path.bonds();
            out[2]+=x*path.coordinate(0,b)/path.bonds();
        }
        // Direct spatial sums independent of FFT and susceptibility helpers.
        for(int q=0;q<nc;++q) {
            Complex d=0.,integral=0.;
            const auto momentum=vbs.momenta()[q];
            for(int r=0;r<nc;++r) {
                const auto coordinate=vbs.lattice().cellCoordinate(r);
                const auto e=std::polar(1.,-2*pi*(momentum.reciprocal_x*coordinate.x+momentum.reciprocal_y*coordinate.y));
                d+=path.coordinate(mid,3*r+2)*e/double(nc);
                for(int j=0;j<=Ntau;++j) {
                    const double score=(j==0||j==Ntau?1.:2.)*std::tanh(omega_ph*dtau/2)/omega_ph;
                    integral+=score*path.coordinate(j,3*r+2)*e/double(nc);
                }
            }
            max_fourier_error=std::max(max_fourier_error,std::abs(values[9*q+8]-double(nc)*d*std::conj(d)));
            max_response_error=std::max(max_response_error,std::abs(response[24*q+17]-double(nc)*d*std::conj(integral)));
            out[3+4*q]=values[9*q+8];out[4+4*q]=response[24*q+17];
            out[5+4*q]=response[24*q+2];out[6+4*q]=response[24*q+5];
        }
        for(auto &value:out)value*=phase;
        out.back()=phase;
        return out;
    }
};
std::vector<double> derive(const std::vector<Complex> &raw,int nc) {
    std::vector<double> out(3+2*nc);
    const auto phase=raw.back();
    if(std::abs(phase)<.5)throw std::runtime_error("unexpected free-limit phase denominator");
    for(int j=0;j<3;++j)out[j]=(raw[j]/phase).real();
    for(int q=0;q<nc;++q) {
        out[3+2*q]=(raw[3+4*q]/phase-double(nc)*raw[5+4*q]/phase*std::conj(raw[5+4*q]/phase)).real();
        out[4+2*q]=(raw[4+4*q]/phase-double(nc)*raw[5+4*q]/phase*std::conj(raw[6+4*q]/phase)).real();
    }
    return out;
}
void report(const std::string &label,const std::vector<std::vector<Complex>> &blocks,int nc,const Eigen::MatrixXd &C,double dtau,double omega_ph) {
    std::vector<Complex> mean(blocks[0].size(),0.);
    for(const auto &block:blocks)for(std::size_t j=0;j<mean.size();++j)mean[j]+=block[j]/double(blocks.size());
    const auto estimate=derive(mean,nc);
    std::vector<std::vector<double>> replicas;
    for(const auto &block:blocks) {
        auto deleted=mean;for(std::size_t j=0;j<mean.size();++j)deleted[j]=(mean[j]*double(blocks.size())-block[j])/double(blocks.size()-1);
        replicas.push_back(derive(deleted,nc));
    }
    std::vector<double> average(estimate.size(),0.),errors(estimate.size(),0.);
    for(const auto &replica:replicas)for(std::size_t j=0;j<average.size();++j)average[j]+=replica[j]/double(replicas.size());
    for(const auto &replica:replicas)for(std::size_t j=0;j<average.size();++j)errors[j]+=(replica[j]-average[j])*(replica[j]-average[j]);
    for(auto &error:errors)error=std::sqrt(error*double(blocks.size()-1)/double(blocks.size()));
    const int Ntau=C.rows()-1,mid=Ntau/2;
    double chi=0.;for(int j=0;j<=Ntau;++j)chi+=C(mid,j)*(j==0||j==Ntau?1.:2.)*std::tanh(omega_ph*dtau/2)/omega_ph;
    for(std::size_t j=0;j<estimate.size();++j) {
        const double reference=j==0?0.:j==1?C(mid,mid):j==2?C(mid,0):(j-3)%2==0?C(mid,mid):chi;
        const std::string name=j==0?"X":j==1?"X2":j==2?"X(tau_mid)X(0)":(j-3)%2==0?"Szz_q"+std::to_string((j-3)/2):"chi_zz_q"+std::to_string((j-3)/2);
        // Fixed six-sigma statistical gate; connected subtraction retains its
        // ordinary finite-sample bias. No tolerance is adapted after failures.
        check(label+"/"+name,reference,estimate[j],errors[j],6*errors[j]+2e-10);
    }
}
void run(int L,double omega_ph,int Ntau,double dtau,bool mc,PhononProposal mode=PhononProposal::LocalConditional) {
    const double mass=1.7,tau=Ntau*dtau,tau_mid=tau/2,k=mass*omega_ph*omega_ph;
    HoneycombLattice lattice(L,L);PhononPath path(lattice.bondCount(),Ntau,dtau,{mass,omega_ph});
    const auto C=covariance(Ntau,dtau,mass,omega_ph);
    const std::string label=std::string(mc?(mode==PhononProposal::LocalConditional?"MC_local":"MC_bond"):"independent_draw")+"/L"+std::to_string(L)+"/w"+std::to_string(omega_ph)+"/Ntau"+std::to_string(Ntau)+"/dtau"+std::to_string(dtau);
    std::cout<<"CASE "<<label<<" tau="<<tau<<" tau_mid="<<tau_mid<<" blocks=40 samples_per_block="<<(mc?20:100)<<'\n';
    double maxCov=0.;for(int i=0;i<=Ntau;++i)for(int j=0;j<=Ntau;++j)maxCov=std::max(maxCov,std::abs(C(i,j)-omega_ph/(2*k)*std::exp(-omega_ph*dtau*std::abs(i-j))));
    check(label+"/Gaussian_vs_ground_state_covariance",0.,maxCov,0.,2e-12);
    double exact=0.,trapezoid=0.;for(int j=0;j<=Ntau;++j){exact+=C(Ntau/2,j)*(j==0||j==Ntau?1.:2.)*std::tanh(omega_ph*dtau/2)/omega_ph;trapezoid+=C(Ntau/2,j)*dtau*(j==0||j==Ntau?.5:1.);}
    check(label+"/finite_tau_score_response",(1-std::exp(-omega_ph*tau_mid))/k,exact,0.,2e-12);
    std::cout<<"INFO "<<label<<" chi_ground_state="<<1/k<<" chi_finite_tau="<<exact<<" chi_trapezoid="<<trapezoid<<" trapezoid_bias="<<trapezoid-exact<<'\n';
    // Check the actual action and single-node conditionals against the same
    // independently inverted precision; this catches endpoint double counting.
    Eigen::VectorXd coordinates(Ntau+1);
    std::vector<double> bond(Ntau+1);
    for(int j=0;j<=Ntau;++j)bond[j]=coordinates[j]=.17*std::sin(.71*j)+.03*j;
    const auto probe=path.withBondCoordinates(0,bond);
    const Eigen::MatrixXd Q=C.ldlt().solve(Eigen::MatrixXd::Identity(Ntau+1,Ntau+1));
    check(label+"/production_action_vs_precision",.5*coordinates.dot(Q*coordinates),probe.action(),0.,2e-11);
    double maxConditional=0.,maxRatio=0.,maxHessian=0.;
    for(int j=0;j<=Ntau;++j) {
        const auto conditional=freeConditional(probe,j,0);
        const double expectedMean=coordinates[j]-(Q*coordinates)[j]/Q(j,j);
        maxConditional=std::max(maxConditional,std::abs(conditional.mean-expectedMean));
        maxConditional=std::max(maxConditional,std::abs(conditional.variance-1/Q(j,j)));
        const double delta=.23;
        const double expectedRatio=-delta*(Q*coordinates)[j]-.5*delta*delta*Q(j,j);
        maxRatio=std::max(maxRatio,std::abs(probe.localLogRatio(j,0,coordinates[j]+delta)-expectedRatio));
        const double hessian=probe.withCoordinate(j,0,coordinates[j]+1).action()+probe.withCoordinate(j,0,coordinates[j]-1).action()-2*probe.action();
        maxHessian=std::max(maxHessian,std::abs(hessian-Q(j,j)));
    }
    check(label+"/production_action_Hessian",0.,maxHessian,0.,2e-11);
    check(label+"/production_free_conditional",0.,maxConditional,0.,2e-11);
    check(label+"/production_local_action_ratio",0.,maxRatio,0.,2e-11);
    std::mt19937_64 rng(20261008+L*100+Ntau+(mc?7000:0)+(mode==PhononProposal::BondRefresh?9000:0));
    std::unique_ptr<CoupledPhononPath> state;std::unique_ptr<PhononSampler> sampler;
    if(mc) {
        const KitaevSSHModel model(lattice,p2Couplings(0.));
        const auto generator=honeycombGenerator(model);
        const auto trial=freeGroundState(model);
        // Match the production magnitude walk; certify the signed Pfaffian
        // at every measurement instead of doing unused signed work per move.
        PhononWeightOptions options;
        options.singular_values_only=true;options.merge_node_half_steps=true;
        options.defer_phase_measurement=true;
        options.reuse_signed_prefix_checkpoints=true;options.reuse_lu_cayley=true;
        options.compute_measurement_green=false;
        state.reset(new CoupledPhononPath(path.freeSample(rng),trial,generator,4,options));
        sampler.reset(new PhononSampler(*state,rng,mode));for(int j=0;j<100;++j)sampler->sweep();
    }
    Measurements measurements(lattice);std::vector<std::vector<Complex>> blocks;
    const KitaevSSHModel measurementModel(lattice,p2Couplings(0.));
    double maxBasicError=0.;
    const int samples=mc?20:100;
    for(int block=0;block<40;++block) {
        std::vector<Complex> sum(4+4*lattice.cells(),0.);
        for(int sample=0;sample<samples;++sample) {
            Complex phase=1.;if(mc){for(int j=0;j<4;++j)sampler->sweep();sampler->ensureMeasured();path=sampler->state().path;phase=sampler->state().weight.matter_phase;}else path=path.freeSample(rng);
            if(mc) {
                const auto basic=basicObservables(measurementModel,sampler->state(),false);
                for(int mu=0;mu<3;++mu) {
                    double mean=0.,square=0.;for(int r=0;r<lattice.cells();++r) {
                        const double x=path.coordinate(path.center(),3*r+mu);
                        mean+=x/double(lattice.cells());square+=x*x/double(lattice.cells());
                    }
                    maxBasicError=std::max(maxBasicError,std::abs(basic[mu]-mean));
                    maxBasicError=std::max(maxBasicError,std::abs(basic[3+mu]-square));
                }
            }
            const auto values=measurements.raw(path,phase);for(std::size_t j=0;j<sum.size();++j)sum[j]+=values[j]/double(samples);
        }
        blocks.push_back(std::move(sum));
    }
    report(label,blocks,lattice.cells(),C,dtau,omega_ph);
    check(label+"/independent_Fourier_Szz",0.,measurements.max_fourier_error,0.,2e-11);
    check(label+"/independent_Fourier_response",0.,measurements.max_response_error,0.,2e-11);
    if(mc)check(label+"/production_basicObservables_X_X2",0.,maxBasicError,0.,2e-12);
    if(mc){const auto &counts=sampler->counts();check(label+"/all_free_proposals_accepted",double(counts.proposed),double(counts.accepted),0.,0.);check(label+"/numerical_failures",0.,double(counts.first_failures),0.,0.);}
}
}
int main(int argc,char **argv) {
    try {
        std::cout<<std::setprecision(12);
        // --direct-only avoids the actual sampler; default also exercises both
        // production proposals with the L=2 Kitaev Hamiltonian at g=0.
        const bool directOnly=argc==2&&std::string(argv[1])=="--direct-only";
        if(argc==3&&std::string(argv[1])=="--mc-case") {
            const std::string argument=argv[2];std::size_t used=0;
            const int index=std::stoi(argument,&used);
            if(used!=argument.size()||index<0||index>5)throw std::invalid_argument("MC case index must be 0..5");
            run(2+2*(index/2),1.,4,.5,true,index%2?PhononProposal::BondRefresh:PhononProposal::LocalConditional);
            std::cout<<(failures?"FAIL":"PASS")<<" free_phonon_suite failures="<<failures<<'\n';return failures?1:0;
        }
        if(argc>1&&!directOnly)throw std::invalid_argument("usage: pfqmc_free_phonon_validate [--direct-only | --mc-case 0..5]");
        for(int L:{2,4,6}){run(L,.5,8,.25,false);run(L,1.,4,.5,false);run(L,2.,16,.25,false);}
        run(2,1.,8,.25,false);run(2,1.,16,.5,false);
        if(!directOnly){run(2,1.,4,.5,true);run(2,1.,4,.5,true,PhononProposal::BondRefresh);}
        std::cout<<(failures?"FAIL":"PASS")<<" free_phonon_suite failures="<<failures<<'\n';return failures?1:0;
    }catch(const std::exception &error){std::cerr<<"FAIL free_phonon_suite exception="<<error.what()<<'\n';return 1;}
}
