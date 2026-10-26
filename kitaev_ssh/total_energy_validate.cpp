// Fixed inputs, analytic Gaussian integration and independent oscillator/Fock
// matrices only. No Monte Carlo sweeps, proposals or chain convergence gates.
#define KITAEV_SSH_TEST_FAILURE_INJECTION
#define main original_fock_validation_main
#include "pfqmc_fock_validate.cpp"
#undef main
#include "phonon_mc.h"
#include "physical_matter_trial.h"

namespace {
void linkSpectral() {
    // Independent oscillator spectral expansion. True projected X^2/P^2
    // ladder elements include intermediate states beyond the retained basis.
    for(double mass:{.7,100./9})for(double omega:{.3,1.})for(double dt:{.05,.2,1.7})
    for(double force:{0.,.17,-.23})for(auto xy:std::vector<std::pair<double,double>>{{.2,-.3},{.7,.7},{-.4,.6}}) {
        // A spectral sum at endpoints many diffusion lengths apart is an
        // exponentially small difference of large alternating terms. Keep
        // this double-precision oracle in its resolved domain; the production
        // formula itself has no such endpoint-distance restriction.
        xy.second=xy.first+(xy.second-xy.first)*std::sqrt(dt/mass);
        const int n=4096;const double a=mass*omega,q=force/(mass*omega*omega);
        auto wave=[&](double x) {
            Eigen::VectorXd v(n);const double z=std::sqrt(a)*(x-q);
            v[0]=std::pow(a/pi,.25)*std::exp(-z*z/2);v[1]=std::sqrt(2.)*z*v[0];
            for(int j=1;j+1<n;++j)v[j+1]=std::sqrt(2./(j+1))*z*v[j]-std::sqrt(double(j)/(j+1))*v[j-1];
            for(int j=0;j<n;++j)v[j]*=std::exp(-dt/2*(omega*(j+.5)-force*force/(2*mass*omega*omega)));
            return v;
        };
        const auto left=wave(xy.first),right=wave(xy.second);const double weight=left.dot(right);
        double first=q*weight,second=0.,kinetic=0.;
        for(int j=0;j<n;++j) {
            const double z=left[j]*right[j];second+=(q*q+(j+.5)/a)*z;kinetic+=omega*(2*j+1)/4*z;
            if(j+1<n) {
                const double adjacent=left[j]*right[j+1]+left[j+1]*right[j];
                first+=std::sqrt((j+1)/(2*a))*adjacent;
                second+=2*q*std::sqrt((j+1)/(2*a))*adjacent;
            }
            if(j+2<n) {
                const double off=std::sqrt(double((j+1)*(j+2)))*(left[j]*right[j+2]+left[j+2]*right[j]);
                second+=off/(2*a);kinetic-=omega*off/4;
            }
        }
        const auto measured=phononLinkEnergy(xy.first,xy.second,{mass,omega},dt,force);
        check("Mehler_midpoint_X_vs_spectral",std::abs(measured.mean-first/weight),2e-10);
        check("Mehler_midpoint_X2_vs_spectral",std::abs(measured.second-second/weight),2e-10);
        check("Mehler_midpoint_kinetic_vs_P2",std::abs(measured.kinetic-kinetic/weight),2e-9);
        check("Mehler_midpoint_source_vs_spectral",std::abs(measured.source+force*first/weight),2e-10);
    }
    for(double dt:{1e-8,1e-5,.05,.5,50.,1000.}) {
        const PhononParameters parameters{1.3,.7};
        const long double m=parameters.mass,w=parameters.omega,sigma=std::sqrt(1/(2*m*w)),rho=std::exp(-w*dt),noise=std::sqrt(-std::expm1(-2*w*dt));
        double kinetic=0.,potential=0.;
        // Four deterministic points reproduce the exact Gaussian mean and
        // covariance, hence integrate every quadratic link estimator exactly.
        for(double u:{-1.,1.})for(double v:{-1.,1.}) {
            auto measured=phononLinkEnergy(sigma*u,sigma*(rho*u+noise*v),parameters,dt);
            kinetic+=measured.kinetic/4;potential+=measured.potential/4;
        }
        check(dt<1e-6?"g0_kinetic_zero_point_extreme_small_dt":"g0_kinetic_zero_point",std::abs(kinetic-w/4),dt<1e-6?2e-7:2e-10);
        check("g0_potential_zero_point",std::abs(potential-w/4),2e-12);
    }
}

struct OscillatorED {
    Eigen::MatrixXd x,x2,kinetic,halfPhonon,halfMatter,step;
    Eigen::VectorXd trial;
    OscillatorED(int n,double mass,double omega,double dt,double force,double displacement,double coupling) {
        x=Eigen::MatrixXd::Zero(n,n);x2=x;kinetic=x;
        for(int j=0;j<n;++j) {
            x2(j,j)=(j+.5)/(mass*omega);kinetic(j,j)=omega*(2*j+1)/4;
            if(j+1<n)x(j,j+1)=x(j+1,j)=std::sqrt((j+1)/(2*mass*omega));
            if(j+2<n){x2(j,j+2)=x2(j+2,j)=std::sqrt(double((j+1)*(j+2)))/(2*mass*omega);kinetic(j,j+2)=kinetic(j+2,j)=-omega*std::sqrt(double((j+1)*(j+2)))/4;}
        }
        auto expm=[](const Eigen::MatrixXd &h,double time) {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> e(h);
            if(e.info()!=Eigen::Success)throw std::runtime_error("oscillator ED eigensolve");
            return Eigen::MatrixXd(e.eigenvectors()*(-time*e.eigenvalues().array()).exp().matrix().asDiagonal()*e.eigenvectors().transpose());
        };
        const Eigen::MatrixXd hph=kinetic+mass*omega*omega/2*x2-force*x;
        halfPhonon=expm(hph,dt/2);halfMatter=expm(coupling*x,dt/2);
        step=halfMatter*halfPhonon*halfPhonon*halfMatter;
        trial.resize(n);const double alpha=std::sqrt(mass*omega/2)*displacement;
        trial[0]=std::exp(-alpha*alpha/2);for(int j=1;j<n;++j)trial[j]=trial[j-1]*alpha/std::sqrt(double(j));
    }
    Eigen::Vector4d measure(int intervals,double dt)const {
        Eigen::Vector4d out=Eigen::Vector4d::Zero();
        for(int link:{intervals/2,intervals/2+1}) {
            auto right=trial,left=trial;
            for(int j=0;j<link-1;++j)right=step*right;
            for(int j=0;j<intervals-link;++j)left=step*left;
            right=halfPhonon*halfMatter*right;left=halfPhonon*halfMatter*left;
            const double z=left.dot(right);
            out+=Eigen::Vector4d(left.dot(x*right),left.dot(x2*right),left.dot(kinetic*right),z)/(2*z);
        }
        return out;
    }
};

void coupledQuantumED() {
    // A two-Majorana matter bilinear and THREE quantum oscillators.
    // In each fixed matter sector this is a forced-oscillator problem, but
    // [H_ph, H_matter(X)] != 0: the exact finite-dt split is tested explicitly.
    for(int M:{2,6})for(double dt:{.05,.2})for(double g:{0.,.6})for(double sector:{-1.,1.}) {
        const double mass=1.3,w=.7,k=mass*w*w;
        const std::vector<double> force={0.,.13,-.08},shift={.2,-.1,.3};
        PhononPath base(3,M,dt,{mass,w},shift,force);
        Eigen::MatrixXd precision=Eigen::MatrixXd::Zero(M+1,M+1);
        const double z=w*dt;
        for(int j=0;j<M;++j){precision(j,j)+=mass*w/std::tanh(z);precision(j+1,j+1)+=mass*w/std::tanh(z);precision(j,j+1)=precision(j+1,j)=-mass*w/std::sinh(z);}
        precision(0,0)+=mass*w;precision(M,M)+=mass*w;
        const Eigen::MatrixXd covariance=precision.inverse();
        const Eigen::MatrixXd root=covariance.llt().matrixL();
        std::vector<Eigen::VectorXd> means;
        for(int b=0;b<3;++b) {
            Eigen::VectorXd source(M+1);
            for(int j=0;j<=M;++j)source[j]=force[b]*std::tanh(z/2)/w*(j==0||j==M?1:2)-sector*g*dt*(j==0||j==M?.5:1.);
            source[0]+=mass*w*shift[b];source[M]+=mass*w*shift[b];means.push_back(precision.ldlt().solve(source));
        }
        std::vector<Complex> integrated(6,0.);
        const int dimension=3*(M+1);
        for(int axis=0;axis<dimension;++axis)for(double sign:{-1.,1.}) {
            auto path=base;
            for(int b=0;b<3;++b)for(int j=0;j<=M;++j) {
                const double offset=b==axis/(M+1)?sign*std::sqrt(double(dimension))*root(j,axis%(M+1)):0.;
                path=path.withCoordinate(j,b,means[b][j]+offset);
            }
            // In a fixed eigenstate i*c0*c1=sector, the exact matter
            // insertions are scalars. Use the production phonon insertion,
            // with normalization/coupling independently assembled here.
            std::vector<Complex> values(6,0.);values[0]=2.5*sector;
            for(int link:{M/2,M/2+1})for(int b=0;b<3;++b) {
                const auto ph=phononLinkEnergy(path.coordinate(link-1,b),path.coordinate(link,b),path.parameters(),dt,force[b]);
                values[1]+=.5*sector*g*ph.mean;values[2]+=.5*ph.kinetic;values[3]+=.5*ph.potential;values[4]+=.5*ph.source;
            }
            for(int j=0;j<5;++j)values[5]+=values[j];
            for(int j=0;j<6;++j)integrated[j]+=values[j]/double(2*dimension);
        }
        std::vector<double> exact(6,0.);exact[0]=2.5*sector;
        for(int b=0;b<3;++b) {
            const auto result=OscillatorED(48,mass,w,dt,force[b],shift[b],sector*g).measure(M,dt);
            const auto larger=OscillatorED(64,mass,w,dt,force[b],shift[b],sector*g).measure(M,dt);
            check("coupled_quantum_ED_cutoff",(result-larger).norm(),2e-10);
            exact[1]+=sector*g*result[0];exact[2]+=result[2];exact[3]+=k/2*result[1];exact[4]-=force[b]*result[0];
        }
        for(int j=0;j<5;++j)exact[5]+=exact[j];
        for(int j=0;j<6;++j)check("coupled_quantum_ED_component_"+std::to_string(j),std::abs(integrated[j]-exact[j]),3e-10);
    }
}

void arbitraryCutFockAndRecovery() {
    const auto c=operators(4);
    for(const std::string policy:{"physical","unrestricted"})for(double g:{0.,.4,.9})for(bool production:{false,true}) {
        KitaevSSHModel model(HoneycombLattice(2,2),{1,1,.5,g});auto chosen=selectFreeMatterTrial(model,policy);
        // Independent Clifford parent selects the same specified Gaussian;
        // all propagation and insertions below use full 16-dimensional Fock.
        MatType parent=MatType::Zero(16,16);
        for(int a=0;a<8;++a)for(int b=a+1;b<8;++b)parent-=chosen.trial.G_T(a,b)*c[a]*c[b];
        Eigen::SelfAdjointEigenSolver<MatType> eigen(parent);const auto psi=eigen.eigenvectors().col(0).eval();
        PhononPath path(12,6,.13,{1.3,.7},std::vector<double>(12,.2),std::vector<double>(12,.1));
        for(int j=0;j<path.nodes();++j)for(int b=0;b<12;++b)path=path.withCoordinate(j,b,.4*std::sin(.7*b+.3*j)+.1*std::cos(.2*b*j));
        std::vector<MatType> gate;for(int j=0;j<6;++j){
            auto h=[&](int node){MatType v=MatType::Zero(16,16);for(auto e:edges())v+=Complex(0,e.j+g*path.coordinate(node,e.index))*c[e.a]*c[e.b];return v;};
            gate.push_back(exponential(h(j),path.dt()/2));gate.push_back(exponential(h(j+1),path.dt()/2));
        }
        PhononWeightOptions opts{production,production,production,production,production,true,production};
        CoupledPhononPath state(path,chosen.trial,honeycombGenerator(model),2,opts);state.ensureMeasured();
        const auto before=state.state().configuration_hash;const auto signedBefore=state.state().weight.matter_phase;
        const auto greens=state.energyLinkGreens();const auto values=totalEnergyObservables(model,path,greens);
        std::vector<Complex> exact(6,0.);
        for(int side=0;side<2;++side) {
            const int cut=6+(side==0?-1:1),link=3+side;
            cVecType ket=psi;Eigen::RowVectorXcd bra=psi.adjoint();
            for(int j=0;j<cut;++j)ket=gate[j]*ket;
            for(int j=11;j>=cut;--j)bra=bra*gate[j];
            const Complex weight=(bra*ket)(0);
            MatType expected(8,8);for(int a=0;a<8;++a)for(int b=0;b<8;++b)expected(a,b)=(a==b?1.:0.)-(bra*c[a]*c[b]*ket)(0)/weight;
            check("noncommuting_Fock_odd_cut_Green",(expected-greens[side]).norm(),2e-9);
            for(auto e:edges()) {
                const auto ph=phononLinkEnergy(path.coordinate(link-1,e.index),path.coordinate(link,e.index),path.parameters(),path.dt(),path.force(e.index));
                const Complex bilinear=(bra*Complex(0,1)*c[e.a]*c[e.b]*ket)(0)/weight;
                exact[0]+=.5*e.j*bilinear;exact[1]+=.5*g*ph.mean*bilinear;exact[2]+=.5*ph.kinetic;exact[3]+=.5*ph.potential;exact[4]+=.5*ph.source;
            }
        }
        for(int j=0;j<5;++j)exact[5]+=exact[j];
        for(int j=0;j<6;++j)check("noncommuting_Fock_energy_component_"+std::to_string(j),std::abs(values[j]-exact[j]),2e-9);
        check("energy_measurement_state_identity",before==state.state().configuration_hash&&signedBefore==state.state().weight.matter_phase?0.:1.,0.);
        PhononSampler sampler(state,std::mt19937_64(8123));std::ostringstream oldSampler;sampler.saveCheckpoint(oldSampler);
        const auto direct=sampler.totalEnergy(model);std::ostringstream newSampler;sampler.saveCheckpoint(newSampler);
        check("energy_measurement_RNG_and_counts_unchanged",oldSampler.str()==newSampler.str()?0.:1.,0.);
        state.failureInjection=[](const char*op,bool reference){if(std::string(op)=="energy-measurement"&&!reference)throw std::runtime_error("injected energy measurement failure");};
        const auto recovered=sampler.totalEnergy(model);state.failureInjection={};
        for(int j=0;j<6;++j)check("energy_reference_once_same_state",std::abs(direct[j]-recovered[j]),2e-9);
        check("energy_recovery_is_not_a_proposal",sampler.counts().proposed==0&&sampler.counts().first_failures==1&&sampler.counts().reference_recovered==1?0.:1.,0.);
        state.failureInjection=[](const char*op,bool){if(std::string(op)=="energy-measurement")throw std::runtime_error("persistent energy failure");};
        bool aborted=false;try{sampler.totalEnergy(model);}catch(const std::runtime_error&){aborted=true;}state.failureInjection={};
        check("energy_failed_recovery_aborts_without_proposal",aborted&&sampler.counts().proposed==0&&sampler.counts().reference_failed==1?0.:1.,0.);
    }
}
}
int main() {try {
    std::cout<<std::setprecision(17);linkSpectral();coupledQuantumED();arbitraryCutFockAndRecovery();
    for(const auto &m:metrics)std::cout<<"PASS "<<m.first<<" count="<<m.second.count<<" max_error="<<m.second.error<<" tolerance="<<m.second.tolerance<<'\n';
    return 0;
}catch(const std::exception &e){std::cerr<<"FAIL total-energy "<<e.what()<<'\n';return 1;}}
