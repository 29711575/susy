// Validation-only generic six-bond/four-Majorana core driver.
// This bypasses HoneycombLattice construction; it does not add CLI 2x1 support.
#include "phonon_mc.h"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace kitaev_ssh;
int main(int argc,char**argv) {try {
 if(argc!=7)throw std::invalid_argument("omega dt intervals warmup sweeps output-prefix");
 const double omega=std::stod(argv[1]),dt=std::stod(argv[2]),mass=1/(omega*omega),g=.3;
 const int m=std::stoi(argv[3]),warmup=std::stoi(argv[4]),sweeps=std::stoi(argv[5]);
 const std::string prefix=argv[6];
 if(std::ifstream(prefix+"_blocks.csv"))throw std::runtime_error("output exists");
 const MajoranaGenerator gen=[g](const std::vector<double>&x){
  MatType k=MatType::Zero(4,4);
  k(0,2)=DataType(0,-2*(1.5+g*(x.at(1)+x.at(2))));
  k(1,3)=DataType(0,-2*(1.5+g*(x.at(4)+x.at(5))));
  k(0,3)=DataType(0,-2*(1+g*x.at(0))); k(1,2)=DataType(0,-2*(1+g*x.at(3)));
  for(int i=0;i<4;++i)for(int j=i+1;j<4;++j)k(j,i)=-k(i,j);return k;
 };
 const auto trial=GaussianTrialState::fromMajoranaHamiltonian(gen(std::vector<double>(6,0.)));
 std::mt19937_64 rng(40001);
 CoupledPhononPath state(PhononPath(6,m,dt,{mass,omega}).freeSample(rng),trial,gen,4,{true,true,true,true,true,true});
 PhononSampler sampler(state,rng,PhononProposal::BondRefresh);
 const std::vector<std::string>names{"E","X0","X1","X2","X20","X21","X22"};
 PhaseStatistics stats(names,100);SamplingTimings times;
 std::ofstream blocks(prefix+"_blocks.csv");blocks<<std::setprecision(17);stats.writeHeader(blocks);
 const auto warm=sample(sampler,{warmup,sweeps,1,100},stats,[&](const auto&s){
  const int center=m/2;std::vector<double>x(6);for(int b=0;b<6;++b)x[b]=s.path.coordinate(center,b);
  const auto&G=s.weight.center_green;
  Complex E=Complex(0,-1)*((1+g*x[0])*G(0,3)+(1+g*x[3])*G(1,2)+(1.5+g*(x[1]+x[2]))*G(0,2)+(1.5+g*(x[4]+x[5]))*G(1,3));
  std::vector<Complex>out{E};for(int mu=0;mu<3;++mu)out.push_back((x[mu]+x[3+mu])/2);
  for(int mu=0;mu<3;++mu)out.push_back((x[mu]*x[mu]+x[3+mu]*x[3+mu])/2);return out;
 },[&](const auto&s){s.writeLastBlock(blocks);},&times);
 const auto total=sampler.counts(),prod=total-warm;
 std::cout<<std::setprecision(17)<<"METADATA generic_validation_only=1 Nc=2 Majoranas=4 phonon_bonds=6 tau="<<m*dt<<" omega="<<omega<<" mass="<<mass<<" k=1 dt="<<dt<<" M="<<m<<" seed=40001 warmup="<<warmup<<" sweeps="<<sweeps<<" block=100 qr=4 combined=1 defer=1 green=1 prefix=1 cayley=1 proposal=bond-refresh\n";
 std::cout<<"COUNTS warmup_proposed="<<warm.proposed<<" warmup_accepted="<<warm.accepted<<" measurement_proposed="<<prod.proposed<<" measurement_accepted="<<prod.accepted<<" first_failures="<<total.first_failures<<" reference_recovered="<<total.reference_recovered<<" reference_failed="<<total.reference_failed<<" untrusted_rejected="<<total.untrusted_rejected<<"\n";
 for(int bs:{100,200}){const auto s=bs==100?stats.summary():stats.doubledBlockSummary();
  std::cout<<"PHASE block="<<bs<<" samples="<<s.samples<<" blocks="<<s.blocks<<" real="<<s.phase.value.real()<<" imag="<<s.phase.value.imag()<<" error_real="<<s.phase.error_real<<" error_imag="<<s.phase.error_imag<<" magnitude_error="<<s.phase_magnitude_error<<" positive="<<s.positive<<" negative="<<s.negative<<" nonreal="<<s.nonreal<<" ratio_defined="<<s.ratio_defined<<" jackknife_defined="<<s.jackknife_defined<<" phase_resolved="<<s.phase_resolved<<" too_few_blocks="<<s.too_few_blocks<<"\n";
  for(size_t j=0;j<s.observables.size();++j){const auto&v=s.observables[j];std::cout<<"OBS block="<<bs<<" name="<<names[j]<<" real="<<v.value.real()<<" imag="<<v.value.imag()<<" error_real="<<v.error_real<<" error_imag="<<v.error_imag<<" covariance="<<v.covariance<<"\n";}
 }
 std::cout<<"TIMING warmup="<<times.warmup_seconds<<" measurement="<<times.measurement_seconds<<"\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<"FAILED "<<e.what()<<'\n';return 1;}}
