#include "phonon_mc.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace kitaev_ssh;
MajoranaGenerator gen(double g){return [g](const std::vector<double>&x){MatType k=MatType::Zero(2,2);k(0,1)=DataType(0,-2*(.5+g*x.at(0)));k(1,0)=-k(0,1);return k;};}
std::vector<double> meanPath(const PhononPath&p,double g){
 std::vector<double> m(p.nodes(),0.);const double v=1/(2*p.parameters().mass*p.parameters().omega);
 for(int j=0;j<p.nodes();j++) for(int k=0;k<p.nodes();k++)m[j]+=g*v*std::exp(-p.parameters().omega*p.dt()*std::abs(j-k))*p.dt()*((k==0||k==p.intervals())?.5:1.);
 return m;
}
double analyticLog(const PhononPath&p,double g){double out=0.;for(int j=0;j<p.nodes();j++)out+=p.dt()*((j==0||j==p.intervals())?.5:1.)*(.5+g*p.coordinate(j,0));return out;}
int main(int argc,char**argv){try{
 std::cout<<std::setprecision(17);
 const std::string mode=argc>1?argv[1]:"oracle";
 const double omega=argc>2?std::stod(argv[2]):.05;
 const int sweeps=argc>3?std::stoi(argv[3]):100;
 const std::string init=argc>4?argv[4]:"free";
 const auto seed=argc>5?std::stoull(argv[5]):40955;
 const double g=std::sqrt(2.4);const int M=int(std::round(80/omega));
 const PhononPath blank(1,M,.1,{1/(omega*omega),omega});
 const auto mean=meanPath(blank,g);auto generator=gen(g);
 const auto trial=GaussianTrialState::fromMajoranaHamiltonian(generator({0.}));
 PhononWeightOptions options;options.singular_values_only=true;options.merge_node_half_steps=true;options.defer_phase_measurement=true;
 std::mt19937_64 rng(seed);
 auto initial=blank.freeSample(rng);
 if(init=="equilibrium")for(int j=0;j<initial.nodes();j++)initial=initial.withCoordinate(j,0,initial.coordinate(j,0)+mean[j]);
 if(mode=="oracle"){
  double logerr=0,phaseerr=0,greenerr=0,energerr=0,ratioerr=0,densityerr=0;
  for(int k=0;k<12;k++){
   auto path=blank.freeSample(rng);
   // Analytically tilted free law, constant deformation, and strongly asymmetric profiles.
   for(int j=0;j<path.nodes();j++)path=path.withCoordinate(j,0,path.coordinate(j,0)+(k%3==0?mean[j]:(k%3==1?.93:.4*std::sin(.017*j)+.2)));
   CoupledPhononPath s(path,trial,generator,4,options);
   logerr=std::max(logerr,std::abs(s.state().weight.log_abs_matter-analyticLog(path,g)));
   const auto&w=s.ensureMeasured().weight;
   phaseerr=std::max(phaseerr,std::abs(w.matter_phase-DataType(1)));
   greenerr=std::max(greenerr,std::abs(w.center_green(0,1)-DataType(0,-1)));
   const auto e=-DataType(0,1)*(.5+g*path.coordinate(path.center(),0))*w.center_green(0,1);
   energerr=std::max(energerr,std::abs(e-DataType(-.5-g*path.coordinate(path.center(),0))));
   auto next=blank.freeSample(rng);const auto before=s.state();auto proposal=s.proposeBondPath(0,[&]{std::vector<double>v(next.nodes());for(int j=0;j<next.nodes();j++)v[j]=next.coordinate(j,0);return v;}());
   const auto&p=proposal.candidate();const double simple=p.weight.log_abs_matter-before.weight.log_abs_matter;
   ratioerr=std::max(ratioerr,std::abs(simple-(analyticLog(next,g)-analyticLog(path,g))));
   const double general=proposal.logAbsRatio()+freeBondLogDensity(before.path,0)-freeBondLogDensity(p.path,0);
   densityerr=std::max(densityerr,std::abs(general-simple));
  }
  std::cout<<"omega="<<omega<<" M="<<M<<" exact_center="<<mean[M/2]<<" log_error="<<logerr<<" phase_error="<<phaseerr<<" green_error="<<greenerr<<" energy_error="<<energerr<<" ratio_error="<<ratioerr<<" MH_error="<<densityerr<<'\n';
  if(logerr>1e-8||phaseerr>1e-10||greenerr>1e-10||energerr>1e-9||ratioerr>1e-8||densityerr>1e-8)return 2;
 }else{
  CoupledPhononPath s(initial,trial,generator,4,options);PhononSampler sampler(s,rng,PhononProposal::BondRefresh);
  std::cout<<"sweep,center_X,time_average_X,accepted,proposed,phase_real,energy\n";
  double avg=0.; std::uint64_t guidedAccepted=0; double guidedRatioMax=0.;
  for(int n=1;n<=sweeps;n++){
   if(mode=="guided") {
    auto draw=blank.freeSample(rng);std::vector<double>v(draw.nodes());for(int j=0;j<draw.nodes();j++)v[j]=draw.coordinate(j,0)+mean[j];
    const auto before=s.state();auto proposal=s.proposeBondPath(0,v);const auto&next=proposal.candidate();
    auto oldResidual=before.path,newResidual=next.path;
    for(int j=0;j<draw.nodes();j++){oldResidual=oldResidual.withCoordinate(j,0,oldResidual.coordinate(j,0)-mean[j]);newResidual=newResidual.withCoordinate(j,0,newResidual.coordinate(j,0)-mean[j]);}
    const double logR=proposal.logAbsRatio()+freeBondLogDensity(oldResidual,0)-freeBondLogDensity(newResidual,0);
    guidedRatioMax=std::max(guidedRatioMax,std::abs(logR));
    if(std::log((double(rng()>>12)+.5)*0x1p-52)<std::min(0.,logR)){s.accept(proposal);guidedAccepted++;}else s.reject(proposal);
    s.ensureMeasured();
   }else{sampler.sweep();sampler.ensureMeasured();}
   const auto&w=sampler.state().weight;
   const auto x=sampler.state().path.coordinate(M/2,0);avg+=x;
   std::cout<<n<<','<<x<<','<<avg/n<<','<<(mode=="guided"?guidedAccepted:sampler.counts().accepted)<<','<<(mode=="guided"?n:sampler.counts().proposed)<<','<<w.matter_phase.real()<<','<<(-DataType(0,1)*(.5+g*x)*w.center_green(0,1)).real()<<'\n';
  }
  std::cerr<<std::setprecision(17)<<"omega="<<omega<<" init="<<init<<" seed="<<seed<<" exact_center="<<mean[M/2]<<" measured_center="<<avg/sweeps<<" accepted="<<sampler.counts().accepted<<" proposed="<<sampler.counts().proposed<<" numerical_failures="<<sampler.counts().first_failures<<" guided_accept="<<guidedAccepted<<" guided_ratio_max="<<guidedRatioMax<<'\n';
 }
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
