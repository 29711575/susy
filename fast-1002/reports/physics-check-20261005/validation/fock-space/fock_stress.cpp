// Independent Fock construction based on the production multi-body validation oracle.
// New cases extend it to the production frequencies, strong g, and long contours.
#define main existing_p2_validation_main
#include "../../../../kitaev_ssh/phonon_validate.cpp"
#undef main
int main(){try{
 std::cout<<std::setprecision(17);
 const double g=std::sqrt(2.4);const KitaevSSHModel model(HoneycombLattice(2,2),p2Couplings(g));
 const auto trial=freeGroundState(model);const auto generator=honeycombGenerator(model);const auto c=fockMajoranas();
 Eigen::SelfAdjointEigenSolver<MatType> init(fockH(c,std::vector<double>(12,0.),g));const cVecType psi=init.eigenvectors().col(0);
 for(double omega:{.4,.2,.1,.05}){
  int M=int(std::round(80/omega));std::mt19937_64 rng(83201+M);
  for(int style=0;style<3;style++){
   auto path=PhononPath(12,M,.1,{1/(omega*omega),omega}).freeSample(rng);
   for(int j=0;j<path.nodes();j++)for(int b=0;b<12;b++){
    double shift=style==0?0.:style==1?(b%3==2?.51:.93):.5+.4*std::sin(.023*j+.91*b);
    path=path.withCoordinate(j,b,path.coordinate(j,b)+shift);
   }
   PhononWeightOptions options;options.singular_values_only=true;options.merge_node_half_steps=true;options.defer_phase_measurement=true;
   CoupledPhononPath state(path,trial,generator,4,options);state.ensureMeasured();const auto&w=state.state().weight;
   std::vector<MatType> half;half.reserve(M+1);
   for(int j=0;j<=M;j++)half.push_back(fockExp(fockH(c,path.node(j),g),path.dt()/2));
   cVecType right=psi;Eigen::RowVectorXcd left=psi.adjoint();double lr=0,ll=0;
   auto advanceR=[&](const MatType&f){right=f*right;double norm=right.norm();right/=norm;lr+=std::log(norm);};
   auto advanceL=[&](const MatType&f){left=left*f;double norm=left.norm();left/=norm;ll+=std::log(norm);};
   for(int link=1;link<=M/2;link++){advanceR(half[link-1]);advanceR(half[link]);}
   for(int link=M;link>M/2;link--){advanceL(half[link]);advanceL(half[link-1]);}
   DataType overlap=(left*right)(0);double log=lr+ll+std::log(std::abs(overlap));DataType phase=overlap/std::abs(overlap);
   MatType green=MatType::Zero(8,8);for(int i=0;i<8;i++)for(int j=0;j<8;j++)if(i!=j)green(i,j)=-(left*c[i]*c[j]*right)(0)/overlap;
   DataType energy=(left*fockH(c,path.node(M/2),g)*right)(0)/overlap;
   const auto estimate=matterEnergyAtDisplacement(model,path.node(M/2),w.center_green);
   double le=std::abs(log-w.log_abs_matter),pe=std::abs(phase-w.matter_phase),ge=(green-w.center_green).cwiseAbs().maxCoeff(),ee=std::abs(energy-estimate);
   std::cout<<"omega="<<omega<<" M="<<M<<" style="<<style<<" log_error="<<le<<" phase_error="<<pe<<" green_error="<<ge<<" energy_error="<<ee<<" matter_energy_real="<<estimate.real()<<'\n'<<std::flush;
   if(le>1e-8||pe>1e-10||ge>1e-8||ee>1e-8)return 2;
  }
 }
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
