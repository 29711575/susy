#define KITAEV_SSH_TEST_FAILURE_INJECTION
#define main original_fixed_path_fock_main
#include "pfqmc_fock_validate.cpp"
#undef main
#include "bond_diagnostics.h"
#include <fstream>
#include <filesystem>
using namespace kitaev_ssh;
namespace {
MatType physicalH(const std::vector<MatType>&c,const std::vector<double>&x,double g) {
 MatType h=MatType::Zero(16,16);for(auto e:edges())h+=Complex(0,e.j+g*x[e.index])*c[e.a]*c[e.b];return h;
}
struct PhysicalFock {
 Complex w;MatType G;std::vector<MatType> c;std::vector<cVecType> right;std::vector<Eigen::RowVectorXcd> left;
 PhysicalFock(const PhononPath&p,double g):c(operators(4)) {
  Eigen::SelfAdjointEigenSolver<MatType> es(physicalH(c,std::vector<double>(12),0));
  std::vector<MatType> u;for(int j=0;j<p.intervals();++j){u.push_back(exponential(physicalH(c,p.node(j),g),p.dt()/2));u.push_back(exponential(physicalH(c,p.node(j+1),g),p.dt()/2));}
  right.resize(u.size()+1);left.resize(u.size()+1);right[0]=es.eigenvectors().col(0);left.back()=right[0].adjoint();
  for(int f=0;f<int(u.size());++f)right[f+1]=u[f]*right[f];
  for(int f=int(u.size())-1;f>=0;--f)left[f]=left[f+1]*u[f];w=(left[0]*right[0])(0);
  G=MatType::Identity(8,8);for(int a=0;a<8;++a)for(int b=0;b<8;++b)G(a,b)-=value(c[a]*c[b],p.intervals());
 }
 Complex value(const MatType&o,int cut)const{return (left[cut]*o*right[cut])(0)/w;}
};
Eigen::MatrixXd precision(const PhononPath&p) {
 int M=p.intervals();const double w=p.parameters().omega,m=p.parameters().mass,z=w*p.dt();Eigen::MatrixXd Q=Eigen::MatrixXd::Zero(M+1,M+1);
 for(int j=0;j<M;++j){Q(j,j)+=m*w/std::tanh(z);Q(j+1,j+1)+=m*w/std::tanh(z);Q(j,j+1)=Q(j+1,j)=-m*w/std::sinh(z);}
 Q(0,0)+=m*w;Q(M,M)+=m*w;return Q;
}
Eigen::VectorXd column(const PhononPath&p,int b){Eigen::VectorXd v(p.nodes());for(int j=0;j<p.nodes();++j)v[j]=p.coordinate(j,b);return v;}
double actionOracle(const PhononPath&p,const Eigen::MatrixXd&Q){double s=0;for(int b=0;b<p.bonds();++b){auto x=column(p,b);s+=.5*x.dot(Q*x);}return s;}
double gaussian(const Eigen::VectorXd&x,const Eigen::VectorXd&mu,const Eigen::MatrixXd&Q,double beta=1.){
 auto d=(x-mu).eval();return -.5*(d.dot(Q*d)/(beta*beta)+x.size()*std::log(2*pi*beta*beta)-Q.ldlt().vectorD().array().log().sum());
}
std::string stateFingerprint(const CoupledPhononPath::Snapshot&s){
 std::ostringstream o;o<<std::setprecision(17);s.path.save(o);o<<s.configuration_hash<<s.weight_order_is_local<<s.weight.logAbsTotal()<<s.weight.matter_phase<<s.weight.phase_valid<<s.weight.green_valid;
 auto mat=[&](const MatType&m){o<<m;};mat(s.weight.center_green);
 for(auto f:s.factors){mat(MatType(f.matrix));o<<f.eta;}for(auto f:s.weight_order){mat(MatType(f.matrix));o<<f.eta;}
 if(s.signed_prefix_cache)for(const auto&p:s.signed_prefix_cache->checkpoints){o<<bool(p);if(p){mat(p->q);o<<p->log_det_r<<p->previous_log_overlap<<p->det_r_phase<<p->previous_overlap_phase<<p->log_abs_weight<<p->complex_phase<<p->prefix_hash;}}
 if(s.lu_cayley_cache)for(const auto&p:*s.lu_cayley_cache){o<<bool(p);if(p)mat(p->value());}return o.str();
}
void savePath(const std::string&file,const PhononPath&p){std::ofstream out(file);p.save(out);if(!out)throw std::runtime_error("path save failed");}
void compareMeasurement(CoupledPhononPath&state,const KitaevSSHModel&model,double g){
 const auto&s=state.ensureMeasured();PhysicalFock ref(s.path,g);const auto&w=s.weight;
 check("A_log_abs_weight",std::abs(w.log_abs_matter-std::log(std::abs(ref.w))));check("A_phase",std::abs(w.matter_phase-ref.w/std::abs(ref.w)));
 check("A_full_complex_amplitude",std::abs(std::exp(w.log_abs_matter)*w.matter_phase/ref.w-Complex(1)));
 check("A_center_Green",(w.center_green-ref.G).cwiseAbs().maxCoeff());
 check("A_matter_energy",std::abs(matterEnergyAtDisplacement(model,s.path.node(s.path.center()),w.center_green)-ref.value(physicalH(ref.c,s.path.node(s.path.center()),g),s.path.intervals())));
 const auto measured=bondDiagnostics(model,s);
 for(auto e:edges()) {
  auto op=(Complex(0,1)*ref.c[e.a]*ref.c[e.b]).eval();const auto expected=ref.value(op,s.path.intervals());
  check("A_production_BondK",std::abs(measured[2*e.index]-expected));check("A_production_BondXK",std::abs(measured[2*e.index+1]-s.path.coordinate(s.path.center(),e.index)*expected));
  for(auto f:edges()) {
   const MatType C=MatType::Identity(8,8)-w.center_green;
   const Complex wick=-(C(e.a,e.b)*C(f.a,f.b)-C(e.a,f.a)*C(e.b,f.b)+C(e.a,f.b)*C(e.b,f.a));
   check("A_Wick_four_point",std::abs(wick-ref.value(-ref.c[e.a]*ref.c[e.b]*ref.c[f.a]*ref.c[f.b],s.path.intervals())));
  }
 }
}
void suite(const std::string&root) {
 std::filesystem::create_directories(root);const HoneycombLattice lat(2,2);
 std::ofstream table(root+"/bond-table.csv");table<<"index,A,B,mu,u,J,cell_x,cell_y\n";
 for(auto e:edges()){const auto&b=lat.bond(e.index);check("A_independent_bond_table",double(e.a!=b.a||e.b!=b.b||e.index%3!=int(b.type)||lat.zeroFluxGauge()[e.index]!=1),0.);table<<e.index<<','<<e.a<<','<<e.b<<','<<int(b.type)<<",1,"<<e.j<<','<<b.cell.x<<','<<b.cell.y<<'\n';}
 std::ofstream logs(root+"/MH.csv");logs<<"g,policy,mode,recovered,decision,log_matter_ratio,phonon_log_ratio,log_q_forward,log_q_reverse,full_log_MH,production_general_MH\n"<<std::setprecision(17);
 for(double g:{.2,.4}) {
  const KitaevSSHModel model(lat,{1,1,.5,g});auto trial=freeGroundState(model);auto gen=honeycombGenerator(model);
  for(int production:{0,1}) {
   PhononWeightOptions opts{bool(production),bool(production),bool(production),bool(production),bool(production),true,bool(production)};
   for(int kind=0;kind<3;++kind) {
    PhononPath p(12,10,.2,{1,1});for(int j=0;j<p.nodes();++j)for(int b=0;b<12;++b)p=p.withCoordinate(j,b,kind==0?0.:kind==1?.3*std::cos(.43*b):.3*std::sin(.43*b+.71*j)+.2*std::cos(.7*b-.11*j));
    savePath(root+"/g"+std::to_string(g)+"-kind"+std::to_string(kind)+".path",p);
    CoupledPhononPath state(p,trial,gen,4,opts);compareMeasurement(state,model,g);
    if(kind==2){auto c=operators(4);auto a=physicalH(c,p.node(1),g),b=physicalH(c,p.node(2),g);check("A_noncommuting_physical_path",(a*b-b*a).norm()>1e-5?0.:1.,0.);}
   }
   for(int mode=0;mode<3;++mode)for(int recover:{0,1})for(int accept:{0,1}) {
    PhononPath old(12,10,.2,{1,1});for(int b=0;b<12;++b)for(int j=0;j<old.nodes();++j)old=old.withCoordinate(j,b,.2*std::sin(.43*b+.31*j));auto next=old;
    const auto Q=precision(old);const auto C=Q.ldlt().solve(Eigen::MatrixXd::Identity(old.nodes(),old.nodes())).eval();
    const int bond=2,node=5;const double beta=.7,rho=std::sqrt(1-beta*beta),force=.19;
    auto qmode=makeQModeOptions(lat,1,1,1,beta);GuidedBondOptions guide;guide.beta=beta;guide.forces.assign(12,force);
    if(mode==0)next=old.withCoordinate(node,bond,old.coordinate(node,bond)+.6);
    if(mode==1){auto x=column(old,bond);for(int j=0;j<old.nodes();++j)x[j]+=.6+.2*std::cos(.2*j);next=old.withBondCoordinates(bond,std::vector<double>(x.data(),x.data()+x.size()));}
    if(mode==2)for(const auto&b:lat.bonds())if(b.type==BondType::Z)for(int j=0;j<old.nodes();++j)next=next.withCoordinate(j,b.index,old.coordinate(j,b.index)+.5*((b.cell.x+b.cell.y)%2?-1.:1.)*(.7+.1*std::cos(.2*j)));
    auto qdensity=[&](const PhononPath&from,const PhononPath&to){
     if(mode==0){auto x=column(from,bond);double mu=x[node]-(Q*x)[node]/Q(node,node),var=1/Q(node,node);return -.5*(std::log(2*pi*var)+std::pow(to.coordinate(node,bond)-mu,2)/var);}
     if(mode==1){Eigen::VectorXd src(old.nodes());for(int j=0;j<old.nodes();++j)src[j]=force*(j==0||j==old.intervals()?1.:2.)*std::tanh(old.dt()/2);auto mu=(C*src).eval();return gaussian(column(to,bond),(mu+rho*(column(from,bond)-mu)).eval(),Q,beta);}
     Eigen::VectorXd a=Eigen::VectorXd::Zero(old.nodes()),b=a;
     for(int j=0;j<old.nodes();++j)for(int x=0;x<2;++x)for(int y=0;y<2;++y){const int k=3*(2*x+y)+2;const double f=((x+y)%2?-1.:1.)/2;a[j]+=f*from.coordinate(j,k);b[j]+=f*to.coordinate(j,k);}return gaussian(b,(rho*a).eval(),Q,beta);
    };
    auto ratio=[&](){return std::log(std::abs(PhysicalFock(next,g).w/PhysicalFock(old,g).w))-actionOracle(next,Q)+actionOracle(old,Q)+qdensity(next,old)-qdensity(old,next);};
    if(ratio()>0)std::swap(old,next);const double expected=ratio();check("A_nontrivial_unfavorable_MH",expected < -1e-7?0.:1.,0.);
    const double u=accept?std::exp(std::max(-600.,expected-1.)):(1+std::exp(expected))/2;
    CoupledPhononPath state(old,trial,gen,4,opts);std::vector<NumericalEvent> events;
    PhononSampler sampler(state,std::mt19937_64(20261017),mode==1?PhononProposal::GuidedBond:PhononProposal::LocalConditional,[&](const auto&e,const auto&){events.push_back(e);},qmode,mode==1?guide:GuidedBondOptions{});
    const auto before=stateFingerprint(state.state());
    if(recover)state.failureInjection=[](const char*,bool reference){if(!reference)throw std::runtime_error("A same-candidate recovery test");};
    MoveResult move;
    if(mode==0)move=sampler.attempt(node,bond,next.coordinate(node,bond),u);
    if(mode==1){auto x=column(next,bond);move=sampler.attemptGuidedBond(bond,std::vector<double>(x.data(),x.data()+x.size()),u);}
    if(mode==2){QModeAmplitudes a(1,std::vector<double>(old.nodes()));for(int j=0;j<old.nodes();++j)for(int x=0;x<2;++x)for(int y=0;y<2;++y)a[0][j]+=((x+y)%2?-1.:1.)/2*next.coordinate(j,3*(2*x+y)+2);move=sampler.attemptQMode(a,u);}
    state.failureInjection={};check("A_full_Fock_Gaussian_MH",std::abs(move.general_log_ratio-expected));check("A_simplified_MH",std::abs(move.log_ratio-expected));check("A_MH_decision",move.accepted==bool(accept)?0.:1.,0.);
    if(!accept)check("A_reject_deep_state",before==stateFingerprint(state.state())?0.:1.,0.);
    else for(int j=0;j<old.nodes();++j)for(int b=0;b<12;++b)check("A_accept_coordinates",std::abs(state.state().path.coordinate(j,b)-next.coordinate(j,b)),2e-14);
    compareMeasurement(state,model,g);
    if(recover)check("A_recovery_same_candidate",events.size()==1&&events[0].recovered?0.:1.,0.);
    // Independently force measurement certification recovery after a lazy move.
    if(production&&accept){auto p=state.state().path;CoupledPhononPath lazy(p,trial,gen,4,opts);lazy.failureInjection=[](const char*op,bool reference){if(std::string(op)=="measurement"&&!reference)throw std::runtime_error("A deferred measurement recovery");};PhononSampler ms(lazy,std::mt19937_64(7));ms.ensureMeasured();lazy.failureInjection={};compareMeasurement(lazy,model,g);check("A_measurement_recovered",ms.counts().reference_recovered==1?0.:1.,0.);}
    std::string stem=root+"/MH-g"+std::to_string(g)+"-p"+std::to_string(production)+"-m"+std::to_string(mode)+"-r"+std::to_string(recover)+"-a"+std::to_string(accept);savePath(stem+"-old.path",old);savePath(stem+"-new.path",next);
    logs<<g<<','<<production<<','<<mode<<','<<recover<<','<<accept<<','<<std::log(std::abs(PhysicalFock(next,g).w/PhysicalFock(old,g).w))<<','<<-actionOracle(next,Q)+actionOracle(old,Q)<<','<<qdensity(old,next)<<','<<qdensity(next,old)<<','<<expected<<','<<move.general_log_ratio<<'\n';
   }
  }
 }
}
}
int main(int argc,char**argv){try{if(argc!=2)throw std::invalid_argument("usage: interacting_fock_validate NEW_OUTPUT_DIRECTORY");if(std::filesystem::exists(argv[1]))throw std::invalid_argument("output exists");std::cout<<std::setprecision(17);suite(argv[1]);for(auto&m:metrics)std::cout<<"PASS "<<m.first<<" count="<<m.second.count<<" max_error="<<m.second.error<<" tolerance="<<m.second.tolerance<<'\n';return 0;}catch(const std::exception&e){std::cerr<<"FAIL A "<<e.what()<<'\n';return 1;}}
