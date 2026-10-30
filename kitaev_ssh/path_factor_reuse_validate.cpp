// Deterministic ownership/certificate checks. No MC, sampler or RNG.
#include "phonon_mc.h"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>
using namespace kitaev_ssh;
namespace {
int checks=0;
void require(bool ok,const char *what){++checks;if(!ok)throw std::runtime_error(what);}
bool same(const MatType &a,const MatType &b){return a.rows()==b.rows()&&a.cols()==b.cols()&&std::memcmp(a.data(),b.data(),sizeof(DataType)*a.size())==0;}
std::uint64_t hash(const PureFactorMatrix &a){std::uint64_t h=1469598103934665603ULL;a.mixHash(h);return h;}
using Pool=PureFactorMatrix::DenseStoragePool;
PureFactorMatrix make(Pool &pool,int n,double value){return PureFactorMatrix::fromValidatedDense(&pool,[&](MatType &m){m=MatType::Constant(n,n,DataType(value,-value/3));});}
template<class F> std::string fallbackFailure(F &&evaluate){
    bool runtime=false,invalid=false;std::string message;
    try{evaluate();}
    catch(const std::invalid_argument &e){invalid=true;message=e.what();}
    catch(const std::runtime_error &e){runtime=true;message=e.what();}
    require(!invalid,"overflow fallback bypassed runtime-error recovery with invalid_argument");
    require(runtime,"overflow fallback did not fail through runtime_error");
    require(message=="invalid merged Majorana node factor","overflow fallback changed original failure message");
    return message;
}
void ownership(){
    static_assert(!std::is_copy_constructible<Pool>::value,"pool copies must be explicit/nonsharing");
    static_assert(!std::is_copy_assignable<Pool>::value,"pool assignment must not alias");
    static_assert(std::is_move_constructible<Pool>::value,"pool move expected");
    Pool pool;pool.beginBatch();auto held=make(pool,4,2.);auto copy=held;
    const MatType saved=*held.denseStorage();const auto oldhash=hash(held);
    const auto *owner=held.denseStorage();
    for(int n:{4,8,2,4}){pool.beginBatch();auto next=make(pool,n,3.);require(next.denseStorage()!=owner,"live owner overwritten");require(same(*held.denseStorage(),saved)&&hash(held)==oldhash,"published factor mutated");}
    require(copy.sharesStorageWith(held),"factor copy must pin same immutable owner");
    auto moved=std::move(copy);require(moved.sharesStorageWith(held),"factor move lost immutable owner");
    const auto reused=pool.slot_reuses;pool.beginBatch();{auto next=make(pool,6,1.);require(next.rows()==6,"pool resize failed");}
    require(pool.slot_reuses>reused,"released slot not reused");
    auto output=held;
    bool threw=false;pool.beginBatch();try{PureFactorMatrix::tryFromValidatedDense(&pool,[](MatType &m){m=MatType::Identity(3,3);throw std::runtime_error("writer");},output);}catch(const std::runtime_error &){threw=true;}
    require(threw&&output.sharesStorageWith(held)&&same(*held.denseStorage(),saved),"throwing writer changed output/live owner");
    pool.beginBatch();const bool finite=PureFactorMatrix::tryFromValidatedDense(&pool,[](MatType &m){m=MatType::Constant(3,3,DataType(std::numeric_limits<double>::infinity(),0));},output);
    require(!finite&&output.sharesStorageWith(held)&&hash(output)==oldhash,"nonfinite writer published output");
    Pool destination=std::move(pool);destination.beginBatch();auto next=make(destination,4,7.);
    require(same(*held.denseStorage(),saved),"moving pool invalidated external owners");
    require(!next.certifiedRealStorage(),"generic pool writer acquired certificate");
}
void certificate(){
    for(int L:{2,4}){
        KitaevSSHModel model(HoneycombLattice(L,L,0),p2Couplings(.8));
        std::vector<double> x(model.lattice().bondCount());for(std::size_t b=0;b<x.size();++b)x[b]=.1*std::sin(.17*b);
        const MatType k=model.frameworkGeneratorWith(x);BipartiteSpectrum spectral;
        require(spectral.compute(k),"fixed bipartite spectrum failed");Pool pool;
        const MatType *certifiedSlot=nullptr;
        for(double a:{0.,.05,.1}){
            pool.beginBatch();PureFactorMatrix factor;
            require(spectral.exponentialFactor(a,factor,&pool,true),"certified assembly failed");
            certifiedSlot=factor.denseStorage();
            MatType old;spectral.exponential(a,old);
            require(same(*factor.denseStorage(),old),"certification changed complex coefficients or signed zeros");
            require(hash(factor)==hash(PureFactorMatrix::fromValidatedDense(old)),"certificate changed coefficient hash");
            const auto *real=factor.certifiedRealStorage();require(real!=nullptr,"certificate absent");
            PureRealMagnitudeWorkspace workspace;
            auto generic=PureFactorMatrix::fromValidatedDense(old);
            const auto *reference=workspace.factor(generic,0,1);require(reference!=nullptr,"generic exact-real conversion failed");
            require(workspace.factor(generic,0,1)==reference&&workspace.cache_hits==1,"generic immutable cache miss");
            require(real->rows()==reference->rows()&&std::memcmp(real->data(),reference->data(),sizeof(double)*real->size())==0,"certificate differs from exact legacy D-dagger-B-D conversion");
            require(workspace.factor(factor,0,1)==real&&workspace.certified_hits==1,"certificate not borrowed directly");
            require(workspace.conversions==1,"certified path performed conversion");
            auto slice=PureProjectorSlice::fromValidatedDense(factor,1.,"certificate");
            require(slice.matrix.sharesStorageWith(factor),"slice copied certified owner");
            require(slice.identity_residual==(old-MatType::Identity(old.rows(),old.cols())).norm(),"identity residual changed");
            MatType extra=old;extra(0,0)+=DataType(0.,1e-30);
            auto unsupported=PureFactorMatrix::fromValidatedDense(extra);
            require(!unsupported.certifiedRealStorage()&&workspace.factor(unsupported,0,1)==nullptr,"tiny extra complex coupling bypassed exact structure fallback");
            workspace.releaseFactorOwners();
            require(workspace.factor(generic,0,1)!=nullptr&&workspace.conversions==2,"released generic owner reused stale conversion");
            auto before=factor;const auto beforehash=hash(before);pool.beginBatch();
            require(!spectral.exponentialFactor(1e308,factor,&pool,true),"overflow factor unexpectedly published");
            require(factor.sharesStorageWith(before)&&hash(factor)==beforehash,"overflow assembly damaged previous result");
            // Overflow recovery itself belongs to prepareOwnedPath: independent
            // reference/cache validators must confirm its original fallback.
        }
        // Every certified owner above has now been released. An ordinary
        // writer must invalidate the companion left in that reusable slot.
        const auto reused=pool.slot_reuses;pool.beginBatch();auto ordinary=make(pool,k.rows(),1.);
        require(pool.slot_reuses>reused&&ordinary.denseStorage()==certifiedSlot,"ordinary writer did not reuse released certified slot");
        require(!ordinary.certifiedRealStorage(),"ordinary writer retained stale real certificate");
        PureRealMagnitudeWorkspace ordinaryWorkspace;
        require(ordinaryWorkspace.factor(ordinary,0,1)==nullptr&&ordinaryWorkspace.certified_hits==0,
                "ordinary writer borrowed stale certified real coefficients");
        MatType incompatible=k;const int n=k.rows()/2;
        incompatible(0,n)+=DataType(1e-30,0.);incompatible(n,0)+=DataType(1e-30,0.);
        BipartiteSpectrum workspace;PhononNodeSpectrum general;
        require(general.compute(incompatible,true,workspace)&&!general.reduced(),"tiny Hermitian real cross coupling entered reduced spectrum");
        Pool referencePool;PureFactorMatrix referenceFactor;
        require(general.exponentialFactor(.05,referenceFactor,&referencePool,true),"general complex assembly failed");
        require(!referenceFactor.certifiedRealStorage(),"general spectrum issued a real certificate");
        MatType expected;general.exponential(.05,expected);
        require(same(*referenceFactor.denseStorage(),expected),"general complex assembly changed reference coefficients");
    }
}
void overflowFallback(){
    // K=[0,-i t I;i t I,0]. At dt*t=800 the full exponential and
    // half*half overflow, while every half-step coefficient remains finite.
    // Endpoints stay mild so the interior recovery is the failing gate.
    MajoranaGenerator::IntoFunction into=[](const std::vector<double> &x,MatType &k){
        k.setZero(4,4);for(int i=0;i<2;++i){k(i,2+i)=DataType(0.,-x[0]);k(2+i,i)=DataType(0.,x[0]);}
    };
    MajoranaGenerator generator([into](const std::vector<double> &x){MatType k;into(x,k);return k;},into,true,true);
    auto path=PhononPath(1,2,1.,{1.,1.}).withCoordinates({{.02},{.1},{.02}});
    const auto extreme=path.withCoordinates({{.02},{800.},{.02}});
    BipartiteSpectrum spectral;require(spectral.compute(generator(extreme.node(1))),"overflow fixture spectrum failed");
    MatType half,full;spectral.exponential(.5,half);spectral.exponential(1.,full);
    require(half.allFinite(),"overflow fixture half-step is already nonfinite");
    require(!full.allFinite()&&!(half*half).eval().allFinite(),"overflow fixture misses full/half-square failure");
    Pool pool;
    const auto prepared=prepareOwnedPath(path,4,generator,false,nullptr,false,&pool);
    const auto held=prepared.weightOrder[1].matrix;const MatType saved=*held.denseStorage();
    const auto oldhash=hash(held);
    const auto referenceMessage=fallbackFailure([&]{(void)prepareOwnedPath(extreme,4,generator,false,nullptr,false);});
    const auto pooledMessage=fallbackFailure([&]{(void)prepareOwnedPath(extreme,4,generator,false,nullptr,false,&pool);});
    require(referenceMessage==pooledMessage,"pooled/reference overflow failure differs");
    require(same(*held.denseStorage(),saved)&&hash(held)==oldhash,"overflow fallback mutated a published factor");
    const auto recovered=prepareOwnedPath(path,4,generator,false,nullptr,false,&pool);
    for(const auto &s:recovered.weightOrder)require(s.matrix.allFinite(),"failed pool slot prevented normal assembly");
    require(same(*held.denseStorage(),saved)&&hash(held)==oldhash,"normal recovery overwrote a held factor");

    MatType phi=MatType::Zero(4,2);const double scale=1./std::sqrt(2.);
    for(int i=0;i<2;++i){phi(i,i)=scale;phi(2+i,i)=DataType(0.,scale);}
    CoupledPhononPath state(path,GaussianTrialState::fromPhi(phi),generator,2,{true,true,true});
    const auto snapshot=state.state();const auto configuration=state.configurationHash();
    std::vector<MatType> snapshotFactors;for(const auto &s:snapshot.weight_order)snapshotFactors.emplace_back(s.matrix);
    fallbackFailure([&]{(void)state.proposeBondPath(0,{.02,800.,.02});});
    require(state.configurationHash()==configuration&&state.state().weight.log_abs_matter==snapshot.weight.log_abs_matter,
            "overflow fallback changed accepted Snapshot");
    for(std::size_t i=0;i<snapshotFactors.size();++i)
        require(same(MatType(snapshot.weight_order[i].matrix),snapshotFactors[i])&&same(MatType(state.state().weight_order[i].matrix),snapshotFactors[i]),
                "overflow fallback overwrote retained Snapshot coefficients");
    auto next=state.proposeBondPath(0,{.02,.11,.02});
    require(std::isfinite(next.candidate().weight.log_abs_matter),"normal proposal failed after overflow fallback");
    state.accept(next);state.ensureMeasured();
    require(state.state().weight.phase_valid&&state.state().weight.green_valid,"normal proposal measurement failed after overflow fallback");
    for(std::size_t i=0;i<snapshotFactors.size();++i)
        require(same(MatType(snapshot.weight_order[i].matrix),snapshotFactors[i]),"normal recovery overwrote retained Snapshot");
}
void proposals(){
    KitaevSSHModel model(HoneycombLattice(4,4,0),p2Couplings(.8));
    PhononPath path(model.lattice().bondCount(),8,.1,{100.,.1});
    std::vector<std::vector<double>> x(path.nodes(),std::vector<double>(path.bonds()));
    for(int j=0;j<path.nodes();++j)for(int b=0;b<path.bonds();++b)x[j][b]=.08*std::sin(.11*j+.07*b);
    path=path.withCoordinates(x);const auto trial=freeGroundState(model);
    CoupledPhononPath state(path,trial,honeycombGenerator(model),4,{true,true,true,true,true,true,true});
    const auto saved=state.state();const auto initial=state.configurationHash();
    auto fingerprints=[](const std::vector<PureProjectorSlice> &s){std::vector<std::uint64_t> out;for(const auto &f:s)out.push_back(hash(f.matrix));return out;};
    const auto initialFactors=fingerprints(saved.weight_order);
    std::vector<double> v(path.nodes());for(int j=0;j<path.nodes();++j)v[j]=x[j][0]+.01*std::cos(.3*j);
    auto a=state.proposeBondPath(0,v);const auto aHashes=fingerprints(a.candidate().weight_order);const auto aLog=a.candidate().weight.log_abs_matter;
    for(int j=0;j<path.nodes();++j)v[j]=x[j][1]-.02*std::sin(.2*j);
    auto b=state.proposeBondPath(1,v);
    require(state.configurationHash()==initial,"uncommitted candidates changed accepted path");
    require(fingerprints(saved.weight_order)==initialFactors,"old Snapshot factors overwritten");
    require(fingerprints(a.candidate().weight_order)==aHashes&&a.candidate().weight.log_abs_matter==aLog,"second proposal overwrote first proposal");
    state.reject(b);require(fingerprints(a.candidate().weight_order)==aHashes,"reject overwrote held candidate");
    state.accept(a);state.ensureMeasured();require(fingerprints(saved.weight_order)==initialFactors,"accept/measurement overwrote retained Snapshot");
    // Numerical accept/reject, same-candidate uncached recovery, signed Green
    // and Fock checks deliberately reuse existing validators and common replay.
}
}
int main(){try{ownership();certificate();overflowFallback();proposals();std::cout<<"{\"checks\":"<<checks<<"}\nPATH_FACTOR_REUSE_VALIDATE_OK\n";}catch(const std::exception &e){std::cerr<<"PATH_FACTOR_REUSE_VALIDATE_FAIL: "<<e.what()<<'\n';return 1;}}
