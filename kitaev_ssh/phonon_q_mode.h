#ifndef KITAEV_SSH_PHONON_Q_MODE_H
#define KITAEV_SSH_PHONON_Q_MODE_H
#include "honeycomb_lattice.h"
#include "phonon_path.h"
#include <istream>
#include <sstream>
#include <set>
namespace kitaev_ssh {
struct QModeOptions {
    int every=0,n1=0,n2=0;
    double beta=1.;
    std::string bond_types="z";
    std::vector<std::vector<double>> basis;
};
inline QModeOptions makeQModeOptions(const HoneycombLattice &lat,int n1,int n2,int every,double beta,const std::string &bondTypes="z") {
    if(every<0 || !std::isfinite(beta) || beta<=0 || beta>1 || n1<0 || n1>=lat.Lx() || n2<0 || n2>=lat.Ly())
        throw std::invalid_argument("invalid q-mode interval, beta or grid index");
    std::vector<BondType> types;
    std::set<char> seen;
    for(char c:bondTypes) {
        if((c!='x'&&c!='y'&&c!='z') || !seen.insert(c).second)throw std::invalid_argument("q-mode bonds must be distinct x/y/z");
        types.push_back(c=='x'?BondType::X:c=='y'?BondType::Y:BondType::Z);
    }
    if(types.empty())throw std::invalid_argument("empty q-mode bond list");
    const double u=double(n1)/lat.Lx(),v=double(n2)/lat.Ly()-double(lat.shear())*n1/(double(lat.Lx())*lat.Ly());
    auto integral=[](double x){return std::abs(x-std::round(x))<1e-12;};
    const bool zero=integral(u)&&integral(v), self=integral(2*u)&&integral(2*v);
    if(every && zero)throw std::invalid_argument("q-mode requires nonzero momentum");
    QModeOptions q; q.every=every;q.n1=n1;q.n2=n2;q.beta=beta;q.bond_types=bondTypes;
    if(zero)return q;
    q.basis.assign(types.size()*(self?1:2),std::vector<double>(lat.bonds().size(),0.));
    const double scale=std::sqrt((self?1.:2.)/lat.cells()),twopi=2*std::acos(-1.);
    for(std::size_t t=0;t<types.size();++t)for(const auto &b:lat.bonds())if(b.type==types[t]) {
        const double angle=twopi*(u*b.cell.x+v*b.cell.y);
        const std::size_t offset=t*(self?1:2);
        q.basis[offset][b.index]=scale*std::cos(angle);
        if(!self)q.basis[offset+1][b.index]=scale*std::sin(angle);
    }
    return q;
}
// File entries are integer n1 n2; q and -q span the same real subspace.
// Duplicate subspaces are rejected instead of silently biasing pool frequencies.
inline std::vector<QModeOptions> readQModePool(std::istream &in,const HoneycombLattice &lat,
                                              int every,double beta,const std::string &bondTypes="z") {
    std::vector<QModeOptions> pool;std::string line;int lineNumber=0;
    while(std::getline(in,line)) {
        ++lineNumber;line=line.substr(0,line.find('#'));std::istringstream row(line);
        row>>std::ws;if(row.eof())continue;
        int n1,n2;std::string extra;
        if(!(row>>n1>>n2) || row>>extra)throw std::invalid_argument("q pool line "+std::to_string(lineNumber)+" must be n1 n2");
        auto q=makeQModeOptions(lat,n1,n2,every,beta,bondTypes);
        if(q.basis.empty())throw std::invalid_argument("Gamma forbidden in q pool");
        auto coords=[&](const QModeOptions &a){return std::pair<double,double>{double(a.n1)/lat.Lx(),
            double(a.n2)/lat.Ly()-double(lat.shear())*a.n1/(double(lat.Lx())*lat.Ly())};};
        const auto uv=coords(q);auto integral=[](double a){return std::abs(a-std::round(a))<1e-12;};
        for(const auto &old:pool) {
            const auto ab=coords(old);
            if((integral(uv.first-ab.first)&&integral(uv.second-ab.second)) ||
               (integral(uv.first+ab.first)&&integral(uv.second+ab.second)))
                throw std::invalid_argument("duplicate q/-q subspace in q pool");
        }
        pool.push_back(std::move(q));
    }
    if(in.bad() || pool.empty())throw std::invalid_argument("unreadable or empty q pool");
    return pool;
}
using QModeAmplitudes=std::vector<std::vector<double>>;
inline QModeAmplitudes qModeAmplitudes(const PhononPath &path,const QModeOptions &q) {
    QModeAmplitudes a(q.basis.size(),std::vector<double>(path.nodes(),0.));
    for(std::size_t k=0;k<q.basis.size();++k) {
        if(q.basis[k].size()!=std::size_t(path.bonds()))throw std::invalid_argument("q-mode basis bond count");
        for(int j=0;j<path.nodes();++j)for(int b=0;b<path.bonds();++b)
            if(q.basis[k][b]!=0.)a[k][j]+=q.basis[k][b]*(path.coordinate(j,b)-path.freeMean(j,b));
    }
    return a;
}
inline PhononPath qModeCandidate(const PhononPath &path,const QModeOptions &q,const QModeAmplitudes &a) {
    const auto old=qModeAmplitudes(path,q);
    if(a.size()!=old.size())throw std::invalid_argument("q-mode amplitude count");
    std::vector<std::vector<double>> x; for(int j=0;j<path.nodes();++j)x.push_back(path.node(j));
    for(std::size_t k=0;k<a.size();++k) {
        if(a[k].size()!=std::size_t(path.nodes()))throw std::invalid_argument("q-mode amplitude node count");
        for(int j=0;j<path.nodes();++j) {
            if(!std::isfinite(a[k][j]))throw std::invalid_argument("nonfinite q-mode amplitude");
            for(int b=0;b<path.bonds();++b)if(q.basis[k][b]!=0.)x[j][b]+=q.basis[k][b]*(a[k][j]-old[k][j]);
        }
    }
    return path.withCoordinates(x);
}
// Normalized density on the selected subspace (the complement is unchanged).
inline double qModeTransitionLogDensity(const PhononPath &path,const QModeAmplitudes &from,const QModeAmplitudes &to,double beta) {
    if(!(beta>0 && beta<=1) || from.size()!=to.size())throw std::invalid_argument("invalid q-mode density");
    const double rho=std::sqrt((1-beta)*(1+beta)),variance=beta*beta/(2*path.alpha()),noise=variance*path.noiseFactor();
    auto logNormal=[](double x,double v){return -.5*(std::log(2*std::acos(-1.)*v)+x*x/v);};
    double log=0.;
    for(std::size_t k=0;k<from.size();++k) {
        if(from[k].size()!=std::size_t(path.nodes()) || to[k].size()!=from[k].size())throw std::invalid_argument("q-mode density shape");
        double previous=to[k][0]-rho*from[k][0];log+=logNormal(previous,variance);
        for(int j=1;j<path.nodes();++j) {const double y=to[k][j]-rho*from[k][j];log+=logNormal(y-path.decay()*previous,noise);previous=y;}
    }
    if(!std::isfinite(log))throw std::runtime_error("nonfinite q-mode density");
    return log;
}
}
#endif
