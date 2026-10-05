#ifndef KITAEV_SSH_FERMION_OBSERVABLES_H
#define KITAEV_SSH_FERMION_OBSERVABLES_H
#include "free_projector.h"
#include "phase_statistics.h"
#include <ostream>

namespace kitaev_ssh {
// Equal-time matter Majoranas, not spin-spin or electron correlations.
// Paths use the existing directed-u convention, without (-i)^length.
class FermionObservables {
public:
    struct Channel {std::string name; int dx,dy,alpha,beta; bool string;};
    explicit FermionObservables(const KitaevSSHModel &model):lattice_(model.lattice()) {
        for(bool string:{false,true}) for(int dx=0;dx<lattice_.Lx();++dx)
            for(int dy=0;dy<lattice_.Ly();++dy) for(int a=0;a<2;++a) for(int b=0;b<2;++b) {
                Channel channel{std::string("F_")+(string?"string":"bare")+"_dx"+std::to_string(dx)+
                    "_dy"+std::to_string(dy)+"_"+(a?"B":"A")+(b?"B":"A"),dx,dy,a,b,string};
                channels_.push_back(channel); pairs_.emplace_back();
                for(int x=0;x<lattice_.Lx();++x) for(int y=0;y<lattice_.Ly();++y) {
                    const auto path=pathFor(x,y,channel);
                    const int end=lattice_.pathEnd(path);
                    const auto endpoint=lattice_.translatedCell({x,y},dx,dy);
                    if(end!=lattice_.siteIndex(endpoint.x,endpoint.y,static_cast<Sublattice>(b)))
                        throw std::logic_error("fermion path endpoint mismatch");
                    pairs_.back().push_back({path.start,end,string?lattice_.gaugeProduct(path,model.gauge()):1});
                }
            }
    }
    const auto &channels() const {return channels_;}
    std::vector<std::string> names() const {
        std::vector<std::string> out;for(const auto &c:channels_)out.push_back(c.name);return out;
    }
    // Translate on the starting sublattice along +a1 then +a2, using the
    // nonnegative periodic cell displacement, then take z if alpha != beta.
    BondPath pathFor(int x,int y,const Channel &c) const {
        BondPath path{lattice_.siteIndex(x,y,static_cast<Sublattice>(c.alpha)),{}};
        for(int axis=0;axis<2;++axis) for(int n=0;n<(axis?c.dy:c.dx);++n) {
            const auto next=lattice_.translatedCell({x,y},axis==0?1:0,axis==1?1:0);
            const int nx=next.x,ny=next.y;
            const auto type=axis?BondType::Y:BondType::X;
            if(c.alpha==0) {
                path.steps.push_back({lattice_.bondIndex(x,y,BondType::Z),1});
                path.steps.push_back({lattice_.bondIndex(nx,ny,type),-1});
            } else {
                path.steps.push_back({lattice_.bondIndex(nx,ny,type),-1});
                path.steps.push_back({lattice_.bondIndex(nx,ny,BondType::Z),1});
            }
            x=nx;y=ny;
        }
        if(c.alpha!=c.beta)path.steps.push_back({lattice_.bondIndex(x,y,BondType::Z),c.alpha==0?1:-1});
        return path;
    }
    std::vector<Complex> measure(const MatType &green) const {
        if(green.rows()!=lattice_.sites() || green.cols()!=lattice_.sites() || !green.allFinite())
            throw std::invalid_argument("invalid fermion measurement Green");
        std::vector<Complex> out(channels_.size(),0.);
        for(std::size_t k=0;k<pairs_.size();++k) {
            for(const auto &p:pairs_[k])out[k]+=double(p.sign)*(p.i==p.j?Complex(1.):Complex(-green(p.i,p.j)));
            out[k]/=double(lattice_.cells());
        }
        return out;
    }
    void writeMetadata(std::ostream &out) const {
        out<<"{\"definition\":\"C_ij=delta_ij-G_ij; string=product_directed_u*C_ij\","
            <<"\"normalization\":\"1/Nc sum over source cells, before phase reweighting\","
            <<"\"time\":\"center_equal_time\",\"scope\":\"matter_majorana_not_spin\","
            <<"\"path_rule\":\"nonnegative dx along +a1, then dy along +a2 on source sublattice, then z if sublattices differ; PBC wrap explicit\","
            <<"\"directed_u\":\"u_BA=-u_AB; no (-i)^length factor\",\"channels\":[";
        for(std::size_t k=0;k<channels_.size();++k) {
            if(k)out<<',';const auto &c=channels_[k];
            const auto a1=lattice_.primitive(0),a2=lattice_.primitive(1);
            out<<"{\"name\":\""<<c.name<<"\",\"kind\":\""<<(c.string?"string":"bare")
                <<"\",\"cell_displacement\":["<<c.dx<<','<<c.dy<<"],\"sublattices\":\""<<(c.alpha?"B":"A")<<(c.beta?"B":"A")
                <<"\",\"unwrapped_site_displacement\":["<<c.dx*a1.x+c.dy*a2.x<<','<<c.dx*a1.y+c.dy*a2.y+c.beta-c.alpha<<']';
            if(c.string) {
                out<<",\"origin_zero_directed_bonds\":[";const auto path=pathFor(0,0,c);
                for(std::size_t n=0;n<path.steps.size();++n) {if(n)out<<',';out<<'['<<path.steps[n].bond<<','<<path.steps[n].direction<<']';}
                out<<']';
            }
            out<<'}';
        }
        out<<"]}";
    }
private:
    struct Pair {int i,j,sign;};
    HoneycombLattice lattice_;
    std::vector<Channel> channels_;
    std::vector<std::vector<Pair>> pairs_;
};
}
#endif
