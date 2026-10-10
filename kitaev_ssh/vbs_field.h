#ifndef KITAEV_SSH_VBS_FIELD_H
#define KITAEV_SSH_VBS_FIELD_H
#include "honeycomb_lattice.h"
#include "phonon_path.h"
#include "phase_statistics.h"
namespace kitaev_ssh {
// Fixed real z-bond projection. h couples to unrescaled X, not g*X.
class VBSField {
public:
    VBSField(const HoneycombLattice &lat,int n1,int n2,double phase,double h)
        : nc_(lat.cells()),n1_(n1),n2_(n2),phase_(phase),h_(h),profile_(lat.bondCount(),0.),cos_(lat.bondCount(),0.),sin_(lat.bondCount(),0.) {
        if(n1<0 || n1>=lat.Lx() || n2<0 || n2>=lat.Ly() || !std::isfinite(phase) || !std::isfinite(h))
            throw std::invalid_argument("invalid allowed VBS source momentum/phase/field");
        u_=double(n1)/lat.Lx();v_=double(n2)/lat.Ly()-double(lat.shear())*n1/(lat.Lx()*lat.Ly());
        for(const auto &b:lat.bonds()) if(b.type==BondType::Z) {
            const double a=2*pi*(u_*b.cell.x+v_*b.cell.y);
            cos_[b.index]=std::cos(a);sin_[b.index]=-std::sin(a);
            profile_[b.index]=std::cos(a+phase);profile2_+=profile_[b.index]*profile_[b.index];
        }
    }
    std::vector<double> forces() const {auto f=profile_;for(auto &x:f)x*=h_;return f;}
    std::vector<std::string> names() const {
        return {"Field_m","Field_score","Field_m_score","Field_m2","Field_m4","Field_absD2","Field_absD4"};
    }
    std::vector<Complex> measure(const PhononPath &p) const {
        double m=0,score=0,dr=0,di=0;
        for(int b=0;b<p.bonds();++b) {
            const double x=p.coordinate(p.center(),b);
            m+=profile_[b]*x;dr+=cos_[b]*x;di+=sin_[b]*x;
            if(profile_[b]!=0.) for(int j=0;j<p.nodes();++j) score+=profile_[b]*p.sourceNodeWeight(j)*p.coordinate(j,b);
        }
        // Path-independent derivative term cancels in Cov(m,score), retained for auditing.
        score+=h_*profile2_/p.parameters().spring()*p.intervals()*(p.dt()-2*p.sourceNodeWeight(0));
        m/=nc_;dr/=nc_;di/=nc_;const double a2=dr*dr+di*di;
        return {m,score,m*score,m*m,m*m*m*m,a2,a2*a2};
    }
    static std::vector<std::string> derivedNames() {return {"m_parallel","chi_parallel","absD2","B4_complex"};}
    static std::vector<Complex> derive(const std::vector<Complex> &r,std::size_t o) {
        return {r.at(o),r.at(o+2)-r.at(o)*r.at(o+1),r.at(o+5),r.at(o+6)/(r.at(o+5)*r.at(o+5))};
    }
    int n1() const{return n1_;} int n2() const{return n2_;}
    double u() const{return u_;} double v() const{return v_;}
    double phase() const{return phase_;} double h() const{return h_;}
    double profileMeanSquare() const{return profile2_/nc_;}
private:
    int nc_,n1_,n2_;double phase_,h_,u_,v_,profile2_=0.;
    std::vector<double> profile_,cos_,sin_;
};
}
#endif
