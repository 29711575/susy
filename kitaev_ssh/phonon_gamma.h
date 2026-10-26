#ifndef KITAEV_SSH_PHONON_GAMMA_H
#define KITAEV_SSH_PHONON_GAMMA_H
#include "phonon_q_mode.h"
namespace kitaev_ssh {
// Separate opt-in kernel. A real normalized constant basis for ONE bond type;
// no cosine/sine pair, no periodic closure in imaginary time.
struct GammaOptions { int every=0; double beta=.5; std::vector<QModeOptions> modes; };
inline GammaOptions makeGammaOptions(const HoneycombLattice &lat,int every,double beta) {
    if(every<0 || !std::isfinite(beta) || beta<=0 || beta>1)
        throw std::invalid_argument("invalid Gamma cadence/beta");
    GammaOptions g;g.every=every;g.beta=beta;
    for(int mu=0;mu<3;++mu) {
        QModeOptions q;q.every=every;q.beta=beta;q.bond_types=std::string(1,"xyz"[mu]);
        q.basis.assign(1,std::vector<double>(lat.bondCount(),0.));
        for(const auto &b:lat.bonds())if(int(b.type)==mu)q.basis[0][b.index]=1/std::sqrt(double(lat.cells()));
        g.modes.push_back(std::move(q));
    }
    return g;
}
// Coefficients are centered on the projection of path.freeMean(j,b).
// pCN innovation C_ij = exp(-omega*dt*|i-j|)/(2*m*omega) includes BOTH
// independent endpoints and the unique center. No trial/source is changed.
template<class Rng,class Normal>
inline QModeAmplitudes drawProjectedPCN(const PhononPath &p,const QModeOptions &q,Rng &rng,Normal &normal) {
    auto a=qModeAmplitudes(p,q);const double rho=std::sqrt((1-q.beta)*(1+q.beta));
    const double sigma=std::sqrt(1/(2*p.alpha())),noise=sigma*std::sqrt(p.noiseFactor());
    for(auto &row:a) {
        double xi=sigma*normal(rng);row[0]=rho*row[0]+q.beta*xi;
        for(int j=1;j<p.nodes();++j){xi=p.decay()*xi+noise*normal(rng);row[j]=rho*row[j]+q.beta*xi;}
    }
    return a;
}
}
#endif
