#include "stage1_dense.h"
#include "../ordered_cdw_20260909/odd_sector_ed.h"
#include <iostream>
using namespace driven_stage1;
void close(DataType a, DataType b, const char *message) {
    require(std::abs(a-b) < 1e-9 * std::max(1.,std::abs(b)),message);
}
int main() { try {
    for(int L:{4,6,18,98}) for(int sub:{0,1}) {
        ModelParameters p;p.L=L;
        const int parity=(L/2)%2 ? -1:1;
        const auto trial=makeCdwInitialState(p,sub,parity);
        require(trial.fermionParity()==parity,"initial parity");
        for(double V:{4.,6.}) {
            const auto o=observe(p,V,trial.G_T);
            close(o.energy/L,-V/4,"initial E/L");
            close(o.S_pi,.25,"initial S_pi");close(o.S_pi_dq,0,"initial S_pi_dq");
            close(o.R_CDW,1,"initial R");
            close(staggeredDensity(trial.G_T),sub?-.5:.5,"initial m_pi");
        }
        bool rejected=false;
        try{makeCdwInitialState(p,sub,-parity);}catch(const std::exception &){rejected=true;}
        require(rejected,"wrong parity not rejected");
        if(L>6)continue;
        const auto gamma=gammaMatrices(L);
        int mask=0;for(int j=sub;j<L;j+=2)mask|=1<<j;
        cVecType fock=cVecType::Zero(1<<L);fock(mask)=1;
        require(relativeError(denseStateGreen(fock,gamma),trial.G_T)<1e-12,"CDW Fock covariance");
        p.dt=.05;p.hs=0;p.analytic_hs_factors=true;
        for(double V:{4.,6.}) {
            const auto protocol=PureImaginaryTimeProtocol::constant(V,.1,p.dt);
            // Sum all HS fields independently of the Gaussian walker.
            close((denseHsSameContourState(p,protocol,gamma,fock)-
                   denseSameContourState(p,protocol,gamma,fock)).norm(),0,"CDW HS identity");
            if(L==6) {
                Model ed(L);Vec x(ed.N);x[mask>>1]=1;
                ed.trotter(x,V,p.dt);ed.trotter(x,V,p.dt);
                auto phase=[&](unsigned n){int count=0;for(int j=0;j<L;j+=2)count+=(n>>j)&1;
                    return std::pow(DataType(0,1),count);};
                cVecType mapped=cVecType::Zero(1<<L);
                for(int i=0;i<ed.N;++i){unsigned n=(i<<1)|((__builtin_popcount(unsigned(i))&1)^1);
                    mapped(n)=x[i]*phase(n)/phase(mask);}
                close((mapped-denseSameContourState(p,protocol,gamma,fock)).norm(),0,
                    "full odd-sector ED/CDW dense evolution mismatch");
            }
            std::mt19937_64 rng(990920);
            auto cfg=makeMirroredDrivenContour(p,protocol,rng).configuration;
            PureFastOptions local;local.stabilized_local_updates=true;
            PureFastOptions strict;
            PureProjectorFastWalker a(trial,cfg,16,PureFastRunMode::FastStrict,local,
                PureFastInitializationPolicy::MirroredTheoremZ2Plus);
            PureProjectorFastWalker b(trial,cfg,16,PureFastRunMode::FastStrict,strict,
                PureFastInitializationPolicy::MirroredTheoremZ2Plus);
            auto before=denseContour(p,protocol,trial,fock,cfg);
            require(before.weight.real()>0,"mirrored initial sign");
            auto indices=proposalIndices(cfg);
            std::uniform_real_distribution<double> uniform(0,1);
            // Each branch and both sweep directions; accepted and rejected flips.
            for(int sweep=0;sweep<2;++sweep) {
                if(sweep)std::reverse(indices.begin(),indices.end());
                for(int index:indices) {
                    const auto proposal=flipProposal(p,protocol,cfg,index,uniform(rng));
                    auto candidate=cfg;
                    candidate.slices[index].matrix=proposal.new_factor;
                    candidate.slices[index].eta=proposal.new_eta;
                    candidate.hs_fields[index]=proposal.new_hs;
                    const auto after=denseContour(p,protocol,trial,fock,candidate);
                    const DataType ratio=after.weight/before.weight;
                    const auto ar=a.propose(proposal),br=b.propose(proposal);
                    require(ar.ratio.ok()&&br.ratio.ok()&&!ar.terminated&&!br.terminated,"CDW proposal failed");
                    close(ar.ratio.ratio,ratio,"local/Fock ratio");
                    close(br.ratio.ratio,ratio,"strict/Fock ratio");
                    const bool accepted=proposal.uniform<std::min(1.,std::abs(ratio));
                    require(ar.accepted==accepted&&br.accepted==accepted,"CDW Metropolis mismatch");
                    if(accepted){cfg=std::move(candidate);before=after;}
                    require(a.z2Sign()==(before.weight.real()>0?1:-1)&&a.z2Sign()==b.z2Sign(),"CDW sign mismatch");
                    const auto ga=a.measurementGreen(),gb=b.measurementGreen();
                    require(ga.ok()&&gb.ok()&&relativeError(ga.green,before.green)<1e-9&&
                        relativeError(gb.green,before.green)<1e-9,"CDW center Green mismatch");
                    const auto oa=observe(p,V,ga.green),od=observe(p,V,before.green);
                    close(oa.energy,od.energy,"CDW energy");close(oa.S_pi,od.S_pi,"CDW S_pi");
                    close(oa.S_pi_dq,od.S_pi_dq,"CDW S_pi_dq");close(oa.R_CDW,od.R_CDW,"CDW R");
                    close(staggeredDensity(ga.green),staggeredDensity(before.green),"CDW m_pi");
                }
            }
        }
    }
    std::cout<<"PASS CDW exact states L4/6/18/98, both translations, parity rejection; L4/6 HS sum, local/strict/Fock ratios, decisions, signs, Green and observables; L6 odd-sector ED/dense evolution\n";
    return 0;
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;} }
