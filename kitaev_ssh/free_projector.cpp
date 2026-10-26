#include "physical_matter_trial.h"
#include <iomanip>
#include <iostream>
#include <map>
#include <set>

using namespace kitaev_ssh;

int main(int argc, char **argv) {
    try {
        std::map<std::string,std::string> args;
        const std::set<std::string> allowed{"Lx","Ly","shear","Jx","Jy","Jz","tau-per-branch","theta","dt","block","matter-trial"};
        for (int n=1; n<argc; ++n) {
            const std::string key=argv[n];
            if (key=="--help") {
                std::cout<<"Fixed zero-flux, Wilson (+1,+1), g=0 matter projector (no MC).\n"
                    <<"--Lx 4 --Ly 4 --shear 0 --Jx 1 --Jy 1 --Jz 0.6 --tau-per-branch 0.4 --dt 0.05 --block 4\n"
                    <<"Legacy --theta is accepted as an alias for --tau-per-branch.\n";
                std::cout<<"--matter-trial physical|unrestricted (default physical); no full spin projection.\n";
                return 0;
            }
            if (key.rfind("--",0)!=0 || !allowed.count(key.substr(2)) || n+1>=argc || args.count(key.substr(2)))
                throw std::invalid_argument("unknown, duplicate or incomplete argument: "+key);
            args[key.substr(2)]=argv[++n];
        }
        auto real=[&](const std::string &key, double fallback) {
            if (!args.count(key)) return fallback;
            std::size_t end=0; const double value=std::stod(args.at(key),&end);
            if (end!=args.at(key).size() || !std::isfinite(value)) throw std::invalid_argument("invalid --"+key);
            return value;
        };
        auto integer=[&](const std::string &key, int fallback) {
            const double value=real(key,fallback);
            if (value!=std::floor(value) || value<std::numeric_limits<int>::min() || value>std::numeric_limits<int>::max())
                throw std::invalid_argument("invalid integer --"+key);
            return int(value);
        };
        const int lx=integer("Lx",4),ly=integer("Ly",4),shear=integer("shear",0);
        const HoneycombLattice lattice(lx,ly,shear);
        Couplings j; j.Jx=real("Jx",1.); j.Jy=real("Jy",1.); j.Jz=real("Jz",.6);
        const KitaevSSHModel model(lattice,j);
        const auto selection=selectFreeMatterTrial(model,args.count("matter-trial")?args.at("matter-trial"):"physical");
        const auto &trial=selection.trial;
        if (args.count("theta") && args.count("tau-per-branch"))
            throw std::invalid_argument("use only --tau-per-branch (legacy --theta is an alias)");
        const double tauPerBranch=args.count("tau-per-branch")?real("tau-per-branch",.4):real("theta",.4);
        const double dt=real("dt",.05);
        const auto result=projectStatic(model,trial,tauPerBranch,dt,integer("block",4));
        DataType trialB[3]={}, projectedB[3]={};
        for (const auto &bond:lattice.bonds()) {
            const int mu=bond.type==BondType::X?0:(bond.type==BondType::Y?1:2);
            const DataType factor=DataType(0,-1)*double(model.gauge()[bond.index])/double(lattice.cells());
            trialB[mu]+=factor*trial.G_T(bond.a,bond.b);
            projectedB[mu]+=factor*result.green(bond.a,bond.b);
        }
        std::cout<<std::setprecision(17)<<"{\"status\":\"complete\",\"scope\":\"fixed_gauge_matter_sector\","
            <<"\"Lx\":"<<lattice.Lx()<<",\"Ly\":"<<lattice.Ly()<<",\"shear\":"<<shear
            <<",\"spatial_boundary\":\"pbc\""
            <<",\"pbc_geometry\":\""<<(shear==0?"axis_aligned":"sheared")<<"\",\"geometry_id\":\""<<lattice.geometryId()<<"\""
            <<",\"torus_matrix\":[["<<lx<<','<<shear<<"],[0,"<<ly<<"]]"
            <<",\"period_vectors_cartesian\":[["<<lattice.periodVector(0).x<<','<<lattice.periodVector(0).y
            <<"],["<<lattice.periodVector(1).x<<','<<lattice.periodVector(1).y<<"]]"
            <<",\"momentum_grid_convention\":\"n1/Lx, n2/Ly-shear*n1/(Lx*Ly), index=n1*Ly+n2\""
            <<",\"sites\":"<<lattice.sites()
            <<",\"majoranas\":"<<model.majoranas()<<",\"bonds\":"<<lattice.bondCount()
            <<",\"Jx\":"<<j.Jx<<",\"Jy\":"<<j.Jy<<",\"Jz\":"<<j.Jz<<",\"g\":0"
            <<",\"plaquette_flux\":1,\"wilson_loops\":[1,1],\"directed_winding_products\":["
            <<lattice.gaugeProduct(lattice.windingPath(0),model.gauge())<<','<<lattice.gaugeProduct(lattice.windingPath(1),model.gauge())<<']'
            <<",\"matter_parity\":"<<result.matter_parity<<",\"matter_trial_sector\":";
        writeMatterTrialSelection(std::cout,selection);
        std::cout<<",\"tau_per_branch\":"<<tauPerBranch<<",\"dt\":"<<dt
            <<",\"energy_total\":"<<result.energy<<",\"energy_per_site\":"<<result.energy/lattice.sites()
            <<",\"log_norm_squared\":"<<result.log_norm_squared<<",\"green_residual\":"<<result.green_residual
            <<",\"solve_residual\":"<<result.solve_residual<<",\"overlap_rcond\":"<<result.overlap_rcond
            <<",\"trial_energy_total\":"<<matterEnergy(model,trial.G_T).real()
            <<",\"trial_to_projected_green_max_error\":"<<(trial.G_T-result.green).cwiseAbs().maxCoeff()
            <<",\"bond_operator\":\"u*i*c_A*c_B=-i*u*G_code(A,B)\",\"bond_order\":[\"x\",\"y\",\"z\"]"
            <<",\"trial_bond_mean\":[";
        for(int mu=0;mu<3;++mu) { if(mu)std::cout<<',';std::cout<<'['<<trialB[mu].real()<<','<<trialB[mu].imag()<<']'; }
        std::cout<<"],\"projected_bond_mean\":[";
        for(int mu=0;mu<3;++mu) { if(mu)std::cout<<',';std::cout<<'['<<projectedB[mu].real()<<','<<projectedB[mu].imag()<<']'; }
        std::cout<<"],\"phonon_sampling\":false,\"physical_spin_projection\":false}\n";
        return 0;
    } catch (const GaussianTrialStateError &e) {
        std::cerr<<"kitaev_ssh_free: "<<gaussianTrialStateStatusName(e.code())<<": "<<e.what()
            <<"; model, gauge and PBC were not altered\n"; return 1;
    } catch (const std::exception &e) { std::cerr<<"kitaev_ssh_free: "<<e.what()<<'\n'; return 1; }
}
