#include "vbs_observables.h"
#include "vbs_candidates.h"
#include "phonon_mc.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>

using namespace kitaev_ssh;
namespace {
std::map<std::string,double> errors;
// O(1000) floating operations on O(1) fields: 5e-12 absolute bound,
// conservative roundoff allowance, not a statistical acceptance threshold.
void check(const std::string &name,double error,double tolerance=5e-12) {
    errors[name]=std::max(errors[name],error);
    if(!std::isfinite(error) || error>tolerance)throw std::runtime_error(name+" failed");
}
void require(bool condition,const std::string &name) {
    if(!condition)throw std::runtime_error(name+" failed");
}
std::vector<VBSQPoint> parseText(const std::string &text,const HoneycombLattice &lattice) {
    std::istringstream in(text);
    return parseVBSQPoints(in,lattice);
}
// Independent direct cell sum at a q point (test-local reimplementation of the
// finite-window estimator; the all-mode FFT provides the second reference).
Complex directF(const HoneycombLattice &lattice,const VBSQPoint &q,const std::vector<double> &x,int mu) {
    Complex sum=0.;
    for(int r=0;r<lattice.cells();++r) {
        const auto c=lattice.cellCoordinate(r);
        sum+=std::polar(1.,-2*pi*(q.phase_u*c.x+q.phase_v*c.y))*x[3*r+mu];
    }
    return sum;
}

// Build only the Szz entries used by the categorical candidate diagnostics.
// The rest of the selected ratio vector stays zero, which makes these tests
// independent of the nonlinear eigensolver and of QMC sampling.
std::vector<Complex> candidateRatios(const VBSObservables &vbs,
                                     const std::vector<Complex> &raw,
                                     const std::vector<Complex> &projected) {
    require(raw.size()==std::size_t(vbs.nq()) && projected.size()==raw.size(),
            "candidate_test_ratio_dimension");
    std::vector<Complex> ratios(vbs.selectedStride(),0.);
    const double nc=vbs.lattice().cells();
    for(int j=0;j<vbs.nq();++j) {
        ratios[VBSObservables::selS(j)+8]=raw[j]*nc;
        ratios[VBSObservables::selPS(j,vbs.nq())+8]=projected[j]*nc;
    }
    return ratios;
}

std::vector<CandidateReplica> repeatedCandidateReplicas(const std::vector<Complex> &ratios,int count) {
    CandidateReplica replica; replica.ratios=ratios; replica.valid=true;
    return std::vector<CandidateReplica>(std::size_t(count),replica);
}

CandidateDiagnostics repeatedCandidateDiagnostics(const VBSObservables &vbs,
                                                  const std::vector<Complex> &ratios,
                                                  int replicas=4) {
    return runCandidateDiagnostics(vbs,ratios,0,
                                   repeatedCandidateReplicas(ratios,replicas),true);
}

const CandidateGroupDiag &candidateGroup(const CandidateDiagnostics &diag,const std::string &name) {
    for(const auto &group:diag.groups) if(group.name==name) return group;
    throw std::runtime_error("missing candidate group "+name);
}

int candidatePoint(const CandidateDiagnostics &diag,const std::string &label) {
    for(std::size_t j=0;j<diag.points.size();++j)
        if(diag.points[j].label==label) return int(j);
    throw std::runtime_error("missing candidate point "+label);
}

std::vector<std::string> sortedLabels(std::vector<std::string> labels) {
    std::sort(labels.begin(),labels.end());
    return labels;
}

std::vector<std::string> classLabels(const CandidateDiagnostics &diag,
                                     const CandidateDiagnostics::PlusMinusClass &klass) {
    std::vector<std::string> labels;
    for(int j:klass.members) labels.push_back(diag.points[j].label);
    return sortedLabels(std::move(labels));
}

std::vector<std::vector<std::string>> allClassLabels(const CandidateDiagnostics &diag) {
    std::vector<std::vector<std::string>> result;
    for(const auto &klass:diag.classes) result.push_back(classLabels(diag,klass));
    std::sort(result.begin(),result.end());
    return result;
}

// Compare by physical labels, not input indices. This is the invariant that
// a reordered q-file must preserve; display indices are intentionally ignored.
void compareCandidateOrder(const CandidateDiagnostics &a,const CandidateDiagnostics &b) {
    require(sortedLabels([&] { std::vector<std::string> out; for(const auto &p:a.points) out.push_back(p.label); return out; }())
                == sortedLabels([&] { std::vector<std::string> out; for(const auto &p:b.points) out.push_back(p.label); return out; }()),
            "candidate_order_point_labels");
    for(const auto &ga:a.groups) {
        const auto &gb=candidateGroup(b,ga.name);
        check("candidate_order_"+ga.name+"_max",std::abs(ga.max_mz2_real-gb.max_mz2_real));
        require(ga.available==gb.available && ga.center_defined==gb.center_defined,
                "candidate_order_"+ga.name+"_availability");
        require(sortedLabels(ga.winner_labels)==sortedLabels(gb.winner_labels),
                "candidate_order_"+ga.name+"_winner_labels");
        for(std::size_t t=0;t<ga.participant_labels.size();++t) {
            const auto label=ga.participant_labels[t];
            auto it=std::find(gb.participant_labels.begin(),gb.participant_labels.end(),label);
            require(it!=gb.participant_labels.end(),"candidate_order_participant_mapping");
            const std::size_t u=std::size_t(it-gb.participant_labels.begin());
            check("candidate_order_unique_frequency",
                  std::abs(ga.unique_win_frequency[t]-gb.unique_win_frequency[u]));
            check("candidate_order_tie_frequency",
                  std::abs(ga.tie_membership_frequency[t]-gb.tie_membership_frequency[u]));
        }
    }
    require(allClassLabels(a)==allClassLabels(b),"candidate_order_pm_classes");
    for(const auto &pa:a.points) {
        const int j=candidatePoint(b,pa.label);
        const auto &pb=b.points[j];
        check("candidate_order_raw_score",std::abs(pa.mz2_raw_real-pb.mz2_raw_real));
        check("candidate_order_projected_score",std::abs(pa.mz2_proj_real-pb.mz2_proj_real));
        require(pa.class_self_conjugate==pb.class_self_conjugate,
                "candidate_order_self_conjugate");
    }
}

PhaseBlock candidateBlock(Complex phase,const std::vector<Complex> &ratios) {
    PhaseBlock block;
    block.count=1; block.phase=phase; block.numerator.resize(ratios.size());
    for(std::size_t j=0;j<ratios.size();++j) block.numerator[j]=phase*ratios[j];
    return block;
}

Complex weightedEntry(const std::vector<PhaseBlock> &blocks,std::size_t index,int omitted=-1) {
    Complex phase=0.,numerator=0.;
    for(std::size_t b=0;b<blocks.size();++b) if(int(b)!=omitted) {
        phase+=blocks[b].phase; numerator+=blocks[b].numerator[index];
    }
    return numerator/phase;
}

double scalarJackknifeRealError(const std::vector<Complex> &replicas) {
    require(replicas.size()>=2,"candidate_test_jackknife_replica_count");
    Complex mean=0.;
    for(const auto z:replicas) mean+=z/double(replicas.size());
    double variance=0.;
    for(const auto z:replicas) variance+=std::pow(z.real()-mean.real(),2);
    return std::sqrt((replicas.size()-1.)/replicas.size()*variance);
}

// 1) PBC grid consistency on non-square lattices: the selected list enumerates
//    every grid momentum; raw F/Delta/S and derived mz2 must agree with the
//    all-mode FFT values and the direct sums; both keeps the all block bitwise.
void gridConsistency() {
    const HoneycombLattice lattice(5,6);
    const int nc=lattice.cells(),nq=nc;
    std::vector<VBSQPoint> qs;
    for(int nx=0;nx<lattice.Lx();++nx)for(int ny=0;ny<lattice.Ly();++ny)
        qs.push_back(makeQPoint("q_"+std::to_string(nx)+"_"+std::to_string(ny),
                                double(nx)/lattice.Lx(),double(ny)/lattice.Ly(),lattice));
    const VBSObservables all(lattice),sel(lattice,VBSMode::Selected,qs),
        both(lattice,VBSMode::Both,qs);
    std::mt19937_64 rng(910);std::normal_distribution<double> normal;
    std::vector<double> x(3*nc);
    for(auto &a:x)a=normal(rng);
    const auto valuesAll=all.measure(x);
    const auto direct=all.directFourier(x);
    const auto valuesSel=sel.measure(x);
    const auto valuesBoth=both.measure(x);
    for(int j=0;j<nq;++j) {
        const int q=qs[j].n1*lattice.Ly()+qs[j].n2;
        for(int mu=0;mu<3;++mu) {
            // all-mode layout: S_q at 9*q, Delta_q at 9*Nc+3*q.
            const Complex fSel=valuesSel[sel.selD(j,nq)+mu]*double(nc);
            check("grid_selected_all_Delta",std::abs(fSel/double(nc)-valuesAll[9*nc+3*q+mu]));
            check("grid_selected_direct_Delta",std::abs(fSel/double(nc)-direct[q][mu]/double(nc)));
            check("grid_direct_all_Delta",std::abs(direct[q][mu]/double(nc)-valuesAll[9*nc+3*q+mu]));
            check("grid_selected_direct_F",std::abs(fSel-direct[q][mu]));
            for(int nu=0;nu<3;++nu)
                check("grid_selected_all_S",std::abs(valuesSel[sel.selS(j)+3*mu+nu]-valuesAll[9*q+3*mu+nu]));
        }
    }
    // Derived mz2/leading eigenvalue from single-configuration ratios: a raw
    // measurement with unit phase IS the ratio vector derive() consumes.
    const auto derivedAll=all.derive(valuesAll);
    const auto derivedSel=sel.derive(valuesSel);
    for(int j=0;j<nq;++j) {
        const int q=qs[j].n1*lattice.Ly()+qs[j].n2;
        check("grid_selected_all_mz2",std::abs(derivedSel[4*j+1]-derivedAll[4*q+1]));
        check("grid_selected_all_leading_eigenvalue",std::abs(derivedSel[4*j+2]-derivedAll[4*q+2]));
    }
    for(std::size_t i=0;i<std::size_t(12*nc);++i)
        if(valuesBoth[i]!=valuesAll[i])throw std::runtime_error("both_all_prefix_bitwise failed");
    errors["both_all_prefix_bitwise"]=0.;
    // Also on a second non-square shape with a uniform selected subset (Nq != Nc).
    const HoneycombLattice lattice2(4,6);
    const VBSObservables all2(lattice2);
    std::vector<VBSQPoint> qs2;
    for(int nx:{0,1,3})for(int ny:{0,2,5})
        qs2.push_back(makeQPoint("p_"+std::to_string(nx)+"_"+std::to_string(ny),
                                 double(nx)/4,double(ny)/6,lattice2));
    const VBSObservables sel2(lattice2,VBSMode::Selected,qs2);
    const int nc2=lattice2.cells();
    std::vector<double> y(3*nc2);
    for(auto &a:y)a=normal(rng);
    const auto valuesAll2=all2.measure(y),valuesSel2=sel2.measure(y);
    const auto direct2=all2.directFourier(y);
    for(int j=0;j<sel2.nq();++j) {
        const int q=qs2[j].n1*lattice2.Ly()+qs2[j].n2;
        for(int mu=0;mu<3;++mu) {
            check("grid46_selected_all_Delta",
                  std::abs(valuesSel2[9*sel2.nq()+3*j+mu]-valuesAll2[9*nc2+3*q+mu]));
            check("grid46_selected_direct_F",
                  std::abs(valuesSel2[9*sel2.nq()+3*j+mu]*double(nc2)-direct2[q][mu]));
        }
    }
}

// 2) Uniform field: raw F(Q)=Nc*v_mu*W_L(Q); the per-configuration uniform
//    projection annihilates every mode including Q=0.
void uniformField() {
    const HoneycombLattice lattice(4,6);
    const int nc=lattice.cells();
    const std::string list="z 0 0\n"
        "g 0.25 0.5\n"
        "g2 0.75 0.1666666666666667\n"
        "off 0.13 0.71\n"
        "off2 -0.31 2.47\n";
    const auto qs=parseText(list,lattice);
    const VBSObservables sel(lattice,VBSMode::Selected,qs);
    const double v[3]={1.,-.4,.7};
    std::vector<double> x(3*nc);
    for(int r=0;r<nc;++r)for(int mu=0;mu<3;++mu)x[3*r+mu]=v[mu];
    const auto out=sel.measure(x);
    const int nq=sel.nq();
    for(int j=0;j<nq;++j) {
        for(int mu=0;mu<3;++mu) {
            const Complex expected=v[mu]*qs[j].window; // F/Nc
            check("uniform_raw_Delta",std::abs(out[sel.selD(j,nq)+mu]-expected));
            check("uniform_projected_Delta",std::abs(out[sel.selPD(j,nq)+mu]));
            for(int nu=0;nu<3;++nu)
                check("uniform_raw_S",std::abs(out[sel.selS(j)+3*mu+nu]
                    -nc*v[mu]*std::conj(v[nu])*std::norm(qs[j].window)));
        }
    }
    // The Q=0 raw block is the uniform background used by the projection.
    int zero=-1;
    for(int j=0;j<nq;++j)if(qs[j].on_grid && qs[j].n1==0 && qs[j].n2==0)zero=j;
    require(zero>=0,"uniform_zero_point_present");
    for(int mu=0;mu<3;++mu)
        check("uniform_background_Delta0",std::abs(out[sel.selBg(nq)+mu]-v[mu]));
}

// 3) Known modulation X=X0+A*cos(2pi(u0*x+v0*y)+phi): the measured F must
//    match both an independent direct sum and the closed-form window sum
//    F(Q)=Nc[X0*W(Q)+(A/2)(e^{i phi} W(Q-Q0)+e^{-i phi} W(Q+Q0))].
//    On-grid Q0 with X0=0 has the exact amplitude Nc*A/2; off-grid momenta are
//    deliberately NOT held to that amplitude (window weight and +/-Q0 leakage).
void knownModulation() {
    const HoneycombLattice lattice(5,6);
    const int nc=lattice.cells();
    const double X0=0.,A=.8,phi=.31;
    struct Case { double u0,v0; bool onGrid; const char *name; };
    const Case cases[]={{1./5,2./6,true,"grid_q0"},{0.31,0.47,false,"offgrid_q0"}};
    const std::string list="zero 0 0\n"
        "on 0.4 0.5\n"
        "off 0.13 0.71\n";
    for(const auto &mc:cases) {
        const auto qs=parseText(list,lattice);
        const VBSObservables sel(lattice,VBSMode::Selected,qs);
        std::vector<double> x(3*nc);
        for(int r=0;r<nc;++r) {
            const auto c=lattice.cellCoordinate(r);
            const double value=X0+A*std::cos(2*pi*(mc.u0*c.x+mc.v0*c.y)+phi);
            for(int mu=0;mu<3;++mu)x[3*r+mu]=value;
        }
        const auto out=sel.measure(x);
        for(int j=0;j<sel.nq();++j) {
            const auto &q=qs[j];
            const Complex wM=windowFactor(lattice,q.phase_u-mc.u0,q.phase_v-mc.v0);
            const Complex wP=windowFactor(lattice,q.phase_u+mc.u0,q.phase_v+mc.v0);
            const Complex analytic=double(nc)*(X0*q.window
                +(A/2.)*(std::polar(1.,phi)*wM+std::polar(1.,-phi)*wP));
            for(int mu=0;mu<3;++mu) {
                const Complex measured=out[sel.selD(j,sel.nq())+mu]*double(nc);
                check(std::string("modulation_")+mc.name+"_direct",
                      std::abs(measured-directF(lattice,q,x,mu)));
                check(std::string("modulation_")+mc.name+"_analytic",
                      std::abs(measured-analytic));
            }
        }
        // On-grid Q0, X0=0: W(Q0)=0 (nonzero grid point) and W(0)=1, so the
        // only surviving term is the +Q0 sideband: F(Q0)=Nc*(A/2)*e^{i phi}.
        if(mc.onGrid) {
            const auto q0=makeQPoint("exact",mc.u0,mc.v0,lattice);
            require(q0.on_grid,"modulation_grid_q0_on_grid");
            const Complex measured=directF(lattice,q0,x,0);
            check("modulation_grid_amplitude",std::abs(measured-double(nc)*(A/2.)*std::polar(1.,phi)));
        }
        // Off-grid Q0: intentional absence of an A/2 amplitude assertion; the
        // window sum above is the only valid reference there.
    }
}

// 4) Real fields: F(-Q)=F(Q)*, hence S(-Q)=F*(F*)^dagger=conj(S(Q))=S(Q)^T
//    (S(Q) itself is Hermitian), and in particular Szz(-Q)=Szz(Q).
//    Covers one on-grid and one off-grid momentum via minusQPoint.
void conjugationSymmetry() {
    const HoneycombLattice lattice(4,5);
    const int nc=lattice.cells();
    const auto plus=parseText("on 0.25 0.5\noff 0.17 0.63\n",lattice);
    std::vector<VBSQPoint> qs=plus;
    for(const auto &q:plus)qs.push_back(minusQPoint(q,lattice));
    const VBSObservables sel(lattice,VBSMode::Selected,qs);
    std::mt19937_64 rng(7712);std::normal_distribution<double> normal;
    std::vector<double> x(3*nc);
    for(auto &a:x)a=normal(rng);
    const auto out=sel.measure(x);
    const int nq=sel.nq();
    for(int j=0;j<int(plus.size());++j) {
        const int m=j+int(plus.size());
        require(qs[m].label==qs[j].label+"_minus","minus_label");
        for(int mu=0;mu<3;++mu) {
            const Complex fP=out[sel.selD(j,nq)+mu],fM=out[sel.selD(m,nq)+mu];
            check("minusQ_F_conjugate",std::abs(fM-std::conj(fP)));
            for(int nu=0;nu<3;++nu)
                check("minusQ_S_transpose",std::abs(out[sel.selS(m)+3*mu+nu]
                    -out[sel.selS(j)+3*nu+mu]));
        }
        check("minusQ_Szz_equal",std::abs(out[sel.selS(m)+3*2+2]-out[sel.selS(j)+3*2+2]));
        check("minusQ_window_conjugate",std::abs(qs[m].window-std::conj(qs[j].window)));
    }
}

// 5) Layout/statistics: name-measure-derive dimension consistency for the
//    three modes with Nq != Nc; complex-phase joint delete-block jackknife on
//    the selected block against direct ratio references; q-file error cases.
void layoutAndStatistics() {
    const HoneycombLattice lattice(2,3);
    const int nc=lattice.cells();
    const auto qs=parseText("a 0 0\nb 0.5 0.5\nc 0.13 0.71\n",lattice);
    const int nq=int(qs.size());
    require(nq!=nc,"Nq_differs_from_Nc");
    const VBSObservables all(lattice),sel(lattice,VBSMode::Selected,qs),both(lattice,VBSMode::Both,qs);
    std::vector<double> x(3*nc);
    std::mt19937_64 rng(2049);std::normal_distribution<double> normal;
    for(auto &a:x)a=normal(rng);
    check("layout_all_names",std::abs(double(all.names().size())-12*nc));
    check("layout_sel_names",std::abs(double(sel.names().size())-(24*nq+3)));
    check("layout_both_names",std::abs(double(both.names().size())-(12*nc+24*nq+3)));
    check("layout_all_measure",std::abs(double(all.measure(x).size())-all.names().size()));
    check("layout_sel_measure",std::abs(double(sel.measure(x).size())-sel.names().size()));
    check("layout_both_measure",std::abs(double(both.measure(x).size())-both.names().size()));
    check("layout_all_derive",std::abs(double(all.derive(all.measure(x)).size())-22*nc));
    check("layout_sel_derive",std::abs(double(sel.derive(sel.measure(x)).size())-(8*nq+3)));
    check("layout_both_derive",std::abs(double(both.derive(both.measure(x)).size())-(22*nc+8*nq+3)));
    const auto bothNames=both.names(),allNames=all.names(),selNames=sel.names();
    for(std::size_t j=0;j<bothNames.size();++j)
        if(bothNames[j]!=(j<allNames.size()?allNames[j]:selNames[j-allNames.size()]))
            throw std::runtime_error("both_name_order failed");
    errors["both_name_order"]=0.;

    // Complex phases through the joint jackknife, selected block.
    PhaseStatistics stats(sel.names(),4);
    std::vector<Complex> phases,delta0zz,s0zz;
    for(int t=0;t<64;++t) {
        for(int j=0;j<3*nc;++j)x[j]=.4*std::sin(.71*j+.13*t)+.02*t;
        const auto phase=std::polar(1.,.8*std::sin(.29*t));
        const auto m=sel.measure(x);
        phases.push_back(phase);
        delta0zz.push_back(m[sel.selD(0,nq)+2]);        // raw Dsel_s0_z = F_z/Nc
        s0zz.push_back(m[sel.selS(0)+3*2+2]);           // raw Ssel_s0_zz = |F_z|^2/Nc
        stats.add(phase,m);
    }
    const auto result=summarizeDerivedBlocks(stats.blocks(),[&](const auto &r){return sel.derive(r);});
    require(result.phase_resolved && result.jackknife_defined,"selected_jackknife_defined");
    // Raw observable: leave-one-block-out ratio of the phase-weighted mean.
    auto rawReference=[&](int omitted) {
        Complex phase=0.,num=0.;
        for(int t=0;t<64;++t)if(t/4!=omitted) { phase+=phases[t];num+=phases[t]*delta0zz[t]; }
        return num/phase;
    };
    const auto rawSummary=summarizeBlocks(stats.blocks());
    const auto &rawEstimate=rawSummary.observables[9*nq+2];
    check("selected_complex_raw_ratio",std::abs(rawEstimate.value-rawReference(-1)));
    std::vector<Complex> replicas;Complex mean=0.;
    for(int b=0;b<16;++b){replicas.push_back(rawReference(b));mean+=replicas.back()/16.;}
    double vr=0.,vi=0.,cov=0.;
    for(auto z:replicas){z-=mean;vr+=z.real()*z.real()*15/16;vi+=z.imag()*z.imag()*15/16;cov+=z.real()*z.imag()*15/16;}
    check("selected_jackknife_real_error",std::abs(rawEstimate.error_real-std::sqrt(vr)));
    check("selected_jackknife_imag_error",std::abs(rawEstimate.error_imag-std::sqrt(vi)));
    check("selected_jackknife_covariance",std::abs(rawEstimate.covariance-cov));
    // Derived mz2_s0 = Ssel_s0_zz/Nc from the same joint replicas.
    auto mz2Reference=[&](int omitted) {
        Complex phase=0.,num=0.;
        for(int t=0;t<64;++t)if(t/4!=omitted) { phase+=phases[t];num+=phases[t]*s0zz[t]; }
        return num/(phase*double(nc));
    };
    check("selected_derived_mz2",std::abs(result.observables[1].value-mz2Reference(-1)));
    const auto coarse=summarizeDerivedBlocks(mergeAdjacentBlocks(stats.blocks()),
                                             [&](const auto &r){return sel.derive(r);});
    for(std::size_t j=0;j<result.observables.size();++j)
        check("selected_blocking_same_total_ratio",
              std::abs(result.observables[j].value-coarse.observables[j].value));
    // Candidate machinery sees exactly the two non-zero participants here.
    check("selected_candidate_count",std::abs(result.observables[8*nq+2].value-2.));

    // q-file parsing: valid forms and every required failure mode.
    require(parseText("# comment\n\n// also comment\np 0.5, 0.5 \nq\t0.25\t0\n",lattice).size()==2,
            "parse_whitespace_and_comments");
    struct ExpectThrow { const char *text; const char *name; };
    const ExpectThrow bad[]={
        {"","parse_empty"},
        {"\n   \n# only comments\n","parse_only_comments"},
        {"a 0.1 0.2\nb oops 0.3\n","parse_non_numeric"},
        {"a 0.1\n","parse_missing_field"},
        {"a 0.1 0.2 9\n","parse_extra_field"},
        {"a nan 0.2\n","parse_nan"},
        {"a inf 0.2\n","parse_inf"},
        {"bad!label 0.1 0.2\n","parse_bad_label"},
        {"a 0.1 0.2\na 0.3 0.4\n","parse_duplicate_label"},
        {"a 0.25 0\nb 1.25 0\n","parse_duplicate_momentum"},
        {"a 0.25 0\nb -0.75 0\n","parse_duplicate_momentum_negative"},
        {"a 0.5 0.5\nb 1.5 -0.5\n","parse_duplicate_momentum_diagonal"},
    };
    for(const auto &b:bad) {
        bool threw=false;
        try { parseText(b.text,lattice); } catch(const std::invalid_argument &) { threw=true; }
        require(threw,b.name);
    }
    errors["parse_error_cases"]=0.;
}

// A large uniform background can make an off-grid raw finite-window value
// larger than a weak genuine grid modulation. The groups must expose that
// distinction instead of silently calling the off-grid raw winner iVBS.
void candidateBackgroundClassification() {
    const HoneycombLattice lattice(4,4);
    const int nc=lattice.cells();
    const auto qs=parseText("zero 0 0\ngrid 0.5 0.5\nleak 0.03 0.04\n",lattice);
    const VBSObservables sel(lattice,VBSMode::Selected,qs);
    std::vector<double> x(3*nc);
    for(int r=0;r<nc;++r) {
        const auto c=lattice.cellCoordinate(r);
        const double value=8.0+0.1*std::cos(2*pi*(.5*c.x+.5*c.y));
        for(int mu=0;mu<3;++mu)x[3*r+mu]=mu==2?value:0.;
    }
    const auto ratios=sel.measure(x);
    const auto diag=repeatedCandidateDiagnostics(sel,ratios);
    const auto &grid=candidateGroup(diag,"grid_raw");
    const auto &raw=candidateGroup(diag,"offgrid_raw");
    const auto &projected=candidateGroup(diag,"offgrid_projected");
    require(grid.participant_labels==std::vector<std::string>{"grid"},
            "candidate_background_grid_partition");
    require(raw.participant_labels==std::vector<std::string>{"leak"}
                && projected.participant_labels==raw.participant_labels,
            "candidate_background_offgrid_partition");
    require(raw.max_mz2_real>grid.max_mz2_real,
            "candidate_background_legacy_mixed_raw_dominance");
    require(grid.winner_labels==std::vector<std::string>{"grid"},
            "candidate_background_grid_winner");
    require(raw.winner_labels==std::vector<std::string>{"leak"},
            "candidate_background_offgrid_raw_winner");
    const int leak=candidatePoint(diag,"leak");
    require(diag.points[leak].group=="offgrid" && !diag.points[leak].on_grid,
            "candidate_background_offgrid_not_grid");
    require(diag.points[leak].mz2_proj_real<diag.points[leak].mz2_raw_real,
            "candidate_background_projection_reduces_leakage");
    require(diag.points[leak].W_abs2>0.1,"candidate_background_window_signal");
}

// Reordering a q-file changes display indices but cannot change physical
// labels, group maxima, ±Q classes, or resampling frequencies.
void candidateOrderInvariance() {
    const HoneycombLattice lattice(4,4);
    const auto qsA=parseText("zero 0 0\nqp 0.25 0.25\nqm -0.25 -0.25\n"
                            "off 0.13 0.41\noffm -0.13 -0.41\n",lattice);
    const auto qsB=parseText("offm -0.13 -0.41\nzero 0 0\nqm -0.25 -0.25\n"
                            "off 0.13 0.41\nqp 0.25 0.25\n",lattice);
    const VBSObservables a(lattice,VBSMode::Selected,qsA),b(lattice,VBSMode::Selected,qsB);
    std::mt19937_64 rng(318);std::normal_distribution<double> normal;
    std::vector<double> x(3*lattice.cells());
    for(auto &value:x)value=normal(rng);
    const auto da=repeatedCandidateDiagnostics(a,a.measure(x));
    const auto db=repeatedCandidateDiagnostics(b,b.measure(x));
    compareCandidateOrder(da,db);
}

// Exact ties, self-conjugate points, and deliberately switching replicas are
// categorical outcomes, not a continuous error bar on an integer index.
void candidateTiesAndClasses() {
    const HoneycombLattice lattice(4,4);
    const auto qs=parseText("qp 0.25 0.25\nqm -0.25 -0.25\nself 0.5 0.5\nzero 0 0\n",lattice);
    const VBSObservables sel(lattice,VBSMode::Selected,qs);
    const std::vector<Complex> center=candidateRatios(
        sel,{Complex(2.),Complex(2.),Complex(.5),Complex(0.)},
        {Complex(2.),Complex(2.),Complex(.5),Complex(0.)});
    const auto tied=repeatedCandidateDiagnostics(sel,center);
    const auto &grid=candidateGroup(tied,"grid_raw");
    require(grid.winner_set==std::vector<int>{0,1},"candidate_exact_tie_winner_set");
    require(grid.unique_win_counts[0]==0 && grid.unique_win_counts[1]==0,
            "candidate_exact_tie_unique_count");
    check("candidate_exact_tie_membership_frequency_qp",std::abs(grid.tie_membership_frequency[0]-1.));
    check("candidate_exact_tie_membership_frequency_qm",std::abs(grid.tie_membership_frequency[1]-1.));
    require(tied.points[0].plus_minus_class==tied.points[1].plus_minus_class,
            "candidate_pm_pair_class");
    require(!tied.classes[tied.points[0].plus_minus_class].self_conjugate
                && tied.classes[tied.points[0].plus_minus_class].members.size()==2,
            "candidate_pm_pair_not_self_conjugate");
    require(tied.points[2].class_self_conjugate
                && tied.classes[tied.points[2].plus_minus_class].members.size()==1,
            "candidate_self_conjugate_class");

    std::vector<CandidateReplica> switching;
    const double values[][4]={{2.1,1.9,.5,0.},{1.9,2.1,.5,0.},
                              {2.1,1.9,.5,0.},{1.9,2.1,.5,0.}};
    for(const auto &row:values) {
        CandidateReplica replica;
        replica.ratios=candidateRatios(sel,
            {Complex(row[0]),Complex(row[1]),Complex(row[2]),Complex(row[3])},
            {Complex(row[0]),Complex(row[1]),Complex(row[2]),Complex(row[3])});
        replica.valid=true; switching.push_back(std::move(replica));
    }
    const auto switched=runCandidateDiagnostics(sel,center,0,switching,true);
    const auto &sg=candidateGroup(switched,"grid_raw");
    require(sg.stability_defined,"candidate_switching_stability_defined");
    check("candidate_switching_qp_unique_frequency",std::abs(sg.unique_win_frequency[0]-.5));
    check("candidate_switching_qm_unique_frequency",std::abs(sg.unique_win_frequency[1]-.5));
    check("candidate_switching_qp_tie_frequency",std::abs(sg.tie_membership_frequency[0]-.5));
    check("candidate_switching_qm_tie_frequency",std::abs(sg.tie_membership_frequency[1]-.5));
}

// Empty groups are unavailable/null semantics, not zero-valued order
// parameters. This covers the three useful q-list shapes explicitly.
void candidateEmptyGroups() {
    const HoneycombLattice lattice(4,4);
    struct Case { const char *text; bool grid,offgrid; };
    const Case cases[]={{"zero 0 0\ngrid 0.5 0.5\n",true,false},
                        {"zero 0 0\noff 0.13 0.41\n",false,true},
                        {"zero 0 0\n",false,false}};
    for(const auto &test:cases) {
        const auto qs=parseText(test.text,lattice);
        const VBSObservables sel(lattice,VBSMode::Selected,qs);
        const std::vector<Complex> ratios(sel.selectedStride(),0.);
        const auto diag=runCandidateDiagnostics(sel,ratios,0,{},true);
        const auto &grid=candidateGroup(diag,"grid_raw");
        const auto &raw=candidateGroup(diag,"offgrid_raw");
        const auto &projected=candidateGroup(diag,"offgrid_projected");
        require(grid.available==test.grid && raw.available==test.offgrid
                    && projected.available==test.offgrid,
                "candidate_empty_group_availability");
        for(const auto *group:{&grid,&raw,&projected}) if(!group->available) {
            require(group->participants.empty() && group->winner_set.empty()
                        && !group->center_defined && std::isnan(group->max_mz2_real),
                    "candidate_empty_group_null_semantics");
        }
        for(const auto &point:diag.points) if(point.group=="zero_excluded")
            require(point.plus_minus_class<0,"candidate_zero_excluded_from_pm_classes");
    }
}

// Hand-check complex phase reweighting, negative real scores, group maxima,
// winner replicas, and one exactly invalid delete-one-block denominator.
void candidateComplexReweighting() {
    const HoneycombLattice lattice(4,4);
    const auto qs=parseText("q1 0.25 0.25\nq2 0.5 0.5\nzero 0 0\n",lattice);
    const VBSObservables sel(lattice,VBSMode::Selected,qs);
    // The last three unit phases sum to zero, so omitting the first block
    // gives one exactly invalid denominator while the other replicas remain
    // valid. None of the valid phases is itself the total phase.
    const Complex phases[]={1.,std::polar(1.,pi/6.),std::polar(1.,5.*pi/6.),
                            std::polar(1.,3.*pi/2.)};
    const std::vector<std::vector<Complex>> values={
        {Complex(-1.),Complex(10.),Complex(0.)},
        {Complex(-1.),Complex(10.),Complex(0.)},
        {Complex(-1.),Complex(10.),Complex(0.)},
        {Complex(0.),Complex(10.),Complex(0.)}};
    std::vector<PhaseBlock> blocks;
    for(std::size_t b=0;b<values.size();++b) {
        const auto ratios=candidateRatios(sel,values[b],
                                          {Complex(0.),Complex(0.),Complex(0.)});
        blocks.push_back(candidateBlock(phases[b],ratios));
    }
    const std::vector<Complex> fullRaw={
        weightedEntry(blocks,VBSObservables::selS(0)+8)/double(lattice.cells()),
        weightedEntry(blocks,VBSObservables::selS(1)+8)/double(lattice.cells()),
        Complex(0.)};
    const auto fullRatios=candidateRatios(sel,fullRaw,
                                          {Complex(0.),Complex(0.),Complex(0.)});
    const auto replicas=candidateReplicas(blocks);
    const auto diag=runCandidateDiagnostics(sel,fullRatios,0,replicas,true);
    const int q1=candidatePoint(diag,"q1"),q2=candidatePoint(diag,"q2");
    const Complex expectedQ1=weightedEntry(blocks,VBSObservables::selS(0)+8)
                             /double(lattice.cells());
    const Complex expectedQ2=weightedEntry(blocks,VBSObservables::selS(1)+8)
                             /double(lattice.cells());
    check("candidate_complex_q1_center",
          std::abs(diag.points[q1].mz2_raw_real-expectedQ1.real()));
    check("candidate_complex_q1_imaginary_diagnostic",
          std::abs(diag.points[q1].mz2_raw_imag-expectedQ1.imag()));
    check("candidate_complex_q2_center",
          std::abs(diag.points[q2].mz2_raw_real-expectedQ2.real()));
    require(expectedQ1.real()<0. && diag.points[q1].mz2_raw_real<0.,
            "candidate_negative_score_preserved");
    const auto &grid=candidateGroup(diag,"grid_raw");
    require(grid.winner_labels==std::vector<std::string>{"q2"},
            "candidate_complex_reweighted_winner");
    require(diag.valid_replicas==3 && diag.invalid_replicas==1
                && grid.valid_replicas==3 && grid.invalid_replicas==1,
            "candidate_invalid_replica_count");
    require(!grid.jackknife_defined && !diag.points[q1].raw_jackknife_defined
                && !diag.points[q2].raw_jackknife_defined,
            "candidate_invalid_replica_jackknife_undefined");
    require(std::isnan(diag.points[q1].mz2_raw_error_real)
                && std::isnan(diag.points[q1].mz2_raw_covariance)
                && std::isnan(diag.points[q2].mz2_raw_error_real),
            "candidate_invalid_replica_errors_null");
    check("candidate_invalid_replica_frequency_denominator",
          std::abs(grid.unique_win_frequency[1]-1.));
    check("candidate_invalid_replica_nonwinner_frequency",
          std::abs(grid.unique_win_frequency[0]));
}

// Candidate errors use the original B delete-one-block replicas. A partial
// replica set can still support a conditional winner diagnostic, but it can
// never be relabeled as a complete ordinary jackknife.
void candidateReplicaValidityStates() {
    const HoneycombLattice lattice(4,4);
    const auto qs=parseText("grid 0.5 0.5\nzero 0 0\n",lattice);
    const VBSObservables sel(lattice,VBSMode::Selected,qs);
    std::vector<CandidateReplica> all;
    for(double value:{1.,2.,3.,4.}) {
        CandidateReplica replica;
        replica.ratios=candidateRatios(sel,{Complex(value),Complex(0.)},
                                       {Complex(value),Complex(0.)});
        replica.valid=true; all.push_back(std::move(replica));
    }
    const auto center=candidateRatios(sel,{Complex(2.5),Complex(0.)},
                                      {Complex(2.5),Complex(0.)});
    const auto complete=runCandidateDiagnostics(sel,center,0,all,true);
    const auto &completeGroup=candidateGroup(complete,"grid_raw");
    require(complete.total_replicas==4 && complete.valid_replicas==4
                && complete.invalid_replicas==0 && completeGroup.jackknife_defined
                && completeGroup.stability_defined,
            "candidate_complete_replica_states");
    require(!completeGroup.conditional_on_valid_replicas
                && std::isfinite(complete.points[0].mz2_raw_error_real)
                && std::isfinite(complete.points[0].mz2_raw_covariance),
            "candidate_complete_replica_errors_defined");
    check("candidate_complete_replica_error",
          std::abs(complete.points[0].mz2_raw_error_real-
                   scalarJackknifeRealError({Complex(1.),Complex(2.),Complex(3.),Complex(4.)})));

    auto partial=all;
    partial[0].valid=false; partial[0].invalid_reason="test_invalid_replica";
    const auto partialDiag=runCandidateDiagnostics(sel,center,0,partial,true);
    const auto &partialGroup=candidateGroup(partialDiag,"grid_raw");
    require(partialDiag.total_replicas==4 && partialDiag.valid_replicas==3
                && partialDiag.invalid_replicas==1 && !partialGroup.jackknife_defined
                && partialGroup.stability_defined && partialGroup.conditional_on_valid_replicas,
            "candidate_partial_replica_states");
    require(std::isnan(partialDiag.points[0].mz2_raw_error_real)
                && std::isnan(partialDiag.points[0].mz2_raw_covariance)
                && std::abs(partialGroup.unique_win_frequency[0]-1.)<1e-12,
            "candidate_partial_replica_null_error_conditional_frequency");

    const auto none=runCandidateDiagnostics(sel,center,0,{},true);
    const auto &noneGroup=candidateGroup(none,"grid_raw");
    require(none.total_replicas==0 && none.valid_replicas==0 && none.invalid_replicas==0
                && !noneGroup.stability_defined && std::isnan(noneGroup.unique_win_frequency[0])
                && std::isnan(none.points[0].mz2_raw_error_real),
            "candidate_zero_replica_unavailable");

    const auto one=runCandidateDiagnostics(sel,center,0,std::vector<CandidateReplica>{all[0]},true);
    require(one.valid_replicas==1 && !candidateGroup(one,"grid_raw").stability_defined
                && std::isnan(candidateGroup(one,"grid_raw").unique_win_frequency[0]),
            "candidate_one_replica_unavailable");
    const auto two=runCandidateDiagnostics(sel,center,0,
        std::vector<CandidateReplica>{all[0],all[1]},true);
    require(two.valid_replicas==2 && candidateGroup(two,"grid_raw").stability_defined
                && candidateGroup(two,"grid_raw").jackknife_defined,
            "candidate_two_replica_states");

    const auto noCenter=runCandidateDiagnostics(sel,center,0,all,false);
    const auto &noCenterGroup=candidateGroup(noCenter,"grid_raw");
    require(!noCenter.center_defined && !noCenterGroup.center_defined
                && noCenterGroup.winner_set.empty() && !noCenterGroup.stability_defined
                && std::isnan(noCenterGroup.unique_win_frequency[0])
                && std::isnan(noCenter.points[0].mz2_raw_error_real),
            "candidate_undefined_center_suppresses_selection");

    auto nonfinite=center;
    nonfinite[VBSObservables::selS(0)+8]=Complex(std::numeric_limits<double>::quiet_NaN(),0.);
    const auto nonfiniteDiag=runCandidateDiagnostics(sel,nonfinite,0,all,true);
    require(!candidateGroup(nonfiniteDiag,"grid_raw").center_defined
                && !nonfiniteDiag.points[0].raw_center_defined
                && std::isnan(nonfiniteDiag.points[0].mz2_raw_error_real),
            "candidate_nonfinite_center_is_unavailable");
    require(candidateReplicas({}).empty(),"candidate_empty_blocks_are_unavailable");
    PhaseBlock badBlock;
    badBlock.count=1; badBlock.phase=1.; badBlock.numerator.assign(sel.selectedStride(),0.);
    badBlock.numerator[0]=Complex(std::numeric_limits<double>::quiet_NaN(),0.);
    const auto badReplicas=candidateReplicas(std::vector<PhaseBlock>{badBlock,badBlock});
    require(badReplicas.size()==2 && !badReplicas[0].valid && !badReplicas[1].valid
                && badReplicas[0].invalid_reason=="nonfinite_block_data",
            "candidate_nonfinite_block_is_invalid");
}

// S(-Q)=S(Q)^T is a linear transpose relation and survives the same complex
// reweighting. The conjugation/Hermitian relation uses conjugated weights and
// can have a visible finite-sample residual; no symmetry-forcing post-process
// is allowed to hide it.
void finiteSampleTransposeAndConjugation() {
    const HoneycombLattice lattice(2,2);
    const auto qs=parseText("q 0.25 0.25\nqm -0.25 -0.25\n",lattice);
    const VBSObservables sel(lattice,VBSMode::Selected,qs);
    std::vector<Complex> first(sel.selectedStride(),0.),second(sel.selectedStride(),0.);
    const int qzx=VBSObservables::selS(0)+3*2+0;
    const int qxz=VBSObservables::selS(0)+3*0+2;
    const int m_zx=VBSObservables::selS(1)+3*2+0;
    const int m_xz=VBSObservables::selS(1)+3*0+2;
    first[qzx]=1.+Complex(0.,2.); first[qxz]=std::conj(first[qzx]);
    first[m_zx]=std::conj(first[qzx]); first[m_xz]=first[qzx];
    second[qzx]=3.-Complex(0.,1.); second[qxz]=std::conj(second[qzx]);
    second[m_zx]=std::conj(second[qzx]); second[m_xz]=second[qzx];
    const auto blocks=std::vector<PhaseBlock>{candidateBlock(1.,first),
                                                candidateBlock(Complex(0.,1.),second)};
    check("finite_sample_config_transpose_q",
          std::abs(first[m_zx]-first[qxz]));
    check("finite_sample_config_transpose_q2",
          std::abs(second[m_zx]-second[qxz]));
    check("finite_sample_config_transpose_reverse_q",
          std::abs(first[m_xz]-first[qzx]));
    check("finite_sample_config_transpose_reverse_q2",
          std::abs(second[m_xz]-second[qzx]));
    const auto summary=summarizeBlocks(blocks);
    const Complex q=summary.observables[qzx].value;
    const Complex minus=summary.observables[m_zx].value;
    require(std::abs(minus-summary.observables[qxz].value)<1e-12,
            "finite_sample_transpose_relation_preserved");
    require(std::abs(summary.observables[m_xz].value-q)<1e-12,
            "finite_sample_transpose_reverse_preserved");
    require(std::abs(minus-std::conj(q))>1e-6,
            "finite_sample_conjugation_residual_visible");
    errors["finite_sample_conjugation_residual_visible"]=0.;
}

// 6) RNG trajectory regression: with the same seed, running the sampler in
//    all/selected/both measurement modes must leave acceptance counts, center
//    node coordinates, log weight and matter phase bitwise identical, because
//    the measurement layer never touches the RNG or the sampler state.
struct Trajectory {
    MoveCounts counts;
    std::vector<std::vector<double>> center;
    std::vector<double> logWeight;
    std::vector<Complex> phase;
    std::vector<Complex> basic;
};
Trajectory runTrajectory(VBSMode mode,const std::vector<VBSQPoint> &qs) {
    const KitaevSSHModel model(HoneycombLattice(2,2),p2Couplings(.3));
    const PhononParameters phonon(1.7,.8);
    const int warmup=5,measurementSweeps=10,stride=2;
    std::mt19937_64 rng(42);
    const auto initial=PhononPath(model.lattice().bondCount(),4,.2,phonon).freeSample(rng);
    CoupledPhononPath state(initial,freeGroundState(model),honeycombGenerator(model),4,
                            PhononWeightOptions{true,true});
    PhononSampler sampler(state,rng,PhononProposal::LocalConditional);
    const VBSObservables vbs(model.lattice(),mode,mode==VBSMode::All?std::vector<VBSQPoint>{}:qs);
    Trajectory tr;
    for(int i=0;i<warmup;++i)sampler.sweep();
    for(int i=0;i<measurementSweeps;++i) {
        sampler.sweep();
        if((i+1)%stride==0) {
            const auto &s=sampler.state();
            tr.center.push_back(s.path.node(s.path.center()));
            tr.logWeight.push_back(s.weight.log_abs_matter);
            tr.phase.push_back(s.weight.matter_phase);
            const auto b=basicObservables(model,s);
            tr.basic.insert(tr.basic.end(),b.begin(),b.end());
            vbs.measure(s.path.node(s.path.center())); // measured, must not disturb the chain
        }
    }
    tr.counts=sampler.counts();
    return tr;
}
void rngTrajectory() {
    const HoneycombLattice lattice(2,2);
    const int nc=lattice.cells();
    const auto qs=parseText("g 0.5 0.5\noff 0.13 0.71\nzero 0 0\n",lattice);
    const auto all=runTrajectory(VBSMode::All,qs);
    const auto sel=runTrajectory(VBSMode::Selected,qs);
    const auto both=runTrajectory(VBSMode::Both,qs);
    require(all.counts.proposed==sel.counts.proposed && all.counts.accepted==sel.counts.accepted
            && both.counts.proposed==sel.counts.proposed && both.counts.accepted==sel.counts.accepted,
            "rng_acceptance_counts");
    require(all.center==sel.center && all.center==both.center,"rng_center_coordinates");
    require(all.logWeight==sel.logWeight && all.logWeight==both.logWeight,"rng_log_weight");
    require(all.phase==sel.phase && all.phase==both.phase,"rng_matter_phase");
    require(all.basic==sel.basic && all.basic==both.basic,"rng_basic_observables");
    errors["rng_trajectory_bitwise"]=0.;
    // The selected all-grid point must reproduce the all-mode FFT value on the
    // same configuration: the measurement layers agree where they overlap.
    const std::vector<double> x=all.center.front();
    const VBSObservables allV(lattice),selV(lattice,VBSMode::Selected,qs);
    const auto va=allV.measure(x),vs=selV.measure(x);
    const int j=0; // "g 0.5 0.5" is on the L=2 grid (n1=n2=1)
    require(qs[j].on_grid && qs[j].n1==1 && qs[j].n2==1,"rng_grid_point_on_grid");
    const int qidx=qs[j].n1*lattice.Ly()+qs[j].n2;
    for(int mu=0;mu<3;++mu)
        check("rng_selected_matches_all",
              std::abs(vs[9*selV.nq()+3*j+mu]-va[9*nc+3*qidx+mu]));
}
}
int main(int argc,char **argv) {
    try {
        // Machine-readable bridge used by the Python generator end-to-end
        // test.  This deliberately calls the production q-file parser and
        // makeQPoint path rather than reimplementing their rules in Python.
        if(argc==6 && std::string(argv[1])=="--inspect-q-file") {
            auto integer=[](const char *text,const char *name) {
                std::size_t end=0;const std::string value(text);long parsed;
                try{parsed=std::stol(value,&end);}catch(...){throw std::invalid_argument(std::string("invalid ")+name);}
                if(end!=value.size() || parsed<std::numeric_limits<int>::min() ||
                   parsed>std::numeric_limits<int>::max())throw std::invalid_argument(std::string("invalid ")+name);
                return int(parsed);
            };
            const HoneycombLattice lattice(integer(argv[2],"Lx"),integer(argv[3],"Ly"),integer(argv[4],"shear"));
            const auto points=readVBSQFile(argv[5],lattice);
            std::cout<<std::setprecision(17)<<std::boolalpha<<"{\"points\":[";
            for(std::size_t j=0;j<points.size();++j) {
                if(j)std::cout<<',';const auto &q=points[j];
                std::cout<<"{\"index\":"<<j<<",\"label\":\""<<q.label
                         <<"\",\"input_fraction\":["<<q.u_input<<','<<q.v_input
                         <<"],\"reciprocal_fraction\":["<<q.u<<','<<q.v
                         <<"],\"phase_fraction\":["<<q.phase_u<<','<<q.phase_v
                         <<"],\"on_grid\":"<<q.on_grid<<",\"n\":";
                if(q.on_grid)std::cout<<'['<<q.n1<<','<<q.n2<<']';else std::cout<<"null";
                std::cout<<",\"grid_snap_residual\":"<<q.snap_residual<<'}';
            }
            std::cout<<"]}\n";
            return 0;
        }
        if(argc!=1)throw std::invalid_argument(
            "usage: kitaev_ssh_vbs_q_validate [--inspect-q-file Lx Ly shear PATH]");
        gridConsistency();uniformField();knownModulation();conjugationSymmetry();
        layoutAndStatistics();
        candidateBackgroundClassification();candidateOrderInvariance();
        candidateTiesAndClasses();candidateEmptyGroups();candidateComplexReweighting();
        candidateReplicaValidityStates();finiteSampleTransposeAndConjugation();rngTrajectory();
        for(const auto &e:errors)std::cout<<e.first<<" max_error="<<std::setprecision(12)<<e.second<<" tolerance=5e-12\n";
        std::cout<<"vbs_q_validate PASS (synthetic fields are tests, not physical plots).\n";
        return 0;
    }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
