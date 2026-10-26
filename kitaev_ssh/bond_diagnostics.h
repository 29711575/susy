#ifndef KITAEV_SSH_BOND_DIAGNOSTICS_H
#define KITAEV_SSH_BOND_DIAGNOSTICS_H
#include "phonon_mc.h"
namespace kitaev_ssh {
// Optional validation output only. K_b=i*c_A*c_B, WITHOUT u_b. The gauge
// belongs to H_M=sum u_b*(J_mu+g*X_b)*K_b, not to the definition of K_b.
inline std::vector<std::string> bondDiagnosticNames(const KitaevSSHModel &m) {
 std::vector<std::string> n;
 for(const auto &b:m.lattice().bonds()) {n.push_back("BondK_b"+std::to_string(b.index));n.push_back("BondXK_b"+std::to_string(b.index));}
 for(const auto *stem:{"BondK_","BondXK_"})for(char mu:std::string("xyz"))n.push_back(std::string(stem)+mu);
 return n;
}
inline std::vector<Complex> bondDiagnostics(const KitaevSSHModel &m,const CoupledPhononPath::Snapshot &s) {
 const auto &G=s.weight.center_green;
 if(!s.weight.green_valid || G.rows()!=m.majoranas())throw std::invalid_argument("bond diagnostics require measured center Green");
 std::vector<Complex> n;std::array<Complex,3> k{},xk{};
 for(const auto &b:m.lattice().bonds()) {
  const Complex z=Complex(0,-1)*G(b.a,b.b),xz=s.path.coordinate(s.path.center(),b.index)*z;
  n.push_back(z);n.push_back(xz);k[int(b.type)]+=z/double(m.lattice().cells());xk[int(b.type)]+=xz/double(m.lattice().cells());
 }
 for(auto z:k)n.push_back(z);for(auto z:xk)n.push_back(z);
 return n;
}
}
#endif
