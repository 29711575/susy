#ifndef KITAEV_SSH_VBS_MOMENTA_H
#define KITAEV_SSH_VBS_MOMENTA_H

#include "honeycomb_lattice.h"
#include <algorithm>
#include <complex>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace kitaev_ssh {

// Measurement mode for the valence-bond-solid (VBS) bond-field Fourier
// analysis. NONE of these modes changes the Hamiltonian, the phonon PBC, the
// QMC sampling weight or the random-number stream; they only select which
// Fourier components are evaluated at measurement time.
enum class VBSMode { All, Selected, Both };

inline VBSMode parseVBSMode(const std::string &name) {
    if (name == "all") return VBSMode::All;
    if (name == "selected") return VBSMode::Selected;
    if (name == "both") return VBSMode::Both;
    throw std::invalid_argument("--vbs-mode must be one of: all, selected, both");
}
inline const char *vbsModeName(VBSMode mode) {
    return mode == VBSMode::All ? "all" : mode == VBSMode::Selected ? "selected" : "both";
}

// A single requested measurement wave vector Q = u*b1 + v*b2, with b1,b2 the
// code's own reciprocal vectors (HoneycombLattice::reciprocal). Fractional
// reciprocal coordinates are the primary input; Cartesian coordinates are
// always derived, never read from the paper. Off-grid points never receive a
// fabricated integer (n1,n2) index.
struct VBSQPoint {
    std::string label;
    double u_input = 0., v_input = 0.; // exactly as entered
    double u = 0., v = 0.;             // normalized into [0,1) modulo reciprocal vectors
    double phase_u = 0., phase_v = 0.; // fractions used in exp(-i Q.R); snapped when on grid
    Position cartesian{0., 0.};
    bool on_grid = false;
    int n1 = -1, n2 = -1;  // dual-grid labels, valid only when on_grid
    double snap_residual = 0.; // fractional distance between normalized input and the grid bin
    std::complex<double> window{0., 0.}; // W_L(Q) = (1/Nc) sum_R exp(-i Q.R)
};

struct VBSGridPoint {
    int n1 = 0, n2 = 0, index = 0;
    double u = 0., v = 0.;
    Position cartesian{0., 0.};
};

inline double normalizeReciprocalFraction(double value) {
    value -= std::floor(value);
    if (value >= 1.) value = 0.;
    return value;
}

// The character labels are dual coordinates of the sheared quotient, not two
// independent rectangular fractions.  The returned order is the stable
// storage order index=n1*Ly+n2 used by all-mode observables and JSON output.
inline std::vector<VBSGridPoint> allowedVBSMomenta(const HoneycombLattice &lattice) {
    const auto b1 = lattice.reciprocal(0), b2 = lattice.reciprocal(1);
    std::vector<VBSGridPoint> result;
    result.reserve(lattice.cells());
    for (int n1 = 0; n1 < lattice.Lx(); ++n1) for (int n2 = 0; n2 < lattice.Ly(); ++n2) {
        const double u = double(n1) / lattice.Lx();
        const double v = normalizeReciprocalFraction(
            double(n2) / lattice.Ly() -
            double(lattice.shear()) * double(n1) /
                (double(lattice.Lx()) * lattice.Ly()));
        result.push_back({n1,n2,n1*lattice.Ly()+n2,u,v,
                          {u*b1.x+v*b2.x,u*b1.y+v*b2.y}});
    }
    return result;
}

// Fractional distance on the unit circle, result in [-0.5,0.5).
inline double fractionalDelta(double a, double b) {
    double d = a - b;
    d = d - std::floor(d);      // [0,1)
    if (d > 0.5) d -= 1.;       // shortest signed separation
    return d;
}

// W_L(Q) = (1/Nc) sum_R exp(-i Q.R); the cell sum factorizes over x and y.
// The closed form is the Dirichlet kernel sum_{n=0}^{L-1} e^{-i2pi f n}
//   = e^{-i pi f (L-1)} sin(pi f L)/sin(pi f),
// keeping the sign of the real ratio so the phase is correct off-grid.
inline std::complex<double> windowFactor(const HoneycombLattice &lattice, double u, double v) {
    using Complex = std::complex<double>;
    auto kernel = [](int L, double f) {
        f -= std::floor(f);
        const double s = std::sin(pi * f);
        if (std::abs(s) < 1e-14) return Complex(double(L), 0.);
        const double ratio = std::sin(pi * f * L) / s; // real, may be negative
        return Complex(ratio, 0.) * std::polar(1., -pi * f * (L - 1));
    };
    return kernel(lattice.Lx(), u) * kernel(lattice.Ly(), v) / double(lattice.cells());
}

// Build a validated, normalized q point from label and raw fractions.
inline VBSQPoint makeQPoint(const std::string &label, double uIn, double vIn,
                            const HoneycombLattice &lattice) {
    if (!std::isfinite(uIn) || !std::isfinite(vIn))
        throw std::invalid_argument("VBS q point " + label + ": non-finite fraction");
    static const std::string allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
                                       "0123456789_.+-";
    if (label.empty() || label.find_first_not_of(allowed) != std::string::npos)
        throw std::invalid_argument("VBS q point label must be non-empty and use only "
                                    "[A-Za-z0-9_.+-]: " + label);
    VBSQPoint q;
    q.label = label;
    q.u_input = uIn; q.v_input = vIn;
    q.u = normalizeReciprocalFraction(uIn);
    q.v = normalizeReciprocalFraction(vIn);
    q.phase_u = q.u; q.phase_v = q.v;
    const double tol = 1e-9;
    for (const auto &grid : allowedVBSMomenta(lattice)) {
        const double du = fractionalDelta(q.u, grid.u);
        const double dv = fractionalDelta(q.v, grid.v);
        if (std::abs(du) <= tol && std::abs(dv) <= tol) {
            if (q.on_grid) throw std::logic_error("duplicate grid bin match");
            q.on_grid = true; q.n1 = grid.n1; q.n2 = grid.n2;
            q.snap_residual = std::hypot(du, dv);
            q.phase_u = grid.u; q.phase_v = grid.v; // exact grid fractions for the summation
        }
    }
    const auto b1 = lattice.reciprocal(0), b2 = lattice.reciprocal(1);
    q.cartesian = {q.u * b1.x + q.v * b2.x, q.u * b1.y + q.v * b2.y};
    q.window = windowFactor(lattice, q.phase_u, q.phase_v);
    return q;
}

// Parse the q-list text format. Each non-comment line is:
//   label  u  v
// with Q = u*b1 + v*b2. Blank lines and comments beginning with '#' or '//'
// are allowed. Throws on malformed lines, non-numeric fractions, empty lists
// and duplicate labels / duplicate normalized momenta.
inline std::vector<VBSQPoint> parseVBSQPoints(std::istream &in, const HoneycombLattice &lattice) {
    std::vector<VBSQPoint> points;
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        for (char &c : line) if (c == '\t' || c == ',') c = ' ';
        auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        auto slash = line.find("//");
        if (slash != std::string::npos) line.erase(slash);
        std::istringstream ls(line);
        std::string label, us, vs, extra;
        if (!(ls >> label)) continue; // blank line
        if (!(ls >> us >> vs) || (ls >> extra))
            throw std::invalid_argument("q-file line " + std::to_string(lineNo) +
                                        ": expected 'label u v'");
        auto number = [&](const std::string &token) {
            std::size_t end = 0; double value;
            try { value = std::stod(token, &end); }
            catch (...) { throw std::invalid_argument("q-file line " + std::to_string(lineNo) +
                                                      ": invalid number '" + token + "'"); }
            if (end != token.size() || !std::isfinite(value))
                throw std::invalid_argument("q-file line " + std::to_string(lineNo) +
                                            ": invalid number '" + token + "'");
            return value;
        };
        const double u = number(us), v = number(vs);
        VBSQPoint q = makeQPoint(label, u, v, lattice);
        for (const auto &existing : points) {
            if (existing.label == q.label)
                throw std::invalid_argument("q-file: duplicate label '" + q.label + "'");
            const double du = fractionalDelta(existing.u, q.u);
            const double dv = fractionalDelta(existing.v, q.v);
            if (std::hypot(du, dv) < 1e-9)
                throw std::invalid_argument("q-file: duplicate momentum for labels '" +
                                            existing.label + "' and '" + q.label + "'");
        }
        points.push_back(std::move(q));
    }
    if (points.empty()) throw std::invalid_argument("--vbs-q-file contains no q points");
    return points;
}

inline std::vector<VBSQPoint> readVBSQFile(const std::string &path,
                                          const HoneycombLattice &lattice) {
    std::ifstream file(path);
    if (!file) throw std::invalid_argument("cannot open --vbs-q-file: " + path);
    return parseVBSQPoints(file, lattice);
}

// Precomputed cell phases for a q point: exp(-i Q.R_cell).
inline std::vector<std::complex<double>> qPointPhases(const HoneycombLattice &lattice,
                                                      const VBSQPoint &q) {
    using Complex = std::complex<double>;
    std::vector<Complex> phases(lattice.cells());
    for (int r = 0; r < lattice.cells(); ++r) {
        const auto c = lattice.cellCoordinate(r);
        phases[r] = std::polar(1., -2. * pi * (q.phase_u * c.x + q.phase_v * c.y));
    }
    return phases;
}

// The negated momentum -Q, normalized the same way (label suffixed _minus).
inline VBSQPoint minusQPoint(const VBSQPoint &q, const HoneycombLattice &lattice) {
    // Re-run the actual grid matcher.  For a shear, negating n2 independently
    // is generally wrong even when Q is on the grid.
    VBSQPoint m = makeQPoint(q.label + "_minus", -q.u, -q.v, lattice);
    return m;
}

} // namespace kitaev_ssh
#endif
