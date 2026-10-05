#ifndef KITAEV_SSH_HONEYCOMB_LATTICE_H
#define KITAEV_SSH_HONEYCOMB_LATTICE_H

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace kitaev_ssh {

constexpr double pi = 3.141592653589793238462643383279502884;
enum class Sublattice { A = 0, B = 1 };
enum class BondType { X = 0, Y = 1, Z = 2 };
struct Position { double x, y; };
struct Cell { int x, y; };
struct CellReduction { Cell representative; Cell winding; };
struct SiteCoordinate { int x, y; Sublattice sublattice; };
struct Bond {
    int index, a, b;
    BondType type;
    Cell cell;
    Cell cell_delta;       // unwrapped B-cell minus A-cell
    Cell boundary_wrap;    // unwrapped B = wrapped B + w1*T1 + w2*T2
    Position displacement;
};
struct DirectedBond { int bond; int direction; }; // +1: stored A -> B
struct BondPath { int start; std::vector<DirectedBond> steps; };

class HoneycombLattice {
public:
    HoneycombLattice(int lx, int ly, int shear = 0) : lx_(lx), ly_(ly), shear_(shear) {
        if (lx < 2 || ly < 2 ||
            shear < 0 || shear >= lx ||
            static_cast<long long>(lx) * ly > std::numeric_limits<int>::max() / 3)
            throw std::invalid_argument("honeycomb requires Lx,Ly >= 2, 0<=shear<Lx and representable counts");
        nc_ = lx * ly;
        for (int x = 0; x < lx_; ++x) for (int y = 0; y < ly_; ++y) {
            for (int mu = 0; mu < 3; ++mu) {
                const int dx = mu == 0 ? -1 : 0, dy = mu == 1 ? -1 : 0;
                const auto reduced = reduceCell(static_cast<std::int64_t>(x) + dx,
                                                static_cast<std::int64_t>(y) + dy);
                const int bx = reduced.representative.x, by = reduced.representative.y;
                const auto a1 = primitive(0), a2 = primitive(1);
                const int id = 3 * cellIndex(x, y) + mu;
                bonds_.push_back({id, siteIndex(x,y,Sublattice::A), siteIndex(bx,by,Sublattice::B),
                    static_cast<BondType>(mu), {x,y}, {dx,dy},
                    reduced.winding,
                    {dx*a1.x+dy*a2.x, 1+dx*a1.y+dy*a2.y}});
                by_type_[mu].push_back(id);
            }
        }
        // A(R) -> B(R) -> A(R+a1) -> B(R+a1-a2)
        //      -> A(R+a1-a2) -> B(R-a2) -> A(R).
        for (int x = 0; x < lx_; ++x) for (int y = 0; y < ly_; ++y)
            plaquettes_.push_back({siteIndex(x,y,Sublattice::A), {
                {bondIndex(x,y,BondType::Z),1}, {bondIndexWrapped(x+1,y,BondType::X),-1},
                {bondIndexWrapped(x+1,y,BondType::Y),1}, {bondIndexWrapped(x+1,y-1,BondType::Z),-1},
                {bondIndexWrapped(x+1,y-1,BondType::X),1}, {bondIndexWrapped(x,y,BondType::Y),-1}}});
    }

    int Lx() const { return lx_; }
    int Ly() const { return ly_; }
    int shear() const { return shear_; }
    int cells() const { return nc_; }
    int sites() const { return 2*nc_; }
    int bondCount() const { return 3*nc_; }
    std::array<std::array<int,2>,2> torusMatrix() const {
        return {{{lx_, shear_}, {0, ly_}}};
    }
    Position periodVector(int direction) const {
        if (direction != 0 && direction != 1) throw std::out_of_range("period direction");
        const auto a1 = primitive(0), a2 = primitive(1);
        if (direction == 0) return {lx_*a1.x, lx_*a1.y};
        return {shear_*a1.x + ly_*a2.x, shear_*a1.y + ly_*a2.y};
    }
    std::string geometryId() const {
        return std::string("honeycomb-pbc-Lx") + std::to_string(lx_) +
               "-Ly" + std::to_string(ly_) + "-shear" + std::to_string(shear_);
    }
    // Reduce (x,y) against the column basis T1=(Lx,0), T2=(shear,Ly).
    // The returned integers obey
    // (x,y) = representative + w1*T1 + w2*T2, including negative inputs.
    CellReduction reduceCell(std::int64_t x, std::int64_t y) const {
        const auto yDivision = floorDivide(static_cast<__int128>(y), ly_);
        const auto xTemporary = static_cast<__int128>(x) -
                                yDivision.quotient * static_cast<__int128>(shear_);
        const auto xDivision = floorDivide(xTemporary, lx_);
        return {{checkedInt(xDivision.remainder, "cell x representative"),
                 checkedInt(yDivision.remainder, "cell y representative")},
                {checkedInt(xDivision.quotient, "cell x winding"),
                 checkedInt(yDivision.quotient, "cell y winding")}};
    }
    Cell translatedCell(Cell cell, std::int64_t dx, std::int64_t dy) const {
        cellIndex(cell.x, cell.y);
        return reduceCell(static_cast<std::int64_t>(cell.x) + dx,
                          static_cast<std::int64_t>(cell.y) + dy).representative;
    }
    int cellIndexWrapped(std::int64_t x, std::int64_t y) const {
        const auto c = reduceCell(x,y).representative;
        return cellIndex(c.x,c.y);
    }
    int cellIndex(int x, int y) const {
        if (x < 0 || x >= lx_ || y < 0 || y >= ly_) throw std::out_of_range("cell coordinate");
        return x*ly_+y;
    }
    Cell cellCoordinate(int cell) const {
        if (cell < 0 || cell >= nc_) throw std::out_of_range("cell index");
        return {cell/ly_,cell%ly_};
    }
    int siteIndex(int x, int y, Sublattice sub) const {
        const int s = static_cast<int>(sub);
        if (s != 0 && s != 1) throw std::invalid_argument("sublattice");
        return s*nc_+cellIndex(x,y);
    }
    int siteIndexWrapped(std::int64_t x, std::int64_t y, Sublattice sub) const {
        const auto c = reduceCell(x,y).representative;
        return siteIndex(c.x,c.y,sub);
    }
    SiteCoordinate siteCoordinate(int site) const {
        if (site < 0 || site >= sites()) throw std::out_of_range("site index");
        const auto c = cellCoordinate(site%nc_);
        return {c.x,c.y,static_cast<Sublattice>(site/nc_)};
    }
    static Position primitive(int direction) {
        if (direction != 0 && direction != 1) throw std::out_of_range("primitive direction");
        return {(direction == 0 ? 1 : -1)*std::sqrt(3.)/2,1.5};
    }
    static Position reciprocal(int direction) {
        if (direction != 0 && direction != 1) throw std::out_of_range("reciprocal direction");
        return {(direction == 0 ? 1 : -1)*2*pi/std::sqrt(3.),2*pi/3};
    }
    Position position(int site) const {
        const auto c = siteCoordinate(site);
        return {std::sqrt(3.)*(c.x-c.y)/2,1.5*(c.x+c.y)+(c.sublattice==Sublattice::B ? 1 : 0)};
    }
    int bondIndex(int x, int y, BondType type) const {
        const int mu = static_cast<int>(type);
        if (mu < 0 || mu >= 3) throw std::invalid_argument("bond type");
        return 3*cellIndex(x,y)+mu;
    }
    int bondIndexWrapped(std::int64_t x, std::int64_t y, BondType type) const {
        const auto c = reduceCell(x,y).representative;
        return bondIndex(c.x,c.y,type);
    }
    const Bond &bond(int id) const { return bonds_.at(id); }
    const std::vector<Bond> &bonds() const { return bonds_; }
    const std::vector<int> &bonds(BondType type) const { return by_type_.at(static_cast<int>(type)); }
    const std::vector<BondPath> &plaquettes() const { return plaquettes_; }
    int pathEnd(const BondPath &path) const {
        siteCoordinate(path.start);
        int current = path.start;
        for (const auto &step : path.steps) {
            const auto &b = bond(step.bond);
            if (step.direction != 1 && step.direction != -1) throw std::invalid_argument("path direction must be +/-1");
            if (current != (step.direction == 1 ? b.a : b.b)) throw std::invalid_argument("disconnected bond path");
            current = step.direction == 1 ? b.b : b.a;
        }
        return current;
    }
    void validateGauge(const std::vector<int> &u) const {
        if (u.size() != bonds_.size()) throw std::invalid_argument("gauge bond count");
        for (int v : u) if (v != 1 && v != -1) throw std::invalid_argument("u must be +/-1");
    }
    int gaugeProduct(const BondPath &path, const std::vector<int> &u) const {
        validateGauge(u); pathEnd(path);
        int product = 1;
        for (const auto &step : path.steps) product *= step.direction*u[step.bond];
        return product;
    }
    int plaquetteFlux(int cell, const std::vector<int> &u) const {
        return -gaugeProduct(plaquettes_.at(cell),u);
    }
    BondPath windingPath(int direction) const {
        if (direction != 0 && direction != 1) throw std::out_of_range("winding direction");
        BondPath path{siteIndex(0,0,Sublattice::A),{}};
        int x = 0, y = 0;
        auto appendNegativeStep = [&](BondType type, int dx, int dy) {
            path.steps.push_back({bondIndex(x,y,type),1});
            const auto next = translatedCell({x,y},dx,dy);
            x = next.x; y = next.y;
            path.steps.push_back({bondIndex(x,y,BondType::Z),-1});
        };
        if (direction == 0) {
            for (int n = 0; n < lx_; ++n) appendNegativeStep(BondType::X,-1,0);
        } else {
            // -T2 is represented by -shear*a1 followed by -Ly*a2.  A
            // boundary bond can carry both winding integers; the quotient
            // reduction, rather than separate axis moduli, keeps that data.
            for (int n = 0; n < shear_; ++n) appendNegativeStep(BondType::X,-1,0);
            for (int n = 0; n < ly_; ++n) appendNegativeStep(BondType::Y,0,-1);
        }
        return path;
    }
    int wilsonLoop(int direction, const std::vector<int> &u) const {
        const auto path = windingPath(direction);
        return ((path.steps.size()/2)%2 ? -1 : 1)*gaugeProduct(path,u);
    }
    std::vector<int> zeroFluxGauge(int wx = 1, int wy = 1) const {
        if ((wx != 1 && wx != -1) || (wy != 1 && wy != -1)) throw std::invalid_argument("Wilson sectors must be +/-1");
        std::vector<int> u(bondCount(),1);
        for (const auto &b : bonds_) {
            if (b.boundary_wrap.x % 2 != 0) u[b.index] *= wx;
            if (b.boundary_wrap.y % 2 != 0) u[b.index] *= wy;
        }
        return u;
    }

private:
    struct FloorDivision { __int128 quotient, remainder; };
    static FloorDivision floorDivide(__int128 value, __int128 divisor) {
        if (divisor <= 0) throw std::invalid_argument("floor division requires a positive divisor");
        __int128 quotient = value / divisor, remainder = value % divisor;
        if (remainder < 0) { --quotient; remainder += divisor; }
        return {quotient,remainder};
    }
    static int checkedInt(__int128 value, const char *what) {
        if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
            throw std::overflow_error(std::string(what) + " overflow");
        return static_cast<int>(value);
    }
    int lx_, ly_, shear_, nc_;
    std::vector<Bond> bonds_;
    std::array<std::vector<int>,3> by_type_;
    std::vector<BondPath> plaquettes_;
};

} // namespace kitaev_ssh
#endif
