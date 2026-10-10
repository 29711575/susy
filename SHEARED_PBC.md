# Integer-sheared periodic torus

This tree implements one discrete finite-size geometry parameter, `--shear s`.
It is an integer quotient of the honeycomb cell coordinates. It is not an
arbitrary continuous U(1) twist, and it is not a way to make an arbitrary
irrational momentum exactly commensurate on a finite lattice.

## Scope and boundary convention

The stored representative cells are

```text
0 <= x < Lx,  0 <= y < Ly
cell = x*Ly + y
bond = 3*cell + mu,  mu = x,y,z
```

The primitive vectors are the existing `a1,a2`; the two period vectors are

```text
T1 = Lx*a1
T2 = s*a1 + Ly*a2
```

The integer period matrix has columns equal to the cell-coordinate periods:

```text
M = [[Lx, s],
     [ 0, Ly]]
```

Thus `det(M)=Lx*Ly=Nc`. The CLI accepts `Lx,Ly >= 2` and
`0 <= s < Lx`; invalid or overflowing inputs are rejected. `boundary_wrap`
stores the two integer winding numbers `(w1,w2)` in this basis, rather than
independent axis-modulo flags. A bond's local `displacement` remains the
physical one-step bond vector; it is never replaced by the long endpoint
separation across the seam.

## Quotient reduction and seams

Every cell translation uses one floor-division reduction. For an integer
coordinate `(x,y)`,

```text
w2   = floor(y/Ly)
ybar = y - w2*Ly
xtmp = x - w2*s
w1   = floor(xtmp/Lx)
xbar = xtmp - w1*Lx
```

The result is

```text
(x,y) = (xbar,ybar) + w1*(Lx,0) + w2*(s,Ly).
```

The implementation uses explicit floor division for negative values and
checked wide intermediate arithmetic. Consequently,

```text
(x+Lx,y) is equivalent to (x,y)
(x,y+Ly) is equivalent to (x-s,y)
```

For example, the y bond from `A(0,0)` ends at `B((0+s) mod Lx,Ly-1)`. For
`(Lx,Ly,s)=(12,3,3)`, this is `B(3,2)`, not `B(0,2)`. The same quotient
interface is used for bonds, plaquettes, winding paths, fermion paths,
gauge strings, and real-space correlation translations. It satisfies
`Tx*Ty=Ty*Tx`, `Tx^Lx=I`, and `Ty^Ly=Tx^(-s)`.

The plaquette construction retains the directed-bond convention and uses the
same wrapped bond lookup for all six steps. The four `zeroFluxGauge(wx,wy)`
choices have `wx,wy` in `{ -1, +1 }`; the boundary sign is computed from the
parity of both winding components. All plaquettes have zero flux and the two
constructed Wilson loops have the requested signs. Production still uses
`(+1,+1)`. The validation binary checks all four sectors, but this does not
prove which sector is the global physical-spin ground state.

The QMC still samples a fixed-gauge matter sector. Gauge variables are not
updated, and no full physical-spin projection or all-sector sign theorem is
claimed. At nonzero shear the finite Hamiltonian and QMC weight generally
change because the finite graph seam changes. The existing optimized A3,
trial-state and spectrum objects are constructed per lattice instance; they
must not be shared between different shears.

## Allowed VBS momenta

For `Q=u*b1+v*b2`, periodic characters obey

```text
Lx*u           in Z
s*u + Ly*v     in Z
```

The complete table is therefore

```text
u = n1/Lx
v = n2/Ly - s*n1/(Lx*Ly)
0 <= n1 < Lx,  0 <= n2 < Ly
index = n1*Ly + n2
```

The labels `n1,n2` index the dual character table. For `s != 0`, they are not
the old pair of independently sampled rectangular fractions. `minusQPoint`
normalizes and searches the actual table, which is necessary because
independently replacing `n2` by `-n2` is generally wrong.

For `shear=0`, the table and the all-mode FFT reduce to the existing rectangular
path. For `shear!=0`, all-mode measurement uses

```text
G(n1,y) = sum_x exp(-2*pi*i*n1*x/Lx) X(x,y)
F(n1,n2) = sum_y exp(-2*pi*i*n2*y/Ly)
                    * exp(+2*pi*i*s*n1*y/(Lx*Ly)) G(n1,y)
```

That is an x FFT, a shear phase, and a y FFT. The phase is a Fourier kernel
from the quotient, not a complex hopping or a U(1) boundary condition.
Selected mode always uses direct sums with the same table. The finite-window
quantity is

```text
W_L(Q) = (1/Nc) sum over representative R exp(-i*Q.R)
F_proj(Q) = F(Q) - W_L(Q)*F(0)
```

Every nonzero allowed character has `W_L=0`, but `W_L=0` is not an on-grid
test. In particular, for the bare target
`Q_bare=(0.4195693767448338,0.5804306232551661)`, the `(12,3,3)` torus
has `12*u=5.034832520938005`, so the bare point is still off-grid. Its
nearest allowed point is a separate quantity,
`q_L=(5/12,-5/12)`, normalized to `(5/12,7/12)`, with
`(n1,n2)=(5,0)` and `index=15`; the reported Cartesian mismatch is
`0.0210597347773431`. The bare target is never silently replaced by this
grid point. The rational point `(5/12,7/12)` is off-grid on `(12,3,0)` but
has `W_L` zero by direct representative-cell summation; it is on-grid on
`(12,3,3)` with the same `(5,0)` label. Thus the bare target's nonzero
window and the rational point's zero window must not be conflated.

The structure-factor convention is

```text
F_mu(Q) = sum_R exp(-i*Q.R) X_mu(R)
S_mu,nu(Q) = <F_mu(Q) conj(F_nu(Q))>/Nc
C_mu,nu(d) = (1/Nc) sum_R X_mu(R) X_nu(R+d)
C(d) = sum_Q exp(-i*Q.d) S(Q)/Nc
```

The last line uses the same quotient translation as the lattice. The sign is
matched to the `R+d` definition above. The VBS JSON preserves normalized
fractions, dual labels, Cartesian momenta, windows, geometry identity and the
raw/projected/candidate validity states.

## Fermions and Wilson sectors

Fermion endpoint paths use the same `translatedCell` operation as the lattice.
The path rule remains: take nonnegative steps along `+a1`, then `+a2`, then a
z bond when the endpoint sublattice differs. The directed gauge string uses
the actual directed bond signs, including seam winding. The phonon momentum
grid is independent of the matter Wilson signs; a non-default Wilson sector
must use its own matter spectrum and is not compared to the production sector
as though the finite spectra were identical.

## Prediction and geometry comparison

`predict_vbs_momenta.py` accepts `--shear` and reports separate `bare`,
`effective`, `kekule`, and `user` rows. It uses the Cartesian reciprocal metric
for nearest points and reports both dimensionless `metric_distance` and
physical Cartesian distance. Its mismatch is

```text
DeltaQ = Q_target - q_L
delta1 = DeltaQ.T1 = 2*pi*Lx*Delta_u
delta2 = DeltaQ.T2 = 2*pi*(s*Delta_u + Ly*Delta_v)
```

where `(Delta_u,Delta_v)` is the metric-minimal representative modulo
reciprocal lattice vectors. Bare and effective rows are search candidates,
not interacting peak locations.

The Kekule points `(1/3,-1/3)` and `(-1/3,1/3)` are on-grid precisely when

```text
Lx % 3 == 0 and (s-Ly) % 3 == 0.
```

Therefore `(Lx,Ly,s)=(6,4,1)` is a useful counterexample: `Ly` is not
divisible by 3, but the Kekule points are exactly allowed. The script also
reports the shortest reduced period, shape ratio and the unmodified
`(+1,+1)` zero-mode risk; it does not choose a physical geometry for the
caller and never runs QMC.

For the target above, the checked 36-cell examples are:

```text
(6,6,0):   nearest distance 0.3497005581116156
(12,3,0):  nearest distance 0.3497005581116157
(12,3,3):  nearest distance 0.0210597347773431
```

The improvement factor of about 16.605 is only a reciprocal-space distance
comparison for this bare target. It is not a physical error, order-parameter
accuracy, runtime, or convergence factor. The `(12,3,3)` torus has equivalent
short period basis `(3,3)` and `(6,-6)`, with physical lengths `9` and
`6*sqrt(3)`. The storage ratio `12/3` is not its physical aspect ratio.

## Metadata and commands

Every new result records at least:

```text
spatial_boundary: "pbc"
pbc_geometry: "axis_aligned" or "sheared"
shear: s
geometry_id: honeycomb-pbc-Lx<Lx>-Ly<Ly>-shear<s>
torus_matrix: [[Lx,s],[0,Ly]]       # columns are T1,T2
period_vectors_cartesian: [T1,T2]
momentum_grid_convention: n1/Lx, n2/Ly-shear*n1/(Lx*Ly), index=n1*Ly+n2
twist_applied: false                 # no continuous U(1) phase
```

Typical runs are:

```bash
./build-vbs/kitaev_ssh_mc --Lx 12 --Ly 3 --shear 3 --Jx 1 --Jy 1 --Jz 0.5 \
  --intervals 4 --dt 0.2 --mass 1.7 --omega-ph 0.8 --g 0.3 --seed 42 \
  --warmup 100 --sweeps 400 --measure-every 2 --block-samples 10 \
  --vbs 1 --vbs-mode all --output PREFIX

./build-vbs/kitaev_ssh_mc --Lx 12 --Ly 3 --shear 3 --fermions 1 \
  --intervals 2 --dt 0.2 --warmup 20 --sweeps 80 --measure-every 1 \
  --block-samples 10 --output PREFIX

python3 kitaev_ssh/predict_vbs_momenta.py --Lx 12 --Ly 3 --shear 3 \
  --Jx 1 --Jy 1 --Jz 0.5 --output-dir REPORTS --prefix 12x3_s3

python3 kitaev_ssh/compare_sheared_geometries.py \
  --geometry 6,6,0 --geometry 12,3,0 --geometry 12,3,3 \
  --target 0.4195693767448338,0.5804306232551661 \
  --Jx 1 --Jy 1 --Jz 0.5 --output-dir REPORTS --prefix three_36_cell
```

Historical JSON without shear is accepted as rectangular only when it has the
known old shape. Files carrying partial new geometry metadata are rejected by
the comparison, plotting and candidate-summary readers. Different geometry
IDs cannot be compared as eager/deferred replicas. No script should infer a
sheared `q_index` by applying the old independent `nx,ny` formula.

## Validation boundary

`kitaev_ssh_sheared_pbc_validate` checks reduction, seam bonds, plaquettes,
all Wilson sectors, free spectra, A3 assembly, staged FFT versus direct sums,
inverse correlations, windows, the off-grid/window-zero counterexample and
fermion gauge strings for `(5,3,1)`, `(5,3,4)`, and `(12,3,3)`. The existing
`kitaev_ssh_vbs_q_validate` remains passing for the rectangular measurement
suite. Local `(3,2,1)` chains cover all/selected/both and eager/deferred, with
an additional `--fermions 1` pair; their four-sample statistical state is
`too_few_blocks` by design.

No production-scale sheared chain, thermodynamic extrapolation, projection or
time-step convergence, global Wilson-sector minimization, sign-problem proof,
phase-transition conclusion, or full physical-spin projection is part of this
delivery. Geometry and numerical consistency tests do not establish those
physics claims.
