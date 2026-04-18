# Design

> **Document hierarchy:** VISION -> ARCHITECTURE -> **DESIGN**
> -> PSEUDOCODE -> Code. For goals and principles, see
> `VISION.md`. For repository layout and module map, see
> `ARCHITECTURE.md`.

---

## 1. LAT K-Point Integration

### 1.1 Overview and Motivation

The current Brillouin-zone integration uses Gaussian broadening:
each eigenvalue is smeared by a Gaussian of width sigma, and the
contributions are summed with k-point weights. This introduces
an arbitrary broadening parameter that affects the shape of the
DOS and requires dense k-point meshes for convergence.

The Linear Analytic Tetrahedral (LAT) method (Bloechl, Jepsen,
& Andersen, PRB 49, 16223, 1994) decomposes the BZ into
tetrahedra and integrates analytically within each one. This
eliminates the broadening parameter and provides better accuracy
at lower k-point densities.

### 1.2 Tetrahedra Generation

The uniform Monkhorst-Pack mesh defines a grid of nA x nB x nC
parallelepipeds. Each parallelepiped has 8 corners and is
decomposed into exactly 6 tetrahedra that tile without overlap.
The standard decomposition (Bloechl 1994) shares the main
diagonal M1-M8:

```
Parallelepiped corners at grid position (a, b, c):
  M1 = (a,   b,   c  )    M5 = (a+1, b+1, c  )
  M2 = (a+1, b,   c  )    M6 = (a+1, b,   c+1)
  M3 = (a,   b+1, c  )    M7 = (a,   b+1, c+1)
  M4 = (a,   b,   c+1)    M8 = (a+1, b+1, c+1)

Six tetrahedra sharing diagonal M1-M8:
  T1: M1, M2, M5, M8      T4: M1, M4, M7, M8
  T2: M1, M3, M5, M8      T5: M1, M4, M6, M8
  T3: M1, M3, M7, M8      T6: M1, M2, M6, M8
```

The mesh is periodic (the BZ is a torus), so indices wrap with
modular arithmetic: `mod(a, nA) + 1` etc. The total count is:

  numTetrahedra = 6 * nA * nB * nC

All tetrahedra span an equal fraction of the BZ:

  tetraVol = 1 / numTetrahedra

The tetrahedra reference the FULL uniform mesh, not the
IBZ-reduced kpoints. Under Option A, eigenvalues are
computed at IBZ points during SCF, then unfolded to the
full mesh for post-processing via a mapping array.

The `generateTetrahedra` call must happen AFTER the axial
kpoint counts are known (after `computeAxialKPoints` for
style code 2, or after reading for style codes 0 and 1).
It is called from `initializeKPoints`, not `readKPoints`.

`tetraVol` is computed in `initializeKPoints` after the
tetrahedra are generated.

**New module-level data in O_KPoints:**
```fortran
integer :: numTetrahedra
integer :: numFullMeshKP
real(kind=double) :: tetraVol
integer, allocatable, dimension(:,:) :: tetrahedra
    ! (4, numTetrahedra) -- indices into the full mesh
integer, allocatable, dimension(:) :: fullKPToIBZKPMap
    ! (numFullMeshKP) -- maps each full mesh kpoint
    ! to its IBZ representative index
```

The `fullKPToIBZKPMap` is produced by the IBZ folding in
`initializeKPointMesh` (from the `kPointTracker` array).
For each full-mesh index i, `fullKPToIBZKPMap(i)` gives
the IBZ kpoint index. Eigenvalue lookup for the TDOS
becomes:
  `eigenValues(band, fullKPToIBZKPMap(fullMeshIdx), spin)`

### 1.3 LAT TDOS (Eigenvalues Only)

The LAT TDOS loop structure inverts relative to Gaussian
broadening:

- **Gaussian (current):** outer = k-points, inner = states,
  innermost = energy bins (Gaussian smear)
- **LAT:** outer = bands, middle = tetrahedra, innermost =
  energy bins (analytic formula)

For each band n and each tetrahedron, sort the 4 corner
eigenvalues: e1 <= e2 <= e3 <= e4. The DOS contribution g(E)
per tetrahedron (Lehmann-Taut / Bloechl):

| Range          | Formula                                    |
|----------------|--------------------------------------------|
| E < e1         | 0                                          |
| e1 <= E < e2   | 3(E-e1)^2 / [(e2-e1)(e3-e1)(e4-e1)]       |
| e2 <= E < e3   | Bloechl eqs. 14-16 (see reference below)   |
| e3 <= E < e4   | 3(e4-E)^2 / [(e4-e1)(e4-e2)(e4-e3)]       |
| E >= e4        | 0                                          |

The middle range (e2 <= E < e3) has a more complex formula
involving cross-terms between all four eigenvalues. The exact
expressions are in Bloechl 1994, equations 14-16. The Bloechl
correction terms (eqs. 22-24) substantially improve accuracy at
lower k-point densities and should be included.

**Degenerate-corner guards:** When two or more corner energies
coincide, denominators in the analytic formulas go to zero. All
formulas must include guards: `if (abs(e2-e1) < eps) then ...`.

**Units:** Eigenvalues are in Hartree. The DOS output is in
states/eV. The Bloechl formulas give DOS in units of
1/(energy units of e), so a Hartree-to-eV conversion factor
is needed in the output, consistent with the current code's
use of `sigmaSqrtPi / hartree`.

**Unified corner DOS subroutine.** The per-tetrahedron TDOS
g(E) is the *sum* of four per-corner density weights
dw_c/dE (the energy derivatives of the Bloechl cumulative
corner weights, eqs. 18-21). These per-corner derivatives
are the more fundamental quantity: the TDOS needs only
their sum, but the PDOS (section 1.4) needs them
individually to weight each corner's Mulliken projection.

A single subroutine
`bloechlCornerDOSWt(E, eps, cornerDOSWt_LAT)` computes
the four dw_c/dE values. It serves both paths:

  TDOS:  dosContrib = sum(cornerDOSWt_LAT(1:4))
  PDOS:  pdosComp(alpha) += cornerDOSWt_LAT(c)
         * V_T * proj(alpha, n, k_c)

This eliminates duplicated case logic. The inline
dosContrib formulas in `computeTDOS_LAT` are replaced by
a call to `bloechlCornerDOSWt` + sum. The subroutine
follows the same case structure (Cases 0-3) as
`bloechlCornerWeights` but returns the energy derivative
of each corner weight (the per-corner spectral density)
rather than the cumulative value.

The identity
`sum(cornerDOSWt_LAT) == dosContrib` provides a built-in
self-consistency check.

**Distinction from `bloechlCornerWeights`.** The
cumulative corner integration weights `cornerIntgWt_LAT`
(returned by `bloechlCornerWeights`) give the fraction of
the tetrahedron's occupation attributed to each corner up
to energy E. They are dimensionless and are the correct
quantity for integrated properties like
`electronPopulation_LAT` (section 1.5), which evaluates
them at a single energy (the Fermi level). The corner DOS
weights `cornerDOSWt_LAT` (returned by
`bloechlCornerDOSWt`) have units of 1/energy and are the
correct quantity for the energy-resolved DOS density.
Both subroutines share the same case structure and
intermediate variables; only the final expressions
differ.

**Diagnostic integral unit fix.** The trapezoidal integral
that computes "Spin States Calculated" (the integrated
area under the TDOS) multiplies the TDOS (in states/eV)
by deltaDOS (stored in Hartree). This produces
states/27.211 instead of states. The integral must include
a hartree conversion factor: `deltaDOS * hartree`. This
same fix applies to the integrated-area diagnostic in
both `computeTDOS_LAT` and `computeDOS`.

**BZ weight normalization.** The Gaussian DOS path uses
`kPointWeight` as its BZ integration weight. By
convention, `kPointWeight` sums to 2.0 (set as
`weightSum` in kpoints.f90). This factor of 2 accounts
for the two electron spin states per band in a
spin-unpolarized calculation: each band contributes
`sum(kPointWeight)/spin` = 2/1 = 2 spin states. The
tetrahedron BZ integration uses `tetraVol`, which sums
to 1.0 (just the geometric BZ fraction). To produce
DOS values on the same scale as the Gaussian path, the
LAT accumulation must include `sum(kPointWeight)` as a
multiplicative factor alongside `tetraVol`. Both
`computeTDOS_LAT` and `integratePDOS_LAT` require this
factor. For spin-polarized calculations (spin=2), the
factor becomes 2/2 = 1 spin state per band, which is
also correct.

### 1.4 LAT PDOS (Energy-Resolved Partial DOS)

The corner DOS weights `cornerDOSWt_LAT` returned by
`bloechlCornerDOSWt` (section 1.3) determine how much of
each corner's partial-DOS projection to include at energy
E. These weights are occupation-independent: they
distribute spectral density among the four tetrahedron
corners at any energy, whether occupied or unoccupied.
They are the LAT replacement for the Gaussian broadening
+ `kPointWeight` mechanism used in the current PDOS path.

The PDOS contribution from tetrahedron T, band n is:

  dPDOS_alpha(E) = (V_T/V_BZ) * sum_{c=1..4}
      cornerDOSWt_LAT(c) * p_alpha(k_sigma(c), n)

where `cornerDOSWt_LAT(c)` is the Bloechl corner DOS
weight from `bloechlCornerDOSWt`, p_alpha(k, n) is the
Mulliken projection of channel alpha onto band n at
k-point k (`oneValeRealAccum` in the current code), and
sigma is the permutation that sorts eigenvalues. The
`cornerDOSWt_LAT` values depend on E, the four sorted
corner eigenvalues, and the corner index c -- they are
recomputed for every (tetrahedron, band, energy point)
combination and never stored as a module-level array.

**The fundamental constraint:** p_alpha is needed at all 4
corner k-points of each tetrahedron simultaneously. The
current code computes projections one k-point at a time and
immediately discards them. This forces a two-pass design:

**Pass 1 -- k-streaming (compute projections):**
  For each k-point, read eigenvectors from HDF5, compute
  Mulliken projections for all bands and channels, store in
  P(alpha, band, kpoint).

**Pass 2 -- tetrahedron integration:**
  For each band and tetrahedron, sort corner eigenvalues,
  compute `bloechlCornerDOSWt`, accumulate weighted
  projections into PDOS.

**Memory strategy (resolved D1):** Store projections
only at IBZ k-points: P(alpha, band, k_IBZ). When
assembling tetrahedron corners, look up the IBZ
representative via fullKPToIBZKPMap and apply atom
permutation on-the-fly to map the channel index:

  P(alpha, n, k_full) =
      P(permuted_alpha, n, fullKPToIBZKPMap(k_full))

where the channel permutation follows from atomPerm
and the operation stored in fullKPToIBZOpMap(k_full).
For mode 0 (per atom-type, per l-shell), R preserves
species, so no channel permutation is needed.

This reduces memory by the IBZ reduction factor
(typically 4-48x depending on point group symmetry):
- Moderate system (200 channels, 500 states, 500 IBZ
  kpts from 2000 full): ~0.4 GB
- Large system (5000 channels, 500 states, 1000 IBZ
  kpts from 5000 full): ~20 GB

The atom permutation infrastructure (atomPerm,
fullKPToIBZOpMap) is shared with the Q* and bond order
fix (section 2.4). One additional array is introduced:
invAtomPerm (section 2.4, item 4), which provides the
inverse mapping R^{-1}(A) for backward channel
permutation during tetrahedron corner assembly.

**Inverse atom permutation.** The channel permutation
for modes 1-2 requires R^{-1}(A) where R is the forward
operation (k_IBZ → k_full) stored in fullKPToIBZOpMap.
Since atomPerm stores R(A), we build invAtomPerm(R, B)
= A where atomPerm(R, A) = B, giving R^{-1}(B)
directly. It is built in O_AtomicSites alongside
atomPerm, with array shape (numPointOps, numAtomSites).

**Per-mode channel permutation rules:**

  Mode  Channel           Permutation rule
  ───────────────────────────────────────────────────
  0     per-type, per-l   None: type-level sum is
                          invariant under R
  1     per-atom total    invAtomPerm(R, atomIdx)
  2     per-atom, per-l   invAtomPerm remaps atom;
                          l-shell offset unchanged
                          (same species ⇒ same
                          orbital structure)
  3     per-atom, per-lm  Not supported: requires
                          D^l(R) rotation matrices
  ───────────────────────────────────────────────────

For modes 1-2 a precomputed channelPermTable(R, alpha)
avoids repeated index decode/encode in the inner loop.
Mode 0 needs no table (identity mapping). Mode 2
decodes alpha into (atom, l-offset), permutes the atom
via invAtomPerm, and re-encodes using the permuted
atom's cumulative offset in cumulNumDOS.

**Mode 3 restriction.** When kPointIntgCode == 1 and
detailCodePDOS == 3, the program stops with a clear
error message. Individual Cartesian Gaussian projections
(px, py, pz separately) mix under rotation via D^l(R)
(section 2.3). Atom permutation alone does not suffice
and the full rotation matrices are not available.

**Refactored computeDOS structure.** computeDOS gains
an internal branch on kPointIntgCode inside the spin
loop. The setup phase (pdosIndex construction,
cumulNumDOS, allocations) and output phase (file
writing, normalization) are shared between Gaussian and
LAT. Only the computation phase -- filling
pdosComplete -- differs:

  computeDOS(inSCF):
    Setup (shared, unchanged)
    do h = 1, spin
      if kPointIntgCode == 1:
        LAT two-pass → fills pdosComplete
      else:
        Gaussian single-pass → fills pdosComplete
      endif
      Output (shared): normalization, file writing
    enddo

The Gaussian single-pass is unchanged, preserving its
memory efficiency (no projection array needed).

**OLCAO.F90 dispatch.** The calling sequence becomes:

  if (kPointIntgCode == 1) then
    call computeTDOS_LAT
  endif
  call computeDOS(inSCF)

computeTDOS_LAT remains the validated eigenvalue-only
TDOS writer (fort.60/61). Inside computeDOS, the LAT
branch writes PDOS (fort.70/71) and localization index
(fort.80/81) but skips TDOS output (already written).
The Gaussian branch writes all three as before.

**Normalization.** The Gaussian path normalizes
pdosComplete by electronFactor (ratio of exact electron
count to Gaussian-broadened count) to correct for
broadening-tail truncation. For LAT, tetrahedron corner
weights provide exact BZ integration, so electronFactor
≈ 1.0. The LAT branch computes and logs this ratio as a
diagnostic (to fort.20) but does not apply it to
pdosComplete.

### 1.5 electronPopulation_LAT for Integrated Properties

For integrated (energy-summed) partial properties -- effective
charge, bond order -- the existing k-point loops barely need
to change if we precompute the LAT analog of
`electronPopulation`.

Define `electronPopulation_LAT(n, k, spin)`: the fractional
electron occupation of state (n, k) as determined by
tetrahedron integration. It answers the same physical
question as `electronPopulation` -- "how occupied is this
state, weighted for BZ integration?" -- but computed via the
LAT method instead of Gaussian broadening + Fermi filling.

  electronPopulation_LAT(n, k) =
      sum over all tetrahedra T containing corner k {
          (V_T / V_BZ) * w_c(E_Fermi)
      }

where w_c(E_Fermi) is the cumulative Bloechl corner
weight from `bloechlCornerWeights` (eqs. 18-21 evaluated
at the Fermi energy). Size: numStates x numKPoints x spin --
for 500 bands, 2000 kpts: ~8 MB.

**Naming convention:** The `_LAT` suffix identifies the
integration method. If future methods fill the same role
(a per-state occupation weight for BZ integration), they
follow the pattern `electronPopulation_XXX`.

**Usage:** In `computeBond` and the effective-charge branch
of `computeDOS`, replace:
  `electronPopulation(stateSpinKPointIndex)`
with:
  `electronPopulation_LAT(j, i, h)` (band j, kpoint i,
  spin h)

This unification means:
- Bond order loop: structure unchanged, swap weight source
- Effective charge loop: same
- Energy-resolved PDOS: `bloechlCornerDOSWt` (two-pass,
  occupation-independent; see section 1.4)
- TDOS: simplest case, eigenvalues only

The integration method (Gaussian vs. LAT) becomes a
pluggable parameter in a common projection-then-weight
framework.

**Distinction from `bloechlCornerDOSWt`.** The corner
DOS weights `cornerDOSWt_LAT` (sections 1.3-1.4) are
energy-resolved, occupation-independent, and transient
(recomputed per tetrahedron per energy point).
`electronPopulation_LAT` uses the cumulative corner
integration weights `cornerIntgWt_LAT` from
`bloechlCornerWeights`, evaluated once at the Fermi
energy; it is occupation-integrated, stored as a
module-level array, and has the same lifecycle as
`electronPopulation`.

**Resolved (D2):** Replacing `electronPopulation` with
`electronPopulation_LAT` provides the correct occupation
weight, but does not by itself cover bond order
accumulation correctly. The accumulation loop must also
apply atom permutation to distribute each IBZ k-point's
contribution correctly. See section 2.5 for the full
analysis.

---

## 2. IBZ Correctness for Eigenvector-Dependent Quantities

### 2.1 Physical Basis for Symmetry in K-Space

Understanding why the IBZ reduction works for some quantities
and fails for others requires tracing the symmetry argument
from its origin in atomic geometry through to its consequences
for computed properties.

**From atomic symmetry to Hamiltonian symmetry.** A crystal's
point group is defined by the nuclear positions: a symmetry
operation R maps every nucleus onto another nucleus of the
same species. The electronic Hamiltonian -- kinetic energy,
electron-nuclear attraction, and electron-electron interaction
-- inherits this symmetry because R leaves the potential
energy landscape unchanged. Formally, [H, R] = 0, so the
eigenstates of H must transform as irreducible
representations of the point group. This is what forces
e_n(k) = e_n(Rk) for any point-group operation R.

The key point: although the quantity of interest is the
electronic wavefunction, we are relying on the atomic
geometry to set our expectations for electronic symmetry.
This reliance is rigorously justified (within the
Born-Oppenheimer approximation) because H is entirely
determined by the nuclear configuration. The one caveat is
symmetry-broken electronic phases (magnetic ordering, charge
ordering, Jahn-Teller distortions) where the electronic
ground state spontaneously adopts lower symmetry than the
nuclear geometry -- but these are special cases outside the
scope of a standard DFT/HF calculation.

**Basis closure under symmetry operations.** Does the choice
of an LCAO basis with Cartesian Gaussian atomic orbitals
undermine the symmetry argument? No -- provided the basis is
properly constructed.

The eigenvalues of H are basis-independent; the symmetry
relation e_n(k) = e_n(Rk) is a statement about the physics,
not the representation. The practical requirement is that the
basis set must be *closed* under the point-group operations.
For atom-centered functions this means: every symmetry-
equivalent atom (related by some operation R) must carry the
same set of basis functions. When R maps atom A to atom B, it
also maps the basis functions on A to the corresponding
functions on B. If the basis is identical on both atoms, the
full basis set is invariant under R, and the Hamiltonian
matrix respects the symmetry.

For Cartesian Gaussians there is one subtlety: the angular
parts (x^a * y^b * z^c) do not individually transform as
irreducible representations. A rotation mixes them. For
example, the six degree-2 Cartesian functions (xx, yy, zz,
xy, xz, yz) span a reducible representation (5 d-type + 1
s-type). Under R, a single Cartesian Gaussian maps to a
linear combination of Cartesian Gaussians of the same degree
on the rotated atom. As long as the complete set of
Cartesians for each degree is included on each equivalent
atom, the basis remains closed and the symmetry is preserved.

**Eigenvalues vs. eigenvectors.** The symmetry argument
guarantees that eigenvalues are invariant: e_n(k) = e_n(Rk).
Eigenvalues are scalar quantities unchanged by unitary
transformation, so the IBZ reduction is safe for any property
that depends only on eigenvalues -- total DOS, total energy,
band structure.

Eigenvectors, however, are *not* invariant. The eigenvectors
at Rk are related to those at k by a unitary transformation
that mixes the orbital expansion coefficients according to R.
They are different vectors in the basis representation. This
distinction is the root cause of the IBZ problem described in
the following subsections: any property that depends on the
wavefunction expansion coefficients (PDOS, bond order,
effective charge) cannot be correctly computed by simply
scaling the IBZ representative's contribution by the star
multiplicity.

**Time-reversal symmetry.** In addition to the crystallographic
point group, time-reversal symmetry (e_n(k) = e_n(-k)) often
doubles the effective symmetry, even when -k is not related
to k by any point-group operation alone. Most codes fold this
in when constructing the IBZ.

### 2.2 The Problem

Observed in KNbO3: bond orders for all O atoms are identical
(within machine precision) when using the full k-point mesh,
but differ for each O atom when using IBZ-reduced k-points.
This is incorrect -- the crystal symmetry requires them to be
identical.

### 2.3 Root Cause: Eigenvector Transformation Law

The root cause of the IBZ bug is that eigenvector-dependent
quantities do not simply scale with the star multiplicity.
This subsection develops the transformation rules needed to
determine which quantities are safe and which are not.

**LCAO wavefunctions and Mulliken projections.** In LCAO
the wavefunction for band n at k-point k is:

  psi_n(k) = sum_mu  c_mu(n,k) * phi_mu(k)

where mu indexes every orbital on every atom in the unit
cell. The coefficients c_mu(n,k) are the eigenvector
components (`valeVale` in OLCAO). The Mulliken population
of orbital mu in state (n,k) partitions the state's
electron count among basis functions:

  p_mu(n,k) = Re[ c_mu*(n,k)
                * sum_nu c_nu(n,k) * S_{mu,nu}(k) ]

where S is the overlap matrix. These satisfy the partition
sum_mu p_mu(n,k) = 1. Every property of interest --
effective charge Q*, bond order, PDOS -- is built from
sums of p_mu over specific index ranges (per atom, per
atom pair, per angular momentum shell). The Mulliken
overlap (bond order contribution) between orbitals mu on
atom A and nu on atom B is similarly:

  b_{mu,nu}(n,k) = Re[ c_mu*(n,k)
                      * c_nu(n,k) * S_{mu,nu}(k) ]

**How a symmetry operation transforms the eigenvector.**
A point-group operation R acts on three things at once:

  (a) Atoms: R maps atom A to atom R(A), preserving
      species and basis function types.
  (b) Orbitals: angular parts transform by the
      representation matrix D^l(R) for each l-shell.
      Example -- 90-degree rotation around z on p-orbs:

              [ 0  -1   0 ]
     D(R)  = [ 1   0   0 ]   px -> -py
              [ 0   0   1 ]   py -> px, pz -> pz

  (c) K-points: k maps to Rk in the Brillouin zone.

The complete eigenvector transformation is:

  C_n(Rk) = D(R) * C_n(k)                          (1)

where D(R) is block-diagonal: one block per atom, each
block being D^l(R) for that atom's orbital set. The
overlap matrix transforms consistently:

  S(Rk) = D(R) * S(k) * D(R)^dagger                (2)

**The Mulliken vector argument.** Define the Mulliken
vector M(n,k) = S(k) * C_n(k), so that the population
is p_mu(n,k) = Re[ c_mu*(n,k) * M_mu(n,k) ].

Under R, both the eigenvector and the Mulliken vector
transform by the same matrix:

  C_n(Rk) = D(R) * C_n(k)                  [from (1)]

  M(n,Rk) = S(Rk) * C_n(Rk)
           = [D(R) S(k) D(R)^dag] [D(R) C_n(k)]
           = D(R) * S(k) * C_n(k)
           = D(R) * M(n,k)                         (3)

The cancellation D(R)^dag D(R) = I in the middle step
is the key identity -- it ensures that eigenvectors and
their overlap-weighted counterparts transform in lockstep.

**Shell-sum invariance (proof).** The Mulliken population
of orbital mu at the symmetry-related point Rk is:

  p_mu(n,Rk) = Re[ sum_{a,b}  D(R)*_{mu,a} c_a*
                             * D(R)_{mu,b}  M_b    ]

(suppressing (n,k) arguments for brevity). Now sum over
all mu in a complete l-shell on atom R(A). These mu are
the images of the l-shell on atom A under D(R), so the
sum invokes the unitarity of D(R) within that subspace:

  sum_mu  D(R)*_{mu,a} * D(R)_{mu,b} = delta_{a,b}

The off-diagonal terms (a != b) vanish, giving:

  sum_{mu in l-shell on R(A)}  p_mu(n, Rk)

      = Re[ sum_a  c_a*(n,k) * M_a(n,k) ]

      = sum_{a in l-shell on A}  p_a(n, k)         (4)

The l-shell-summed Mulliken population on atom R(A) at
Rk equals the l-shell-summed population on atom A at k.
Atom totals obey the same relation by summing all shells.
Bond order between atom pairs follows by the same proof
applied to the two-atom coefficient product:

  B(R(A), R(B), n, Rk) = B(A, B, n, k)             (5)

In every case, no explicit rotation matrices are needed
-- only the atom relabeling A -> R(A).

**Why individual orbitals break the pattern.** For a
single orbital mu (not summed over a shell), the Mulliken
population at Rk is:

  p_mu(n,Rk) = Re[ sum_{a,b}  D(R)*_{mu,a} D(R)_{mu,b}
                             * c_a* * M_b              ]

Without a sum over mu, unitarity cannot be invoked and
the cross terms (a != b) persist. The individual-orbital
projection depends on the full D(R) matrix, not simply
on which atom mu belongs to.

Example: p-orbitals under 90-degree rotation around z.
Suppose the coefficients at k for one band are c_px=0.5,
c_py=0.3, c_pz=0.1. Applying D(R) gives at Rk:
c_px = -0.3, c_py = 0.5, c_pz = 0.1.

Simplified Mulliken projections (|c|^2):

  At  k:  |c_px|^2 = 0.25   |c_py|^2 = 0.09
  At Rk:  |c_px|^2 = 0.09   |c_py|^2 = 0.25

The px and py projections swap -- they are NOT related
by a simple atom permutation. But the p-shell sum is
0.35 at both k and Rk, preserved by ||c||^2 = ||Dc||^2.
This is why PDOS modes that sum over complete l-shells
(detailCodePDOS 0-2) work with atom permutation, while
individual-orbital PDOS (detailCodePDOS 3) does not.


### 2.4 The Atom Permutation Fix

From equation (4), the fix for any quantity that sums
Mulliken projections over a complete l-shell or over an
entire atom does not require rotating eigenvectors -- it
requires only knowing which atom maps to which under each
symmetry operation: the atom permutation table.

**Corrected accumulation.** For each IBZ k-point k_i,
instead of multiplying the projection by the star
multiplicity, loop over each operation R_s in the star
of k_i and accumulate into the permuted atom indices:

  Effective charge:
    chargeContrib(R_s(A)) += p_A(n,k_i) * f(n,k_i)

  Bond order:
    bondContrib(R_s(A), R_s(B)) += b(A,B,n,k_i)
                                 * f(n,k_i)

where f(n,k_i) is the occupation weight.

**Example 1: charge in a mirror-symmetric 1D chain.**
Two atoms per cell (A, B), related by mirror m that
swaps A with B. Each has one s-orbital. The IBZ contains
one k-point; the star is {k, mk} with size 2.

Eigenvectors at k for the bonding band:

  c_A = 0.8,  c_B = 0.6

At the mirror image mk, D(m) swaps coefficients:

  c_A = 0.6,  c_B = 0.8   (same eigenvalue)

Simplified Mulliken projections (|c|^2):

  At  k:   p_A = 0.64    p_B = 0.36
  At mk:   p_A = 0.36    p_B = 0.64

Note the atom permutation: p_A(mk) = p_B(k).

Naive IBZ weighting (WRONG -- this is the KNbO3 bug):

  charge(A) = 2 * 0.64 = 1.28
  charge(B) = 2 * 0.36 = 0.72    symmetry violated

Atom permutation (CORRECT):

  From identity: charge(A) += 0.64  charge(B) += 0.36
  From mirror:   charge(B) += 0.64  charge(A) += 0.36

  Result:  charge(A) = 1.00,  charge(B) = 1.00

**Example 2: bond order with C3 rotation symmetry.**
Three atoms A, B, C with 120-degree rotation symmetry.
R1 (120 deg): A->B, B->C, C->A. R2 (240 deg): A->C,
B->A, C->B. The IBZ k-point has star size 3.

Mulliken overlaps at k_1 for one occupied band:

  b(A,B) = 0.15    b(A,C) = 0.10    b(B,C) = 0.20

Naive star-weight multiplication (WRONG):

  BO(A,B) = 3*0.15 = 0.45
  BO(A,C) = 3*0.10 = 0.30
  BO(B,C) = 3*0.20 = 0.60          C3 violated

Atom permutation (CORRECT):

  R0: BO(A,B)+=0.15  BO(A,C)+=0.10  BO(B,C)+=0.20
  R1: BO(B,C)+=0.15  BO(B,A)+=0.10  BO(C,A)+=0.20
  R2: BO(C,A)+=0.15  BO(C,B)+=0.10  BO(A,B)+=0.20

Collecting (BO is symmetric):

  BO(A,B) = 0.15 + 0.10 + 0.20 = 0.45
  BO(B,C) = 0.20 + 0.15 + 0.10 = 0.45
  BO(A,C) = 0.10 + 0.20 + 0.15 = 0.45     C3 OK

**Required infrastructure:**

1. `atomPerm(numOps, numAtomSites)`: atom permutation
   table. atomPerm(R, A) = B means operation R maps atom
   A to atom B.
2. `fullKPToIBZOpMap(numFullMeshKP)`: for each full-
   mesh k-point i, the index of the symmetry operation
   R such that R(k_IBZ) = k_full(i) -- i.e., the
   operation that maps the IBZ representative to the
   full-mesh k-point. Currently `fullKPToIBZKPMap`
   gives only the IBZ index, not which operation did
   the mapping. For an IBZ k-point itself, the stored
   operation is the identity.
3. Star decomposition: for each IBZ k-point, the set
   of operations in its star. This falls naturally out
   of `fullKPToIBZOpMap` by collecting all full-mesh
   k-points that share the same IBZ representative.
4. `invAtomPerm(numOps, numAtomSites)`: inverse atom
   permutation. invAtomPerm(R, B) = A where
   atomPerm(R, A) = B, i.e., R^{-1}(B) = A. Built
   in O_AtomicSites alongside atomPerm. Required for
   LAT PDOS channel unfolding (section 1.4): when
   assembling projections at full-mesh corner k_f
   mapped from IBZ point k_i by forward operation R,
   the channel index must be transformed by R^{-1}
   to reference the stored IBZ-kpoint projection.


### 2.5 Per-Quantity Implications

The table below summarizes what each computed quantity
requires for correct IBZ unfolding, based on the shell-
sum invariance (4) and bond order invariance (5):

  Quantity                  Needed            Why
  ────────────────────────────────────────────────────
  TDOS (eigenvalues)        fullKPToIBZKPMap  e(Rk)=e(k)
  Q* (effective charge)     atom perm         eq. (4)
  Bond order                atom perm         eq. (5)
  PDOS mode 0 (type, l)    nothing extra      *
  PDOS mode 1 (atom total) atom perm         eq. (4)
  PDOS mode 2 (atom, l)    atom perm         eq. (4)
  PDOS mode 3 (atom, lm)   D^l(R) matrices   **
  ────────────────────────────────────────────────────

*Mode 0 sums projections over all atoms of the same
type. Since point-group operations permute atoms only
within the same species, the type-level sum is
automatically invariant -- no correction needed.

**Mode 3 resolves individual Cartesian Gaussian
components (px, py, pz separately). Under rotation these
mix via D^l(R), so atom permutation alone does not
suffice. Correct unfolding requires the full
representation matrices for each l-shell. This is
deferred -- mode 3 is rarely used in practice.

**Resolution of D2.** The open question asked whether
replacing electronPopulation with electronPopulation_LAT
covers bond order accumulation correctly. The answer is
no. The occupation weight from electronPopulation_LAT is
correct (it properly sums tetrahedron contributions for
each IBZ k-point), but the Mulliken projection is
computed only at the IBZ representative. Multiplying a
correct weight by a non-permuted projection distributes
charge incorrectly among symmetry-equivalent atoms.

The fix has two parts that work together:

  (a) electronPopulation_LAT provides the correct
      occupation weight per (band, kpoint, spin).
  (b) The accumulation loop distributes each IBZ
      k-point's contribution across atom pairs using
      the atom permutation table (section 2.4).

Together, (a) and (b) give correct per-atom Q* and per-
pair bond order. Either alone is insufficient. Note that
total charge summed over all atoms IS correct with (a)
alone -- the error is only in the per-atom distribution.

### 2.6 Relation to LAT and Pragmatic Options

Both the LAT and Gaussian integration paths share the same
IBZ symmetry issue. LAT does not bypass it -- eigenvectors
still come only from IBZ k-points regardless of how the
occupation weights are computed. The distinction between
the two paths is purely in the weight calculation, not in
how projections must be unfolded.

**LAT-specific note on PDOS.** For LAT PDOS (section
1.4), the tetrahedron integration needs Mulliken
projections at all four corners of each tetrahedron
simultaneously. The corners are full-mesh k-point indices.
If we only diagonalize at IBZ k-points, each full-mesh
corner's projection must be "unfolded" from its IBZ
representative via atom permutation. Concretely, if
full-mesh corner k_f maps to IBZ point k_i via operation
R (stored in fullKPToIBZOpMap). Because R maps k_IBZ
to k_f (forward direction: R(k_IBZ) = k_f), the
inverse R^{-1} maps atoms at k_f back to those at
k_IBZ:

  p_{atom A, l-shell}(k_f, n) =
      p_{atom R^{-1}(A), l-shell}(k_i, n)

This is equation (4) applied to corner unfolding.

**Decided approach (Option A):** Use IBZ for SCF
diagonalization, unfold eigenvalues (and eventually
eigenvectors) to the full BZ for post-processing. The SCF
benefits from fewer diagonalizations, while tetrahedra
reference the full mesh for integration.

**Workflow scenarios:**

OLCAO reads separate kpoint files for SCF (`fort.15` =
`kp-scf.dat`) and PSCF (`fort.16` = `kp-pscf.dat`). Each
kpoint set goes through `readKPoints` + `initializeKPoints`
independently. The SCF and PSCF meshes can differ in
density, style code, and integration code.

Supported combinations involving LAT:
1. SCF with IBZ → DOS/bond with LAT (within SCF phase,
   same kpoints): The SCF diagonalizes at IBZ points. When
   `doDOS_SCF=1` and `kPointIntgCode=1`, the DOS routine
   unfolds eigenvalues to the full mesh via `fullKPToIBZKPMap`
   and uses tetrahedra for integration.
2. SCF with IBZ → PSCF with LAT (different, typically
   denser mesh): The PSCF reads its own kpoint file, builds
   its own mesh, IBZ reduces it, builds its own tetrahedra
   and `fullKPToIBZKPMap`. The PSCF diagonalizes at its IBZ
   points, then LAT DOS unfolds to its full mesh.
3. SCF with IBZ → PSCF with Gaussian (standard current
   behavior): No tetrahedra needed. IBZ kpoints and weights
   used directly.

The `fullKPToIBZKPMap`, `tetrahedra`, `numTetrahedra`, and
`tetraVol` are per-kpoint-set data stored in `O_KPoints`.
They are rebuilt each time `initializeKPoints` runs (once
for SCF, once for PSCF). When `kPointIntgCode == 0`
(Gaussian), they are not allocated.
`electronPopulation_LAT` lives in `O_Populate` alongside
its sibling `electronPopulation`; tetrahedra data is
passed in from O_KPoints as arguments.

**Implementation strategy:**

- The IBZ folding in `initializeKPointMesh` saves the
  full-to-IBZ mapping as
  `fullKPToIBZKPMap(numFullMeshKP)`. For each full-mesh
  index i, `fullKPToIBZKPMap(i)` gives the IBZ kpoint
  index whose eigenvalues are identical.
- `fullKPToIBZOpMap(numFullMeshKP)` (new): for each
  full-mesh index i, the index of the point group
  operation R such that R(k_IBZ) = k_full(i). The
  folding loop applies each operation to the IBZ
  representative and checks for matches among full-mesh
  k-points; the matching operation is stored. For an
  IBZ k-point itself, the identity operation is stored.
  Saved alongside `fullKPToIBZKPMap` during IBZ folding
  in `initializeKPointMesh`.
- `atomPerm(numPointOps, numAtomSites)` (new): for each
  point-group operation and atom, the index of the image
  atom. Built once during `initializeKPoints` from the
  point-group operations and atomic positions.
- `generateTetrahedra` uses `numAxialKPoints` to build
  tetrahedra referencing full-mesh indices (1 to
  `numFullMeshKP`). Called from `initializeKPoints`
  after the mesh is constructed.
- `tetraVol = 1 / numTetrahedra`, computed in
  `initializeKPoints` after tetrahedra are generated.
- For TDOS, eigenvalue lookup at full-mesh corner k is:
    `eigenValues(band, fullKPToIBZKPMap(k), spin)`
- For PDOS, projection at full-mesh corner k for an
  l-shell channel on atom A is:
    `p(R^{-1}(A), l, band, fullKPToIBZKPMap(k))`
  where R = operation(fullKPToIBZOpMap(k)) maps the
  IBZ representative to k (forward direction). See
  section 2.3 equation (4).
- **Style code 0 warning (resolved D3):** When
  `kPointStyleCode == 0` (explicit k-point list), OLCAO
  does not build the full mesh internally and therefore
  cannot construct `fullKPToIBZKPMap`, `fullKPToIBZOpMap`,
  or `atomPerm`. Emit a prominent warning at initialization
  that decomposition properties (effective charge, bond
  order, PDOS) will not be correct unless the user has
  taken extreme care to provide a symmetric k-point mesh.
  For style codes 1 and 2, OLCAO builds the full mesh
  and all symmetry maps internally, so the atom permutation
  fix works for both Gaussian and LAT integration paths.
  makeinput.py no longer produces style code 0 files;
  mesh mode (`-kp`) now writes style code 1.

Note: `kPointWeight` is irrelevant for LAT integration --
the DOS contribution is determined entirely by the
analytic formula over each tetrahedron. The weight array
still matters for SCF electron counting (charge density
integration).

---

## 3. Density-Based K-Point Input Pipeline

### 3.1 Motivation

The existing workflow requires the user to specify explicit
axial k-point counts (e.g. `-kp 4 4 4`). This forces the user
to know the lattice geometry in advance and manually choose
counts that give adequate sampling density. A density-based
option (`-kpd D`) lets the user specify a single number -- the
minimum k-point volume density (kpoints per unit reciprocal-
space volume, Bohr^-3) -- and the program determines the per-
axis counts automatically. The total kpoint count is at least
`D * recipCellVolume`, distributed as uniformly as possible
across the three axes. This gives users a geometry-independent
knob: the same density value produces the same sampling
quality regardless of cell size.

This mechanism already exists in `uolcao` (`kPointStyleCode=2`,
`computeAxialKPoints`). The goal here is to expose it through
`makeinput.py` so that the normal user workflow supports it.

### 3.2 Pipeline When Using Density Mode

When the user passes `-kpd`, `-scfkpd`, or `-pscfkpd`:

1. `makeinput.py` extracts the point group operations from
   the space group database (`_extract_point_ops`) and
   writes kpoint files (`kp-scf.dat` and/or `kp-pscf.dat`)
   directly, bypassing the `makeKPoints` executable.
2. Each file uses `kPointStyleCode=2` and contains:
   - `KPOINT_STYLE_CODE` = 2
   - `KPOINT_INTG_CODE` = 0 (default; histogram)
   - `MIN_KP_LINE_DENSITY` = the user's volume density
     (label is historical; the value is a volume density)
   - `KP_SHIFT_A_B_C` = the shift (auto or user-specified)
   - `NUM_POINT_OPS` = number of point group operations
   - `POINT_OPS` = the 3x3 rotation matrices (abc coords)
3. At runtime, `uolcao` reads this file in `readKPoints`
   (including the point group operations), then in
   `initializeKPoints` calls:
   - `computeAxialKPoints` — density to axial counts
   - `computeRecipPointOps` — convert abc point ops to
     reciprocal-space abc operations
   - `initializeKPointMesh(1)` — build the uniform mesh
     and fold it to the IBZ using the reciprocal-space
     point group operations
   - `convertKPointsToXYZ` — transform to Cartesian

### 3.3 Interaction with Existing Options

- `-kpd D` sets both SCF and PSCF density to D.
- `-scfkpd D` sets only the SCF density.
- `-pscfkpd D` sets only the PSCF density.
- Density mode is all-or-nothing: if any density option is
  given, both SCF and PSCF use density mode. An unset group
  defaults to a density of 1 (equivalent to a single
  k-point per axis). Explicit-mesh options (`-kp`, `-scfkp`,
  `-pscfkp`) are mutually exclusive with density options;
  if both are given, density wins and a warning is printed.
- `-kpshift` applies to both pipelines.
- `-printbz` is skipped when density mode is active. A note
  is printed informing the user that BZ visualization only
  works with an explicit k-point mesh.

### 3.4 File Format for Style Code 2

The file written by `makeinput.py` must match what `readKPoints`
in `kpoints.f90` expects. Specifically:

```
KPOINT_STYLE_CODE
2
KPOINT_INTG_CODE
0
MIN_KP_LINE_DENSITY
<density_value>
KP_SHIFT_A_B_C
<shift_a> <shift_b> <shift_c>
NUM_POINT_OPS
<n>
POINT_OPS
<3x3 matrix for op 1, one row per line>

<3x3 matrix for op 2, one row per line>
...
```

The point group operations are the rotational parts of the
space group symmetry operations with translations stripped.
They are given in fractional (abc) coordinates. Blank lines
between operations are optional (for readability) and are
skipped by the Fortran reader.

### 3.5 Impact on Memory Estimation and Summary

When using density mode, the exact number of k-points and the
axial mesh dimensions are not known at `makeinput` time (they
depend on the reciprocal lattice geometry, which is computed
inside `uolcao`). Two consequences:

- **Memory estimation:** Skipped entirely in density mode.
  To be revisited in the future.
- **Summary output:** The k-point count and mesh dimensions
  cannot be printed. The summary should indicate that
  density mode is active and print the density value instead
  of the count and mesh array.

---

## 4. UFF Bond Parameter Database

### 4.1 Motivation

The current `bonds.dat` enumerates specific element pairs
with hand-tuned force constants and rest bond lengths.  It
covers only six elements (H, B, C, N, O, Si) across 14 bond
types.  Every force constant is a uniform 5000.0 kcal/mol/A^2
-- roughly 15 times stiffer than physically realistic values
for typical covalent bonds.  Expanding coverage to 36 or more
elements in the pair-listing format is impractical: 36
elements yield up to 666 unique pairs, and 54 elements
yield 1,485.

A scalable alternative stores per-element parameters and
computes the equilibrium bond length and harmonic force
constant for any pair on the fly, using the Universal Force
Field (UFF) formulas (Rappe et al. 1992).  The new file is
named `bond_parameters.dat` to distinguish it from the legacy
`bonds.dat` format.

### 4.2 UFF Bond Stretching Model

The UFF describes bond stretching with a harmonic potential:

  E = K_ij * (r - r_ij)^2

where K_ij is the force constant (kcal/mol/A^2) and r_ij
is the equilibrium bond length (Angstroms).  The LAMMPS
`bond_style harmonic` uses this same convention: the K
parameter absorbs the factor of 1/2 that appears in the
physics textbook form E = (1/2) k x^2.

For any element pair (i, j), UFF defines:

**Equilibrium bond length:**

  r_ij = r_i + r_j - r_EN

  r_EN = r_i * r_j * (sqrt(chi_i) - sqrt(chi_j))^2
         / (chi_i * r_i + chi_j * r_j)

where r_i and r_j are single-bond covalent radii
(Angstroms) and r_EN is an electronegativity correction
that shortens bonds between elements of unequal
electronegativity.  For homonuclear bonds (same element),
chi_i = chi_j so r_EN = 0 and the equilibrium length
reduces to r_ij = 2 * r_i.

**Force constant (LAMMPS harmonic convention):**

  K_ij = 332.06 * Zstar_i * Zstar_j / r_ij^3

The prefactor 332.06 = 664.12 / 2 absorbs the 1/2 that
converts from the UFF convention E = (1/2) k (r-r0)^2 to
the LAMMPS convention E = K (r-r0)^2.  Zstar_i and Zstar_j
are the UFF effective charges (dimensionless).  The 664.12
constant carries units of kcal*A/mol and encodes the
fundamental relationship between bond stiffness, effective
nuclear charges, and bond length.

Each element requires only three tabulated parameters:

  Parameter  Meaning                       Units
  --------------------------------------------------
  r_i        Single-bond covalent radius   Angstroms
  Zstar_i    Effective charge              (none)
  chi_i      GMP electronegativity         eV

The formula is symmetric in i and j, so the computed values
are independent of element ordering.  However, the calling
code continues to enforce the Z_1 <= Z_2 convention used
throughout the codebase (bond analysis output, tag
construction, pair matching).  The `get_bond_params` method
accepts Z arguments in any order but the callers in
`create_lammps_files` canonicalize to Z_1 <= Z_2 before
constructing bond tags, exactly as the current code does.

**Validation against established force fields.**  UFF-derived
values agree with AMBER within 10-20 %, which is typical
inter-force-field variation for single bonds:

  Bond   UFF K       UFF r0   AMBER K     AMBER r0
         (kcal/      (A)      (kcal/      (A)
         mol/A^2)             mol/A^2)
  ----------------------------------------------------
  C-H    ~331        1.11     ~340        1.09
  C-C    ~350        1.51     ~310        1.53
  C-N    ~360        1.44     ~337        1.47
  O-H    ~540        0.98     ~553        0.96
  Si-O   ~285        1.63     --          1.61

### 4.3 Element Coverage

The UFF provides parameters for every element from Z = 1
through Z = 103.  The initial database covers Z = 1 through
Z = 54 (hydrogen through xenon), spanning:

- All main-group elements through the 5th period
- All 3d transition metals (Sc through Zn)
- All 4d transition metals (Y through Cd)
- Halogens, chalcogens, and pnictogens
- Noble gases (He, Ne, Ar, Kr, Xe)

Noble gases are included for table completeness.  Their very
small effective charges yield negligible force constants, so
they will not produce meaningful bonds in practice.

**Contiguity requirement.**  The table must contain every
element from Z = 1 through `NUM_UFF_ELEMENTS` with no gaps.
The validation check in `get_bond_params` tests
`z > num_uff_elements`, so a gap would leave an uninitialized
slot that could silently produce wrong results.  Extension
beyond Z = 54 requires appending rows for every Z up to the
new maximum -- no code changes are needed.

### 4.4 New bond_parameters.dat File Format

The new file `bond_parameters.dat` replaces `bonds.dat`.  Its
format changes from pair-enumeration to a per-element parameter
table.  Comment lines (beginning with `#`) are permitted and
skipped by the reader.  The tagged-section structure is
preserved for consistency with other OLCAO data files:

```
# UFF bond stretching parameters.
#
# Source: Rappe, A. K.; Casewit, C. J.; Colwell, K. S.;
#   Goddard, W. A., III; Skiff, W. M.
#   J. Am. Chem. Soc. 1992, 114, 10024-10035.
#   DOI: 10.1021/ja00051a040
#
# For any element pair (i, j):
#   r_ij = r_i + r_j - r_EN
#   r_EN = r_i r_j (sqrt(chi_i) - sqrt(chi_j))^2
#          / (chi_i r_i + chi_j r_j)
#   K_ij = 332.06 Zstar_i Zstar_j / r_ij^3
#          (kcal/mol/A^2, LAMMPS harmonic convention)
NUM_UFF_ELEMENTS
54
UFF_BOND_PARAMS
#  Z    r_i     Zstar_i  chi_i    Element
   1   0.3540  0.7120   4.5280  # H
   2   0.8490  0.0980   9.6600  # He
   3   1.3360  1.0260   3.0060  # Li
  ...
  54   ...                       # Xe
```

Each data line provides: the atomic number Z, the covalent
radius r_i (Angstroms), the effective charge Zstar_i, and the
GMP electronegativity chi_i (eV).  The element comment at the
end of each line is optional but aids readability.

**Reader semantics.**  The Z column on each row is used as
the array index: the reader stores the parameters at
position Z in the arrays, not at the sequential row number.
This makes the file order-independent and robust against
accidental reordering.  If a Z value appears that is outside
the range 1..`NUM_UFF_ELEMENTS`, the reader exits with an
error.

### 4.5 Bond Scale Factor

A new parameter `bond_parameter_scale` provides a global
multiplier for all bond force constants.  Its default value
(1.0) is defined in `condenserc.py`, loaded into the
`Condense` object by `assign_rc_defaults()`, and can be
overridden by a `bond_parameter_scale` keyword in the
condense.in input file (read by `parse_input_file()`):

  bond_parameter_scale 0.9

The value is dimensionless and multiplies every computed K_ij
before writing the LAMMPS Bond Coeffs section.  Values below
1.0 loosen all bonds; values above 1.0 stiffen them.

Only K_ij is scaled -- the equilibrium bond length r_ij is
left unchanged.  This lets the user tune overall bond rigidity
while preserving the equilibrium geometry of the system.  The
motivation is that UFF force constants are approximate by
nature (10-20 % inter-force-field variation is typical), so a
global rescaling provides a simple empirical knob for tuning
the dynamic behavior of the condensation simulation without
modifying the underlying database.

### 4.6 Impact on condense.py

The `BondData` class is restructured:

1. **Reading.**  `init_bond_data()` reads the per-element
   parameter table (r_i, Zstar_i, chi_i) and stores it in
   arrays indexed by atomic number Z.  This replaces the
   pair-enumerated `hooke_bond_coeffs` list.

2. **Querying.**  A new method `get_bond_params(z1, z2)`
   computes and returns (K_ij, r_ij) for any element pair
   using the UFF formulas from section 4.2.  Argument order
   does not matter (the formula is symmetric).

3. **Bond lookup in create_lammps_files.**  The existing
   linear scan over `hooke_bond_coeffs` is replaced by a
   single call to `get_bond_params(z1, z2)` -- O(1) per
   lookup instead of O(n).

4. **Bond lookup in normalize_types.**  The same linear
   scan appears a second time in `normalize_types()`,
   where unique bond types are matched to coefficients
   for rewriting LAMMPS files with unified type indices.
   This scan is replaced by `get_bond_params(z1, z2)` in
   exactly the same way as item 3.

5. **Scale factor plumbing.**  The default value of
   `bond_parameter_scale` (1.0) is defined in
   `condenserc.py` and loaded by `assign_rc_defaults()`.
   It can then be overridden by the `bond_parameter_scale`
   keyword in condense.in (read by `parse_input_file()`).
   The multiplier is applied to K_ij in **both** output
   paths: `create_lammps_files` (initial LAMMPS data file)
   and `normalize_types` (rewritten LAMMPS data file with
   unified type indices).  Both paths must produce the
   same scaled force constants.

6. **Error handling.**  If either Z exceeds the range of
   the parameter table, print the element symbol and Z
   number and exit with a clear message directing the user
   to extend the bond_parameters.dat table.

### 4.7 Backward Compatibility

The new `bond_parameters.dat` is not readable by the Perl
`BondData.pm` module, which expects the old `bonds.dat`
pair-listing format.  Since `condense.py` is the active
development path and `BondData.pm` belongs to the deprecated
Perl toolchain, this is accepted.  Users who still need the
Perl `condense` script can retain a local copy of the old
`bonds.dat` file.

**Build system.**  `src/data/CMakeLists.txt` must be updated
to install `bond_parameters.dat` instead of `bonds.dat` in
the DATABASES list.  The old `bonds.dat` is removed from the
install set.

### 4.8 Geometry-Derived Angle Parameters

#### 4.8.1 Motivation

The old `angles.dat` file listed 56 explicit triplet entries
covering only seven elements (H, B, C, N, O, Si).  Every
entry used a uniform spring constant k = 500.0 kcal/mol/rad^2
regardless of the element triplet.  Adding a new element
required manually enumerating every triplet and rest angle
it participates in -- an unsustainable maintenance burden and
a frequent source of "Cannot find angle in the database"
failures when the system contains elements outside the seven.

Section 4.8 of the previous design (preserved below in
section 4.8.2) analyzed why the per-element UFF strategy
that worked for bonds does not transfer directly to angles:
the same element triplet can have multiple physically
distinct rest angles (e.g., C-C-C at 60, 109.5, 120, and
180 degrees), and the UFF angle potential is a cosine
Fourier series rather than a simple harmonic.

The key insight is that the OLCAO bond analysis already
computes the actual bond angles for every atom in every
molecule.  These observed angles encode the real electronic
structure -- hybridization, strain, ring membership, and
neighbor effects -- for the specific system at hand.  They
are more accurate than any generic lookup table, and they
are already available at runtime.  The design below uses
these angles directly as equilibrium values, eliminating
the need for an external angle database entirely.

#### 4.8.2 Prior Analysis (retained for reference)

**Why per-element UFF does not transfer directly to angles.**

1. **Multiple rest angles per triplet.**  The same element
   triplet (e.g., C-C-C) appears in angles.dat with several
   distinct equilibrium angles (60 deg for cyclopropane,
   108 deg for cyclopentane, 180 deg for linear chains).
   A per-element UFF lookup gives one natural angle per
   vertex atom type (e.g., C_3 = 109.47 deg), which cannot
   distinguish these chemical environments.

2. **UFF angle potential form mismatch.**  The UFF angle
   bending potential is a cosine Fourier series (Rappe
   eq. 8), not a simple harmonic.  LAMMPS `angle_style
   harmonic` uses E = K (theta - theta_0)^2.  Adopting
   UFF angle parameters would require either approximating
   the cosine form as harmonic near theta_0 or switching
   to a different LAMMPS angle style -- both are significant
   changes beyond the data file.

3. **Complex force constant formula.**  The UFF angle force
   constant K_IJK (Rappe eq. 13) depends on the bond lengths
   of both arms (r_IJ and r_JK), all three effective charges,
   and the vertex atom's natural angle.  This is considerably
   more involved than the clean two-element bond formula
   K_ij = 332.06 * Zstar_i * Zstar_j / r_ij^3.

The geometry-derived approach (sections 4.8.3-4.8.9) resolves
all three issues: it uses observed angles (solving 1), feeds
them into a LAMMPS harmonic potential (avoiding 2), and uses
a simplified force constant formula based on bond stiffnesses
already computed by `get_bond_params()` (simplifying 3).

#### 4.8.3 Approach: Cluster Observed Angles by Triplet

For each molecule in the system, the OLCAO bond analysis
(via `bond_analysis.py`) already computes every bond angle.
The `create_lammps_files` method already iterates over these
angles and constructs triplet tags of the form
(Z_end1, Z_vertex, Z_end2).  Currently it searches
`angles.dat` for a matching entry.  The new approach replaces
that database lookup with the following procedure:

1. **Collect.**  For each angle instance, extract the
   full triplet (Z1, Zv, Z2) with Z1 <= Z2, and the
   observed angle theta_obs.

2. **Cluster by triplet.**  Group all angle instances
   that share the same (Z1, Zv, Z2) triplet.  Within
   each triplet group, sort the observed angles and
   greedy-merge values into a cluster while two
   conditions hold: (a) the candidate is within +/-
   `angle_cluster_tolerance` degrees of the running
   mean, and (b) the resulting cluster span (max - min)
   remains within `2 * angle_cluster_tolerance`.  The
   spread cap prevents a long chain of closely-spaced
   observations from silently sweeping values from
   opposite ends of a wide distribution into a single
   cluster.  When either condition fails, the current
   cluster is finalized and a new one begins at the
   candidate.  The cluster's rest angle theta_0 is the
   mean of its members.  The same spread cap is applied
   when `normalize_types()` re-clusters across sources
   (see 4.8.8 item 4a), so local and cross-source
   clustering use consistent semantics.

3. **Assign types.**  Each cluster becomes one LAMMPS
   angle type.  Every angle instance is assigned to the
   cluster whose mean it contributed to.

**Example.**  Suppose carbon vertex atoms yield observed
angles of 108.3, 109.1, 109.8, 120.2, 119.7, and 60.1
degrees, all for the C-C-C triplet, with
`angle_cluster_tolerance = 5.0`:

- Cluster 1: {60.1} -> theta_0 = 60.1 (ring)
- Cluster 2: {108.3, 109.1, 109.8} -> theta_0 = 109.1 (sp3)
- Cluster 3: {119.7, 120.2} -> theta_0 = 120.0 (sp2)

This produces three angle types instead of six individual
entries.

#### 4.8.4 Force Constant Formula

The angular spring constant K is computed from the UFF
per-element parameters already stored in
`bond_parameters.dat` (section 4.4).  The formula uses
the bond stiffnesses of the two arms:

  K_angle = C_angle * sqrt(K_bond_IJ * K_bond_JK)

where K_bond_IJ is the UFF harmonic bond force constant
for the (Z1, Zv) pair and K_bond_JK is for the (Zv, Z2)
pair, both obtained from `get_bond_params()`.  The
geometric mean captures the essential physics: stiffer
bonds produce stiffer angles.  The calibration constant
C_angle is dimensionless and converts bond stiffness
(kcal/mol/A^2) into an angular stiffness scale
(kcal/mol/rad^2).

Unlike the UFF bond constant (332.06, well-established
from the Rappe paper), C_angle is a project-specific
heuristic with no published source (see Provenance below).
It is therefore exposed as a user-tunable keyword in
condense.in:

  angle_stiffness_coeff 0.15

The default value (0.15) is defined in `condenserc.py`.
Together with `angle_parameter_scale` (section 4.8.5),
the user has two complementary controls:
`angle_stiffness_coeff` sets the base conversion from
bond stiffness to angle stiffness, while
`angle_parameter_scale` applies a uniform global
multiplier on top.  The final force constant written to
LAMMPS is:

  K_final = angle_stiffness_coeff
            * sqrt(K_bond_IJ * K_bond_JK)
            * angle_parameter_scale

**Provenance.**  This formula is a project-specific
heuristic, not drawn from a published force field.  The
full UFF angle bending force constant K_IJK (Rappe et al.
eq. 13) depends on the bond lengths of both arms, all
three effective charges, the equilibrium angle, and uses
a cosine Fourier expansion rather than a harmonic
potential.  Adopting the full UFF angle treatment would
require either switching LAMMPS to a cosine angle style
or performing a non-trivial harmonic approximation of the
Fourier series near each equilibrium angle.  The geometric
mean heuristic sidesteps both issues by staying within
the LAMMPS `angle_style harmonic` framework (Thompson et
al. 2022; E = K (theta - theta_0)^2) while still
producing element-dependent K values that track the
underlying bond stiffnesses.

**Calibration.**  Published harmonic angle force constants
for small organic molecules typically fall in the range
30-100 kcal/mol/rad^2.  For reference, the AMBER ff94
force field (Cornell et al. 1995) assigns C-C-C angles
K ~ 40 kcal/mol/rad^2 and H-C-H angles K ~ 35
kcal/mol/rad^2; the OPLS-AA force field (Jorgensen et al.
1996) gives similar values.  These are considerably softer
than the uniform k = 500 used in the old `angles.dat`.

Typical UFF bond force constants from `get_bond_params()`
are 200-700 kcal/mol/A^2.  For a C-C-C angle, both arms
give K_bond ~ 470 kcal/mol/A^2, so sqrt(470 * 470) = 470.
The default `angle_stiffness_coeff` of 0.15 yields
K_angle ~ 70 kcal/mol/rad^2, which is within the range
of published values.  Users should calibrate this value
against a known system (e.g., a small organic molecule
with published force field parameters).

**Note on the uniform k = 500 in the old database.**  The
old `angles.dat` used k = 500 for every entry.  This is
extremely stiff -- roughly 5-10x typical literature values.
It is not physically motivated; it appears to have been
chosen as a "rigid enough" default.  The computed K values
from the formula above will be significantly softer and
more physically realistic.  If the user needs the old
stiff behavior, `angle_parameter_scale` can be set to a
large value (e.g., 5.0-7.0).

#### 4.8.5 Angle Scale Factor

A new parameter `angle_parameter_scale` provides a global
multiplier for all computed angle force constants, following
the same pattern as `bond_parameter_scale` (section 4.5).
Its default value (1.0) is defined in `condenserc.py`,
loaded into the `Condense` object by `assign_rc_defaults()`,
and can be overridden by a keyword in condense.in:

  angle_parameter_scale 0.8

The value is dimensionless and multiplies every computed
K_angle before writing the LAMMPS Angle Coeffs section.
Values below 1.0 loosen all angular springs; values above
1.0 stiffen them.  Only K_angle is scaled -- the rest angle
theta_0 is left unchanged.

#### 4.8.6 Angle Cluster Tolerance

A new parameter `angle_cluster_tolerance` controls how
aggressively observed angles are merged into shared types.
Its default value (5.0 degrees) is defined in
`condenserc.py` and can be overridden in condense.in:

  angle_cluster_tolerance 3.0

**Clustering algorithm.**  Within each (Z1, Zv, Z2)
triplet group, the observed angles are sorted in
ascending order and then merged greedily: the first angle
starts a new cluster; each subsequent angle is added to
the current cluster if it falls within
`angle_cluster_tolerance` of the cluster's running mean,
otherwise it starts a new cluster.  This greedy approach
is simple, deterministic, and keeps the type count low.

Note that a chain of angles spaced just under the
tolerance apart (e.g., 105, 108, 111, 114 with a 5-degree
tolerance) will merge into one cluster because each new
member is compared to the running mean, not to the first
member.  This is the intended behavior: it favors fewer,
broader clusters, which reduces the angle type count and
lowers the risk of bond/react type-mismatch failures.

**Interaction with bond/react type count.**  LAMMPS
bond/react requires that the atom, bond, and angle types
in the pre- and post-reaction templates match the types
in the main data file.  Every distinct angle type in the
system increases the combinatorial space that must be
consistent across all files.  A larger tolerance produces
fewer, coarser angle types, which reduces the risk of
type-mismatch failures in bond/react.  A smaller tolerance
preserves finer geometric detail but creates more types.

The default of 5.0 degrees is a practical compromise.
For systems with many distinct molecular species or
complex reaction networks, increasing the tolerance to
8-10 degrees may be necessary to keep the type count
manageable.

#### 4.8.7 Look-Ahead Angles for Bond/React Products (deferred)

**The problem.**  The clustering procedure in section 4.8.3
discovers angle types from the initial molecular geometries.
But LAMMPS bond/react creates new bonds between molecules,
and those new bonds produce new angles that did not exist
in any isolated molecule.

Consider B12H12 and CH4.  In the isolated molecules, no
C-B bond exists, so no C-B-H or C-B-B angle is ever
observed.  After bond/react fires and creates a C-B bond,
the post-reaction template would need angle types for
every triplet that includes the new bond.  If those types
are not present in the LAMMPS data file's Angle Coeffs
section, bond/react would fail.

**Current state of the code.**  The Perl `makeReactions`
script (and its Python port `make_reactions.py`) adds one
new *bond* between the trigger atoms in the post-reaction
template but does **not** add any new *angles*.  The Perl
source (lines 2508-2513) states this explicitly:

  "In the future it might be necessary to *add* bond
   angles through the bonding atoms (after the S are
   removed), but presently we do not do that.  (Hence
   the bonded molecules may be too floppy.)"

A commented-out `addBondAngle` subroutine (Perl lines
2990-3141) shows that an attempt was started: it computes
angles via the law of cosines from post-reaction
coordinates, builds angle tags, and registers new angle
types.  But the subroutine was never activated.

Because the post-reaction templates do not currently
contain any new angles, there are no novel angle types
for condense.py to "look ahead" to.  The look-ahead
mechanism in condense.py and the angle-creation logic
in makeReactions are two halves of the same problem.

**Phased approach.**  This work is split into two steps
to keep each change testable:

1. **Step 1 (this design, sections 4.8.3-4.8.6 and
   4.8.8):**  Replace the angles.dat database with
   geometry-derived clustering and computed force
   constants.  The system works for all intra-molecular
   angles that already exist in the pre- and post-
   reaction templates.  Post-reaction bonding sites
   remain "floppy" (same as today) because no new
   angles are created by makeReactions.

2. **Step 2 (future work):**  Activate angle creation in
   `make_reactions.py` (porting and completing the
   commented-out Perl subroutine).  Once post-reaction
   templates carry the new angles, condense.py's
   `normalize_types()` will automatically pick them up
   during its template-scanning pass.  The rest angle
   theta_0 can be computed from the post-reaction atom
   coordinates (law of cosines, as the Perl prototype
   does), and K_angle can be computed from the same
   formula (section 4.8.4).  The clustering tolerance
   should be applied when deciding whether a novel
   post-reaction angle merges with an existing type or
   creates a new one.

Step 2 will also need to address whether new *bond* types
(not just angles) can appear in post-reaction templates.
The bond case is simpler because `get_bond_params()` can
always compute K and r0 for any element pair, but the
type must still be registered in the unified type list.

  >> DESIGN QUESTION D5a (deferred to step 2): When
  >> post-reaction angles are added to the templates,
  >> should novel angles always be added as distinct
  >> types to ensure the post-reaction geometry is
  >> exactly preserved, or should they be merged into
  >> the existing cluster list if they fall within
  >> `angle_cluster_tolerance` of an existing type?
  >> Merging keeps the type count down (reducing
  >> bond/react fragility), but exact preservation may
  >> matter for the product geometry.

#### 4.8.8 Impact on condense.py (step 1)

The `AngleData` class is eliminated.  No external data
file is read for angles.  The changes are:

1. **AngleData class: remove.**  The class and its
   `init_angle_data()` method are deleted from every
   source file that currently carries them -- both the
   copy in `condense.py` and the duplicated copy in
   `make_reactions.py` (lines ~113-193) along with its
   `self.angle_data` instantiation (line ~568).  All
   references to `self.angle_data` and
   `ad.hooke_angle_coeffs` are removed from
   `create_lammps_files()`, `normalize_types()`, and
   `make_reactions.py`'s template-emission path (the
   hooke_angle_coeffs scan near line 2518).

2. **Clustering step: add (shared helper).**  A new
   helper routine implements the cluster-by-triplet
   algorithm from section 4.8.3.  Input: the list of
   all angle instances with their (Z1, Zv, Z2,
   theta_obs) tuples from one producer's scope.
   Output: a list of local angle types, each with
   (Z1, Zv, Z2, theta_0) and an observation count,
   plus a mapping from each instance to its local
   type index.  This helper is shared by both
   producers -- `create_lammps_files()` in
   condense.py and the template-emission path in
   `make_reactions.py` -- so that local clustering
   semantics are byte-identical across sources.  The
   cross-source unification described in item 4 then
   operates on the per-source outputs.

3. **Force constant computation: add.**  For each
   angle type, compute K_angle using the formula from
   section 4.8.4, calling `get_bond_params()` for the
   two arm bond stiffnesses.  Apply
   `angle_stiffness_coeff` and `angle_parameter_scale`
   (both factor into the K_angle formula).  This runs
   locally in every producer; `normalize_types()`
   recomputes the same K values in item 4b, which is
   safe because K_angle depends only on the triplet
   and not on theta_0.

4. **normalize_types(): cross-source unification.**
   `normalize_types()` is the authoritative
   cross-source clusterer for angle types, not a
   passive consumer of producer-emitted theta_0
   values.  Each producer (`create_lammps_files()`
   and `make_reactions.py`) locally clusters its
   own observations and emits one type per local
   cluster, with the cluster-mean theta_0 carried
   in the tag tail "{theta_0_local} {t}".  Because
   the two producers see different observation
   populations of the same physical angle, their
   local theta_0 for a chemically identical triplet
   can differ by a few tenths of a degree.  A plain
   string comparison on the tag tail would split
   such cases into distinct types and break
   bond/react type-ID matching across lammps.dat
   and the reaction templates.  `normalize_types()`
   absorbs this drift in four steps:

   a. **Cross-source clustering.**  Collect every
      angle type emitted by every source -- the
      lammps.dat produced by `create_lammps_files()`
      and every reaction template produced by
      `make_reactions.py` -- each carrying
      (z1, zv, z2, theta_0_local, obs_count,
      source, local_type_id).  Group by canonical
      triplet (z1 <= z2).  Within each group, sort
      by theta_0_local and greedy-merge while the
      candidate is within `angle_cluster_tolerance`
      of the running cluster mean.  Apply a spread
      cap (max cluster span <= 2 *
      `angle_cluster_tolerance`) to prevent greedy
      chaining across a wide distribution.  The
      running mean, weighted by obs_count, is the
      final canonical theta_0 for the merged
      cluster.

   b. **Force constant computation.**  For each
      final cluster, compute K_angle via
      `get_bond_params()` using the formula from
      section 4.8.4, with `angle_stiffness_coeff`
      and `angle_parameter_scale` applied (same
      formula as item 3).  Because K_angle depends
      only on the arm bond stiffnesses, which are a
      function of the triplet (z1, zv, z2) alone
      and not of theta_0, the K_angle computed here
      matches whatever any producer locally
      computed for a member cluster -- cross-source
      merging does not alter K_angle, only theta_0.

   c. **Tag rewrite and type-ID remap.**  Each
      (source, local_type_id) pair maps to exactly
      one final cluster id.  Walk the Angles
      section of lammps.dat and every angle
      reference in every reaction template, and
      rewrite the per-angle type id to the global
      cluster id.  Rewrite the tag tail
      "{theta_0_local} {t}" to
      "{theta_0_final} {global_t}" so any tool that
      later inspects the tag sees a consistent
      value.  The rewrite is deterministic given
      the cluster map, so repeated runs on
      identical inputs produce byte-identical
      output.

   d. **Cluster-map diagnostic.**  Emit a log file
      or log section listing, for every final
      cluster: global id, canonical theta_0,
      (z1, zv, z2), and every contributing
      (source, local_theta_0, obs_count) tuple.
      Students debugging a bond/react type mismatch
      should be able to open this file and see at a
      glance which observations were merged, which
      were split, and why.  This diagnostic is the
      main debuggability payback for routing all
      clustering through `normalize_types()` rather
      than accepting per-producer string tags.

5. **Parameter plumbing: add.**  Introduce three new
   force-field parameters on the `Condense` class with
   hardcoded defaults in `Condense.__init__`:
   `angle_stiffness_coeff` (default 0.15),
   `angle_parameter_scale` (default 1.0), and
   `angle_cluster_tolerance` (default 5.0).  These
   follow the `bond_parameter_scale` precedent:
   force-field parameters are deliberately kept out of
   `condenserc.py` and `ScriptSettings` to avoid
   cluttering the CLI with rarely-touched knobs.  User
   overrides are accepted through matching keywords in
   `condense.in`, parsed by `parse_input_file()`.  The
   values are applied wherever K_angle is computed --
   in `create_lammps_files()`, in the template-
   emission path of `make_reactions.py`, and in
   `normalize_types()`.  `angle_cluster_tolerance`
   additionally governs wherever local or cross-
   source clustering happens (the same three sites,
   plus the spread-cap policy `2 * tolerance`).

6. **angles.dat: retire.**  Remove from
   `src/data/CMakeLists.txt` DATABASES list.  Remove
   from `share/` install target.  The file may be kept
   in the repository for historical reference but is no
   longer read by any code path.

Note: the look-ahead pass described in section 4.8.7 is
deferred to step 2.  Step 1 handles all angles that
already exist in the pre- and post-reaction templates
(which is the same set that the old angles.dat handled).
The bonding-site "floppiness" is unchanged from the
current behavior.

#### 4.8.9 Backward Compatibility

The new approach is not backward-compatible with the old
`angles.dat` format.  Since `condense.py` is the active
development path and the Perl toolchain is deprecated,
this is accepted (same reasoning as section 4.7 for bonds).

The behavioral difference is that angle rest values will
now come from the system's own geometry rather than a
curated database.  For well-prepared input structures
(which is the expected use case), this produces identical
or better rest angles.  For poorly prepared structures,
the rest angles will reflect the input geometry -- which
is arguably more honest than imposing idealized angles
that the structure does not actually have.

Users who relied on the old uniform k = 500 behavior can
approximate it by setting `angle_parameter_scale` to a
large value, though the per-triplet variation in K will
still be present.

---

## References

P. E. Bloechl, O. Jepsen, O. K. Andersen, "Improved
tetrahedron method for Brillouin-zone integrations,"
Phys. Rev. B 49, 16223 (1994). Key equations:
- Eqs. 14-16: analytic DOS formulas (total per
  tetrahedron)
- Eqs. 18-21: corner weights. Cumulative form
  `cornerIntgWt_LAT` from `bloechlCornerWeights` for
  integrated properties; energy derivatives
  `cornerDOSWt_LAT` from `bloechlCornerDOSWt` for
  energy-resolved DOS/PDOS
- Eqs. 22-24: correction terms for improved accuracy

A. K. Rappe, C. J. Casewit, K. S. Colwell, W. A.
Goddard III, W. M. Skiff, "UFF, a Full Periodic Table
Force Field for Molecular Mechanics and Molecular
Dynamics Simulations," J. Am. Chem. Soc. 1992, 114,
10024-10035.  DOI: 10.1021/ja00051a040
- Table 1: per-element parameters (r_i, Zstar_i,
  chi_i) used for bond stretching (section 4) and
  as inputs to the angle K heuristic (section 4.8.4)
- Eq. 2: natural bond length with electronegativity
  correction
- Eq. 3: bond stretching force constant formula
- Eq. 8: angle bending potential (cosine Fourier
  expansion -- not used directly; section 4.8.2
  explains why the harmonic approximation is preferred)
- Eq. 13: full UFF angle bending force constant K_IJK
  (not adopted; the geometric-mean heuristic in
  section 4.8.4 is used instead)

W. D. Cornell, P. Cieplak, C. I. Bayly, I. R. Gould,
K. M. Merz Jr., D. M. Ferguson, D. C. Spellmeyer,
T. Fox, J. W. Caldwell, P. A. Kollman, "A Second
Generation Force Field for the Simulation of Proteins,
Nucleic Acids, and Organic Molecules," J. Am. Chem.
Soc. 1995, 117, 5179-5197.  DOI: 10.1021/ja00124a002
- Referenced in section 4.8.4 for calibration context:
  typical harmonic angle force constants for organic
  molecules (C-C-C ~ 40, H-C-H ~ 35 kcal/mol/rad^2)

W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives,
"Development and Testing of the OPLS All-Atom Force
Field on Conformational Energetics and Properties of
Organic Liquids," J. Am. Chem. Soc. 1996, 118,
11225-11236.  DOI: 10.1021/ja9621760
- Referenced in section 4.8.4 for calibration context:
  independent confirmation that organic angle force
  constants fall in the 30-100 kcal/mol/rad^2 range

A. P. Thompson, H. M. Aktulga, R. Berger, et al.,
"LAMMPS - a flexible simulation tool for particle-based
materials modeling at the atomic, meso, and continuum
scales," Comp. Phys. Comm. 2022, 271, 108171.
DOI: 10.1016/j.cpc.2021.108171
- `angle_style harmonic`: E = K (theta - theta_0)^2
  convention used throughout section 4.8
