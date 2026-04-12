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
