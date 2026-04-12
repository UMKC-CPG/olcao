# Architecture

> **Document hierarchy:** VISION -> **ARCHITECTURE** -> DESIGN
> -> PSEUDOCODE -> Code. For goals and principles, see
> `VISION.md`.

---

## 1. Repository Layout

```
src/
  uolcao/              Primary production program
    OLCAO.F90          Top-level dispatcher
    kpoints.f90        O_KPoints module (mesh, weights, phases)
    dos.F90            O_DOS module (TDOS, PDOS)
    bond.F90           Bond order (Mulliken overlap population)
    populate.F90       Electron population / occupancy
  makeKPoints/
    makekpoints.F90    Standalone k-point mesh generator
  kinds.f90            Shared precision kinds
  constants.f90        Shared physical/mathematical constants
dev/
  VISION.md            Goals and principles
  ARCHITECTURE.md      This document
  DESIGN.md            Algorithmic design
  PSEUDOCODE.md        Algorithm specifications
  TODO.md              Task list by level
```

---

## 2. Module Map

Modules directly affected by the current development:

- **O_KPoints** (`kpoints.f90`): Stores the k-point mesh,
  weights, and phase factors. Declares `kPointIntgCode`
  (0=Gaussian, 1=LAT). LAT data already implemented:
  `numTetrahedra`, `tetraVol`, `tetrahedra(:,:)`,
  `fullKPToIBZKPMap(:)`, `generateTetrahedra`,
  `computeTetraVol`, `initializeKPointMesh` (with IBZ
  folding). Phase F additions for correct IBZ unfolding
  of eigenvector-dependent quantities:
  `fullKPToIBZOpMap(:)` (which point group operation
  mapped each full-mesh k-point to its IBZ
  representative), `atomPerm(:,:)` (atom permutation
  table under each point group operation), and
  `buildAtomPerm` (builds the table from point ops and
  fractional atom positions).
- **O_DOS** (`dos.F90`): Two DOS paths: `computeIterationTDOS`
  (in-SCF convergence monitoring) and `computeDOS` (full
  TDOS/PDOS post-processing). Will gain a LAT branch in
  `computeDOS` dispatched on `kPointIntgCode`.
- **bond.F90**: Computes bond order and effective charge
  (Q*). Already uses `electronPopulation_LAT` when
  `kPointIntgCode == 1`. Phase F additions: the
  accumulation loop will buffer per-IBZ-kpoint
  projections and distribute them across the star of
  each IBZ k-point using `atomPerm` and
  `fullKPToIBZOpMap` from O_KPoints, giving correct
  per-atom Q* and per-pair bond order under IBZ
  reduction.
- **populate.F90** (`O_Populate`): Computes
  `electronPopulation` (Gaussian/Fermi-filling path)
  which folds in `kPointWeight`.
  `computeElectronPopulation_LAT` (already implemented)
  produces `electronPopulation_LAT` using tetrahedra
  data from O_KPoints. Both occupation arrays serve the
  same downstream consumers (bond order, effective
  charge); the integration method (`kPointIntgCode`)
  selects which one is used. Upstream of bond order and
  effective charge.
- **makeKPoints** (`makekpoints.F90`): Legacy standalone
  program that generates explicit k-point lists with IBZ
  reduction (`foldMesh`). Separate executable, not linked
  into uolcao. Being retired -- its mesh-building and IBZ
  reduction capabilities have been ported into uolcao's
  `initializeKPointMesh`. Retained only for backward
  compatibility; new workflows should not depend on it.
- **makeinput.py** (`src/scripts/makeinput.py`): Top-level
  orchestrator that prepares all input files for OLCAO. For
  k-points it supports three pipelines:
  1. *Mesh mode* (`-kp`, `-scfkp`, `-pscfkp`): writes a
     style-code-1 k-point file with axial counts and shift.
     OLCAO builds the full mesh internally and reduces to
     the IBZ, giving it full access to the mesh topology
     and symmetry maps needed for correct decomposition
     properties.
  2. *Density mode* (`-kpd`, `-scfkpd`, `-pscfkpd`):
     writes a style-code-2 k-point file with density and
     shift. OLCAO computes axial counts from the density
     and reciprocal cell geometry, then builds and reduces
     the mesh identically to mesh mode.
  3. *Explicit mode* (style code 0 in the k-point file):
     a pre-built list of k-points with weights, read
     directly by OLCAO with no internal mesh construction.
     Supported for special cases (e.g., hand-crafted
     k-point sets), but OLCAO emits a prominent warning
     that decomposition properties (effective charge,
     bond order, PDOS) will not be correct unless the
     user has taken extreme care to provide a symmetric
     mesh. Not produced by makeinput.

---

## 3. Dependency Graph

```
makeinput.py (top-level orchestrator)
  +-- (mesh mode: -kp)
  |     writes kp-scf.dat / kp-pscf.dat (style 1)
  +-- (density mode: -kpd)
  |     writes kp-scf.dat / kp-pscf.dat (style 2)
  +-- makeKPoints (legacy, no longer called by makeinput)

OLCAO.F90 (top-level dispatcher)
  +-- O_KPoints (kpoints.f90)
  |     +-- O_Lattice (recipCellVolume, invRealVectors)
  |     +-- O_AtomicSites (atom positions for buildAtomPerm)
  +-- O_DOS (dos.F90)
  |     +-- O_KPoints (eigenvalues, tetrahedra, weights)
  |     +-- O_PSCFIntg / O_SCFIntg (eigenvectors via HDF5)
  +-- populate.F90 (O_Populate)
  |     +-- O_KPoints (kPointWeight, tetrahedra,
  |     |     eigenvalues, fullKPToIBZKPMap for LAT)
  +-- bond.F90
  |     +-- O_Populate (electronPopulation or
  |     |     electronPopulation_LAT)
  |     +-- O_KPoints (fullKPToIBZKPMap,
  |     |     fullKPToIBZOpMap, atomPerm,
  |     |     numFullMeshKP for star distribution)
  |     +-- O_PSCFIntg (eigenvectors via HDF5)
```

---

## 4. Build System

- Fortran 90; CMake out-of-source build in `build/`
- Compiler set via `$FC` (must match HDF5-compiled compiler)
- Install prefix via `$OLCAO_DIR`
- HDF5 required for primary I/O in uolcao
- Supported compilers: gfortran, ifort

---

## 5. Key Existing Infrastructure

These elements exist in the code and are leveraged by
the current design:

- `kPointStyleCode`: integer controlling how k-points
  are specified (kpoints.f90)
  - 0 = explicit list (legacy; not produced by
    makeinput; OLCAO emits a warning that
    decomposition properties may be incorrect)
  - 1 = axial counts + shift (mesh built internally;
    primary mode for `-kp` / `-scfkp` / `-pscfkp`)
  - 2 = minimum density + shift (mesh built internally
    via `computeAxialKPoints`; used by `-kpd` etc.)
- `kPointIntgCode`: integer read from input
  (0=Gaussian, 1=LAT) (kpoints.f90)
- `readKPoints`: reads the k-point input file; branches
  on `kPointStyleCode`. For style codes 1 and 2, reads
  point group operations into `abcPointOps`
  (kpoints.f90)
- `computeAxialKPoints`: converts `minKPointDensity`
  into `numAxialKPoints(:)` using `recipMag` and
  `recipCellVolume` (kpoints.f90)
- `initializeKPointMesh`: builds the uniform mesh from
  `numAxialKPoints` and `kPointShift`; when
  `applySymmetry=1`, folds the mesh to the IBZ using
  `abcRecipPointOps` and saves `fullKPToIBZKPMap`
  (kpoints.f90)
- `computeRecipPointOps`: converts the abc point group
  ops from the kpoint file into reciprocal-space abc
  ops using `realVectors` and `recipVectors` from
  O_Lattice (kpoints.f90)
- `generateTetrahedra`: tiles the full uniform mesh
  with tetrahedra (6 per sub-cube) using
  `getIndexFromIndices` for periodic wrapping.
  Called from `initializeKPoints` when
  `kPointIntgCode == 1` (kpoints.f90)
- `getIndexFromIndices(a,b,c)`: converts mesh indices
  to the linear k-point index (kpoints.f90)
- `energyEigenValues(n, i, spin)`: fully in memory
  after SCF
- `recipCellVolume`, `invRealVectors`: available from
  O_Lattice after `initializeLattice`
- `atomSites(:)%cartPos`, `atomSites(:)%atomTypeAssn`:
  atom positions and type assignments from
  O_AtomicSites, available after input parsing

---

## 6. Compute Architecture Direction

### 6.1 Configurable Precision

The program will support compilation in both double
precision (64-bit) and single precision (32-bit) via a
compile-time kind parameter (e.g., `wp` for "working
precision") defined in `kinds.f90`. All floating-point
declarations, literal constants, and MPI/HDF5 type tags
must use this parameter so that a single preprocessor
flag or CMake option switches the entire build.

**Motivation.** Single precision doubles SIMD throughput
on CPU (8 floats vs. 4 doubles in AVX-256) and is
essential for GPU performance, where double-precision
throughput is 2x slower on data-center GPUs and up to
32x slower on consumer hardware.

**Numerical stability.** Single precision provides only
~7 significant digits, which raises concerns in
accumulation loops where small contributions are added
to large sums. The lattice loop in the integral engine
(`gaussOverlapOL` and siblings in `integrals.F90`) sums
contributions from nearest to farthest replicated cells;
far-away terms are individually small but collectively
significant. At double precision this is safe (~15
digits of headroom). At single precision, compensated
(Kahan) summation or smallest-first accumulation order
may be required for specific accumulation sites.
Stability analysis must accompany the precision switch,
targeting at minimum the lattice-sum accumulators and
the eigenvalue solver interface.

### 6.2 Inner-Loop Vectorization

The alpha-pair iteration inside the lattice loop of
`gaussOverlapOL` (and the analogous nuclear and
three-center routines) is currently control-flow-
dominated: a `do while (.true.)` loop with mode
switching, early exits, and per-pair negligibility
tests. This structure prevents SIMD vectorization.

**Target restructure.** Separate the selection phase
(which alpha pairs survive the `alphaDist` threshold for
a given atom pair and lattice distance) from the compute
phase (evaluate the Gaussian overlap integral). The
selection phase produces a packed list of surviving alpha
pair indices; the compute phase processes that list in a
tight, branchless loop amenable to SIMD.

This "gather, filter, compute in bulk" pattern is the
same restructure needed for GPU offload, so the two
goals reinforce each other.

### 6.3 GPU Offload

Long-term, the integral engine and eigenvalue solver are
candidates for GPU acceleration. The restructured inner
loops from 6.2 translate almost directly to GPU kernels:
the packed alpha-pair list becomes a kernel launch over
surviving pairs. Single precision from 6.1 is required
to achieve full GPU throughput.

**Practical path:**
1. Introduce `wp` kind parameter (6.1) and validate
   numerical stability at single precision.
2. Restructure alpha-pair loops into gather/compute
   phases (6.2); verify SIMD vectorization on CPU.
3. Offload the compute phase to GPU (OpenACC, CUDA
   Fortran, or OpenMP target) as a later step.

Each phase is independently useful: phase 1 halves
memory footprint and improves cache behavior; phase 2
speeds up CPU execution; phase 3 adds GPU capability.

### 6.4 Reference: Prior Vectorized Integrals

An earlier vectorized implementation exists outside this
repository at:

```
/home/rulisp/lewis/CPG/cpg-repo/v34/src/olcao/
  integralsSCF.vec.F90    (SCF integrals)
  integralsPSCF.vec.F90   (post-SCF integrals)
  gaussIntegrals.vec.f90  (vectorized primitives)
```

The vectorized `gaussIntegrals.vec.f90` primitives are
also tracked in this repo under `src/olcao/`,
`src/uolcao/`, and `src/upolcao/`.

These files demonstrate the gather/compute separation
described in 6.2. The alpha-pair selection loop
(branchy, scalar) collects surviving pairs into
`lmAlphaPairs` organized by angular momentum type. A
packing step linearizes them into `orderedAlphaPairs`
with `segIndices` marking contiguous same-type regions.
A single batched call (`overlap2CIntgVec`) then
processes all pairs, enabling SIMD over the contiguous
segments. This pattern is the starting point for the
restructure in A4.
