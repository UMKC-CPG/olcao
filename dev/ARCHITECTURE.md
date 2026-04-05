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
  weights, and phase factors. Already declares
  `kPointIntgCode` (0=Gaussian, 1=LAT) and has a
  `generateTetrahedra` stub. Will gain: `numTetrahedra`,
  `tetraVol`, `tetrahedra(:,:)`.
- **O_DOS** (`dos.F90`): Two DOS paths: `computeIterationTDOS`
  (in-SCF convergence monitoring) and `computeDOS` (full
  TDOS/PDOS post-processing). Will gain a LAT branch in
  `computeDOS` dispatched on `kPointIntgCode`.
- **bond.F90**: Computes bond order using
  `electronPopulation`. Will use
  `electronPopulation_LAT` when
  `kPointIntgCode == 1`.
- **populate.F90** (`O_Populate`): Computes
  `electronPopulation` (Gaussian/Fermi-filling path) which
  folds in `kPointWeight`. Will gain
  `computeElectronPopulation_LAT`, which produces
  `electronPopulation_LAT` using tetrahedra data passed
  from O_KPoints. Both occupation arrays serve the same
  downstream consumers (bond order, effective charge);
  the integration method (`kPointIntgCode`) selects which
  one is used. Upstream of bond order and effective charge.
- **makeKPoints** (`makekpoints.F90`): Standalone program that
  generates the k-point file. Has a complete IBZ reduction
  (`foldMesh`). This is a separate executable, not a module
  linked into uolcao.
- **makeinput.py** (`src/scripts/makeinput.py`): Top-level
  orchestrator that prepares all input files for OLCAO. For
  k-points it supports two pipelines:
  1. *Mesh mode* (`-kp`, `-scfkp`, `-pscfkp`): writes
     `kpSpecs.dat` and calls `makeKPoints` to produce an
     explicit k-point list (style code 0).
  2. *Density mode* (`-kpd`, `-scfkpd`, `-pscfkpd`): bypasses
     `makeKPoints` entirely and writes k-point files directly
     with style code 2 (density + shift). The mesh is computed
     later inside `uolcao` by `computeAxialKPoints`.

---

## 3. Dependency Graph

```
makeinput.py (top-level orchestrator)
  +-- makeKPoints (mesh mode: -kp)
  |     writes kpSpecs.dat -> reads kpSpecs.out (style 0)
  +-- (density mode: -kpd)
        writes kp-scf.dat / kp-pscf.dat directly (style 2)

OLCAO.F90 (top-level dispatcher)
  +-- O_KPoints (kpoints.f90)
  |     +-- O_Lattice (recipCellVolume for tetraVol)
  +-- O_DOS (dos.F90)
  |     +-- O_KPoints (eigenvalues, tetrahedra, weights)
  |     +-- O_PSCFIntg / O_SCFIntg (eigenvectors via HDF5)
  +-- populate.F90 (O_Populate)
  |     +-- O_KPoints (kPointWeight, tetrahedra,
  |     |     eigenvalues, fullToIBZMap for LAT path)
  +-- bond.F90
  |     +-- O_Populate (electronPopulation or
  |     |     electronPopulation_LAT)
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

These elements already exist in the code and are leveraged by
the current design:

- `kPointStyleCode`: integer controlling how k-points are
  specified -- kpoints.f90 lines ~17-27
  - 0 = explicit list (written by makeKPoints)
  - 1 = axial counts + shift (mesh built internally)
  - 2 = minimum density + shift (mesh built internally
    via `computeAxialKPoints`)
- `kPointIntgCode`: integer read from input (0=Gaussian, 1=LAT)
  -- kpoints.f90 lines ~28-41
- `readKPoints`: reads the k-point input file; branches on
  `kPointStyleCode` -- kpoints.f90 lines ~95-196
- `computeAxialKPoints`: converts `minKPointDensity` into
  `numAxialKPoints(:)` using `recipMag` and
  `recipCellVolume` -- kpoints.f90 lines ~816-891
- `initializeKPointMesh`: builds the uniform mesh from
  `numAxialKPoints` and `kPointShift`; when
  `applySymmetry=1`, folds the mesh to the IBZ using
  `abcRecipPointOps` -- kpoints.f90
- `computeRecipPointOps`: converts the abc point group ops
  read from the kpoint file into reciprocal-space abc ops
  using `realVectors` and `recipVectors` from O_Lattice
  -- kpoints.f90
- `generateTetrahedra`: stub with conceptual comments but empty
  body -- kpoints.f90 lines ~986-1019
- `getIndexFromIndices(a,b,c)`: converts mesh indices to the
  linear k-point index -- kpoints.f90 lines ~965-983
- `readKPoints` already calls `generateTetrahedra` when
  `kPointIntgCode == 1` -- kpoints.f90 lines ~192-194
- `energyEigenValues(n, i, spin)`: fully in memory after SCF
- `recipCellVolume`: available from O_Lattice after
  `initializeLattice`
