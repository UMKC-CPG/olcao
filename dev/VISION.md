# Vision

## Purpose

OLCAO (Orthogonalized Linear Combination of Atomic Orbitals) is an
all-electron electronic structure code that uses periodic boundary
conditions and an LCAO basis to compute electronic properties of a
broad range of material systems: crystals, amorphous solids,
nanoparticles, molecules, interfaces, and grain boundaries. It is
an academic code used for research and student training.

The current development focus is on two closely related issues in
Brillouin-zone integration: implementing the Linear Analytic
Tetrahedral (LAT) method and correcting the treatment of
eigenvector-dependent quantities under IBZ symmetry reduction.

## Goals

1. **Implement LAT k-point integration.** Replace the current
   Gaussian-broadening approach for Brillouin-zone integration with
   the Linear Analytic Tetrahedral method (Bloechl, Jepsen, &
   Andersen, PRB 49, 16223, 1994). This provides exact analytic
   integration within each tetrahedron, eliminating the arbitrary
   broadening parameter and improving accuracy at lower k-point
   densities.
2. **Correct partial properties under IBZ reduction.** Ensure that
   eigenvector-dependent quantities (PDOS, bond order, effective
   charge) are computed correctly when the k-point mesh is reduced
   to the irreducible Brillouin zone. The current approach of
   multiplying IBZ contributions by star multiplicity is incorrect
   for these quantities.

## Design Principles

1. **Eigenvalue symmetry != eigenvector symmetry.** Point-group
   operations guarantee e_n(k) = e_n(Rk), but eigenvectors at Rk
   are related to those at k by a basis transformation that mixes
   orbital coefficients. Any integration scheme that weights
   eigenvector-dependent quantities by IBZ star multiplicity alone
   is incorrect. This was confirmed empirically: bond orders for
   symmetry-equivalent O atoms in KNbO3 disagreed under IBZ
   reduction but agreed with the full mesh.
2. **Full mesh for post-processing.** Use the IBZ for the expensive
   SCF diagonalization, but reconstruct the full BZ mesh before
   computing partial properties. For eigenvalue-only quantities
   (TDOS), unfolding eigenvalues by symmetry suffices. For
   eigenvector-dependent quantities, either use the full mesh
   directly or apply atom-permutation corrections.
3. **Phased implementation.** Implement LAT TDOS first (eigenvalues
   only), validate against Gaussian-broadened results at high
   k-point density, then extend to integrated partial properties
   via electronPopulation_LAT, and finally to energy-resolved
   PDOS with cornerIntgWt_LAT. Each
   phase must be validated before proceeding to the next.
