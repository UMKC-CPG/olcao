# Task List

> **Document hierarchy:** Tasks are organized by the level of
> the design chain they affect. Each item cites the relevant
> document section.

---

## VISION

- [ ] V1. Confirm design principles against additional test
  cases beyond KNbO3 (VISION Principles.1)

---

## ARCHITECTURE

- [x] A1. Add numTetrahedra, tetraVol, tetrahedra(:,:),
  numFullMeshKP, fullKPToIBZKPMap(:) to O_KPoints module
  (ARCHITECTURE 2, DESIGN 2.5)
- [x] A2. electronPopulation_LAT lives in populate.F90
  (O_Populate), alongside electronPopulation. Tetrahedra
  data passed from O_KPoints as arguments.
  (ARCHITECTURE 2)
- [ ] A3. Introduce compile-time working precision kind
  (`wp`) in kinds.f90; propagate to all floating-point
  declarations, constants, and I/O type tags so a single
  build flag switches between double and single precision
  (ARCHITECTURE 6.1)
- [ ] A4. Restructure alpha-pair inner loops in
  integrals.F90 (gaussOverlapOL and siblings): separate
  selection phase (which pairs survive alphaDist test)
  from compute phase (evaluate integral), producing a
  packed list for branchless SIMD execution
  (ARCHITECTURE 6.2)
- [ ] A5. Identify accumulation sites that require
  compensated (Kahan) summation or reordering for
  numerical stability at single precision; implement
  guards and validate (ARCHITECTURE 6.1)
- [ ] A6. GPU offload of restructured integral compute
  phase via OpenACC, CUDA Fortran, or OpenMP target
  (ARCHITECTURE 6.3)

---

## DESIGN

- [x] D1. Resolve memory strategy for PDOS projection array
  P(alpha, n, k): store at IBZ k-points only, apply
  atom permutation on-the-fly during tetrahedron
  corner assembly (DESIGN 1.4)
- [x] D2. Verify that electronPopulation_LAT correctly
  replaces electronPopulation for bond order in all
  cases (DESIGN 1.5, 2.5) -- resolved: weight is
  correct but accumulation needs atom permutation
- [x] D3. Design IBZ warning/detection for Gaussian PDOS
  path -- resolved: makeinput mesh mode now writes
  style code 1 (axial counts), so OLCAO builds the
  full mesh internally for both Gaussian and LAT.
  Style code 0 (explicit list) retained for special
  cases with a prominent warning (DESIGN 2.6)
- [x] D4. Write PSEUDOCODE for atom permutation fix:
  atomPerm table (section 4), fullKPToIBZOpMap save
  (section 5), corrected Q* (section 6), corrected
  bond order (section 7) (DESIGN 2.4, 2.5)

---

## PSEUDOCODE

- [x] P1. Transcribe exact Bloechl middle-range DOS formula
  for e2 <= E < e3 (PSEUDOCODE 2, Bloechl eqs. 14-16)
- [x] P2. Transcribe cornerIntgWt_LAT formulas for partial
  properties (PSEUDOCODE 3a, derived from first principles)

---

## CODE

### Phase A -- LAT TDOS (eigenvalues only)

- [x] C1a. Save fullKPToIBZKPMap from initializeKPointMesh
  IBZ folding (preserve kPointTracker as module data)
- [x] C1b. Move generateTetrahedra call from readKPoints
  to initializeKPoints (after mesh is built)
- [x] C1c. Complete generateTetrahedra in kpoints.f90
  (PSEUDOCODE 1)
- [x] C2. Compute tetraVol in initializeKPoints after
  lattice initialization (DESIGN 1.2)
- [x] C3. Implement computeTDOS_LAT in dos.F90 using
  fullKPToIBZKPMap for eigenvalue unfolding (PSEUDOCODE 2)
- [x] C4. Validate LAT TDOS against Gaussian broadening at
  high k-point density
- [ ] C4a. Add Bloechl correction terms (eqs. 22-24) to
  computeTDOS_LAT for improved accuracy at lower
  k-point densities (DESIGN 1.3)

### Phase B -- electronPopulation_LAT (integrated properties)

- [x] C5. Implement computeElectronPopulation_LAT
  (PSEUDOCODE 3) -- in populate.F90
- [x] C6. Modify computeBond to use
  electronPopulation_LAT when kPointIntgCode == 1
  (DESIGN 1.5) -- dispatch via statePopulation local
- [x] C7. Modify effective charge in computeBond to use
  electronPopulation_LAT (DESIGN 1.5) -- handled by
  C6 statePopulation dispatch with 2/spin factor

### Phase C -- LAT PDOS (energy-resolved)

- [x] C8. Implement bloechlCornerDOSWt subroutine:
  per-corner DOS density weights, used by both TDOS
  and PDOS (DESIGN 1.3, PSEUDOCODE 2a)
- [x] C8a. Refactor computeTDOS_LAT to call
  bloechlCornerDOSWt + sum, replacing inline
  dosContrib formulas (DESIGN 1.3, PSEUDOCODE 2)
- [x] C8b. Fix deltaDOS * hartree unit in integrated-
  area diagnostic in computeTDOS_LAT and computeDOS
  (DESIGN 1.3)
- [x] C9. Fix integratePDOS_LAT to call
  bloechlCornerDOSWt instead of bloechlCornerWeights
  (DESIGN 1.4, PSEUDOCODE 8.3)
- [x] C9a. Validate: rerun KNbO3 LAT PDOS and verify
  Spin States Calculated ≈ Spin States Expected and
  LAT PDOS matches Gaussian PDOS shape/magnitude --
  validated: Spin States 76.06 vs 76.86 expected
  (1% gap from band-edge effects), per-atom electron
  counts match Gaussian, PDOS shape correct. Minor
  O-atom symmetry spread (<0.1% integrated) noted.

### Phase D -- Density-based k-point input (DESIGN 3)

- [x] C11. Add `-kpd`, `-scfkpd`, `-pscfkpd` CLI options to
  makeinput.py argument parser (DESIGN 3.1)
- [x] C12. Add `_write_density_kp_file()` to makeinput.py that
  writes style-code-2 k-point files directly (DESIGN 3.4)
- [x] C13. Modify `_make_kp()` in makeinput.py to dispatch
  between mesh mode (makeKPoints) and density mode
  (direct write); all-or-nothing (DESIGN 3.2, 3.3)
- [x] C14. Handle summary output for density mode (print
  density value instead of kp count) (DESIGN 3.5)
- [x] C15. Validate end-to-end: makeinput.py -kpd produces
  a file that uolcao readKPoints correctly parses as
  style code 2
- [x] C16. Embed point group operations in style-code-2
  kpoint file (_extract_point_ops + updated writer)
- [x] C17. Read point group ops in readKPoints style-2
  branch (NUM_POINT_OPS + POINT_OPS labels)
- [x] C18. Port computeRecipPointOps into kpoints.f90
  (convert abc point ops to reciprocal-space)
- [x] C19. Port IBZ mesh folding into initializeKPointMesh
  (foldMesh algorithm from makeKPoints)
- [x] C20. Wire computeRecipPointOps call into
  initializeKPoints style-2 path
- [x] C21. Validate IBZ reduction: compare density-mode
  kpoint count against makeKPoints for same mesh

### Phase E -- Mesh-mode conversion and style-code-0 warning

- [x] C22. Convert makeinput.py mesh mode (`-kp`, `-scfkp`,
  `-pscfkp`) to write style-code-1 k-point files
  (axial counts + shift + point ops) instead of
  calling makeKPoints (ARCHITECTURE 2, DESIGN 2.6)
- [x] C23. Add style-code-1 branch to readKPoints in
  kpoints.f90: read axial counts, shift, and point
  ops; wire computeRecipPointOps into
  initializeKPoints style-1 path (ARCHITECTURE 5)
- [x] C24. Add prominent warning in initializeKPoints
  when kPointStyleCode == 0 that decomposition
  properties may be incorrect (DESIGN 2.6)
- [x] C25. Validate: makeinput.py -kp produces a style-
  code-1 file that uolcao parses and reduces
  correctly -- validated by user test run

### Phase F -- Atom permutation fix (DESIGN 2.4)

- [x] C26. Save fullKPToIBZOpMap in initializeKPointMesh
  folding loop: store operation index m when a match
  is found, store 1 (identity) for IBZ representatives
  (PSEUDOCODE 5) -- also renamed fullToIBZMap to
  fullKPToIBZKPMap across kpoints.f90, dos.F90, and
  populate.F90 for consistency with design docs
- [x] C27. Implement buildAtomPerm in atomicSites.f90:
  build atomPerm(numPointOps, numAtomSites) from point
  group operations and fractional atom positions
  (PSEUDOCODE 4) -- data and builder both in
  O_AtomicSites; imports abcPointOps from O_KPoints
- [x] C28. Wire buildAtomPerm call into OLCAO.F90 (not
  initializeKPoints, to avoid circular dependency
  between O_AtomicSites and O_KPoints) for style
  codes 1 and 2, in both SCF and PSCF paths
- [x] C29. Corrected Q* accumulation in computeBond:
  buffer per-atom projection per IBZ kpoint, then
  distribute across star via atomPerm (PSEUDOCODE 6).
  Also applies same star distribution to
  atomOrbitalCharge (per-atom, per-QN_l). Unified
  code path for all style codes: style code 0 now
  sets up trivial identity IBZ maps in kpoints.f90
  (numPointOps=1, identity abcPointOps, identity
  fullKPToIBZKPMap/OpMap) so buildAtomPerm and star
  distribution work without special-casing
- [x] C30. Corrected bond order accumulation in
  computeBond: buffer per-pair overlap per IBZ kpoint
  in ibzBondProj, then distribute across star via
  atomPerm with both atom indices permuted
- [x] C31. Validate: run KNbO3 with IBZ and verify that
  symmetry-equivalent atoms produce identical Q* and
  bond order -- confirmed correct

---

## ARCHIVE

<!-- Resolved items go here. -->
