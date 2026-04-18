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
- [x] P3. Write pseudocode for angle clustering and force
  constant computation: greedy cluster-by-triplet, K from
  geometric mean of arm bond stiffnesses, integration
  into create_lammps_files and normalize_types
  (PSEUDOCODE 10, DESIGN 4.8.3-4.8.8)

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

### Phase G -- UFF bond parameter database (DESIGN 4)

- [x] C32. Create bond_parameters.dat with UFF per-element
  parameters for Z=1-54 (DESIGN 4.4)
- [x] C33. Rewrite BondData class in condense.py to read
  bond_parameters.dat and compute K_ij, r_ij via
  get_bond_params() (DESIGN 4.6 items 1-3)
- [x] C34. Add bond_parameter_scale to condenserc.py,
  ScriptSettings, Condense.__init__, and
  parse_input_file (DESIGN 4.5, 4.6 item 5)
- [x] C35. Replace linear bond scans in create_lammps_files
  and normalize_types with get_bond_params() calls,
  applying bond_parameter_scale (DESIGN 4.6 items 3-5)
- [x] C36. Update CMakeLists.txt to install
  bond_parameters.dat instead of bonds.dat
  (DESIGN 4.7)
- [ ] C37. Validate: run a condense.py job end-to-end and
  verify that LAMMPS Bond Coeffs contain realistic
  UFF force constants and bond lengths

### Phase H -- Geometry-derived angle parameters (DESIGN 4.8)

#### Step 1: Replace angles.dat with clustering + computed K

- [x] D5. Design geometry-derived angle parameter system:
  cluster observed angles by triplet, compute K from UFF
  bond stiffnesses, add angle_parameter_scale and
  angle_cluster_tolerance keywords (DESIGN 4.8.1-4.8.6)
- [x] C38. Add angle_stiffness_coeff (default 0.15),
  angle_parameter_scale (default 1.0), and
  angle_cluster_tolerance (default 5.0) to
  Condense.__init__ and parse_input_file (DESIGN
  4.8.4-4.8.6, 4.8.8 item 5).  Follows the
  bond_parameter_scale precedent: force-field params
  live in Condense.__init__ with condense.in override,
  not in condenserc.py or ScriptSettings.  Landed in a
  prior session; verified present at lines 767-780
  (defaults) and 917-947 (parse_input_file).
- [ ] C39. Implement angle clustering in
  create_lammps_files: collect (Z1, Zv, Z2, theta_obs)
  tuples, cluster within tolerance, compute K_angle from
  get_bond_params(), replace angles.dat lookup
  (DESIGN 4.8.3, 4.8.4, 4.8.8 items 1-3)

  Plan (captured 2026-04-16 — pending before writing code):

  - [ ] C39.1. Add two helper methods on Condense, inserted
    just before create_lammps_files:
    - _cluster_angles(observations): PSEUDOCODE 10a.  Group
      observations by canonical (z1, zv, z2) triplet (with
      z1 <= z2), sort each group by observed theta, greedy-
      merge while BOTH |theta - running_mean| <=
      self.angle_cluster_tolerance AND the resulting cluster
      span stays within 2 * self.angle_cluster_tolerance
      (spread cap, consistent with 10e's cross-source step).
      Return (angle_types, angle_type_map), where each
      angle_type carries z1, zv, z2, theta_0, and a
      representative base_tag copied from the first
      observation in the cluster.
    - _compute_angle_k(z1, zv, z2): PSEUDOCODE 10b.  Return
      self.angle_stiffness_coeff
      * sqrt(K_arm1 * K_arm2)
      * self.angle_parameter_scale,
      pulling K_arm1/K_arm2 from
      self.bond_data.get_bond_params().

  - [ ] C39.2. Refactor the angle section of the atom loop
    in create_lammps_files (current lines ~1521-1637) into
    a collect-only pass.  Build each observation dict
    {'z_trip', 'theta_obs', 'base_tag', 'vertex_atom',
    'end_atom_1', 'end_atom_2'} and append to a new
    angle_observations list declared alongside bond_count
    and angle_count near line 1390.  Remove the
    angles.dat lookup, the tag-uniqueness dedup, and the
    in-loop type_id assignment.

  - [ ] C39.3. After the atom loop completes, add a
    clustering block that:
    - Calls self._cluster_angles(angle_observations).
    - Sets angle_count = len(angle_observations) and
      num_local_angle_types = len(angle_types).
    - Builds local_angle_tags[t] as the string
      "{base_tag} {theta_0:.4f} {t}" so the final two
      tokens ({theta_0} {t}) match the tag-tail slot
      that normalize_types expects (10c Phase 3, 10d
      Phase 3, 10f Phase B all share this format).
    - Builds local_angle_coeffs[t] = [None, K_angle,
      theta_0] via _compute_angle_k.
    - Walks angle_observations in collection order and
      populates angle_bonded_atoms and
      ordered_angle_type with local type ids (which
      normalize_types remaps to global ids in 10f
      Phase A) so the per-atom ordering the downstream
      LAMMPS writer depends on is preserved exactly.
    - Exports per-local-type obs_count (slot 5 of
      each angle_types entry from 10a) so 10e can
      weight cross-source clustering by observation
      population.

  - [ ] C39.4. Delete the now-unused local alias
    ``ad = self.angle_data`` inside create_lammps_files
    (around line 1338).  The AngleData class itself, and
    the identical alias in normalize_types, stay until
    C41 — normalize_types is updated in C40.

  - [ ] C39.5. Sequencing constraint (not a side effect
    to live with): C39 and C39a must land together, and
    the precursor reaction-template DB must be rebuilt
    before any end-to-end test.  Until both producers
    (create_lammps_files in condense.py and the angle
    tag builder in make_reactions.py) emit the new
    "{theta_0} {t}" tag tail, normalize_types will see
    the same physical angle as two distinct types (at4
    string mismatch across lammps.dat vs templates) and
    bond/react type IDs will disagree, even with C40
    already applied.  The fix is to migrate both
    producers, not to paper over it with a hooke-id
    shim in create_lammps_files.
- [ ] C39a. Port make_reactions.py's existing angle
  handling off angles.dat, in lockstep with C39.
  Remove the duplicated AngleData class (lines ~113-
  193) and the self.angle_data instantiation (line
  ~568).  In the angle tag construction (around line
  2518), replace the hooke_angle_coeffs scan with the
  new approach so the emitted tag tail matches
  "{theta_0} {t}" exactly as C39 writes it in
  create_lammps_files.  Cross-source unification of
  equivalent angle types is done centrally in
  normalize_types (C40), not by string comparison on
  the producer-emitted tag tail.

  Decision resolved (2026-04-17): hybrid approach.
  Each producer (make_reactions.py and
  create_lammps_files) locally clusters its own
  observations within angle_cluster_tolerance and
  emits one type per local cluster with the cluster-
  mean theta_0 in the tag.  This keeps template files
  compact and human-readable -- one type per
  chemically distinct angle environment, not one per
  floating-point observation.  normalize_types (C40)
  then re-clusters across all sources using the same
  tolerance, picks a final canonical theta_0, and
  rewrites tag tails plus remaps type IDs uniformly.
  The cross-source step compares theta_0 as a float
  within tolerance, not as a string, so small local-
  cluster drift between producers is absorbed
  instead of splitting a physical angle in two.

  Plan (captured 2026-04-17 -- pending before
  writing code):

  - [ ] C39a.1. Promote the _cluster_angles helper
    from C39.1 to a shared call site both producers
    use.  Options: (a) keep it on Condense and have
    make_reactions.py import and call it, (b) move
    it to a standalone module-level utility.  Either
    way, both producers must invoke the identical
    routine so local clustering semantics are byte-
    identical.  No new algorithm, just a shared
    entry point.

  - [ ] C39a.2. In make_reactions.py, refactor the
    angle construction loop around line 2518 into a
    collect-only pass.  Build observation dicts
    {'z_trip', 'theta_obs', 'base_tag',
    'vertex_atom', 'end_atom_1', 'end_atom_2'} per
    triplet encountered during template
    construction.  Remove the inline
    hooke_angle_coeffs scan and the in-loop tag-
    tail construction.

  - [ ] C39a.3. After the template's angle loop
    completes, call the shared _cluster_angles on
    the collected observations.  For each returned
    local cluster, build the tag tail
    "{theta_0_local_mean} {t}" where t is the type
    id within the template's local numbering.  Walk
    the observations in collection order and
    populate the per-atom angle_bonded[] and
    angle_tag_id[] structures so downstream template
    writes see the ordering they expect.

  - [ ] C39a.4. Delete the self.angle_data
    instantiation at line ~568 and the duplicated
    AngleData class at lines ~113-193.  The
    AngleData class itself lives on in condense.py
    until C41, which removes both copies and
    retires angles.dat.

  - [ ] C39a.5. Sequencing: C39 and C39a land
    together (same commit or back-to-back).  Until
    both producers emit the new "{theta_0} {t}" tag
    tail AND C40's cross-source re-clustering is
    live in normalize_types, the pipeline is
    inconsistent.  Ordering is C39 + C39a, then
    C40, then C41 (teardown), then C42 (validate).
- [ ] C40. Update normalize_types to own cross-
  source angle unification (DESIGN 4.8.8 item 4).
  Expands beyond the original "replace
  hooke_angle_coeffs scan with get_bond_params()"
  into the authoritative cross-file clustering step
  the hybrid approach requires.

  Plan (captured 2026-04-17 -- pending before
  writing code):

  - [ ] C40.1. Collect all angle types emitted by
    every source: the lammps.dat produced by
    create_lammps_files and every reaction template
    produced by make_reactions.py.  Each incoming
    type carries (z1, zv, z2, theta_0_local, source,
    source_type_id).

  - [ ] C40.2. Group by canonical triplet (z1, zv,
    z2 with z1 <= z2).  Within each group, sort by
    theta_0_local and greedy-merge any pair within
    angle_cluster_tolerance of the running cluster
    mean.  The running mean is the final canonical
    theta_0 for the merged cluster.  Cap the total
    cluster spread (e.g., max theta span <=
    2 * tolerance) to avoid greedy chaining across
    a wide distribution.

  - [ ] C40.3. For each final cluster, compute
    K_angle using get_bond_params() for the two arm
    bonds (same formula as _compute_angle_k from
    C39.1).  Apply angle_stiffness_coeff and
    angle_parameter_scale.

  - [ ] C40.4. Rewrite tag tails and remap type IDs
    across every source file.  Each (source,
    source_type_id) pair maps to one final cluster
    id; the lammps.dat Angles section and every
    template angle reference is rewritten to the
    unified id.  Rewrite the tag tail to the final
    canonical theta_0 so any downstream tool that
    inspects the tag sees a consistent value.

  - [ ] C40.5. Diagnostic (recommended): write a
    cluster-map file or log section listing, for
    every final cluster: final id, canonical
    theta_0, (z1, zv, z2), and the contributing
    (source, local_theta_0) pairs.  This is the
    main payback for the hybrid approach -- someone
    debugging a bond/react mismatch can open the
    file and see exactly which observations got
    merged where.
- [ ] C41. Remove AngleData class and angles.dat from
  build (DESIGN 4.8.8 items 1, 6).  Also rebuild the
  precursor reaction-template DB after this step: any
  on-disk templates carrying the old "{rest_angle}
  {hooke_id}" tail are no longer valid and must be
  regenerated by the updated make_reactions.py (C39a)
- [ ] C42. Validate: run a condense.py job end-to-end
  and verify that LAMMPS Angle Coeffs contain
  physically reasonable force constants and rest angles
  derived from the input geometry

#### Step 2: Look-ahead angles in makeReactions (deferred)

- [ ] D6. Design angle creation in make_reactions.py
  post-reaction templates: port commented-out Perl
  addBondAngle, compute theta_0 from post-reaction
  coordinates, register new angle types
  (DESIGN 4.8.7)

---

## TOOLING (lint helpers)

- [ ] T1. Improve `.claude/commands/scripts/rewrap_prose.py`
  heuristics so it stops over-aggressively breaking
  paragraphs.  Three known false positives observed
  while landing C39.1/C39a.1 on 2026-04-18:
    1. Lines ending in `:` always terminate a paragraph
       (lines 657-658, 683-684), which orphans the prose
       that flows into a colon-introduced equation or
       list ("...via the formula:" -> previous line
       ends mid-sentence with no following words).
    2. The `_LABEL_LINE_RE` heuristic
       (`^\s*\w[\w_\. \t]{2,25}:(\s|$)`) catches
       sentence-final phrases like "Default upstream:
       0.15", treating them as field labels and
       refusing to merge them with the surrounding
       prose.
    3. Hyphenated compounds split across lines
       (e.g. "re-clustering", "geometric-mean") are
       sometimes rejoined with an extra space ("re-
       clustering" -> "re- clustering") instead of
       being un-split.
  Each issue has a workaround at the source level
  (rephrase to avoid the trigger pattern), but the
  heuristics should ideally handle these cases so the
  source can read naturally.  See angle_utils.py and
  test_angle_utils.py for the prose patterns that
  initially tripped the script.

---

## ARCHIVE

<!-- Resolved items go here. -->
