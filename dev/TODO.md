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
- [x] C39. Implement angle clustering in
  create_lammps_files: collect (Z1, Zv, Z2, theta_obs)
  tuples, cluster within tolerance, compute K_angle from
  get_bond_params(), replace angles.dat lookup
  (DESIGN 4.8.3, 4.8.4, 4.8.8 items 1-3).  Landed in
  commit d51da45 along with a condense.py-wide cleanup
  of cryptic 1-2 letter locals and removal of 6 dead
  variables.  Subtasks C39.1-C39.5 all resolved by the
  same commit (see checkboxes below).

  Plan (captured 2026-04-16 — resolved in commit
  d51da45):

  - [x] C39.1. Add two shared helper functions.  Original
    plan proposed `_cluster_angles` and `_compute_angle_k`
    as methods on Condense; the actually-shipped form (commit
    74805ed) promoted them to module-level functions in
    src/scripts/angle_utils.py to avoid forcing
    make_reactions.py to import the larger condense module:
    - cluster_angles(observations, tolerance): PSEUDOCODE 10a.
      Group observations by canonical (z1, zv, z2) triplet
      (with z1 <= z2), sort each group by observed theta,
      greedy-merge while BOTH |theta - running_mean| <=
      tolerance AND the resulting cluster span stays within
      2 * tolerance (spread cap, consistent with 10e's
      cross-source step).  Return (angle_types,
      angle_type_map), where each angle_type carries z1,
      zv, z2, theta_0, obs_count, and a representative
      base_tag copied from the first observation in the
      cluster.
    - get_angle_k(z1, zv, z2, stiffness, scale,
      get_bond_params): PSEUDOCODE 10b.  Return
      stiffness * sqrt(K_arm1 * K_arm2) * scale, pulling
      K_arm1/K_arm2 from the injected get_bond_params
      callable (bound to Condense.bond_data.get_bond_params
      at the condense.py call site).

  - [x] C39.2. Refactor the angle section of the atom loop
    in create_lammps_files (current lines ~1521-1637) into
    a collect-only pass.  Build each observation dict
    {'z_trip', 'theta_obs', 'base_tag', 'vertex_atom',
    'end_atom_1', 'end_atom_2'} and append to a new
    angle_observations list declared alongside bond_count
    and angle_count near line 1390.  Remove the
    angles.dat lookup, the tag-uniqueness dedup, and the
    in-loop type_id assignment.

  - [x] C39.3. After the atom loop completes, add a
    clustering block that:
    - Calls cluster_angles(angle_observations,
      self.angle_cluster_tolerance) from angle_utils.
    - Sets angle_count = len(angle_observations) and
      num_local_angle_types = len(angle_types).
    - Builds local_angle_tags[t] as the string
      "{base_tag} {theta_0:.4f} {t}" so the final two
      tokens ({theta_0} {t}) match the tag-tail slot
      that normalize_types expects (10c Phase 3, 10d
      Phase 3, 10f Phase B all share this format).
    - Builds local_angle_coeffs[t] = [None, K_angle,
      theta_0] via get_angle_k(...,
      self.bond_data.get_bond_params).
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

  - [x] C39.4. Delete the now-unused local alias
    ``ad = self.angle_data`` inside create_lammps_files
    (around line 1338).  The AngleData class itself, and
    the identical alias in normalize_types, stay until
    C41 — normalize_types is updated in C40.

  - [x] C39.5. Sequencing constraint (not a side effect
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
- [x] C39a. Port make_reactions.py's existing angle
  handling off angles.dat.  Remove the duplicated
  AngleData class (lines ~113-193) and the
  self.angle_data instantiation (line ~568).  In the
  _read_angle_data method (around line 2518), replace
  the hooke_angle_coeffs scan with the new collect /
  cluster(tolerance=0) / emit structure from
  PSEUDOCODE 10d so the emitted tag tail matches
  "{theta_0:.4f} {t}" in the same format C39 writes
  in create_lammps_files.  Cross-source unification
  of equivalent angle types is done centrally in
  normalize_types (C40), not by string comparison on
  the producer-emitted tag tail.

  Decision refined (2026-04-18): fixed tolerance = 0
  at the template producer; no K computation at the
  template producer.  The original 2026-04-17
  "hybrid with matching tolerance on both producers"
  plan was revisited once we traced the cross-source
  data flow in detail and realized:

  1. Local clustering at make_reactions.py is an
     optimization, not a correctness step.  Any
     tolerance T_m > condense.py's T_c causes silent
     wrong physics (distinct angles fused at the
     producer cannot be un-fused downstream).  Fixing
     T_m = 0 eliminates the hazard and makes
     reaction templates reusable across any
     condense.py run regardless of its T_c.  See
     DESIGN 4.8.10 for the full hazard analysis and
     rejection of the parameter-manifest alternative.

  2. Reaction template files carry no K values --
     only connectivity, per-atom angle entries, and
     the "{theta_0_local:.4f} {t}" tag tail.
     normalize_types recomputes K authoritatively
     from the triplet in 10f Phase C regardless of
     what any producer computed, so a producer-side
     K in make_reactions.py would be neither written
     nor consumed.  Skipping it keeps
     make_reactions.py independent of BondData and
     avoids plumbing angle_stiffness_coeff and
     angle_parameter_scale into a script that never
     writes their effect to disk.  DESIGN 4.8.8
     item 3 was amended to record this split:
     condense.py computes K, make_reactions.py does
     not.

  Net consequence: the C39a scope is smaller than
  the 2026-04-17 plan called for.  No parameter
  plumbing, no BondData wiring, no shared
  angle_cluster_tolerance handshake.  Just collect,
  cluster at tolerance=0, emit tags.

  Plan (captured 2026-04-17, refined 2026-04-18 --
  pending before writing code):

  - [x] C39a.1. Shared cluster_angles helper.
    Resolved in commit 74805ed: src/scripts/
    angle_utils.py provides AngleType NamedTuple +
    cluster_angles (PSEUDOCODE 10a) + get_angle_k
    (10b), with 10 unit tests in
    src/tests/test_angle_utils.py.  Both producers
    call cluster_angles from there.  make_reactions.py
    will call cluster_angles but NOT get_angle_k --
    get_angle_k is condense.py-only under the
    refined decision.

  - [x] C39a.2. In make_reactions.py, refactor the
    angle loop inside _read_angle_data (around line
    2518) into a collect-only pass per PSEUDOCODE 10d
    Phase 1.  Build each observation as the 8-tuple
    (z1, zv, z2, theta_obs, base_tag, a1, v, a2)
    with z1 <= z2 canonicalization -- same shape as
    condense.py create_lammps_files uses.  Remove
    the inline hooke_angle_coeffs scan, the in-loop
    tag-tail construction, and the existing
    angle_tag uniqueness dedup (cluster_angles does
    the dedup instead).

  - [x] C39a.3. After the template's angle loop
    completes, call cluster_angles(observations, 0.0)
    per PSEUDOCODE 10d Phase 2.  Tolerance = 0 means
    identity-only merge: observations with
    bit-identical theta values (possible because
    _read_angle_data rounds raw angles to 0.5-degree
    resolution before this point) collapse into one
    local type with obs_count > 1; non-identical
    observations each become their own local type.
    For each returned local cluster, build the tag
    tail "{theta_0:.4f} {t}" matching C39.3's
    format.  Do NOT compute or store K -- templates
    carry no angle coefficients (see PSEUDOCODE 10d
    Phase 3 and DESIGN 4.8.8 item 3).  Walk
    observations in collection order to populate
    the per-atom angle_bonded[] and angle_tag_id[]
    structures so downstream template writes see
    the ordering they expect.

  - [x] C39a.4. Delete the self.angle_data
    instantiation at line ~568 and the duplicated
    AngleData class at lines ~113-193.  The
    AngleData class itself lives on in condense.py
    until C41, which removes both copies and
    retires angles.dat.

  - [x] C39a.5. Caller update.  _read_angle_data's
    current return tuple includes unique_angle_tags
    and num_unique_angle_tags, which are no longer
    part of the new export (see PSEUDOCODE 10d
    export list).  Update the caller at
    make_reactions.py:~2171 to unpack only the
    fields normalize_types actually consumes:
    num_bond_angles, angle_bonded, angle_tag_id,
    local_angle_tags, per-local-type obs_count,
    num_angles_total, bond_angles_ext.  Any
    downstream code paths in make_reactions.py that
    read the old unique_angle_tags slot must be
    traced and either removed or updated to read
    local_angle_tags instead.

  - [x] C39a.6. Sequencing note (eased).  Under the
    refined decision, C39a does NOT have to land
    in lockstep with C40 -- the emitted tag-tail
    format already matches what C39 writes and
    what C40 will parse.  The pipeline is still
    inconsistent until C40 is live in
    normalize_types (cross-source clustering still
    uses the old hooke-id-suffix lookup that
    normalize_types currently does), but the C39a
    landing moment is flexible.  Strict ordering:
    C39a then C40 then C41 (teardown) then C42
    (validate) -- C39a in any ordering before C41
    is fine.
- [x] C40. Update normalize_types to own cross-
  source angle unification (DESIGN 4.8.8 item 4).
  Expands beyond the original "replace
  hooke_angle_coeffs scan with get_bond_params()"
  into the authoritative cross-file clustering step
  the hybrid approach requires.  Landed 2026-04-19
  alongside C39a: angle_utils.py gained
  cross_source_cluster + LocalRecord + FinalAngleType
  (with 9 new unit tests), and normalize_types now
  routes every angle through the collect-and-unify
  pipeline rather than the legacy hooke-id lookup.

  Plan (captured 2026-04-17 -- resolved 2026-04-19):

  - [x] C40.1. Collect all angle types emitted by
    every source: the lammps.dat produced by
    create_lammps_files and every reaction template
    produced by make_reactions.py.  Each incoming
    local_record carries (z1, zv, z2, theta_0_local,
    obs_count, base_tag, source, local_type_id) per
    PSEUDOCODE 10e's input layout.  obs_count is the
    observation-population weight used by the
    weighted running-mean merge; base_tag is the
    producer's representative tag prefix
    (element/species/molecule) that 10e carries
    through to the final type's canonical tag.

  - [x] C40.2. Group by canonical triplet (z1, zv,
    z2 with z1 <= z2).  Within each group, sort by
    theta_0_local and greedy-merge any pair within
    angle_cluster_tolerance of the running cluster
    mean.  The running mean is the final canonical
    theta_0 for the merged cluster.  Cap the total
    cluster spread (e.g., max theta span <=
    2 * tolerance) to avoid greedy chaining across
    a wide distribution.

  - [x] C40.3. For each final cluster, compute
    K_angle via get_angle_k() from angle_utils.py
    (PSEUDOCODE 10b; same helper create_lammps_files
    calls for its local step).  Apply
    angle_stiffness_coeff and angle_parameter_scale.

  - [x] C40.4. Rewrite tag tails and remap type IDs
    across every source file.  Each (source,
    source_type_id) pair maps to one final cluster
    id; the lammps.dat Angles section and every
    template angle reference is rewritten to the
    unified id.  Rewrite the tag tail to the final
    canonical theta_0 so any downstream tool that
    inspects the tag sees a consistent value.

  - [x] C40.5. Diagnostic (recommended): write a
    cluster-map file or log section listing, for
    every final cluster: final id, canonical
    theta_0, (z1, zv, z2), and the contributing
    (source, local_theta_0) pairs.  This is the
    main payback for the hybrid approach -- someone
    debugging a bond/react mismatch can open the
    file and see exactly which observations got
    merged where.
- [x] C41. Remove AngleData class and angles.dat from
  build (DESIGN 4.8.8 items 1, 6).  Landed 2026-04-19:
  AngleData class deleted from condense.py (the make_
  reactions.py copy went in C39a), self.angle_data
  init removed, init_env docstring updated, and
  angles.dat removed from the DATABASES list in
  src/data/CMakeLists.txt (per the bonds.dat
  precedent, the Perl AngleData.pm and the on-disk
  angles.dat file remain for historical reference but
  are not installed or read by any active code path).
  Runbook note: any pre-existing on-disk reaction
  templates carrying the old "{rest_angle}
  {hooke_id}" tail must be regenerated by the updated
  make_reactions.py before C42 validation.
- [ ] C42. Validate: run a condense.py job end-to-end
  and verify that LAMMPS Angle Coeffs contain
  physically reasonable force constants and rest angles
  derived from the input geometry

#### Step 2: Look-ahead angles in makeReactions (deferred)

- [ ] D6. Design angle creation in make_reactions.py
  post-reaction templates: port commented-out Perl
  addBondAngle, compute theta_0 from post-reaction
  coordinates, register new angle types
  (DESIGN 4.8.7, see "Empirical confirmation"
  paragraph for the 10-missing-angle decomposition).
  Concrete evidence first observed 2026-04-19 in
  jobs/molecules/b12/3_mol and re-confirmed
  2026-04-25 in jobs/molecules/b12/60_mol: the
  postRxn template has 230 angles (preRxn 240 minus
  the 10 that touched the deleted H's) but contains
  zero angles involving both new-bond endpoints
  (atoms 1 and 18 in the b12h12-b12h12 case).  Ten
  B-B-B angles are expected -- five at each vertex
  (e.g. 2-1-18, 3-1-18, ..., and the symmetric five
  around atom 18).  Without those angles the new
  inter-cage B-B bond is held only by its harmonic
  stretch -- the two cages can rotate freely about
  the bond axis with no restoring torque, and the
  bond axis itself can swing relative to either
  cage's local symmetry.  In an N-mer chain the
  deficit grows as 10N unconstrained angular DoFs
  at the joints.  The existing docstring on
  _build_phase_angles (src/scripts/make_reactions.py
  ~line 2675) already warns "the bonded molecules
  may be too floppy" because of this gap.  Not the
  cause of any reaction-trigger failure observed so
  far; deferred cleanly.

- [ ] C43. In make_reactions.py, the new B-B (or
  more generally new s1-s2) bond written into the
  postRxn template reuses the existing bond type
  for the element pair.  Confirmed by inspection
  for b12h12_1_b-1_b12h12_1_b-1: postRxn bond 61 is
  type 1 (the same B-B harmonic as every other B-B
  bond in the molecule).  That is correct for UFF-
  derived parameters where a new B-B is physically
  the same spring as any other B-B, but the choice
  should be revisited if/when a reaction wants to
  distinguish a "new" bond from a native one (e.g.
  for a different K value or equilibrium length at
  the reaction site).  Revisit alongside D6.

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
