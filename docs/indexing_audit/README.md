# 1-Indexing Audit: Perl → Python Conversion

## Purpose

The Perl scripts in `src/scripts/` (notably `StructureControl.pm`) use a
1-indexed array convention throughout: slot `[0]` of every atom-, axis-,
element-, or species-keyed array is either unused or reserved as a
sentinel, and meaningful data lives in slots `[1..N]`. Loops walk
`foreach (1..N)`, and subroutines that accept atom/axis/element indices
expect 1-based values.

When these scripts are ported to Python, the native 0-indexed list
semantics create pressure to "fix" this to 0-indexed, but doing so
silently changes the interface convention at subroutine boundaries.

## The rule (user-confirmed)

**Preserve Perl conventions in both directions**, not "standardize on
1-indexing throughout Python."

- Where Perl uses 1-indexed arrays (the common case — atom/axis/
  element/species arrays, vector refs, coordinate tuples, etc.),
  Python must also be 1-indexed with a `[None]` or `[""]` sentinel in
  slot 0.
- Where Perl genuinely uses 0-indexed arrays (rare but real — e.g.
  `$absoluteSysQn[0..8]` where `n` is literally the index, or RGBA
  quadruples stored as `[r,g,b,a]` at `[0..3]`), Python must also be
  0-indexed without a sentinel.

The rule applies to all subroutines including throwaway 3-vector
helpers — Perl's `$product[0] = 0.0` sentinels are deliberate and
must be preserved.

## Severity levels

- **DRIFT**: Perl is 1-indexed; Python was converted to 0-indexed.
  End-to-end behavior may be correct because all callers were updated
  consistently, but the interface convention violates the rule.
- **REVERSE DRIFT**: Perl is 0-indexed; Python was forced to
  1-indexed. Violates the rule in the opposite direction.
- **MIXED**: A single method uses one convention for some of its
  parameters/returns and a different convention for others.
- **BUG**: Interface convention differs AND the caller/callee are
  mismatched, producing an off-by-one or wrong-slot access at runtime.
- **OK**: Python preserves the Perl convention exactly.

## Audit scope and order

Priority order (one file at a time):

1. `structure_control.py` ← *in progress*
2. `makeinput.py`
3. `pdb2skl.py`
4. `condense.py`
5. `bond_analysis.py`
6. `make_reactions.py`
7. `mod_struct.py`
8. Remaining converted scripts

For each file the audit proceeds in passes:

- **Pass 1**: Public-API methods that accept atom, axis, element, or
  species indices as parameters.
- **Pass 2**: Internal cross-method calls that pass such indices.
- **Pass 3**: Data-structure slot-0 conventions (are 1-indexed arrays
  initialized with a sentinel at index 0?).
- **Pass 4**: Loop bounds and range expressions.

Each finding is recorded as: *file:line*, Perl convention, Python
convention, severity, recommended fix.

## Per-section files

### `structure_control.py` — Pass 1

| # | Topic | Findings |
|---|-------|----------|
| [1.1](sc_1.1_coord_conv.md) | Coordinate conversion helpers | 4 DRIFT |
| [1.2](sc_1.2_atom_idx.md) | Atom-index accepting methods | 11 OK |
| [1.3](sc_1.3_vector_math.md) | Vector math helpers | 6 DRIFT |
| [1.4](sc_1.4_rotation.md) | Rotation helpers | 1 MIXED, 1 DRIFT, 1 refactor |
| [1.5](sc_1.5_item_pair.md) | Item/pair bounds helpers | 3 OK |
| [1.6](sc_1.6_multigroup.md) | Multi-group builders | 3 OK |
| [1.7](sc_1.7_remaining_pass1.md) | Remaining Pass-1 methods | 1 REVERSE DRIFT, 1 DRIFT, 7 OK |
| [1.8](sc_1.8_element_species.md) | Element/species-index methods | 7 OK |
| [1.9](sc_1.9_callsite_bridges.md) | Call-site bridge scaffolding | 6 sites |
| [1.10](sc_1.10_io_and_setters.md) | I/O methods + scan-point setters | 20 OK (17 spot-checked), 1 DRIFT |
| [1.11](sc_1.11_spacegroup_supercell_prep.md) | Space-group, supercell, prep-surface | 2 OK, 1 consumer of upstream drifts (~25 bridge sites) |
| [1.12](sc_1.12_small_helpers.md) | Small-helper utilities | 2 DRIFT, 3 signature change, 1 behavior fix, 1 OK |
| [1.13](sc_1.13_ring_coord_boo.md) | Ring / coordination / BOO family | 1 DRIFT (inner `[real,imag]` pair), 4 OK |
| [1.14](sc_1.14_bonding_list_family.md) | Bonding list / distance / pore / q-list | 8 OK (downstream consumers only) |
| 1.15 | Getters (`get*Ref` methods) | OK — API removed, direct attribute access |
| [Passes 2–4](sc_pass2_3_4.md) | Cross-method + slot-0 + loop-bounds sweeps | 0 new findings, 1 additional bridge site (1.9.e) |

## Resource-control (rc) file audit summary

Each priority-1–7 script has a companion `*rc.py` resource-
control file.  Checked individually for indexing-sensitive
defaults:

| rc file | Findings | Status |
|---|---|---|
| `makeinputrc.py` | `kp_mesh_scf`/`kp_mesh_pscf` default was `[1, 1, 1]` (0-indexed inner); Perl `@kpMesh[$grp][1..3]` is 1-indexed | **Fixed** during makeinput Cluster MI-B → `[None, 1, 1, 1]` |
| `pdb2sklrc.py` | Scalars only (file names, float buffer, bool) | OK |
| `condenserc.py` | Scalars only (str, bool, float) | OK |
| `bond_analysisrc.py` | Two flat element-name lists (`bridges`, `ions`) consumed as unordered sets by `sc.compute_qn`; box borders exposed as scalar keys (`box_min1..3`, `box_max1..3`) and reassembled into the 1-indexed `box_borders` at construction time | OK |
| `make_reactionsrc.py` | Pure scalars (int, str, float, bool) | OK |
| `mod_structrc.py` | `abc_order`/`xyz_order` default was `[1, 2, 3]` (0-indexed, mismatched with `sc.set_abc_xyz_assignment_order` which reads `[1..3]`); critical crash on every invocation | **Fixed** during mod_struct Cluster MS-D → `[None, 1, 2, 3]` |

Two rc files had list defaults that silently crossed a
`StructureControl` call boundary expecting 1-indexed input.
Both are now fixed.  The other four rc files either store
scalars only or assemble their list-valued settings into
1-indexed layouts at construction time.

---

### `makeinput.py`

- [Pass 1 + gap coverage](mi_pass1.md) — 1 BUG, 1 DRIFT, all other
  subsystems OK.  **Both clusters fixed** (see the "Fix status"
  section of the audit file).  Summary of the findings:
  - **BUG**: `group_block` reads slot 0 of 1-indexed
    `sc.direct_abc[atom]` → runtime `TypeError` when `-block` is
    used.  Latent since the initial port (no tests exercise
    `-block`).
  - **DRIFT**: `kp_mesh_scf` / `kp_mesh_pscf` inner axis is
    0-indexed; Perl `@kpMesh[$group][$axis]` was fully 1-indexed on
    both dimensions.  ~8 consumer sites plus the rc defaults.

### `pdb2skl.py`

- [Pass 1](pdb2skl_pass1.md) — zero findings.  The script is a
  thin CLI→`StructureControl` wrapper (~274 lines) with no
  per-atom, per-axis, per-element, or per-species array access of
  its own.  No fix phase required.  Two pre-existing non-indexing
  observations flagged for a separate follow-up: `-b`/`--buffer`
  is parsed but never applied, and there is minor docstring drift.

### `mod_struct.py`

- [Pass 1](mod_struct_pass1.md) — **4 BUGs, 1 DRIFT, all
  fixed**.
  The findings are:
  - **MS-A** (BUG): `hkl` passed to `sc.prep_surface` as a
    0-indexed 3-element list; prep_surface expects
    1-indexed `[None, h, k, l]`.
  - **MS-B** (BUG): `block_border` passed to
    `sc.cut_block` as a doubly-0-indexed 3×2 list;
    cut_block expects doubly-1-indexed.
  - **MS-C** (BUG): `sphere_loc` passed to
    `sc.cut_sphere` as a 0-indexed 3-element list;
    cut_sphere expects 1-indexed.
  - **MS-D** (BUG): `abc_order` / `xyz_order` passed to
    `sc.set_abc_xyz_assignment_order` as 0-indexed
    `[1, 2, 3]` from the rc defaults.  **Critical**: this
    call happens unconditionally at the start of `main()`,
    so every invocation of `mod_struct.py` crashes with
    `IndexError` *before any operation runs*.  Confirmed
    with a live smoke test.
  - **MS-E** (DRIFT): `op["sc"]` / `op["mirror"]` inner
    layout is 0-indexed; `do_supercell` bridges via
    `i + 1`.  Functional but diverges from Perl's
    1-indexed convention.
  No tests exercise `mod_struct.py`; the critical MS-D bug
  has been latent since the initial port.
  `make_reactions.py` invokes `mod_struct.py` as a
  subprocess in its condense pipeline, so this crash would
  break that pipeline too.

### `make_reactions.py`

- [Pass 1](make_reactions_pass1.md) — **zero findings**.  The
  ~3130-line file is noticeably more disciplined than the
  earlier ports: every nested-list data structure (per-atom,
  per-S-atom, per-bond, per-angle) uses the `[None]`-sentinel
  1-indexed layout on every dimension; every
  `StructureControl` call-boundary matches the post-Cluster-A
  interfaces; and `AngleData.hooke_angle_coeffs` is already
  built with the 1-indexed inner layout that `condense.py`
  needed Cluster CD-B to retrofit.  Four transient
  0-indexed scratch buffers (`trans_vector1`/`2`,
  `plane_coords`, `origin`, `angle_set`) were verified
  against the Perl original and are legitimate 0-indexed
  subprocess-argv / local scratch structures — matches Perl.
  No fix phase required.

### `bond_analysis.py`

- [Pass 1](bond_analysis_pass1.md) — 0 BUGs, 1 DRIFT, **fixed**.
  The single finding (Cluster BA-A) was the `bond_info_temp`
  inner 3-element layout in `_get_bonds_to_show` (Perl
  1-indexed `[1..3]`, Python 0-indexed).  Three touch points
  updated in one function; every other subsystem in the
  ~3429-line file preserves Perl 1-indexed conventions.  No
  test coverage exists for `bond_analysis.py`; post-fix
  verification by inspection + full `structure_control`
  suite (295 passed, 23 skipped — unchanged).

### `condense.py`

- [Pass 1](condense_pass1.md) — 0 BUGs, 9 DRIFT findings, all on
  the *inner axis* of nested-list bookkeeping structures.  Every
  drift is self-consistent between its producer and consumer,
  so the script functions correctly as-is; the fixes are pure
  rule compliance.  **Eight of nine clusters fixed**; CD-E
  (`ordered_species_pair_coeffs`) is intentionally deferred
  until the upstream `element_data.lj_pair_coeffs` tuple drift
  is repaired during the `element_data.py` audit.  Summary of
  the clusters (CD-A through CD-I):
  - CD-A/CD-B: `BondData.hooke_bond_coeffs` (4-elem) and
    `AngleData.hooke_angle_coeffs` (6-elem) inner entries are
    0-indexed; Perl `$hookeBondCoeffs[$hb][1..4]` /
    `$hookeAngleCoeffs[$ha][1..6]` are fully 1-indexed.
  - CD-C/CD-D: `ordered_bonded_atoms` (pair) and
    `angle_bonded_atoms` (triple) local bookkeeping arrays are
    0-indexed inner.
  - CD-E: `ordered_species_pair_coeffs` inherits the
    `element_data.lj_pair_coeffs` tuple drift (pending the
    separate `element_data.py` audit).
  - CD-F/CD-G: `unique_bond_coeffs` and `unique_angle_coeffs`
    in both `create_lammps_files` and `normalize_types`.
  - CD-H/CD-I: `unique_bond_types` (pair) and
    `unique_angle_types` (quadruple) in `normalize_types`.
- No test coverage for `condense.py` in `tests/`; regression
  verification after the fixes is by inspection plus the full
  `structure_control` suite that covers the shared call
  boundaries.

### Other files

- *pending (bond_analysis.py, make_reactions.py, mod_struct.py, remaining scripts)*

## Pass 1 gap verification (all unaudited methods spot-checked)

After Sections 1.1–1.14 were written, a subsequent pass re-read every
method in `structure_control.py` not explicitly covered by a section
file. The following methods were verified clean:

- **Lattice construction / inverse**: `abc_alpha_beta_gamma`,
  `get_angle_sine`, `get_abc_vectors`, `make_inv_or_recip_lattice` —
  all 1-indexed throughout.
- **Reorder helpers**: `reorder_lattice_parameters`,
  `reorder_atom_coordinates` — persistent arrays 1-indexed; the
  internal `indices` temporary is a 0-indexed local implementation
  detail that does not cross a public boundary.
- **Setup / compute**: `compute_implicit_info`,
  `compute_crystal_parameters` — 1-indexed throughout. `temp_array`
  from `element_data.get_vale_orbitals` is legitimately 0-indexed
  (external API, matches Perl).
- **Interaction family**: `finalize_interaction_data` — clean
  except for the already-cataloged `gaussian_broaden` bridge
  (Section 1.14.e).
- **Ring/coordination presumed-clean**: `compute_ring_distribution`,
  `create_coordination_list`, `create_coordination_summary`,
  `compute_min_network_dist_matrix` — confirmed clean. `curr_ring`
  inside ring traversal is a 0-indexed stack buffer matching Perl.
- **Bonding/extension presumed-clean**: `create_ext_bonding_list`,
  `map_ext_to_central`, `compute_atom_mesh_dist`,
  `compute_bond_mesh_dist` — confirmed clean orchestrators.
- **Add connection atoms**: `add_connection_atoms` — atom 1 path
  clean; atoms 2–4 raise NotImplementedError.
- **Initialization**: `_init_state` — every persistent array is
  allocated with the `[None]` or `[None, ...]` sentinel. No slot-0
  hazards.
- **Setters**: `set_num_atoms`, `set_limit_dist`, `set_buffer`,
  `set_ylm_l`, `set_xyz_mesh_points`, `set_scan_element`,
  `setup_data_base_ref`, `set_lattice_from_mag_angle` (1-indexed
  mags/angles), `set_abc_xyz_assignment_order` (1-indexed),
  `set_rpdf_select_atoms` (1-indexed pair storage) — all clean.

Two NEW findings surfaced during gap verification and are cataloged
below (not in a separate section file because they are standalone):

### NEW 1.16.a `set_border(zone, low, high, coord_system)` — DRIFT

**Perl** `setBorder` at `StructureControl.pm:547`: positional args
`$_[2..7]` for the six bound values; storage in `$borderLow[1..3]`
and `$borderHigh[1..3]` is 1-indexed.

**Python** `set_border(self, zone, low, high, coord_system='xyz')`
at `structure_control.py:799`: takes `low` and `high` as 0-indexed
length-3 lists and internally prepends `None` to build
`self.border_low = [None] + list(low)`. Storage is 1-indexed (OK)
but the parameter-convention is inconsistent with every other
Python method in the file that accepts coordinate-shaped arguments.

**External caller**: `bond_analysis.py:1154` currently passes
0-indexed lists; will need to update.

**Severity**: mild DRIFT. Fix: accept `low` and `high` as
1-indexed `[None, v1, v2, v3]` and store directly.

### NEW 1.16.b `get_unique_array(values)` — DRIFT + signature change (dead in this file)

**Perl** `getUniqueArray` at `StructureControl.pm:8403`: takes
`(array_ref, dataForm, doSort)` where `array_ref` is 1-indexed
(slot 0 unused). Returns two 1-indexed arrays: `@uniqueData` and
`@uniqueCount` (occurrence count per unique element).

**Python** `get_unique_array(self, values)` at
`structure_control.py:8726`: takes a 0-indexed input list, drops
`dataForm` and `doSort` flags, returns a single flat sorted list.
Count tracking is gone.

**Severity**: DRIFT + signature change. No in-file callers
(coordination methods use inline dicts instead); no external
callers across `makeinput.py`, `pdb2skl.py`, `condense.py`,
`bond_analysis.py`, `make_reactions.py`, `mod_struct.py`.

**Recommended fix**: delete the method entirely (no callers exist
across all converted Python scripts). The coordination methods
already replicate the Perl behavior inline.

## Running tally (structure_control.py Pass 1 complete — Sections 1.1–1.14 + gap pass)

- **Methods audited**: all methods in `structure_control.py`
- **DRIFT** (Perl 1-indexed → Python 0-indexed): 17
  - Section 1.1: 4 (coord conversion helpers)
  - Section 1.3: 6 (vector math helpers)
  - Section 1.4: 1 (rotate_one_point)
  - Section 1.7: 1 (spherical_angles)
  - Section 1.10: 1 (set_scan_points line-mode)
  - Section 1.12: 2 (stable_sort, gaussian_broaden)
  - Section 1.13: 1 (accumulate_q_bar inner `[real, imag]`,
    propagating to q_bar_per_bond/normalize/correlation)
  - Section 1.16.a: 1 (set_border parameter convention)
  - Section 1.16.b: 1 (get_unique_array — delete)
- **REVERSE DRIFT** (Perl 0-indexed → Python 1-indexed): 1
  - Section 1.7: 1 (`compute_qn` → `absolute_sys_qn`,
    `fractional_sys_qn`)
- **MIXED**: 1
  - Section 1.4: 1 (`define_rot_matrix`)
- **OK**: ~70 (Sections 1.2, 1.5, 1.6, 1.7 most, 1.8, 1.10 most,
  1.11 most, 1.12 most, 1.13 most, 1.14, 1.15, plus gap-pass set)
- **Signature/refactor only** (not indexing): 7
- **Call-site bridges** to be removed after helpers fixed: ~32
  (6 cataloged in Section 1.9 + ~25 inside `prep_surface` + 1 in
  `bond_analysis.py` for `set_border`)

## Fix-order recommendation

The drifts form dependency clusters. Fix in this order:

1. **Cluster A — 3-vector helpers** (Sections 1.1, 1.3, 1.4.b,
   1.7.i). Restore all helpers to 1-indexed input and output.
   Single unified sweep because they share call-sites. This alone
   eliminates ~31 of the bridge sites.
2. **Cluster B — `define_rot_matrix`** (Section 1.4.a). Fix input
   axis to 1-indexed, keep matrix output as-is (already correct).
3. **Cluster C — `compute_qn` reverse drift** (Section 1.7.c).
   Standalone — remove `+1` offsets from `absolute_sys_qn` /
   `fractional_sys_qn`.
4. **Cluster D — BOO inner pair** (Section 1.13.d-e). Edit all
   five related methods together
   (`initialize_interaction_data` type-4 branch,
   `accumulate_q_bar`, `q_bar_per_bond`, `q_bar_normalize`,
   `q_bar_correlation`).
5. **Cluster E — `set_scan_points` line-mode** (Section 1.10.e).
   Standalone; update the one or two callers.
6. **Cluster F — `gaussian_broaden`** (Section 1.12.g). Restore
   1-indexed layout, update the single caller in
   `finalize_interaction_data`.
7. **Cluster G — `stable_sort`** (Section 1.12.f). First confirm
   whether any external script calls it; if not, delete. If so,
   restore Perl's in-place 1-indexed signature.
8. **Cluster H — `set_border` parameter convention**
   (Section 1.16.a). Accept `low` / `high` as `[None, v1, v2, v3]`.
   Update the one external caller in `bond_analysis.py`.
9. **Cluster I — `get_unique_array`** (Section 1.16.b). Delete the
   method: no in-file or cross-script callers. The coordination
   family already inlines the Perl behavior via Python dicts.

After each cluster, run the full `tests/` suite to confirm no
regressions.

## Fix status (structure_control.py)

All nine clusters listed above are implemented as of this
revision.  Summary of what changed:

- **Cluster A** — `direct_xyz2fract_abc`, `direct_xyz2direct_abc`,
  `fract_abc2direct_xyz`, `get_direct_xyz_point`, `get_direct_xyz`,
  `dot_product`, `cross_product`, `cross_product_mag`,
  `normalized_cross_product`, `get_vector_angle`, `get_plane_normal`,
  `spherical_angles`, `rotate_one_point`, `rotate_all_atoms`, and
  the full `prep_surface` body rewritten for 1-indexed vectors.
- **Cluster B** — `define_rot_matrix` now takes a 1-indexed axis
  vector; `rotate_arb_axis` updated accordingly.
- **Cluster C** — `compute_qn` reverts to Perl's 0-indexed
  `absolute_sys_qn` / `fractional_sys_qn` (n is the live index).
- **Cluster D** — `q_bar[atom][k]` inner pair restored to
  `[None, real, imag]` across `initialize_interaction_data`
  (type 4), `accumulate_q_bar`, `q_bar_per_bond`,
  `q_bar_normalize`, and `q_bar_correlation`.
- **Cluster E** — `set_scan_points` line-mode now accepts
  1-indexed `start_xyz` / `end_xyz` vectors.
- **Cluster F** — `gaussian_broaden` now consumes and returns
  1-indexed data; `finalize_interaction_data` passes `self.rpdf`
  directly without the former slice-and-shift bridge.
- **Cluster G** — `stable_sort` **deleted** (no callers anywhere
  in the converted scripts; `sort_atoms` uses Python's built-in
  stable `sorted`).
- **Cluster H** — `set_border` accepts 1-indexed `low` / `high`
  vectors; the single external caller in `bond_analysis.py`
  updated to pass `[None, v1, v2, v3]`.
- **Cluster I** — `get_unique_array` **deleted** (no callers;
  `create_coordination_*` inline their own dict-based counting).

External callers updated: `makeinput.py` (target coord building),
`bond_analysis.py` (`set_border`, `get_direct_xyz_point` cell-rod
loop), `make_reactions.py` (`normalized_cross_product`,
`get_vector_angle`, diff vector plumbing), `mod_struct.py`
(`define_rot_matrix` rot axis, local `_get_plane_normal`).

Test outcomes: `tests/` full suite 295 passed / 23 skipped (down
from 302 passed / 23 skipped before fixes — the 7 deleted tests
belonged to `TestStableSort`, which was removed together with the
`stable_sort` method).
