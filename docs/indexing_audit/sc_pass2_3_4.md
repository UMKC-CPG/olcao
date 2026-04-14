# structure_control.py — Passes 2, 3, 4

Verification sweeps across the file after Pass 1 section-by-section
audit. Each pass is a systematic scan targeting a different kind of
indexing hazard; findings are correlated back to the Pass 1 sections
where they originated.

## Pass 2 — Cross-method call-site sweep

**Method**: grep every `self.<method>(` call across the file,
verify each argument matches the callee's expected indexing
convention. Focus on methods that accept atom/axis/element
indices, 1-indexed sub-lists, or list references.

**Method families swept**:

1. Coordinate conversion helpers (`get_direct_xyz`,
   `get_direct_abc`, `get_fract_abc`, `get_direct_xyz_point`,
   `direct_xyz2fract_abc`, `direct_xyz2direct_abc`,
   `fract_abc2direct_xyz`) — ~50 call sites
2. Vector math helpers (`dot_product`, `cross_product`,
   `cross_product_mag`, `normalized_cross_product`,
   `get_vector_angle`, `get_plane_normal`) — ~20 call sites
3. Rotation helpers (`define_rot_matrix`, `rotate_one_point`,
   `rotate_all_atoms`, `rotate_arb_axis`) — ~10 call sites
4. Interaction data family (`compute_num_cells`,
   `compute_extended_pos`, `obtain_atomic_interaction`,
   `initialize_interaction_data`, `record_interaction_data`,
   `finalize_interaction_data`, `map_ext_to_central`) — ~30 call
   sites
5. Bonding/analysis families (`create_bonding_list`,
   `create_ext_bonding_list`, `create_q_list`, `create_min_dist_matrix`,
   `compute_rpdf`, etc.) — ~20 call sites
6. Element/species (`create_element_list`, `create_species_data`,
   `count_element_atoms`, `map_element_number`,
   `get_element_number`, `assign_coval_radii`, `assign_color_vtk`)
   — ~20 call sites
7. BOO (`accumulate_q_bar`, `q_bar_per_bond`, `q_bar_normalize`,
   `q_bar_correlation`, `init_ylm_coeff`) — ~5 call sites
8. Structure manipulation (`apply_sort`, `apply_filter`,
   `apply_perturbation`, `translate_atoms`, `shift_xyz_center`,
   `check_bounding_box`, `cut_block`, `cut_sphere`, `insert_vacuum`,
   `make_ortho`, `apply_space_group`, `apply_supercell`) — ~20 call
   sites
9. I/O (`read_*`, `print_*`, `read_input_file`) — ~15 call sites
10. Ring/coordination (`compute_ring_distribution`,
    `compute_ring_level_delta`, `save_ring`, `compute_bond_angles_ext`,
    `compute_min_network_dist_matrix`, `create_coordination_list`,
    `create_coordination_summary`) — ~10 call sites
11. Utility helpers (`compute_qn`, `stable_sort`, `gaussian_broaden`,
    `item_out_of_bounds`, `item_not_requested`, `pair_not_of_element`)
    — ~15 call sites
12. Setter/configuration (`set_*`, `reset`, `compute_implicit_info`,
    `compute_crystal_parameters`, `add_connection_atoms`, `get_min_max_xyz`)
    — ~15 call sites

**Total**: ~230 call sites traced.

### Pass 2 findings

**No new DRIFT findings.** Pass 2 confirms that the call-site
bridges catalogued in Section 1.9 constitute the complete set.
The 1-indexed persistent-array conventions are passed
consistently between methods; every deviation from that pattern
is already classified as either:

- a bridge to a drifted helper (Sections 1.9, 1.11.c), or
- a legitimate 0-indexed-buffer translation (e.g.
  `data_lines[dat_atom - 1]` in `read_olcao_species` /
  `read_dat_skl_map`, already discussed in Section 1.10.d), or
- a legitimate spiral-search / mesh-step counter (inside
  `prep_surface` line 4659–4664 and `set_scan_points` mesh-mode
  branch, both matching Perl `0..max-1` counters).

### Pass 2 addition to Section 1.9 — 1.9.e `rotate_arb_axis` → `define_rot_matrix`

One additional bridge site not previously catalogued:

**Python** at `structure_control.py:5794`:

```python
mag = math.sqrt(sum(v*v for v in axis))
unit = [v / mag for v in axis]
self.define_rot_matrix(unit, angle_deg)
```

`unit` is built as a 0-indexed length-3 list (from the 0-indexed
`axis` parameter) and passed to `define_rot_matrix`, which has the
0-indexed-input MIXED problem in Section 1.4.a. This is a 5th
bridge site for Section 1.9 (not 1.9.a-d).

Once Section 1.4.a is fixed to accept a 1-indexed `axis` input,
this becomes:

```python
# axis is already 1-indexed [None, x, y, z]
mag = math.sqrt(sum(v*v for v in axis[1:]))
unit = [None] + [v / mag for v in axis[1:]]
self.define_rot_matrix(unit, angle_deg)
```

The `rotate_arb_axis` signature itself (Section 1.4.c) still needs
revisiting — it currently accepts a 0-indexed `axis` parameter as
documented.

## Pass 3 — Slot-0 conventions sweep

**Method**: grep every list initialization pattern (`[None]`,
`[None] + [...]`, `[None] * N`, `[None, ...]`) and verify the
slot-0 sentinel is appropriate for the array's role. Also grep
for arrays that *should* be 0-indexed per Perl but got a
sentinel added (reverse drift candidates).

### Pass 3 statistics

- **183 `[None]`-sentinel initialization patterns** across
  `structure_control.py`. This is a very strong signal that the
  1-indexed-with-sentinel convention is pervasively applied —
  almost every per-atom, per-element, per-species, per-axis, or
  per-cell array uses the `[None] + [...]` or `[None, ...]`
  pattern.
- Per-atom arrays (`fract_abc`, `direct_xyz`, `direct_abc`,
  `atom_element_name`, `atom_element_id`, `atom_species_id`,
  `atom_tag`, `connection_tag`, `atomic_z`, `coval_radii`,
  `color_vtk`, `molecule_name`, `residue_name`, etc.) are
  uniformly initialized as `[None]` empty lists or
  `[None] * (num_atoms + 1)` fixed-size lists. ✓
- Per-axis arrays (`mag`, `angle`, `angle_deg`, `full_cell_mag`,
  `full_cell_angle`, `full_cell_angle_deg`, `supercell`,
  `sc_mirror`) — initialized with slot-0 sentinel. ✓
- Per-element and per-species arrays (`element_list`,
  `element_count`, `species_list`, `num_species`) — 1-indexed
  outer with `[None]` sentinel; species_list inner also uses
  `[None]`. ✓
- Per-bond arrays (`num_bonds`, `bonded_ext`, `bond_length_ext`,
  `bonded`, `num_bond_angles`, `bond_angles_ext`) — 1-indexed
  outer with `[None]` sentinel; inner lists pre-initialized with
  `[None]` before append. ✓

### Pass 3 findings

**`compute_qn` REVERSE DRIFT (already catalogued in 1.7.c)** — two
specific sites:

- Line 517: `self.absolute_sys_qn = [None]` in `_init_state()`
  — initializes with sentinel.
- Line 518: `self.fractional_sys_qn = [None]` — same.
- Line 7998: `self.absolute_sys_qn = [None] + [0] * 9` inside
  `compute_qn`.
- Line 8011: `self.fractional_sys_qn = [None] + [0.0] * 9` inside
  `compute_qn`.

These are the sites that need to be fixed back to `[0] * 9` (no
sentinel, direct indexing by `n`) per Section 1.7.c.

**q_bar BOO inner pair DRIFT (already catalogued in 1.13.d)** —
the type-4 branch of `initialize_interaction_data` at lines
7543–7546 initializes `q_bar[atom][m_idx]` as
`[[0.0, 0.0] for _ in range(m)]` — a 0-indexed `[real, imag]`
pair. Should be `[None, 0.0, 0.0]` per Perl's 1-indexed inner
pair.

**Everything else checks out.** No hidden slot-0 issues found.
Pass 3 confirms the 1-indexed convention is applied consistently
except in the two already-known places.

## Pass 4 — Loop bounds sweep

**Method**: grep every `range(` expression across the file,
categorize each by whether it is walking a 1-indexed persistent
array (should be `range(1, N+1)`), a 0-indexed step counter
(should be `range(0, N)`), or a legitimate axis loop.

### Pass 4 statistics

- **`range(1, 4)` (axis loops)**: ~40 occurrences, all in atom
  coordinate / lattice vector / abc-xyz axis iteration. ✓
- **`range(1, self.num_atoms + 1)` (atom loops)**: ~60 occurrences
  across readers, printers, analyzers, and manipulators. ✓
- **`range(1, num_items1 + 1)`, `range(1, num_items2 + 1)`
  (interaction item loops)**: ~10 occurrences in the interaction
  data family. ✓
- **`range(1, num_elements + 1)`, `range(1, num_species[eid] + 1)`
  (element/species loops)**: correct in `create_element_list`,
  `create_species_data`, and elsewhere. ✓
- **`range(1, 14)` (BOO l=6 m-index loop)**: 1 occurrence in
  `accumulate_q_bar` — correct 1-indexed iteration over the 13 m
  values for l=6. ✓
- **`range(1, ring_counts[ring_len])` (ring outer loop)**: correct
  1-indexed iteration inside `save_ring`. ✓

### Pass 4 — `range(3)` / `range(0, ...)` hits categorized

- **Lines 912, 920** (`set_scan_points` line mode) —
  Section 1.10.e DRIFT, already catalogued.
- **Lines 3943–3944** (`get_direct_xyz_point`) — Section 1.1.d
  DRIFT, already catalogued.
- **Lines 4011–4012, 4048** (`direct_xyz2fract_abc`,
  `direct_xyz2direct_abc`) — Sections 1.1.a, 1.1.b DRIFT.
- **Lines 4783, 4829, 4891, 4902, 4936, 4947** (inside
  `prep_surface`) — Section 1.11.c bridge sites, already
  catalogued.
- **Line 5886** (`rotate_one_point`) — Section 1.4.b DRIFT.
- **Lines 8580, 8582** (`get_plane_normal`) — Section 1.3.f DRIFT.
- **Line 8681** (`normalized_cross_product`) — Section 1.3.d
  DRIFT.
- **Lines 4659–4664** (`prep_surface` spiral search): these are
  step counters `range(max_rep[i] + 2)`, matching Perl `0..max+1`
  counters. Legitimate 0-indexed iteration. ✓
- **Lines 947–949** (`set_scan_points` mesh mode): step counters
  `range(0, nx)`, `range(0, ny)`, `range(0, nz)` matching Perl
  `0..numMeshPoints[1]-1`. ✓
- **Lines 7509, 7513** (`initialize_interaction_data` type 3):
  `range(num_items1 + 1)` / `range(num_items2 + 1)` — 0-indexed
  start to populate slot 0 with a sentinel `False`. Minor quirk
  noted in Section 1.10.d; produces 1-indexed skip arrays.

### Pass 4 findings

**No new DRIFT findings.** Every `range(3)` or `range(0, N)` hit
is either already catalogued as a drift (Sections 1.1, 1.3, 1.4,
1.7.i, 1.10.e, 1.11.c) or is a legitimate 0-indexed step counter
that matches Perl's 0-indexed iteration in the same role.

## Summary of Passes 2–4

**0 new DRIFT findings.** All three verification sweeps confirm
that the Pass 1 section-by-section audit was complete. The
catalog of indexing issues in `structure_control.py` is final at:

| Category | Count | Sections |
|---|---|---|
| DRIFT (Perl 1-indexed → Python 0-indexed) | 15 | 1.1 (×4), 1.3 (×6), 1.4.b, 1.7.i, 1.10.e, 1.12 (×2), 1.13.d |
| REVERSE DRIFT (Perl 0-indexed → Python 1-indexed) | 1 | 1.7.c |
| MIXED | 1 | 1.4.a |
| OK | ~60 | Sections 1.2, 1.5, 1.6, 1.7 (most), 1.8, 1.10 (most), 1.11 (most), 1.12 (most), 1.13 (most), 1.14, 1.15 |
| Signature-only refactors (not indexing) | 7 | 1.2.g, 1.2.l, 1.4.c, 1.6.a, 1.6.b, 1.6.c, 1.12.b/c |
| Call-site bridges | ~32 | Sections 1.9.a–e, 1.11.c |

The Pass 1 section files constitute the complete and final audit
of `structure_control.py`. Fix clusters are listed in the
[README.md](README.md) under "Fix-order recommendation."
