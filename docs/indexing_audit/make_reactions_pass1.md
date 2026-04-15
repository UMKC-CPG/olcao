# make_reactions.py — Pass 1 audit

## Scope and method

`src/scripts/make_reactions.py` (~3130 lines) is a Python port
of the Perl `makeReactions` script (~3200 lines).  It builds
LAMMPS `bond/react` pre-reaction / post-reaction template
files plus a mapping file from two input molecule directories
that each contain an OLCAO skeleton (.skl) file.  The pipeline
is:

1. `finalize_env` — arrange molecule directories, create the
   reaction directory, copy the .skl files in.
2. `prep_molecules` — centre each molecule in its own cell,
   unify cell sizes so both fit.
3. `catalog_s` — identify surface (S) atoms, trigger atoms,
   binding atoms, and the per-S-atom delete tree.
4. `generate_merged_mols` — translate/rotate molecule 2 to
   align with molecule 1 at each S-atom pair, using
   `mod_struct.py` as a subprocess for the geometry
   transformations.
5. `catalog_elem_spec_mol_id` — assign element names, species
   IDs, and molecule tags to every atom in each merged file.
6. `prune_merged_molecule` — BFS-prune atoms beyond the chain
   length from either S atom, record edge atoms.
7. `make_reaction_templates` — read bond / bond-angle analyses
   for every pruned merged file and emit pre- and
   post-reaction templates plus the reaction map file.

The file has four top-level classes (`AngleData`,
`ScriptSettings`, `MakeReactions`) with ~25 methods total.

Pass 1 approach:

- `Grep` sweep for drift signals: `range(3)`, `range(0, …)`,
  `[atom - 1]`, `[0]`-prefixed writes, `.append([…])` with
  non-`None`-prefixed inner lists.
- Read every method end-to-end, cross-referencing with the
  Perl source at `src/scripts/makeReactions` at every
  non-trivial data-structure decision.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | Angle coefficient database | `AngleData.init_angle_data` | OK — already 1-indexed inner from the start |
| B | CLI / settings | `ScriptSettings.*` | OK |
| C | Pipeline orchestration | `MakeReactions.__init__`, `run`, `finalize_env`, `prep_molecules` | OK |
| D | S-atom cataloguing | `catalog_s`, `_add_atoms_to_delete` | OK |
| E | Merged molecule generation | `generate_merged_mols`, `_merge_skl_files` | OK |
| F | Element/species/mol-id catalog | `catalog_elem_spec_mol_id` | OK |
| G | Pruning | `prune_merged_molecule`, `_mark_atoms_to_keep`, `_write_pruned_skl` | OK |
| H | Reaction template builders | `make_reaction_templates`, `_read_bond_data`, `_read_angle_data`, `_compute_ordered_species`, `_build_phase_angles`, `_build_phase_bonds`, `_check_bonding`, `_check_angle_bonding`, `_write_template`, `_write_map_file` | OK |

**Zero findings.**  No BUGs, no DRIFT, no reverse drift, no
mixed.

## Why zero findings

The coding pattern here is noticeably more disciplined than
the earlier ports in the priority list:

- Every nested-list data structure (`del_tree`, `bonded`,
  `bond_tag_id`, `angle_tag_id`, `bond_angles_ext`,
  `angle_bonded`, `phase_bond_list`, `phase_angle_list`,
  `unique_bond_tags`, `unique_angle_tags`, etc.) is
  initialised with `[None]` sentinels and extended via
  `.append([None, …])` so that both dimensions are 1-indexed.
- `AngleData.hooke_angle_coeffs` is built as
  `[None, Z1, Zv, Z2, k, angle, tol]` on the inner level —
  this is the same shape that `condense.py` **needed** to be
  fixed up to (Cluster CD-B), but `make_reactions.py` got it
  right from the initial port.  The readers at
  `hc[1..6]` in `_read_angle_data` match.
- All atom / S-atom / element / species loops use
  `range(1, N + 1)`.
- Merged and pruned molecule data are keyed by tuples
  (`(s1, s2, atom)`, `(s1, s2, atom, phase)`, etc.) rather
  than nested lists — there is no inner-dim choice to
  make.
- Every call boundary into `StructureControl`
  (`sc.read_input_file`, `sc.set_buffer`,
  `sc.compute_crystal_parameters`, `sc.print_olcao`,
  `sc.read_bond_analysis_bl`, `sc.read_bond_analysis_ba`,
  `sc.bonded[atom][bond]`, `sc.direct_xyz[atom][1..3]`, etc.)
  is used with 1-indexed arguments matching the post-Cluster-A
  interfaces.

## Non-findings (verified consistent)

A few items look like drift at first glance but aren't:

### N.1 `trans_vector1` / `trans_vector2` — legitimate 0-indexed

`generate_merged_mols` builds two local 3-element 0-indexed
lists at lines 1294–1297 and 1340–1344:

```python
trans_vector2 = [
    diff_vector1[ax] * self.s.ss_dist
    for ax in range(1, 4)
]
...
trans_vector1 = [
    self.s_atom_coords1[s1][ax + 1]
    - self.s_atom_coords2[s2][ax + 1]
    for ax in range(3)
]
```

These look like inner-axis drift, but Perl's
`@transVector1` / `@transVector2` are **also 0-indexed**
(see `makeReactions:1042, 1082`):

```perl
$transVector2[$axis-1] = $diffVector1[$axis] * $ssDist;
...
$transVector1[$axis-1] = $sAtomCoords1[$sAtom1][$axis]
```

The reason is that Perl's string interpolation
`"@transVector1"` produces a space-separated concatenation of
all elements, which is exactly the argument format needed by
the subprocess call to `mod_struct.py`.  Both arrays are used
purely as flat buffers of subprocess arg strings
(`str(trans_vector1[0])`, `str(trans_vector1[1])`,
`str(trans_vector1[2])`).  They never cross a
`StructureControl` boundary.

**Verdict**: legitimate 0-indexed.  Matches Perl.  Not a
finding.

### N.2 `plane_coords`, `origin` — legitimate 0-indexed

Same story at lines 1419–1429 and 1441–1445:

```python
plane_coords = [
    self.s_atom_coords1[s1][1],
    self.s_atom_coords1[s1][2],
    self.s_atom_coords1[s1][3],
    self.trigger_coords1[s1][1],
    ...
    final_coord[3],
]

origin = [
    self.s_atom_coords1[s1][1],
    self.s_atom_coords1[s1][2],
    self.s_atom_coords1[s1][3],
]
```

Perl uses `@{$sAtomCoords1[$sAtom1]}[1..3]` (slice syntax) to
produce a 0-indexed flat array for the same purpose
(`makeReactions:1151, 1162`).  Consumers are subprocess args
only — `str(plane_coords[0])`, `str(origin[0])`, etc.

**Verdict**: legitimate 0-indexed.  Matches Perl.  Not a
finding.

### N.3 `diff_vector1`, `diff_vector2`, `temp_vector`, `final_coord` — 1-indexed

The 1-indexed counterparts were fixed during the earlier
Cluster A sweep of `make_reactions.py` alongside
`structure_control.py` — they feed `sc.get_vector_angle` and
`sc.normalized_cross_product` directly, so they need the
1-indexed `[None, x, y, z]` layout.  `generate_merged_mols`
currently uses the correct shape for all four:

```python
diff_vector1 = [None, 0.0, 0.0, 0.0]
for axis in range(1, 4):
    diff_vector1[axis] = ...
...
final_coord = [None] + [
    trans_vector1[ax]
    + self.trigger_coords2[s2][ax + 1]
    for ax in range(3)
]
...
diff_vector2 = [None] + [... for axis in range(1, 4)]
...
temp_vector = [None, diff_vector2[1]*2.0, ..., diff_vector2[3]*4.0]
```

All fed straight into the 1-indexed SC helpers.  OK.

### N.4 `angle_set` — legitimate 0-indexed local

`_build_phase_angles` line 2669–2673:

```python
angle_set = [
    angle_bonded[atom][angle][1],
    atom,
    angle_bonded[atom][angle][2],
]
```

Consumed two lines later as `angle_set[0]` / `angle_set[2]`.
Perl's `@angleSet[0..2]` at `makeReactions:2520–2522` uses
the exact same 0-indexed layout.

**Verdict**: legitimate 0-indexed.  Matches Perl.  Not a
finding.

### N.5 `angle_bonded[curr_atom][angle]` entries `[None, None, None]`

`_read_angle_data` line 2484 initialises each new entry as
`[None, None, None]` before overwriting with the real
`[None, a1, a2]` / `[None, a2, a1]` at lines 2499–2516.  The
3-element placeholder is slightly asymmetric with the
2-element 1-indexed content that follows (slot 3 is left as
`None` and never read), but the layout is compatible and
matches the `angle_bonded[atom][angle][1..2]` read pattern
used in `_build_phase_angles` and `_check_angle_bonding`.

**Verdict**: cosmetic asymmetry only.  Not a finding.

## Non-indexing observations (out of audit scope)

Flagged for a separate follow-up backlog:

- **`angle_bonded` placeholder initialisation** at
  `_read_angle_data:2484` initialises each entry to
  `[None, None, None]` then overwrites with
  `[None, a1, a2]` or `[None, a2, a1]`.  The three-element
  placeholder is a vestige of an earlier draft; since the
  consumer code (in `_build_phase_angles`,
  `_check_angle_bonding`, `_write_template`) only ever reads
  slots 1 and 2, the placeholder could be tightened to
  `[None, None, None]` → `[None, 0, 0]` (or simply removed by
  pre-sizing the list to `angle + 1`).  Cosmetic only.

- **`sc.set_buffer(total_buffer)` is wired correctly** at
  `prep_molecules:844` / `866`, in contrast to the pending
  `pdb2skl.py` follow-up where the same setter is accepted
  on the CLI but never applied.  Mentioned here purely for
  comparison — this is the good pattern.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every subsystem.

**No fix phase required** for `make_reactions.py`.
