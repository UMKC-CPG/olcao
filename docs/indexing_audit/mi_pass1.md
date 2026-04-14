# makeinput.py — Pass 1 audit

## Scope and method

`src/scripts/makeinput.py` (~5855 lines, ~70 functions plus the
`ScriptSettings` class) is a straight port of the Perl `makeinput`
script (~5171 lines).  The audit strategy mirrors the
`structure_control.py` pass:

- **Pass 1** — every function that stores or consumes
  atom/axis/element/species/type indices, lattice vectors, target
  coordinates, or block borders is read and compared against the Perl
  original.  Call sites into `StructureControl` are checked for
  convention match after the Cluster A–I fixes in that module.
- **Pass 2** — grep sweep for `range(3)`, `[... - 1]`, `[... + 1]`,
  `[0]`, `range(0, ...)` and categorise each hit.
- **Pass 3** — slot-0 sentinel sweep across `ScriptSettings`
  attributes.

The findings below are grouped by subsystem instead of by a
section-per-file because the file is already structured function by
function with clear subsystem boundaries; splitting each one into its
own `mi_1.x.md` would mostly produce very short files.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | ScriptSettings data layout | `assign_rc_defaults`, `_parse_*` | 2 mild DRIFT + several OK |
| B | Cell initialisation | `initialize_cell`, `_auto_sybd_path`, `_auto_kp_shift` | OK |
| C | Target coordinate prep | `prepare_targets`, `make_min_dist_matrices` | OK (fixed in SC Cluster A) |
| D | Grouping methods — `group_block` | `group_block` | **BUG** (slot 0 vs None) |
| E | Grouping methods — `group_target`, `group_reduce` | `group_target`, `group_reduce`, helpers | OK |
| F | Relaxation | `relax_group`, `_relax_species`, `_relax_types` | OK |
| G | XANES bookkeeping | `assign_xanes_types`, `_prepare_xanes_list` | OK |
| H | Unit conversion | `_convert_a_to_au` | OK |
| I | K-point file writing | `_make_kp`, `_write_{mesh,density}_kp_file`, `_write_point_ops_block`, `_extract_point_ops`, `_print_kp_in_file` | 1 DRIFT (inner axis) + several OK |
| J | Atom sorting | `_sort_atoms`, `_get_cumulative_types`, `_get_num_atoms_of_type` | OK |
| K | Structure / olcao.dat output | `_print_structure`, `_print_olcao_input`, `_print_basis_set`, `_print_scf_pot`, `_get_elec_and_states`, `_process_sybd_path`, `_get_orbital_states`, `_count_orbital_lines` | OK |
| L | Memory estimate | `_compute_and_print_mem` | OK (pre-existing `kp_num` stays zero — latent bug, not indexing) |
| M | Skeleton / visualisation / submission | `_make_olcao_mi`, `_make_pdb`, `_make_cif`, `_make_sub_file` | OK |
| N | Summary | `print_summary` | OK |

Two findings total: one DRIFT (K-point mesh inner axis convention)
and one runtime BUG (`group_block` reads slot 0 of
`sc.direct_abc[atom]`).  Details below.

---

## A. ScriptSettings data layout — `assign_rc_defaults`, `_parse_*`

Most per-group settings arrays correctly use the 1-indexed
`[0 unused, 1=SCF, 2=PSCF]` layout:

- `kp_note = ["", "(General)", "(General)"]` — 1-indexed with
  sentinel `""` at slot 0. OK.
- `kp_density = [0.0, 0.0, 0.0]` — documented `[0]=unused`, used at
  `kp_density[1]` / `kp_density[2]`. OK.
- `kp_intg_code = [0, 0, 0]` — same layout, same comment. OK.
- `kp_num = [0, 0, 0]` — same layout, used at `kp_num[1]`/`kp_num[2]`.
  OK on indexing (separately flagged below — the variable is *never
  written* after initialisation, so the memory estimator reads zero
  for both SCF and PSCF k-point counts; this is a latent bug
  predating the Python port and orthogonal to the indexing audit).

### A.1 `kp_mesh_scf` / `kp_mesh_pscf` inner axis — DRIFT

**Perl** `@kpMesh` is a doubly 1-indexed array:
`$kpMesh[$kpGroup][$axis]` for `$kpGroup in 1..2` and
`$axis in 1..3` (declaration line 439 of `src/scripts/makeinput`,
first assignments at lines 944, 1034–1059).

**Python** `settings.kp_mesh_scf` and `settings.kp_mesh_pscf` are
flat three-element **0-indexed** Python lists.  The default values
come from `makeinputrc.py:103–104`:

```python
"kp_mesh_scf":      [1, 1, 1],
"kp_mesh_pscf":     [1, 1, 1],
```

Consumers read slots 0..2 directly, without a sentinel:

- `_write_mesh_kp_file` at `makeinput.py:3658` —
  `kp_mesh[0] kp_mesh[1] kp_mesh[2]`.
- `_print_olcao_input` at `makeinput.py:4434` —
  `settings.kp_mesh_pscf[0]`, `[1]`, `[2]`.
- `print_summary` at `makeinput.py:5671–5675` —
  `for axis in range(3): if settings.kp_mesh_scf[axis] >= 10`.
- `print_summary` at `makeinput.py:5740–5748` — same `[0]/[1]/[2]`.

Inside `_make_kp` there is a local 1-indexed outer view:

```python
kp_mesh = [None, settings.kp_mesh_scf, settings.kp_mesh_pscf]
for kp_group in range(1, 3):
    ...
```

So the **outer** group dimension (SCF vs PSCF) is correctly restored
to 1-indexed locally, but the **inner** axis dimension is 0-indexed
throughout.

**Severity**: DRIFT (inner axis only).  Per the audit rule the inner
list should be `[None, nx, ny, nz]`.

**Recommended fix**: change `makeinputrc.py` defaults to
`[None, 1, 1, 1]`, update the CLI parser (`_parse_kp`, `_parse_scfkp`
etc.) to write slots 1..3, and update every reader to use
`[1]/[2]/[3]` or `range(1, 4)`.  This is a systematic sweep across
roughly eight call sites and one default dictionary.

### A.2 `_parse_target` target `loc` format — mild convention choice

**Perl** stores `-atxyz`/`-atabc` target coordinates as a single
string `"x y z"` in `$targetLoc[$target]` (see `makeinput:1164`).
When needed, the consumer splits the string and reads
`@values[1..3]`.

**Python** `_parse_target` at `makeinput.py:1123–1129` stores the
coordinates as a three-element 0-indexed Python list
`t["loc"] = [x, y, z]`.

This is a **Python-side storage change**, not preservation of a Perl
convention.  `prepare_targets` then promotes the 0-indexed `loc` into
a 1-indexed `[None, x, y, z]` for `t["fract_abc"]` / `t["direct_xyz"]`
before passing to `StructureControl`, so the eventual consumer sees
the correct convention.

**Severity**: consistency note only — the `loc` slot only exists as
a transient buffer between command-line parsing and
`prepare_targets`, and it never crosses a `StructureControl` call
boundary.  **Not a finding**.  Leave as-is unless a larger
harmonisation pass is performed.

---

## B. Cell initialisation — OK

`initialize_cell` (`makeinput.py:1331`) copies per-atom and
per-element arrays from `StructureControl` into `settings` with the
correct 1-indexed layout (`[None]` sentinel at slot 0,
`.append(...)` for each entry).  `atom_type_id[0]` and
`num_types[0]` are built with `[None]` placeholders at both outer
and inner levels.

`_auto_sybd_path` (`makeinput.py:1475`) reads
`sc.angle_recip[1..3]`, `sc.full_cell_mag[1..3]`,
`sc.full_cell_angle[1..3]`, `sc.angle[1..3]` — all consistent with
`StructureControl`'s 1-indexed axis convention.

`_auto_kp_shift` (`makeinput.py:1650`) is string-based; no indexed
arrays involved.

---

## C. Target coordinate prep — OK (already fixed)

`prepare_targets` (`makeinput.py:1811`) and `make_min_dist_matrices`
(`makeinput.py:1876`) were updated during the `structure_control`
Cluster A fix sweep.  Both now build and pass 1-indexed
`[None, v1, v2, v3]` vectors into `direct_xyz2fract_abc`,
`fract_abc2direct_xyz`, and `create_min_dist_matrix`.  No drift
remains.

---

## D. `group_block` — BUG (runtime crash)

**File:line**: `makeinput.py:2148–2158`.

```python
for atom in range(1, settings.num_atoms + 1):
    abc = sc.direct_abc[atom]   # [a, b, c]

    if (abc[0] >= borders[0][0] - epsilon and
        abc[0] <= borders[0][1] + epsilon and
        abc[1] >= borders[1][0] - epsilon and
        abc[1] <= borders[1][1] + epsilon and
        abc[2] >= borders[2][0] - epsilon and
        abc[2] <= borders[2][1] + epsilon):
        status = 1
    else:
        status = 0
```

**The bug**: `sc.direct_abc[atom]` is a 1-indexed list
`[None, a, b, c]` (inherited from `StructureControl`, where every
per-atom coordinate row uses the `[None]`-sentinel layout).
Reading `abc[0]` returns `None`; Python 3 raises
`TypeError: '>=' not supported between instances of 'NoneType' and
'float'` on the very first comparison.  The remaining slot
references `abc[1]/abc[2]` correspond to a-axis/b-axis, not b/c as
the code intends, so even if the `None` compared cleanly the
per-axis tests would be shifted by one.

Meanwhile `borders` is constructed at `_parse_block:1158–1162` as a
3×2 **0-indexed** list:

```python
b["borders"] = [
    [tokens[i+1], tokens[i+2]],
    [tokens[i+3], tokens[i+4]],
    [tokens[i+5], tokens[i+6]],
]
```

**Perl** `groupBlock` (line 1918 of `src/scripts/makeinput`) uses
`$directABC_ref->[$atom][1..3]` and
`$blockBorders[axis][from/to][currentBlock]` where axis=1..3,
from/to=1..2 — fully 1-indexed on both axes.

The inner `from/to` is a genuine Perl 1-indexed pair (index 1 = from,
index 2 = to) that the Python port reasonably collapsed to a plain
2-element list — this is OK as an implementation detail.  The
**atom-axis** read is the problem.

**Severity**: **BUG**.  Any invocation with `-block` on the command
line would crash before the first atom is processed.  Since no test
in `tests/` exercises `-block`, the bug has been latent since the
initial Python port.

**Recommended fix**: read `sc.direct_abc[atom][1..3]` instead of
`abc[0..2]`.  The borders array remains 0-indexed (it is a local
bookkeeping structure that does not cross a `StructureControl`
boundary).  Either keep the borders as `[3][2]` 0-indexed and
translate on each read, or promote borders to
`[None] + [[None, from, to] ...]` for full symmetry with the Perl
original.  The minimal fix is the former:

```python
for atom in range(1, settings.num_atoms + 1):
    abc = sc.direct_abc[atom]  # [None, a, b, c]
    if (abc[1] >= borders[0][0] - epsilon and
        abc[1] <= borders[0][1] + epsilon and
        abc[2] >= borders[1][0] - epsilon and
        abc[2] <= borders[1][1] + epsilon and
        abc[3] >= borders[2][0] - epsilon and
        abc[3] <= borders[2][1] + epsilon):
        status = 1
    ...
```

A test that parses `-block` with real atoms should be added at the
same time.

---

## E. `group_target` / `group_reduce` / helpers — OK

- `group_target` at `makeinput.py:2170` reads
  `sc.min_dist[target_row][atom]` — a scalar distance, no axis
  access.  `target_row = settings.num_atoms + target_idx + 1` is a
  correct 0-based → 1-based bridge into the matrix.
- `group_reduce` at `makeinput.py:2221` is internally consistent:
  `atom_level`, `reduced_atoms[atom][level]`, `ra_elements[level]`,
  `ra_species[level]`, `temp_atom_species_id`, `temp_num_species`
  all use the standard `[None]`-sentinel or `[0]`-sized-`(N+1)`
  layouts.  Nested lists are 1-indexed with per-level sentinels.
  Phase 2's CA-matching loop uses a genuinely 0-indexed working
  buffer (`ca_species`, `ca_elements`), which is the same swap-and-
  pop idiom Perl uses (`@caSpecies` is 0-indexed in the Perl
  source).  OK.
- `compare_status_and_request`, `update_from_track_flag`,
  `extend_group`, `init_track_flag` all walk 1-indexed atom/element/
  species ranges correctly and write to 1-indexed settings arrays.
  OK.

---

## F. Relaxation — OK

`_relax_species` (`makeinput.py:2725`) and `_relax_types`
(`makeinput.py:2792`) build `species_flagger` / `type_flagger` and
`species_map` / `type_map` arrays using the `[None]`-sentinel
pattern on every nesting level.  Walks are `range(1, N+1)`.
Assignments honour the `atom_type_id[file_set][atom]` and
`num_types[file_set][element][species]` layouts.  OK.

---

## G. XANES — OK

`assign_xanes_types` (`makeinput.py:2867`) and
`_prepare_xanes_list` (`makeinput.py:3013`):

- `settings.xanes_atoms` is a flat **0-indexed** Python list of
  1-based atom numbers, documented at
  `makeinput.py:3027–3029`.  Translation to the 1-based file-set
  index `file_set = xa_idx + 1` is explicit at lines 2961–2965.
- `xanes_temp` (auto-select scratch buffer) is a 2-D
  `[None][None]`-sentinel structure.  OK.
- `sorted_xanes_atoms` is 1-indexed: index 0 is sentinel, indices
  `1..num_xanes_atoms` hold the sorted position of each XANES
  target.
- The `xanes_atoms[h - 1]` / `xanes_atoms[file_set - 1]` bridges at
  `makeinput.py:3330`, `4136`, `5781` all correctly translate from
  the 1-based file-set index to the 0-based `xanes_atoms` list
  index.  These are the documented semantics, not drift.

---

## H. Unit conversion — OK

`_convert_a_to_au` (`makeinput.py:3395`) iterates
`for abc in range(1, 4)` and `for xyz in range(1, 4)` and writes
`sc.mag[abc]`, `sc.real_lattice[abc][xyz]`,
`sc.direct_abc[atom][axis]`, `sc.direct_xyz[atom][axis]` —
every access uses the 1-indexed axis convention.  OK.

---

## I. K-point file writing

Inside `_make_kp` (`makeinput.py:3432`) the local view
`kp_mesh = [None, settings.kp_mesh_scf, settings.kp_mesh_pscf]`
uses a 1-indexed outer dim.  The group iteration
`for kp_group in range(1, 3)` indexes it correctly.  OK on the
outer dim.

The inner `kp_mesh[kp_group]` is still 0-indexed (see finding A.1).

`_extract_point_ops` (`makeinput.py:3517`) parses the space-group
database file and returns `point_ops` (list of 3×3 matrices) and
`frac_trans` (list of 3-vectors).  Both are plain 0-indexed lists
built from `split()` tokens.  This is a **legitimate 0-indexed
file-parse output**; the consumers (`_write_point_ops_block`,
`_write_mesh_kp_file`) read them as 0-indexed throughout.  OK.

`_write_mesh_kp_file` (`makeinput.py:3595`),
`_write_density_kp_file` (`makeinput.py:3666`), and
`_write_point_ops_block` (`makeinput.py:3733`) only drift on
`kp_mesh[0]/[1]/[2]` (already cataloged as A.1); `row[0]/[1]/[2]`
and `t[0]/[1]/[2]` are legitimate 0-indexed file-buffer reads.  OK
apart from A.1.

`_print_kp_in_file` (`makeinput.py:3775`) is the legacy writer for
the external `makeKPoints` program.  It reads
`sc.real_lattice[abc][1..3]` correctly and prints
`kp_mesh_request[0]/[1]/[2]` — same 0-indexed inner as A.1, same
fix applies.

---

## J. Atom sorting — OK

`_get_cumulative_types` (`makeinput.py:4014`) and `_sort_atoms`
(`makeinput.py:4051`) build all sorted arrays (`sorted_elem_id`,
`sorted_spec_id`, `sorted_type_id`, `sorted_fract_abc`,
`sorted_direct_abc`, `sorted_direct_xyz`) in the 1-indexed layout
with matching `[None, x, y, z]` inner rows.  `index_map` starts with
`[0]` as a sentinel position 0.  `datSkl.map` is written with
`new_pos` 1..num_atoms.  OK.

`_get_num_atoms_of_type` (`makeinput.py:4154`) allocates
`nat = [0] * (total_num_types + 1)` (1-indexed) and accumulates by
`seq_type`.  OK.

---

## K. Structure / olcao.dat output — OK

`_print_structure` (`makeinput.py:4183`) writes
`sc.real_lattice[abc][1..3]` and `sorted_direct_xyz[atom][1..3]`
directly.  The inner element name lookup uses
`settings.element_list[sorted_elem_id[atom]]` — all 1-indexed.  OK.

`_print_olcao_input` (`makeinput.py:4250`) iterates element,
species, and type loops with `range(1, N+1)`; writes
`num_states_used[1..3]`, `num_total_core_states[1..3]`,
`num_total_vale_states[1..3]`, `num_atoms_of_type[seq_type]` — all
1-indexed.  OK.

`_print_basis_set` (`makeinput.py:4525`) — same pattern.  The inner
`num_terms = [int(x) for x in basis_lines[line_idx].split()]` is a
0-indexed file-parsed token list; consumers read `num_terms[0]`,
`core_vals[0]`, `vale_vals[2]` — legitimate 0-indexed parse
buffers, matching the Perl `@values = split(...)` idiom.  OK.

`_print_scf_pot` (`makeinput.py:4825`) — same pattern as
`_print_basis_set`.  Per-file-set element/species/type loops use
1-indexed ranges; the potential file lines are consumed with
0-indexed `pot_lines[0..7]` and `nt_vals[0]` / `charge_vals[0]`
indices — legitimate file-parse buffer access.  OK.

`_get_elec_and_states` (`makeinput.py:4963`) uses
`num_total_core_states[1]` and `num_total_vale_states[basis]` for
`basis in range(1, 4)` — 1-indexed.  OK.

`_process_sybd_path` (`makeinput.py:4999`):

- `angle_name = ["", "alpha", "beta", "gamma"]` — 1-indexed with
  sentinel `""` at slot 0, iterated `for angle in range(1, 4)`.  OK.
- `variable_name = [None]` / `variable_eqn = [None]` with
  `.append(...)` per variable — 1-indexed.  OK.
- Inside the high-symmetry k-point loop (lines 5125–5138),
  `for axis in range(3): expr = vals[axis]; ...` is a **legitimate
  0-indexed file-token access**: `vals` is the tokenised line from
  `prep(sybd_lines[line_idx])`, and the three coordinate columns are
  tokens 0/1/2.  The Perl source uses `@values = split(...)` the
  same way.  OK.

`_get_orbital_states` (`makeinput.py:4736`) and
`_count_orbital_lines` (`makeinput.py:4788`) parse the
contracted-basis file and return a 4-element
`num_states = [0, 0, 0, 0]` — 1-indexed with slot 0 unused for
MB/FB/EB counts.  The inner `vals[0]`, `vals[1]`, `vals[3]` reads
are 0-indexed file-parse tokens, matching Perl.  OK.

---

## L. Memory estimate — OK on indexing (latent `kp_num` bug)

`_compute_and_print_mem` (`makeinput.py:5281`):

- `mem_setup`, `mem_main`, `mem_intg`, `mem_band`, `mem_pdos`,
  `mem_bond`, `mem_optc` are all 4-element 1-indexed lists with
  slot 0 unused, iterated `for bt in range(1, 4)`.  OK.
- `recip = [[None, 0.0, 0.0, 0.0] for _ in range(4)]` with
  `recip[0] = None` — fully 1-indexed 4×4.  OK.
- `prim_reps = [0, 0, 0, 0]` — 4-element, `prim_reps[1..3]`.  OK.
- `num_total_core_states[bt]` / `num_total_vale_states[bt]` —
  1-indexed, OK.

**Latent bug (not indexing)**: `settings.kp_num` is initialised to
`[0, 0, 0]` in `_make_kp:3462` and **never updated**; consequently
`kp1 = settings.kp_num[1]` and `kp2 = settings.kp_num[2]` are always
zero.  Every memory estimate that multiplies by `kp1` or `kp2`
returns zero for those terms, which under-estimates `mem_setup`,
`mem_main`, `mem_band`, and `mem_optc`.  Mesh mode is the affected
path; density mode skips memory estimation entirely.  This is a
pre-existing bug from the density-mode conversion — call it out in
a follow-up task but keep it out of the indexing-audit scope.

---

## M. Skeleton / visualisation / submission — OK

`_make_olcao_mi` (`makeinput.py:5175`) writes
`sc.mag[1..3]`, `sc.angle[1..3]`, and
`sorted_fract_abc[atom][1..3]` / `sorted_direct_xyz[atom][1..3]` —
all 1-indexed.  OK.

`_make_pdb` and `_make_cif` are thin wrappers around
`sc.print_pdb` / `sc.print_cif` — those methods were audited clean
in Section 1.10 of the `structure_control` audit.  OK.

`_make_sub_file` (`makeinput.py:5496`) is string-based; no indexed
arrays beyond dict lookups.  OK.

---

## N. `print_summary` — OK (inherits A.1 for kp_mesh display)

`print_summary` (`makeinput.py:5618`) reads `sc.mag[1..3]`,
`sc.angle[1..3]`, `sc.real_lattice[abc][xyz]` (fully 1-indexed),
`num_total_vale_states[1..3]`, `num_total_core_states[1..3]`,
`num_states_used[1..3]`, `settings.kp_note[1..2]`,
`settings.kp_intg_code[1..2]`, `settings.kp_density[1..2]` — all
correctly 1-indexed.

The one exception is `settings.kp_mesh_scf[0..2]` /
`settings.kp_mesh_pscf[0..2]` at lines 5671–5675 and 5740–5748 —
covered by finding A.1.  OK.

---

## Running tally

- **BUG** (interface convention mismatch causing runtime error): 1
  - Section D: `group_block` reads slot 0 of 1-indexed
    `sc.direct_abc[atom]`.
- **DRIFT** (Perl 1-indexed → Python 0-indexed): 1
  - Section A.1: `kp_mesh_scf` / `kp_mesh_pscf` inner axis; affects
    ~8 consumer sites plus the rc defaults.
- **OK**: every other subsystem.
- **Latent bug (not indexing)**: 1
  - Section L: `kp_num` never written; memory estimator reads 0.

## Fix-order recommendation

Once the user approves, the fix phase for `makeinput.py` should be:

1. **Cluster MI-A — `group_block` bug**.  Single-file change in
   `makeinput.py`; add a smoke test that exercises `-block` on a
   small structure.  Highest priority because it is a live runtime
   crash.
2. **Cluster MI-B — `kp_mesh` inner-axis 1-indexing**.  Update
   `makeinputrc.py` defaults, the argparse parsing, the
   `_make_kp` / `_write_*_kp_file` / `_print_kp_in_file` writers,
   and `print_summary`.  Roughly eight sites.
3. *(Out of scope: Cluster MI-Latent — `kp_num` never written.)*
   Flag as a separate TODO item for a follow-up quality pass.

## Fix status

Both clusters implemented.  Summary of the edits:

- **Cluster MI-A — `group_block`**: the border comparison now
  reads `sc.direct_abc[atom][1..3]` via an axis loop and bridges
  into the 0-indexed ``borders`` local with ``axis - 1``.  This
  matches Perl's 1-indexed row access while keeping ``borders``
  as a local 3×2 bookkeeping structure.  No test currently
  exercises ``-block``; a follow-up smoke test is still
  recommended but is out of the indexing-audit scope.
- **Cluster MI-B — `kp_mesh` inner axis**:
  - ``makeinputrc.py`` defaults updated to
    ``[None, 1, 1, 1]`` for both ``kp_mesh_scf`` and
    ``kp_mesh_pscf``.
  - ``ScriptSettings.reconcile`` now promotes argparse's flat
    ``[nx, ny, nz]`` into ``[None, nx, ny, nz]`` before storing
    in ``self.kp_mesh_scf`` / ``self.kp_mesh_pscf``.  The
    ``(Gamma)`` detection now compares ``list(args.kp) ==
    [1, 1, 1]`` before the promotion.
  - ``_write_mesh_kp_file`` reads ``kp_mesh[1..3]``.
  - ``_print_kp_in_file`` reads ``kp_mesh_request[1..3]`` and
    performs the ``[1, 1, 1]`` gamma check slot by slot.
  - ``_print_olcao_input`` MTOP block reads
    ``settings.kp_mesh_pscf[1..3]``.
  - ``print_summary`` loops ``range(1, 4)`` for the
    formatting-width probe and reads ``[1..3]`` for the final
    mesh print.

Tests: full ``tests/`` suite re-run after both fixes —
**295 passed, 23 skipped** (same as pre-fix baseline — no
``makeinput.py``-specific tests exist, so the regression
coverage is entirely through the ``structure_control`` suite
via the shared helpers that ``makeinput`` relies on).

Follow-up: the ``kp_num`` latent bug in ``_compute_and_print_mem``
(memory estimator reads zero for SCF/PSCF k-point counts in mesh
mode) remains open and should be cataloged as a separate TODO
item.
