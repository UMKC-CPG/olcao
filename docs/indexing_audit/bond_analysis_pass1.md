# bond_analysis.py — Pass 1 audit

## Scope and method

`src/scripts/bond_analysis.py` (~3429 lines) is a Python port of
the Perl `bondAnalysis` script (~3378 lines).  It reads a
structure file, runs bond / bond-angle / coordination /
statistics / visualization computations, and emits one of about
a dozen output formats (ball-and-stick, VTK, OpenDX, OpenSCAD
full, OpenSCAD parts, ring distribution, Q^n distribution,
etc.).

The file is organised around a single `BondAnalysis` class
(plus a `ScriptSettings` companion).  The class has ~30 methods,
most of which are output writers that iterate over the
`StructureControl` bonding arrays (`bonded`, `bond_length_ext`,
`bonded_ext`, `bond_angles_ext`, `angle_bonded`, etc.) and the
per-atom data (`direct_xyz`, `atom_element_name`, etc.) to
produce text files.

Pass 1 approach:

- `Grep` sweep for drift signals: `range(3)`, `range(0, …)`,
  `[atom - 1]`, `[axis - 1]`, `[0] +=`, `.append([` with inner
  nested data.
- Read every hot spot the sweep flagged plus the main
  structural setup methods (`_setup_coords`, `_set_box_border`,
  `run`, `__init__`).
- Compare the remaining non-trivial blocks
  (`_get_bonds_to_show`, `_compute_statistics`,
  `_print_openscad`, `_write_bav_module`, `_write_rod`,
  `_print_bond_angles`, `_print_bond_lengths`, `_print_bond_dx`,
  `_print_vtk_ball_and_stick`) against the Perl original via
  `src/scripts/bondAnalysis`.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings | `ScriptSettings.*`, `_parse_mode` | OK |
| B | Class setup | `__init__`, `run`, `_setup_coords`, `_set_box_border` | OK |
| C | Ball-and-stick output | `_print_ball_and_stick`, `_print_vtk_ball_and_stick` | OK |
| D | Bond angle output | `_print_bond_angles`, `_write_angle_histogram` | OK |
| E | Bond length output | `_print_bond_lengths`, `_write_length_histogram` | OK |
| F | Coordination / Q^n / ring output | `_print_coordination_data`, `_print_qn_distribution`, `_print_ring_distribution`, `_print_bond_oo` | OK |
| G | OpenDX output | `_print_bond_dx` | OK |
| H | OpenSCAD output | `_print_openscad`, `_write_bav_module`, `_write_rod_module`, `_write_rod_label_module`, `_write_sphere_label_module`, `_write_trans_sphere_module`, `_write_trans_key_module`, `_write_trans_sub_cube_module`, `_write_bav_module`, `_write_union_atoms_module`, `_write_rod`, `_write_sphere`, `_add_atom_to_sheet` | OK (one stylistic note below) |
| I | Bond list assembly | `_get_bonds_to_show` | **1 DRIFT** (`bond_info_temp` inner) |
| J | Statistics | `_compute_statistics`, `_print_bond_stats` | OK |
| K | Warning helpers | `_element_warning`, `_species_warning`, `_atom_warning` | OK |

One DRIFT finding.  No runtime bugs.  Every other subsystem
preserves Perl 1-indexed conventions correctly (or uses a
legitimate 0-indexed buffer at a file-format boundary).

## Finding BA-A — `bond_info_temp` inner 3-element entry

**Perl** `src/scripts/bondAnalysis` lines 2694–2783
(`sub getBondsToShow`):

```perl
my @bondInfoExtTemp1;
...
$numBondsExt = -1;     # starts at -1 so the first bond lands at slot 0
foreach $atom (1..$numAtoms)
{
   $bondedAtom1 = $central2ExtItemMap_ref->[$atom];
   foreach $bond (1..$numBonds_ref->[$atom])
   {
      $numBondsExt++;
      if ($bondedAtom1 < $bondedAtom2)
      {
         $bondInfoExtTemp1[$numBondsExt][1] = $bondedAtom1;
         $bondInfoExtTemp1[$numBondsExt][2] = $bondedAtom2;
      }
      ...
      $bondInfoExtTemp1[$numBondsExt][3] =
            $bondLengthExt_ref->[$atom][$bond];
   }
}

@bondInfoExtTemp2 = sort{$a->[2] <=> $b->[2]} @bondInfoExtTemp1;
@bondInfoExt = sort{$a->[1] <=> $b->[1]} @bondInfoExtTemp2;
```

The **outer** dimension of `@bondInfoExtTemp1` is 0-indexed on
purpose (the variable `$numBondsExt` starts at `-1` specifically
so the first increment lands the first real entry at slot 0,
which matches the Perl `sort` comparator's expectations — sort
output is itself 0-indexed).  The **inner** dimension, however,
is 1-indexed: slot 1 is the first bonded atom, slot 2 is the
second bonded atom, and slot 3 is the bond length.

**Python** `_get_bonds_to_show` at `bond_analysis.py:3001`:

```python
bond_info_temp = []
for atom in range(1, sc.num_atoms + 1):
    ...
    for bond in range(1, sc.num_bonds[atom] + 1):
        ...
        bond_info_temp.append([a1, a2, bl])

bond_info_temp.sort(key=lambda x: x[1])  # a2
bond_info_temp.sort(key=lambda x: x[0])  # a1

for entry in bond_info_temp:
    a1, a2, bl = entry
    ...
```

The Python `.append()` builds the same 0-indexed outer
dimension — aligned with Perl.  The **inner** 3-element list,
however, is flattened to 0-indexed `[a1, a2, bl]` instead of
the 1-indexed `[None, a1, a2, bl]` that would mirror Perl's
`[1..3]`.

**Severity**: DRIFT (inner axis only).  This is a **local
temporary** scratch buffer whose only consumers are the two
`sort` calls in the same function and a single unpack
iteration below.  The eventual output (`reduced_bond_info`)
is correctly 1-indexed `[None, ext_atom1, ext_atom2,
bond_length]`, so the drift does not leak past
`_get_bonds_to_show`.  No tests cover this method either, and
producer/consumer pairs are self-consistent — so the fix is
purely rule compliance.

**Recommended fix**: store `[None, a1, a2, bl]` at the append
site, change the two `sort` keys to `lambda x: x[2]` (a2) and
`lambda x: x[1]` (a1), and replace `a1, a2, bl = entry` with
`a1, a2, bl = entry[1], entry[2], entry[3]` (or
`_, a1, a2, bl = entry`).

## Non-findings (verified consistent)

- **`box_borders`** — doubly 1-indexed `[None, [None, min, max],
  ...]` at `__init__` (line 465).  Every read in
  `_set_box_border` uses slots `[axis][1]` / `[axis][2]` with
  `axis in 1..3`.  Matches Perl `$boxBorders[$axis][1..2]`.
  OK.

- **`reduced_bond_info`** — **1-indexed** `[None]` sentinel
  outer, each entry is `[None, ext_atom1, ext_atom2,
  bond_length]` (see `_get_bonds_to_show:3077`).  Every
  consumer reads slots `[1]`, `[2]`, `[3]` — in
  `_print_bond_dx` (lines 1815–1829), `_print_openscad` and
  `_write_bav_module` (lines 2097–2125, 2346–2402, 2645–2695).
  OK.

- **Sort-by-slot idiom** — when Python sorts a list of
  1-indexed lists (e.g. `sorted(..., key=lambda x: x[1])`
  would sort by slot 1 = first meaningful axis), it works fine
  as long as slot 0 is a consistent sentinel.  `bond_info_temp`
  is the only place where the inner is currently 0-indexed
  instead of 1-indexed.

- **`histogram = [0] * (num_points + 1)`** in
  `_write_angle_histogram` and `_write_length_histogram` —
  1-indexed with slot 0 unused, `bucket` always ≥ 1.  OK.

- **Statistics arrays (`avg_elem_bl`, `std_elem_ba`,
  `avg_spec_bl`, etc.)** — allocated as `[-1] * (ne + 1)` or
  `[-1] * (na + 1)` at `_compute_statistics`
  (lines 3106–3132, 3228–3252).  Outer 1-indexed; slot 0 holds
  a harmless `-1` sentinel that is never read.  The 2D species
  arrays (`self.avg_spec_bl[element][species]` etc.) are
  initialised to `[-1] * (ns + 1)` in the element loop, so
  both dimensions are 1-indexed.  Pass-1 and Pass-2
  accumulators (`num_elem_bonds`, `acc_elem_bl`, etc.)
  follow the same layout.  Every read and write uses
  `[elem]` / `[elem][spec]` / `[atom]` with 1-indexed values.
  OK.

- **`_get_bonds_to_show` BO list** — `bond_bo_list` is a
  0-indexed list of `(int, int)` tuples built from
  `values[2]`, `values[3]` of the bond-order file.  Perl
  `$bondBOList[$i][0]` / `[1]` is also 0-indexed inner (see
  `bondAnalysis:2732–2733`).  Matches Perl, legitimate
  0-indexed parse buffer.  OK.

- **VTK atom index translation** — `a1 = atom - 1` /
  `a2 = bonded_atom - 1` at `_print_vtk_ball_and_stick:1281–1282`
  and `atom_c_index = atom - 1` at line 1326 are
  legitimate 0-indexed output-format bridges (VTK arrays start
  at 0).  Matches Perl.  OK.

- **OpenDX atom index translation** — `reduced_bond_info[bond]
  [1] -= 1` / `[2] -= 1` at `_print_bond_dx:1815–1817` is the
  equivalent bridge for OpenDX (also 0-indexed).  Matches
  Perl.  OK.

- **`sc.color_vtk[atom][0..2]` RGBA read** at
  `_print_vtk_ball_and_stick:1356–1358` — legitimate
  0-indexed inner pair, matching the exception noted in SC
  audit Section 1.8.g.  OK.

- **`_write_rod` / `_write_sphere` / `_write_bav_module`** —
  all Cartesian vectors passed around internally are 1-indexed
  `[None, x, y, z]` or `[0, x, y, z]`; every read uses slots
  `[1..3]`.  The single exception is the style quirk at
  `_write_rod` (lines 2793, 2795) where sentinel values are
  `0` instead of `None`:

  ```python
  normal = sc.get_plane_normal(
      [0, 0, 0, 0],
      shifted_p2,
      [0, 0, 0, rod_length],
  )
  ```

  The slot-0 sentinel is used as a harmless 0 instead of a
  `None` in the other SC call sites.  Both work — the
  1-indexed vector helpers in SC only read slots 1..3 — but the
  inconsistency is cosmetic.  Not an audit finding.

- **`sheet_pos1 = [0, 0.0, 0.0, bond_radius]`** in
  `_print_openscad:2299–2300` — same slot-0=0 convention, same
  harmless cosmetic inconsistency.  OK.

- **Per-atom `num_*` and `acc_*` arrays** in
  `_compute_statistics` — 1-indexed with appropriate
  initialisation.  OK.

- **Q^n reader** at `_print_qn_distribution:1709–1720` reads
  `sc.absolute_sys_qn[n]` and `sc.fractional_sys_qn[n]` for
  `n in range(0, 9)`.  This is the **0-indexed** layout after
  the structure_control Cluster C reverse-drift fix (see
  `sc_1.7_remaining_pass1.md`).  Matches the fixed SC
  convention where `n` is literally the index.  OK.

- **`_parse_mode` sub-option loops** — `range(2)`, `range(4)`,
  `range(6)` at lines 728, 769, 799 are simple counter loops
  bounding the number of sub-options greedily consumed by each
  mode flag.  No array indexing concerns.  OK.

## Non-indexing observations (out of audit scope)

- **Slot-0 sentinel style inconsistency.**  Several vector
  locals in the OpenSCAD / rod writers use `[0, …]` as the
  slot-0 sentinel (e.g. `sheet_pos1`, the `[0, 0, 0, 0]` passed
  to `get_plane_normal` in `_write_rod`), while the rest of
  the codebase uses `[None, …]`.  Both are functional; the
  inconsistency is cosmetic.  Flag for a future harmonisation
  pass.

- **`self.avg_*` / `self.std_*` initialised to `[]` in
  `__init__`**.  `_compute_statistics` reinitialises them to
  `[-1] * (ne + 1)` etc., so by the time `_print_bond_stats`
  reads them they are 1-indexed.  However, if the operation is
  not `OP_STATISTICS` the arrays remain empty; any future
  caller that reads them without first invoking
  `_compute_statistics` would hit an IndexError.  Not a
  current issue — `_print_bond_stats` is the only reader and
  it is only dispatched from the stats path.

## Running tally

- **BUG**: 0
- **DRIFT**: 1 — finding BA-A (`bond_info_temp` inner
  3-element layout in `_get_bonds_to_show`).
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every other subsystem.

## Fix-order recommendation

Single cluster, single file, minimal edit:

1. **Cluster BA-A** — `bond_info_temp` inner axis.  Three
   touch points in `_get_bonds_to_show` (one builder, two sort
   keys, one unpack); no external callers or test coverage.
   Trivial to apply.

No test coverage exists for `bond_analysis.py`.  Regression
verification after the fix is by inspection plus the full
`structure_control` test suite (which covers the shared
helpers `bond_analysis.py` calls into: `create_bonding_list`,
`compute_bond_angles_ext`, `create_q_list`, `compute_qn`,
`compute_ring_distribution`, and the bonding / coordinate data
structures they populate).

## Fix status

**Cluster BA-A applied.**  `_get_bonds_to_show` now builds
`bond_info_temp` with 1-indexed inner entries
``[None, a1, a2, bl]``, the two sort keys read slots ``x[2]``
(a2) then ``x[1]`` (a1), and the single unpack iteration uses
``a1, a2, bl = entry[1], entry[2], entry[3]``.

Verification: `ast.parse` on `bond_analysis.py` passes; full
`tests/` suite re-run is **295 passed, 23 skipped** — unchanged
from the pre-fix baseline (as expected, since no
`bond_analysis.py`-specific tests exist and the shared
`structure_control` helpers are unaffected by the local
bookkeeping change).
