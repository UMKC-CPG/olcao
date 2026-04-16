# makeSGDB.py — Pass 1 Audit

**Status: 1 DRIFT (MSG-A).** Script is self-consistent and
produces correct output; the drift is pure rule compliance on the
op-axis of the internal `symmetry_ops` / `symmetry_shifts` tables.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 1 (MSG-A) |
| REVERSE DRIFT | 0 |
| MIXED | 0 |

**No `structure_control` calls at all** — Pass 1 item 2 trivially
satisfied. Perl `StructureControl::prepLine` calls were replaced
with `str.split` / `re.split`.

## Finding MSG-A — `symmetry_ops` / `symmetry_shifts` op-axis — DRIFT

**Perl** `src/scripts/makeSGDB`:

- Declaration lines 36 (`@symmetryOps`), 43 (`@symmetryShifts`).
- `saveSymmetryOperation` lines 263–272: outer axes `(0..2)`
  0-indexed; stores at `[axis][axisContrib][$numSymmetryOps]`
  where `$numSymmetryOps` is pre-incremented to **1** before first
  store (lines 237–240).
- `applyShiftedRepetition` lines 368–384: reads `[$symmetryOp]`
  walking `1..$currentNumSymmetryOps`.
- `writeSymmetryComponents` lines 402–418: reads `[$op]` walking
  `1..$numSymmetryOps`.

Outer two axes (matrix row/column) are 0-indexed; **op-axis is
1-indexed with slot 0 as unused sentinel**.

**Python** `src/scripts/makeSGDB.py`:

- Allocation `read_symmetry_components:352-354` —
  `self.symmetry_ops = [[[], [], []] for _ in range(3)]`,
  `self.symmetry_shifts = [[] for _ in range(3)]`. Op dimension
  starts empty (no `[None]` sentinel).
- Builder `save_symmetry_operation:400-413` — four `.append()`
  calls, so first op lands at slot 0.
- Copy `apply_shifted_repetition:527-540` — `for sym_op in
  range(current_num_ops)` reading `[sym_op]` (0..N-1).
- Writer `write_symmetry_components:575-593` — `for op in
  range(self.num_symmetry_ops)` reading `[op]` (0..N-1).

All three consumers self-consistently walk 0-indexed, so output
file is byte-identical to Perl. **DRIFT only, not a bug.**

**Severity**: DRIFT (inner op-axis only). Internal to
`SpaceGroupDB`, no external callers, no test coverage.

**Recommended fix** — four touch points in one file:

1. `read_symmetry_components:352-354`: initialise each inner op
   list with `[None]` sentinel.
2. `save_symmetry_operation:400-413`: `.append()` calls unchanged
   (each op appends once, landing at correct 1-based slot after
   the sentinel).
3. `apply_shifted_repetition:527`: change to
   `range(1, current_num_ops + 1)`.
4. `write_symmetry_components:575`: change to
   `range(1, self.num_symmetry_ops + 1)`.

## Verified-OK subsystems

- **`ScriptSettings`, `main`, `make_db` driver** — OK.
- **`get_sub_group_name`** — string indexing / split buffers are
  genuinely 0-indexed in both languages.
- **`get_axis_op_components`** — `component_list = [0,0,0,0]` and
  the `chars` buffer are legitimate 0-indexed locals matching
  Perl exactly.
- **`shift_vector`** in `apply_shifted_repetition` — `[0,0,0]`,
  0-indexed in both languages (matches Perl `@shiftVector[0..2]`).
- **Axis loops** — every `range(3)` over the matrix row/column
  axis maps directly onto Perl `foreach $axis (0..2)`.
- **`sub_group` loop** — `range(1, num_sub_groups + 1)` matches
  Perl `foreach $subGroup (1..$numSubGroups)`; `chr(sub_group - 1
  + ord_a)` matches Perl.
- **`make_soft_links`** — pure string/OS operations.
- **Space-group data tables (audit item 5)** — the script is a
  streaming parser; no persistent keyed table of space groups in
  either version.

## Running tally

- **BUG**: 0
- **DRIFT**: 1 (MSG-A)
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every other subsystem

## Fix-order recommendation

Single cluster, single file, four touch points. No external
callers, no test coverage. Fix deferred per the coordinated-sweep
policy.
