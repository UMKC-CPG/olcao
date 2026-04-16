# porosity.py — Pass 1 Audit

**Status: 1 DRIFT (PO-A).** Self-consistent internally, but
`main()` has to explicitly bridge to `sc.compute_pore_map` because
that method still expects the Perl 1-indexed layout.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 1 (PO-A) |
| REVERSE DRIFT | 0 |
| MIXED | 0 |

## Scope and method

- Perl original: `src/scripts/porosity` (~309 lines)
- Python port:   `src/scripts/porosity.py` (~551 lines)
- rc companion:  `src/scripts/porosityrc.py`

Reads a structure file, builds a 3-D sampling mesh over the unit
cell, marks every mesh point that falls inside the scaled covalent
radius of any atom (including periodic images), and reports the
fraction of unmarked mesh points as the system's porosity.

The `StructureControl` call surface is small: `read_input_file`,
`set_limit_dist`, `set_xyz_mesh_points`, `compute_pore_map`, plus
direct attribute reads on `sc.mag` and `sc.ext_pore_map`.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings | `ScriptSettings.*` | **1 DRIFT (PO-A)** |
| B | Resolution helpers | `compute_resolution`, `compute_mesh_and_refine_res` | Consumer of PO-A |
| C | Void counter | `calc_void_points` | Consumer of PO-A |
| D | Porosity compute | `compute_porosity` | OK |
| E | Main driver | `main` | Consumer of PO-A (bridge at line 532) |
| F | rc companion | `porosityrc.parameters_and_defaults` | Consumer of PO-A (rc side) |

## Finding PO-A — `settings.num_mesh_points` is 0-indexed

**Perl** `src/scripts/porosity` lines 50-57, 160-163, 198-229,
262-266, 296:
```perl
my @numMeshPoints;
...
$numMeshPoints[1] = $ARGV[++$number];
$numMeshPoints[2] = $ARGV[++$number];
$numMeshPoints[3] = $ARGV[++$number];
...
foreach $axis (1..3)
   {$axisResolution[$axis] = $mag_ref->[$axis]
                             /$numMeshPoints[$axis];}
...
foreach $aPoint (1..$numMeshPoints[1]) { ... }
...
StructureControl::setXYZMeshPoints(@numMeshPoints[1..3]);
StructureControl::computePoreMap(\@axisResolution,
                                 \@numMeshPoints);
```
`@numMeshPoints` is autovivified at slots `[1..3]`; slot 0 never
referenced. Passed by reference to `computePoreMap`, which also
reads slots `[1..3]`.

**Python** `src/scripts/porosity.py` touch points:

- **rc default** (`porosityrc.py:46`): `"num_mesh_points": [0, 0, 0]`
  — 0-indexed 3-element list.
- **Storage in ScriptSettings** (`porosity.py:128`):
  `self.num_mesh_points = list(rc["num_mesh_points"])`.
- **Storage after `-m` parse** (`porosity.py:240-242`):
  `self.num_mesh_points = list(args.num_mesh_points)` (argparse
  `nargs=3`).
- **Read in `reconcile`** (lines 253-257):
  `nmp[0] * nmp[1] * nmp[2]`.
- **Read in `compute_resolution`** (lines 301-309):
  ```python
  for axis in range(1, 4):
      axis_resolution[axis] = sc.mag[axis] / nmp[axis - 1]
  ```
  Explicit `axis - 1` bridge.
- **Write in `compute_mesh_and_refine_res`** (lines 359-362): same
  `axis - 1` bridge on writes.
- **Read in `calc_void_points`** (lines 400, 409-411): range walks
  `1..N` (matching the 1-indexed `sc.ext_pore_map`), but upper
  bound pulled from 0-indexed slots.
- **Bridge site in `main`** (lines 530-533):
  ```python
  # compute_pore_map expects 1-indexed lists:
  # [None, val_a, val_b, val_c].
  nmp_1idx = [None, nmp[0], nmp[1], nmp[2]]
  sc.compute_pore_map(axis_resolution, nmp_1idx)
  ```
  Exactly the bridge pattern the audit rule discourages.
  `compute_pore_map` is 1-indexed on both arguments (verified at
  `structure_control.py:6720-6830`).

**Severity**: DRIFT (self-consistent, no runtime bug). Directly
analogous to the `kp_mesh_scf` / `kp_mesh_pscf` DRIFT in
`mi_pass1.md` and the `abc_order` / `xyz_order` DRIFT in
`mod_struct_pass1.md`. Functionally closer to the
`makeinputrc.py` case (latent, not crashing) because
`compute_pore_map` currently receives a freshly built 1-indexed
bridge rather than the raw rc list.

**Recommended fix** — single cluster, 10 touch points:

1. `porosityrc.py:46` — `"num_mesh_points": [None, 0, 0, 0]`.
2. `porosity.py:128` — `self.num_mesh_points =
   list(rc["num_mesh_points"])` (copies 1-indexed layout unchanged).
3. `porosity.py:240-242` — prepend `None` after argparse:
   `self.num_mesh_points = [None] + list(args.num_mesh_points)`.
4. `porosity.py:253-257` — `[0], [1], [2]` → `[1], [2], [3]`.
5. `porosity.py:301-309` — drop `axis - 1` bridge; use `nmp[axis]`.
6. `porosity.py:359-362` — drop `axis - 1` bridge; write to
   `settings.num_mesh_points[axis]`.
7. `porosity.py:364-368` — `[0], [1], [2]` → `[1], [2], [3]`.
8. `porosity.py:400, 409-411` — `nmp[0..2]` → `nmp[1..3]`.
9. `porosity.py:524-525` —
   `sc.set_xyz_mesh_points(nmp[1], nmp[2], nmp[3])`.
10. `porosity.py:530-533` — delete `nmp_1idx` scratch; pass
    `settings.num_mesh_points` directly to `sc.compute_pore_map`.

After the fix, the chain from rc defaults through helpers and into
`compute_pore_map` is uniformly 1-indexed and matches the Perl
original line-for-line.

## Non-findings (verified consistent)

- **`axis_resolution`** — allocated `[None, 0.0, 0.0, 0.0]` (lines
  305, 340), written to slots `[1..3]` via `range(1, 4)` loops.
  Matches Perl `@axisResolution[1..3]`. OK.
- **`sc.mag`** reads at lines 308, 345, 355, 361, 450 — all
  1-indexed. Matches Perl `$mag_ref->[1..3]`. OK.
- **`sc.ext_pore_map` triple loop** (lines 409-414). Walks
  `1..N` matching Perl; map allocated as nested dict over keys
  `[1..nx]`, `[1..ny]`, `[1..nz]` at
  `structure_control.py:879-883`. OK.
- **Scalar `sc.*` calls** (`set_limit_dist`, `read_input_file`) —
  scalar arguments, no indexing concern. OK.
- **`record_clp`** — `sys.argv` iteration is legitimate 0-indexed.

## Non-indexing observations (out of scope)

- Python port drops `setupDataBaseRef` as a dedicated sub and
  inlines element-data init in `main()` (lines 494-495).
- Python `record_clp` adds a `Date:` timestamp line not in the
  Perl original.
- No `tests/` coverage for `porosity.py`.

## Running tally

- **BUG**: 0
- **DRIFT**: 1 (PO-A)
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every other subsystem

## Fix-order recommendation

Single cluster (PO-A), two files, 10 touch points, no external
callers, no test coverage. Regression verification after the fix
by inspection plus the full `structure_control` test suite.
