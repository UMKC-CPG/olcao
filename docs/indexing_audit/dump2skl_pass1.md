# dump2skl.py — Pass 1 Audit

**Status: 1 DRIFT (DS-A).** Producer and every consumer are
self-consistent; no runtime hazard.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 1 (DS-A) |
| REVERSE DRIFT | 0 |
| MIXED | 0 |

## Scope

- Perl original: `src/scripts/dump2skl`
- Python port:   `src/scripts/dump2skl.py`
- rc companion:  `src/scripts/dump2sklrc.py`

**Scope note**: `dump2skl.py` imports `element_data.ElementData`
but does **not** import `structure_control`. There are zero
`sc.*` call sites to audit.

## Finding DS-A — `atom_coords` / `atom_types` outer atom axis — DRIFT

**Perl** `src/scripts/dump2skl`:

- Compaction loop at lines 426-437 explicitly uses
  `my $currentID = 1;` and walks `foreach $atomID (1..$highestID)`,
  producing 1-indexed `@atomCoords` / `@atomTypes`.
- Consumers at 508-513 (shift), 527-531 (fractional), 559-565 (skl
  write), and the `inBounds` / `maxCoord` / `minCoord` helpers all
  walk `foreach (1..$totAtoms)` with slot 0 unused.
- Inner axis of `@atomCoords[*][0..2]` is 0-indexed in Perl as well
  — legitimate, matches `@xyz_idx[0..2]`.

**Python** `src/scripts/dump2skl.py`:

- Builder at `read_frame:815-820` produces `atom_types = []` /
  `atom_coords = []` with no `None` sentinel; first valid atom at
  slot 0.
- Docstrings at lines 393, 436, 469, 673, 675 explicitly label the
  lists "0-indexed".
- Six consumers walk `range(tot_atoms)` or the seed-and-walk
  `range(1, tot_atoms)` pattern:
  - `in_bounds` — line 405
  - `max_coord` — lines 448-451
  - `min_coord` — lines 482-485
  - `process_coordinates` — lines 895-898 and 908-911
  - `write_skl` — lines 979-985

**Severity**: DRIFT (outer atom axis only; inner `[x,y,z]` axis
matches Perl's legitimate 0-indexed axis convention).

**Runtime impact**: None. Producer and every consumer are
self-consistent within the file; the lists never cross a
`structure_control` boundary.

**Recommended fix**: Seed `atom_types = [None]` and
`atom_coords = [None]` in `read_frame`, rewrite the six consumer
loops to walk `range(1, tot_atoms + 1)` (and seed `max_coord` /
`min_coord` at `atom_coords[1][axis]`, loop
`range(2, tot_atoms + 1)`), and update the six docstrings to say
"1-indexed with `None` at slot 0".

## Non-findings (verified OK against Perl)

- **`box_bounds[axis][0..1]`** (`read_frame:737-744`, `in_bounds`,
  `process_coordinates`) — Perl `@boxBounds[0..2][0..1]` is doubly
  0-indexed; Python matches exactly. Axis `0..2`, inner `[0]=lo,
  [1]=hi`. OK.
- **`lattice[i][j]`** (`read_frame:718-754`, `main:1076-1085`) —
  Perl `@lattice[0..2][0..2]` is doubly 0-indexed in both the
  triclinic and orthogonal branches; Python matches. OK.
- **`xyz_idx[0..2]`** (`read_frame:771`) — Perl 0-indexed, Python
  `[None, None, None]` 3-element list, 0-indexed access. Matches
  Perl. OK.
- **`elements` (type→name map)** — Perl
  `@elements[$atomType]` 1-indexed via `foreach $atomType
  (1..$numAtomTypes)`. Python `get_lammps_types:317-360` uses a
  **dict** keyed by 1-indexed `atom_type` integer, read by
  `write_skl:980` as `elements[atom_types[i]]`. Dict semantics
  preserve the 1-indexed convention exactly. OK.
- **`timesteps` (frame→timestep)** (`count_frames:518-534`,
  `resolve_frame_and_timestep:599, 604-605`) — Python seeded with
  `[None]`, appended to, read with 1-indexed frame numbers. Matches
  Perl exactly. OK.
- **`ElementData` boundary** (`get_lammps_types:302-356`) —
  `ed.atomic_masses` / `ed.element_names` both 1-indexed; Python
  loop `for z in range(1, num_elements + 1):` mirrors Perl. OK.
- **`num_atom_types` loop** at `dump2skl.py:341` —
  `range(1, num_atom_types + 1)` matches Perl `1..$numAtomTypes`.
- **`raw_types` / `raw_coords` scratch dicts**
  (`read_frame:791-808`) — keyed by original LAMMPS atom IDs
  (already 1-indexed); local compaction buffers. OK.
- **`dump2sklrc.py`** — defaults dict at lines 18-44 stores six
  **scalars only**. No list-valued default crosses an indexing
  boundary. OK.

## Non-indexing observations (out of scope)

- **Triclinic box is not fully processed** — `main:1060-1090`
  falls back to lattice-vector magnitudes and emits a stderr
  warning. Perl had the same limitation.
- **`--np` flag rename** — Python uses long-form `--np` to avoid
  argparse conflict with `-n`/`--name`.
- **`process_coordinates` both mutates in place and returns** —
  documented.

## Fix-order recommendation

Single cluster (DS-A), single file, ~15 touch points: seed
`atom_coords` / `atom_types` with `[None]` in `read_frame`, rewrite
six consumer loops and two seed expressions, update six docstring
slot-0 labels. No external callers, no test coverage.

## Running tally

- **BUG**: 0
- **DRIFT**: 1 (DS-A)
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every other subsystem
