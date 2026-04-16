# insert3cbo.py — Pass 1 Audit

**Status: 3 DRIFT, 0 BUG.** All three drifts are locally self-
consistent (producer/consumer agree), so the script produces
byte-identical output to Perl; the drifts are pure rule violations.

## Scope and method

`src/scripts/insert3cbo.py` (~1018 lines) is a Python port of the
Perl `insert3cbo` script (~611 lines). It rewrites a 2-center bond
order raw file to subsume bonds that belong to a 3-center bond into
fake centroid "atoms", and emits a modified `structure.dat` with
those atoms appended. Both Perl and Python carry a FIX-header note
that the program is broken until re-merged with the current makeBOND
pipeline, and no tests exercise it.

Key observation: `insert3cbo.py` makes **no** `sc.*` /
`StructureControl` calls. The Perl original uses
`StructureControl::prepLine` as a line tokenizer only; the port
inlines that as a local `prep_line` helper. No cross-script
1-indexed call boundary exists in this file — every drift is purely
local producer/consumer bookkeeping.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 3 (IC-A, IC-B, IC-C) |
| REVERSE DRIFT | 0 |
| MIXED | 0 |

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings | `ScriptSettings.*` | OK |
| B | rc companion defaults | `insert3cborc.parameters_and_defaults` | OK |
| C | 2CBO reader — per-atom scalar fields | `read_2cbo` | OK |
| D | 2CBO reader — orbital charges inner | `read_2cbo` | **DRIFT (IC-A)** |
| E | 2CBO reader — bonded-atoms inner | `read_2cbo` | **DRIFT (IC-B)** |
| F | 2CBO reader — bond-angles list | `read_2cbo` | OK (outer-only pass-through) |
| G | Bond-removal helper | `_remove_bond` | consumer of IC-B |
| H | 3CBO reader — lattice vectors 4×4 | `read_3cbo` | OK |
| I | 3CBO reader — lattice indices 4×4 | `read_3cbo` | OK |
| J | 3CBO reader — atom coordinates 4×4 | `read_3cbo` | OK |
| K | 3CBO reader — displacement / bond_len | `read_3cbo` | OK |
| L | 3CBO reader — centroid positions | `read_3cbo` | **DRIFT (IC-C)** |
| M | 3CBO reader — fake-atom record | `read_3cbo` | OK |
| N | 2CBO writer | `write_2cbo` | consumer of IC-A / IC-B |
| O | Structure writer — header + atom rows | `add_to_struct` | OK |
| P | Structure writer — centroid rows | `add_to_struct` | consumer of IC-C |
| Q | Loop bounds (all) | whole file | OK |

## Finding IC-A — `rec['orbital_charges']` inner pair

**Perl**: `src/scripts/insert3cbo:225-230` (`sub read2CBO`)
```perl
foreach $orbital (1..$boData[$atom][8])
{
   @values = StructureControl::prepLine(\*BO2C,$line,'\s+');
   $boData[$atom][9][$orbital][1] = $values[0];
   $boData[$atom][9][$orbital][2] = $values[1];
}
```
Outer orbital index 1-indexed; inner 2-element slot layout
(`[1]` = label, `[2]` = charge) also 1-indexed with slot 0 unused.
Write-back at `insert3cbo:497-501` reads `[1]` / `[2]`.

**Python**: `insert3cbo.py:387-392`
```python
rec['orbital_charges'] = []
for _ in range(num_orb):
    values = prep_line(f)
    rec['orbital_charges'].append(
        [values[0], values[1]]
    )
```
Outer container is a 0-indexed `list`; inner 2-element pair is
0-indexed `[label, charge]`. Consumer at `write_2cbo:764-767` reads
`orb[0]`, `orb[1]`.

**Severity**: DRIFT (both outer and inner axes). Producer/consumer
self-consistent, so output is byte-identical.

**Fix**: store inner as `[None, label, charge]`; outer may remain a
plain Python list (matches the `bond_info_temp` precedent in BA-A).
Update the writer's `orb[0]`/`orb[1]` to `orb[1]`/`orb[2]`.

## Finding IC-B — `rec['bonded_atoms']` inner triple

**Perl**: `src/scripts/insert3cbo:233-239` (`sub read2CBO`)
```perl
foreach $bondedAtom (1..$boData[$atom][10])
{
   @values = StructureControl::prepLine(\*BO2C,$line,'\s+');
   $boData[$atom][11][$bondedAtom][1] = $values[0]; # atom ID
   $boData[$atom][11][$bondedAtom][2] = $values[1]; # length
   $boData[$atom][11][$bondedAtom][3] = $values[2]; # order
}
```
Outer bond index 1-indexed; inner 3-tuple slots `[1]` = atom ID,
`[2]` = bond length, `[3]` = bond order. Modification sites at
`insert3cbo:326, 337, 353, 441-449` splice/assign using `[1]`, `[2]`,
`[3]`; writer at `insert3cbo:503-507` reads the same slots.

**Python**: `insert3cbo.py:401-408`
```python
rec['bonded_atoms'] = []
for _ in range(num_bonded):
    values = prep_line(f)
    rec['bonded_atoms'].append([
        int(values[0]), float(values[1]), float(values[2]),
    ])
```
Outer `list` 0-indexed; inner 3-tuple stored at `[0]`, `[1]`, `[2]`.

Three consumers follow suit:
- `_remove_bond` at `insert3cbo.py:446-452` uses `enumerate` +
  `del atom_rec['bonded_atoms'][i]` and reads `bond[0]`.
- New-3C-bond append at `insert3cbo.py:656-668` stores
  `[next_3c_atom, bl, current_3cbo]`.
- Writer at `insert3cbo.py:774-779` reads `bond[0]`, `bond[1]`,
  `bond[2]`.

**Severity**: DRIFT (both axes). Fully self-consistent end to end.

**Fix**: store inner as `[None, atom_id, length, order]`; update
`_remove_bond` to read `bond[1]`; update the read_3cbo append to
`[None, next_3c_atom, bl, current_3cbo]`; update the writer to read
`[1]`, `[2]`, `[3]`. Outer axis may remain a plain Python list.

## Finding IC-C — `bo3c_positions[bo3c]` centroid coordinates

**Perl**: `src/scripts/insert3cbo:420-423` (`sub read3CBO`)
```perl
@values = StructureControl::prepLine(\*BO3C,$line,'\s+');
$bo3cPositions[$bo3c][1] = $values[1];
$bo3cPositions[$bo3c][2] = $values[2];
$bo3cPositions[$bo3c][3] = $values[3];
```
Outer `$bo3c` runs `1..$num3CBonds`; inner 3-element coordinate
tuple is 1-indexed with slot 0 unused. Writer at `insert3cbo:574-
575` and `603-604` reads `[1..3]`.

**Python**: `insert3cbo.py:632-637`
```python
values = prep_line(f)
bo3c_positions[bo3c] = [
    float(values[1]),
    float(values[2]),
    float(values[3]),
]
```
Outer container is a `dict` keyed by the 1-indexed `bo3c` value (OK
— dict keys faithfully carry 1-indexing), but the inner 3-element
list is stored 0-indexed at `[0]`, `[1]`, `[2]`.

Consumer at `add_to_struct` lines 860-869 and 891-900:
```python
for a in range(1, num_3c_bonds + 1):
    pos = bo3c_positions[a]
    fout.write(
        f" {atom_number}  {type_number}"
        f" {pos[0]} {pos[1]} {pos[2]} 3c\n"
    )
```

**Severity**: DRIFT (inner 3-vector only; outer dict keys are
already 1-indexed). Same class as Cluster A in structure_control,
re-introduced here.

**Fix**: store as `[None, x, y, z]`; update both writers in
`add_to_struct` to read `pos[1]`, `pos[2]`, `pos[3]`.

## Non-findings (verified)

- **`lattice_v` / `lat_idx` / `coords` / `disp` 4×4 scratch
  matrices** (`read_3cbo:521-526, 573-578, 584-598, 611-617`):
  allocated `[[0.0]*4 for _ in range(4)]`, accessed `[1..3][1..3]`,
  matching Perl with slot-0 sentinel rows and columns. OK.
- **`bond_len` 1D scratch** (`read_3cbo:621-627`): `[0.0]*4` seeded,
  filled at `[1..3]`, matching Perl `$bondLength[1..3]`. OK.
- **`bo_data` outer keying**: Python `dict` keyed by 1-indexed atom
  number carries the Perl `@boData[1..numAtoms]` layout faithfully.
- **`values[1..3]` reads throughout `read_3cbo`**: `prep_line` output
  is 0-indexed token list; `[0]` is the first token. The
  "index 1 = second token" pattern reflects the on-disk format
  (each line starts with a tag token). OK.
- **Loop bounds**: every Perl `foreach (1..N)` maps to Python
  `range(1, N+1)` — verified across `read_2cbo:349`, `read_3cbo:531,
  522, 574, 587, 589, 614, 622`, and `add_to_struct:847, 860, 886,
  891`.
- **`bond_angles` outer list** (`read_2cbo:417-420`,
  `write_2cbo:786-787`): opaque string pass-through; not indexed.

## rc companion (`insert3cborc.py`)

Five string defaults (file paths). No list or dict layouts cross a
1-indexed boundary. Joins the scalar-only rc class with
`pdb2sklrc.py`, `condenserc.py`, `make_reactionsrc.py`.

**Verdict**: OK. No rc-level fix needed.

## Running tally

- **BUG**: 0
- **DRIFT**: 3 (IC-A, IC-B, IC-C)
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: everything else

## Fix order

Leaf utility, no test coverage, no cross-script consumers, no
`sc.*` calls, program already FIX-flagged as broken at the pipeline
level. All three clusters can be fixed in a single edit pass.
Suggested order (smallest footprint first):

1. **IC-C** — 1 producer + 2 consumer sites.
2. **IC-B** — 4 touch points across `read_2cbo`, `_remove_bond`,
   `read_3cbo` append, `write_2cbo`.
3. **IC-A** — 2 touch points in `read_2cbo` and `write_2cbo`.

Verification: `ast.parse` + inspection. Live smoke-test is
impractical given the broken-pipeline status noted in the FIX
header.
