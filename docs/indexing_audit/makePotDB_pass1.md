# makePotDB.py — Pass 1 Audit

**Status: zero findings.** Clean port.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | all nine subsystems |

## Scope and method

- Perl original: `src/scripts/makePotDB` (353 lines)
- Python port:   `src/scripts/makePotDB.py` (545 lines)
- rc companion:  `src/scripts/makePotDBrc.py`

The program creates the default atomic potential database
(`atomicPDB`): for each requested element it creates a directory,
writes a `pot1` definition file, then runs `gaussFit` on the
numerical potential from `atomicBDB` to produce an initial `coeff1`.
Optional forking parallelises the fit stage.

**Critical scope observation**: the script does **not** import
`structure_control.py` at all. All per-element data comes from
`element_data.py`, whose arrays (`element_names`, `num_terms_pot`,
`min_term_pot`, `max_term_pot`) are 1-indexed with `[None]`
sentinels — allocated at `element_data.py:140,154-156` and
populated via `.append(...)` at lines 281, 404-413.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings plumbing | `ScriptSettings.*`, `recordCLP` | OK |
| B | Element-range helper | `get_element_range` | OK |
| C | Element-data lookup | `get_element_data` | OK |
| D | Database path lookup | `get_database_paths` | OK |
| E | Directory creation | `make_dirs` | OK |
| F | Potential-definition writer | `make_pot_definition` | OK |
| G | Coefficient fitter | `make_pot_coeffs`, `_fit_one_element` | OK |
| H | Main driver | `main` | OK |
| I | rc companion defaults | `makePotDBrc.py` | OK |

## Verification highlights

### MPDB-B — `get_element_range`

Perl `makePotDB:142-145,173-176,228-231` (three identical blocks):
```perl
if ($targetZ == 0)
   {$elementInit = 1; $elementFin = $numElements;}
else
   {$elementInit = $targetZ; $elementFin = $targetZ;}
foreach $element ($elementInit..$elementFin)
```
Python `makePotDB.py:194-211`:
```python
if target_z == 0:
    return range(1, num_elements + 1)
else:
    return range(target_z, target_z + 1)
```
Walks `1..N`, slot 0 skipped. Matches Perl.

### MPDB-F — `make_pot_definition`

Perl `makePotDB:180-211`: loop variable `$element` used both as a
directory index (`$elementNames_ref->[$element]`) and as the
literal nuclear-charge value (`printf POT "%f %f\n",$element,20.0;`).
The 1-indexed element number **is** the atomic Z.

Python `makePotDB.py:316-343`: same dual use —
`_ed.element_names[element]` for the directory name and
`pot.write(f"{element:f} {20.0:f}\n")` at line 337 for the nuclear
charge. `get_element_data(element)` called with 1-indexed value.

### MPDB-G — coefficient fitter

`coeff_data = [line.rstrip() for line in f]` at line 440 is a
0-indexed file-read buffer. The rewrite at lines 446-448 uses
`coeff_data[0]` for the header and `for line in coeff_data[1:]`
for body rows. This matches the Perl `$coeffData[0]` +
`@coeffData[1..$#coeffData]` split — both legitimately 0-indexed
(file slurps are never 1-indexed with a sentinel).

### MPDB-H — `main`

`num_elements = len(_ed.element_names) - 1` (line 513, with
`# 1-indexed` comment). Because `element_names` is allocated with
`[None]` slot-0 sentinel and populated via `.append`, list length
is `N + 1` and `len - 1 == N`, matching Perl
`ElementData::getNumElements()`.

### MPDB-I — rc companion

`makePotDBrc.py:16-28`:
```python
param_dict = {
    "element": 0,
    "fork": 1,
}
```
Two scalar ints. No list-valued settings. Matches the "scalars
only" rc class.

## Non-findings

- **`num_elements` derivation**: `len(_ed.element_names) - 1`
  correctly yields Perl's `getNumElements()` because of the slot-0
  sentinel.
- **`element` as printf argument**: loop variable serves as both
  a 1-indexed array key and the literal nuclear charge Z — no
  conversion needed.
- **`coeff_data[0]` + `coeff_data[1:]`**: legitimate 0-indexed
  file-read buffer matching Perl `<COEFF>` slurp.
- **No `structure_control` imports** confirmed via grep for
  `sc\.`, `StructureControl`, `structure_control` — zero matches.
- **`range()` bounds**: `get_element_range` is the sole source of
  element-iteration ranges and walks `1..N+1` or `[target_z]`.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every subsystem

## Fix order

None required.
