# vasp2skl.py — Pass 1 Audit

**Status: zero findings.**

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | six clusters (V2-A through V2-F) |

## Scope

- Perl original: `src/scripts/vasp2skl`
- Python port:   `src/scripts/vasp2skl.py`
- rc companion:  `src/scripts/vasp2sklrc.py`

## Cluster verdicts

- **V2-A** `ScriptSettings` (CLI/rc) — OK. Scalar strings only;
  `recordCLP` walks `sys.argv` 0-indexed matching Perl's
  `foreach (0..$#ARGV)`.
- **V2-B** `read_element_list` (vasp2skl.py:227–267) — OK.
  `element_list = [None]` seeded 1-indexed, `.append()` walks from
  slot 1. Matches Perl `$elementList[$numElements]` layout. The
  `tokens[3]` vs Perl `$values[4]` difference is a Perl-vs-Python
  whitespace-split quirk (Perl retains empty slot 0 on leading
  whitespace; Python's `str.split()` drops it), not an indexing
  drift — both resolve to the same real 4th token.
- **V2-C** `read_vasp` / `_read_scaled_vector`
  (vasp2skl.py:274–443) — OK. `atom_counts` is a legitimate
  0-indexed parse buffer (Perl `@atomCounts` is also 0-indexed;
  Perl comment "Note the -1. It starts at 0." preserved verbatim
  at line 370). `atom_element`, `atom_a`, `atom_b`, `atom_c` all
  seeded with `[None]` and built by `.append()` so they are
  1-indexed. Outer loop `range(1, num_elements+1)` at line 369
  mirrors Perl `for ($i=1;$i<=$numElements;$i++)`. `line_values =
  line.split()` at line 397 is a 0-indexed parse buffer at the
  "Direct"/selective-dynamics detector — matches Perl
  `@lineValues`.
- **V2-D** `print_skl` (vasp2skl.py:450–513) — OK. Loop
  `range(1, data['num_atoms']+1)` at line 502, reads
  `data['atom_element'][i]` / `atom_a[i]` / `atom_b[i]` /
  `atom_c[i]` with 1-indexed `i`. Matches Perl
  `for ($i=1;$i<=$numAtoms;$i++)` at line 271.
- **V2-E** `main` — OK. Four function calls pass scalars or the
  already-1-indexed data dict; no indexed containers.
- **V2-F** `vasp2sklrc.py` — OK. Three scalar strings
  (`input_file`, `atom_file`, `output_file`). Joins the
  "scalar-only rc" class.

## Key non-findings

- **No `sc.*` calls**: the script never touches
  `structure_control.py` (the Perl original imports
  `StructureControl` but calls nothing from it).
- **No `element_data.*` calls**: element names come directly from
  POTCAR TITEL lines.
- Both 1-indexed loops (`range(1, num_elements+1)` and
  `range(1, num_atoms+1)`) match their Perl `1..$n` counterparts
  exactly.
- `atom_counts[i-1]` bridge inside the 1-indexed element loop
  exactly replicates Perl's mixed-indexing bridge at the counts/
  element boundary, comment included.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every subsystem

## Fix order

None required.
