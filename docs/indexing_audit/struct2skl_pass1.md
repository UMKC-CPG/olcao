# struct2skl.py — Pass 1 Audit

**Status: zero findings.** Joins the "thin wrapper — no issues" tier
with `pdb2skl.py`, `make_reactions.py`, `skl2isaacs.py`,
`runIsoAtoms.py`.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | all subsystems |

## Scope

- Perl original: `src/scripts/struct2skl` (140 lines)
- Python port:   `src/scripts/struct2skl.py` (296 lines)
- rc companion:  `src/scripts/struct2sklrc.py` (47 lines)

## Cluster-by-cluster verdict

### S2-A — CLI / settings (`ScriptSettings`)

All five instance attributes (`input_file`, `output_file`,
`map_file`, `struct_types`, `make_mol`) are scalars (3 strings, 2
booleans). Matches Perl `setDefaults` / `parseCommandLine` at
`struct2skl:80-132`. No 1-indexed data crosses any interface.
**OK.**

### S2-B — rc companion

`struct2sklrc.py:16-42` returns a dict with 5 scalar keys, zero
list-valued defaults. Joins the scalar-only rc class. **OK.**

### S2-C — `sc.*` call boundary in `main()`

`struct2skl.py:263-285`. Three calls:

- `sc.read_input_file(settings.input_file,
  use_file_species=settings.struct_types)` — signature at
  `structure_control.py:1495`, two scalars.
- `sc.print_olcao(filename=..., title=..., style="cartesian")` —
  signature at `structure_control.py:2819`, three scalars.
- `sc.print_olcao_map(filename=...)` — signature at
  `structure_control.py:2936`, one scalar.

No per-atom / per-axis / per-element / per-species data crosses
the boundary. **OK.**

### S2-D — Local bookkeeping lists

None exist in the port. `main` creates the `sc` object and one
scalar `title` string. **OK.**

### S2-E — Loop bounds

Only one loop in the file (`struct2skl.py:236-238`,
`for argument in sys.argv`). Iterates `sys.argv` which is
legitimately 0-indexed in both Perl (`@ARGV`) and Python. **OK.**

## Non-indexing observations (out of scope)

- Perl help text at `struct2skl:21` claims default output is
  `skl.dat`, but Perl `setDefaults` at `struct2skl:84` actually
  sets `olcao.skl`. Python correctly uses `olcao.skl`. Doc drift
  on the Perl side only.
- `--makemol` is documented but unimplemented in both Perl and
  Python (deliberate carry-forward placeholder).

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: all five clusters

## Fix order

None required.
