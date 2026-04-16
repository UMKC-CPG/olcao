# skl2isaacs.py — Pass 1 Audit

**Status: zero findings.** Joins `pdb2skl.py` and `make_reactions.py`
on the clean list.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | all subsystems |

## Scope

- Perl original: `src/scripts/skl2isaacs`
- Python port:   `src/scripts/skl2isaacs.py`
- rc companion:  `src/scripts/skl2isaacsrc.py`

## Findings

`skl2isaacs.py` is a thin CLI wrapper that only calls:

- `sc.read_input_file(filename, use_file_species=False)`
- `sc.print_isaacs(filename=output_root)`

Both take scalar arguments only — no per-atom, per-axis, per-
element, or per-species data crosses the `StructureControl`
boundary. No local indexed bookkeeping lists, no `range(1, N+1)`
loops over indexed containers. The only loop is a
`for argument in sys.argv` scan that mirrors Perl's 0-indexed argv
walk.

The rc companion `skl2isaacsrc.py` stores two scalar strings
(`input_file`, `output_root`), joining the "scalars only" rc
category.

No fix phase required.

## Running tally

- Methods / subsystems audited: full file
- Findings: 0
- Status: **CLEAN**
