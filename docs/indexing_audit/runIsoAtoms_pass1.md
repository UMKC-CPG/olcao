# runIsoAtoms.py — Pass 1 Audit

**Status: zero findings.** Joins `pdb2skl.py`, `make_reactions.py`,
and `skl2isaacs.py` in the "nothing to fix" tier.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | 6 clusters (RIA-A through RIA-F) |

## Scope

- Perl original: `src/scripts/runIsoAtoms` (219 lines,
  `sub makeAtomPot` @119-154, `sub makePDB` @156-218)
- Python port:   `src/scripts/runIsoAtoms.py` (457 lines)
- rc companion:  `src/scripts/runIsoAtomsrc.py` (41 lines)

The script runs isolated-atom SCF calcs for every element (or a
single target element), then harvests converged potentials into
`atomPDB/`. It imports **only** `element_data.ElementData` — no
`structure_control` calls at all.

## Summary table

| # | Cluster | Function(s) | Verdict |
|---|---------|-------------|---------|
| RIA-A | CLI / rc / settings | `ScriptSettings.*`, `parameters_and_defaults` | OK |
| RIA-B | Element range helper | `get_element_range` | OK |
| RIA-C | Isolated-atom SCF driver | `make_atom_pot` | OK |
| RIA-D | Potential DB builder | `make_pdb` | OK |
| RIA-E | Main driver | `main` | OK |
| RIA-F | rc companion | `runIsoAtomsrc.py` | OK |

## Key verification points

- **RIA-B** `get_element_range` at `runIsoAtoms.py:206-225` returns
  1-indexed inclusive bounds `(1, num_elements)` or
  `(target_z, target_z)`, matching Perl
  `$elementInit=1; $elementFin=$numElements` at
  `runIsoAtoms:127,170`. Callers walk
  `range(element_init, element_fin + 1)` — canonical 1-indexed
  idiom.
- **RIA-C** `make_atom_pot` at `runIsoAtoms.py:232-307`: loop @266
  walks `1..num_elements`; `elem_name = element_names[element]`
  @267 indexes directly with the 1-based counter, matching Perl
  `$elementNames_ref->[$element]` at `runIsoAtoms:134-135,139,144`.
  Docstring @252-254 explicitly documents `element_names` as
  1-indexed.
- **RIA-D** `make_pdb` at `runIsoAtoms.py:314-401`: loop @356 walks
  `1..num_elements`; `elem_name = element_names[element]` @357 and
  all path constructions use the live 1-based counter, matching
  Perl `runIsoAtoms:181-182,195-196,209`.
- **RIA-E** `main` @408-446 pulls `element_names = ed.element_names`
  unchanged; `element_data.py:62-64` docstring + `element_data.py:140`
  confirm `element_names` is 1-indexed with `[None]` sentinel. Clean
  boundary.
- **RIA-F** `runIsoAtomsrc.py:16-36` returns
  `{"target_z": 0, "nocore": False}` — **scalars only**, no list
  defaults that could cross a 1-indexed boundary.

## Non-findings (verified)

- No `range(3)`, `range(0, 3)`, `[atom - 1]`, `[axis - 1]`, or
  `[element - 1]` bridges anywhere in the file.
- No nested-list local bookkeeping — only non-scalar local is
  `remaining` (flat string buffer from the SCF file).
- `recordCLP` argv loop @197-198 iterates `sys.argv` including slot 0
  (program name); matches Perl `foreach $argument (0..$#ARGV)` @101.
- No `sc.*` calls at all.

## Non-indexing observations (out of scope)

- Perl `-help` is hand-rolled (`runIsoAtoms:82-83`); Python relies
  on argparse `-h/--help`. Functionally equivalent.
- `make_pdb` missing-file failures surface as `FileNotFoundError`
  from `shutil.copy2` rather than Perl's formatted `die`.
  Diagnostic-only, not indexing.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: all six clusters (RIA-A through RIA-F)

## Fix order

None required. `runIsoAtoms.py` is clean.
