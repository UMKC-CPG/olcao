# uolcao.py — Pass 1 Audit

**Status: zero findings.** The second converted script (after
`pdb2skl.py`) to require no fix phase.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | every function |

## Scope

- Perl original: `src/scripts/uolcao`
- Python port:   `src/scripts/uolcao.py` (~2385 lines)
- rc companion:  `src/scripts/uolcaorc.py`

`uolcao.py` is a pure orchestration / file-staging driver for the
Fortran `uOLCAO`/`guOLCAO` executables. It parses CLI flags,
stages input files, shells out to Fortran, and renames `fort.*`
outputs.

**Does not call `structure_control` at all.** The single reference
is a comment next to a hard-coded `HARTREE_EV = 27.21138386` at
line 2132, replacing the Perl `StructureControl::getHartree()`.
The Perl original's only `StructureControl::` uses were `prepLine`
(7 call sites) and `getHartree` (1 site), all of which the Python
port correctly inlines as `line.split()` / scalar literal.

## Five Pass-1 coverage categories — all clean

1. **Public entry-point params/returns**: `main()` takes no args;
   `ScriptSettings.__init__` reads argparse + rc only. No atom /
   axis / element / species indices ever cross the API.
2. **`sc.*` calls**: none exist. Grep confirmed.
3. **Local bookkeeping lists**:
   - `alt_basis_list` — unindexed iteration, matches Perl
     `@altBasisList`.
   - `lines = f.readlines()` in `check_gamma_kp` — literal text-
     line offsets 5, 7 (not 1-indexed array slots; Perl reads the
     same physical lines sequentially).
   - `parts` / `values` from `line.split()` — 0-indexed on both
     sides, matching Perl `@values` from `prepLine` (README lists
     as a 0-indexed exception).
   - `sys.argv` iteration, path-component iteration — legitimate
     0-indexed.
4. **Loop bounds**: only two `range()` calls in the whole file,
   both in `update_pacs`: `range(4)` for four fixed header lines,
   and `range(num_core_orbitals)` whose loop variable is `_`
   (unused as an index). Perl uses `1..$numCoreOrbitals` but
   neither side indexes into an array with the loop variable.
5. **rc companion `uolcaorc.py`**: two scalar booleans
   (`serialxyz`, `valgrind`). No lists.

## Verified non-findings

- `EDGE_MAP[edge]` returns `(QN_n, QN_l)` 2-tuple → legitimately
  0-indexed (values are literal quantum numbers, matching the
  README's `absoluteSysQn[0..8]`-style exception).
- `JOB_DEFS` entries are positional 4-tuples `(job_id, job_name,
  def_scf, def_pscf)`, not indexed slot data.
- `check_gamma_kp` reads `lines[5]` / `lines[7]` — text-file line
  offsets, not 1-indexed arrays.
- `parts[-1]` total-energy extraction mirrors Perl
  `$values[$#values]`.

## Non-indexing observations (out of scope)

- `ScriptSettings.recordCLP` is the only non-snake_case method
  name — carry-over from Perl `recordCLP`. Cosmetic.
- `EDGE_MAP.get(fin_scf_edge, (0, 0))` is looked up twice inside
  the inner loop of `update_pacs` (once for `[0]`, once for `[1]`);
  could be hoisted.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every subsystem

## Fix order

None required.
