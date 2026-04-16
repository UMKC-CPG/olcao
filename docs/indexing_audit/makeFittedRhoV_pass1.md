# makeFittedRhoV.py — Pass 1 Audit

**Status: zero findings.** Among the cleanest Pass 1 audit targets
so far.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | all six clusters (A-F) |

## Scope

`makeFittedRhoV.py` (~826 lines) ports Perl `makeFittedRhoV`
(~406 lines). The tool builds a six-column input deck for
`OLCAOrhoV` Fortran90 by reading `structure.dat` for type→element
mapping, subtracting isolated-atom Gaussian coefficients from SCF
coefficients, and invoking the Fortran driver.

**No `structure_control.py` calls exist.** The only Perl
`StructureControl` reference is to the `prepLine` tokenizer, which
the port replaces with a local `prep_line` helper at
`makeFittedRhoV.py:101`.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings | `ScriptSettings.*` | OK |
| B | Structure reader | `read_struct` | OK |
| C | SCF vs iso difference | `compute_diff` | OK |
| D | Temp-file writer | `print_temp_data` | OK |
| E | Fortran driver | `run_olcao_rho_v` | OK |
| F | rc companion | `makeFittedRhoVrc.parameters_and_defaults` | OK |

## Non-findings (verified consistent)

### `num_mesh_points` — 0-indexed, matches Perl

Perl `makeFittedRhoV:141` initializes `@numMeshPoints = (0,0,0)`;
writes at `:160-162` use `$numMeshPoints[0..2]`; validation loop
at `:184-197` uses `foreach $number (1..3)` but explicitly
subtracts 1 on every read `$numMeshPoints[$number-1]`; final
command-line consumption at `:398-400` interpolates all three
0-indexed slots. The Perl array is **genuinely 0-indexed** — no
sentinel, no `[1..3]` access.

Python `makeFittedRhoV.py:197-199, 337` stores `list(...)`
directly; validates via `enumerate(self.num_mesh_points)` at
`:348`; consumes via `mesh[0]`, `mesh[1]`, `mesh[2]` at
`:744-746`. Cleaner than the Perl `[$number-1]` dance but
semantically identical. rc default `[0, 0, 0]` at
`makeFittedRhoVrc.py:90` is 0-indexed to match. OK — this is the
correct preservation of a legitimate 0-indexed Perl structure
(like the `absoluteSysQn[0..8]` case noted in README §1.7). Unlike
`makeinputrc.py`/`mod_structrc.py`, no fix is needed because the
Perl side is 0-indexed too.

### `type_elements` — 1-indexed, matches Perl

Perl `makeFittedRhoV:243-247` pre-increments `$numTypes` before
assigning `$typeElements[$numTypes]`, so slot 0 is unused. Consumer
`foreach $type (1..$numTypes)` at `:295, 298` is 1-indexed.

Python `makeFittedRhoV.py:439-440` initializes
`type_elements = [None]`; `:454-456` does `num_types += 1` then
`type_elements.append(...)`, landing first payload at index 1.
Consumer at `:567-571` uses `for type_idx in range(1, num_types +
1)` with `type_elements[type_idx]`. Fully 1-indexed, matching Perl
exactly.

### `out_values` outer 1-indexed, inner 0-indexed — matches Perl

Perl `makeFittedRhoV:346-351` writes `$outValues[$numTerms][0..5]`
after pre-incrementing `$numTerms` at `:317`, producing an
**asymmetric** structure: outer 1-indexed (slot 0 unused), inner
6-element `[0..5]` genuinely 0-indexed. Writer at `:370-374`
confirms: `foreach $term (1..$numTerms)` outer, `foreach $value
(0..5)` inner.

Python `makeFittedRhoV.py:532-533` initializes
`out_values = [None]`; `:642-649` appends a flat 6-element inner
list. Writer at `:694-698` walks outer `range(1, num_terms + 1)`,
inner via `for value in out_values[term]`. Both dimensions
preserved.

### `vals_np` / `vals_ip` / token lists — 0-indexed, matches Perl

Perl `prepLine` is a whitespace tokenizer returning 0-indexed
token arrays (not 1-indexed data structures); Perl reads tokens
via `$valuesIP[0..4]` etc. Python local `prep_line` at `:101-128`
returns `line.strip().split()` (0-indexed list); consumers read
`vals_ip[0..4]` at `:619, 631-632, 642-649`. The column layout of
the SCF potential file is 0-indexed on both sides.

### Loop bounds

All eight loops in the file verified. Outer element / type /
term loops walk `range(1, N+1)` matching Perl `(1..$n)`; mesh /
skip / coefficient loops walk 0-indexed containers where Perl
also walks 0-indexed.

### rc companion (`makeFittedRhoVrc.py`)

Six scalars plus one list (`num_mesh_points = [0, 0, 0]`). The
list is 0-indexed on both sides of the `ScriptSettings` boundary,
matching Perl's 0-indexed `@numMeshPoints = (0,0,0)`. Unlike
`makeinputrc.py` and `mod_structrc.py`, no fix is required here.

## Non-indexing observations (out of scope)

- Python promotes `num_cells` and `do_open_dx` to rc-configurable
  settings (Perl hard-coded `1 0`), adding `-nc` / `-noDX` CLI
  switches. Behavior enhancement.
- `-noDX` inverts the True default via
  `if args.no_dx: self.do_open_dx = False`.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every subsystem

## Fix order

None required.
