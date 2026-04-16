# skl2vasp.py -- Pass 1 audit

## Scope and method

`src/scripts/skl2vasp.py` (~574 lines) is a Python port of the
Perl `skl2vasp` script (~346 lines).  Both programs take an
OLCAO skeleton file, reorder its atom list so that atoms of
the same element are grouped together (a VASP POSCAR
requirement), parse the space-group designation to determine
a Bravais lattice cell type, and then hand the sorted file
off to `StructureControl` / `structure_control.py` to emit
VASP input files (POSCAR only in the Python port; POTCAR,
KPOINTS, and INCAR remain unported).

The file is structurally a thin CLI-plus-preprocessing
wrapper around two `StructureControl` calls.  It has no
per-atom, per-axis, per-element, or per-species array
bookkeeping of its own -- every atom line is handled as a
single opaque string during the pre-sort stage, and all
structural computation is delegated to `sc.read_input_file`
and `sc.print_vasp`.

Pass 1 approach:

- Read the Perl original end-to-end (all 346 lines).
- Read the Python port end-to-end (all 574 lines).
- Read `skl2vasprc.py` end-to-end.
- Grep sweep in `skl2vasp.py` for drift signals:
  `range(`, `[atom`, `[axis`, `[element`, `[species`,
  `sc.`, `.append([`.
- Cross-check signatures of `sc.read_input_file` and
  `sc.print_vasp` against the call sites in `main()`.
- Compare loop bounds in `pre_sort_skl` against the Perl
  `preSortSkl` subroutine.

## Summary table

| # | Topic                            | Functions                              | Verdict |
|---|----------------------------------|----------------------------------------|---------|
| A | CLI / settings                   | `ScriptSettings.*`                     | OK      |
| B | rc defaults                      | `skl2vasprc.parameters_and_defaults`   | OK      |
| C | Pot sub-type resolver            | `resolve_sub_pot_type`, `SUB_POT_MAP`  | OK      |
| D | Space-group -> cell-type helper  | `get_cell_type_from_sg`                | OK      |
| E | Atom-line pre-sort / header pass | `pre_sort_skl`                         | OK      |
| F | Main orchestration               | `main`                                 | OK      |

**Zero findings.**  No BUGs, no DRIFTs, no REVERSE DRIFTs,
no MIXED clusters.  The script is well-behaved with respect
to the 1-indexing convention because it simply does not
operate on any 1-indexed data of its own -- the only numeric
loop walks line strings, not atoms-as-data.

## Pass 1 coverage

### 1. Public entry point params / returns

The only public entry point is `main()`, which takes no
arguments and returns nothing.  `ScriptSettings.__init__`
is the other "top level" constructor; it also takes no
atom / axis / element / species arguments.  All seven
stored settings (`job_type`, `pot_type`, `sub_pot_type`,
`gw`, `gamma_kpoint`, `cell_type`, `input_file`) are pure
scalars (strings and booleans).  No indexing concerns.
**OK.**

### 2. `sc.*` calls -- do they pass 1-indexed data?

Two `StructureControl` call sites, both in `main()`:

- **`skl2vasp.py:554`** -- `sc.read_input_file(sorted_skl,
  use_file_species=False)`.  Perl original at
  `skl2vasp:100` calls
  `StructureControl::readInputFile($olcaoSklSorted,0)`.
  Both arguments are scalars (string filename + boolean
  flag); Python's `use_file_species=False` matches the
  Perl `0` literal.  **OK.**

- **`skl2vasp.py:559`** -- `sc.print_vasp(filename='POSCAR')`.
  The Perl original at `skl2vasp:103` calls
  `StructureControl::printVASP($cellType, $potType,
  $subPotType, $GKPoint, $jobType)` -- five scalar
  arguments used to drive POSCAR / POTCAR / KPOINTS /
  INCAR generation.  The Python port of `sc.print_vasp`
  only emits POSCAR (documented in the `print_vasp`
  docstring at `structure_control.py:2988`) and takes
  only the output filename.  The dropped arguments are a
  *feature-scope* gap (POTCAR / KPOINTS / INCAR not yet
  ported), not an indexing concern.  All arguments that
  *are* passed are still scalars.  **OK.**

No 1-indexed container ever crosses the `sc.*` call
boundary from this script.

### 3. Local bookkeeping lists -- slot-0 sentinels?

The script creates exactly one list in `pre_sort_skl`:

```python
atom_lines = []
for _ in range(num_atoms):
    atom_lines.append(fin.readline())
...
sorted_lines = sorted(atom_lines, key=str.casefold)
for line in sorted_lines:
    fout.write(line)
```

Compare to the Perl original at `skl2vasp:233-242`:

```perl
foreach $atom (0..$numAtoms-1)
   {$atomLines[$atom] = <SKLIN>;}
@sortedAtomLines = sort {lc($a) cmp lc($b)} @atomLines;
foreach $atom (0..$numAtoms-1)
   {print SKLOUT "$sortedAtomLines[$atom]";}
```

Perl iterates `0..$numAtoms-1` and assigns to
`$atomLines[$atom]`, i.e. **the Perl script itself uses a
0-indexed list here** -- this is one of the rare
legitimate 0-indexed cases in the Perl codebase (a line-
buffer scratch buffer used only for a sort, with no slot-0
sentinel).  The Python `atom_lines.append(...)` produces
the same 0-indexed layout.  Both ends match.  **OK.**

`SUB_POT_MAP` in `skl2vasp.py:88` is a dict keyed by
user-supplied abbreviation strings -- not an index-keyed
container, so the 1-indexing convention does not apply.
**OK.**

### 4. Loop bounds -- 1..N+1 for 1-indexed containers?

Only one numeric loop in the entire file:
`skl2vasp.py:477` -- `for _ in range(num_atoms):`.  This
walks the number of skeleton atom *lines* being consumed
from the input file; the loop variable is explicitly
discarded (`_`).  The body only calls `fin.readline()`
and appends the result to a 0-indexed scratch list.
Perl's equivalent at `skl2vasp:234` also uses a 0-based
range (`0..$numAtoms-1`).  **OK.**

No other `for ... in range(...)` or `while` loops exist
in the file.

### 5. rc companion -- list defaults crossing 1-indexed boundary?

`skl2vasprc.py:16-52` -- every default value in
`parameters_and_defaults` is a pure scalar:

| key            | type | value           |
|----------------|------|-----------------|
| `job_type`     | str  | `"relaxfull"`   |
| `pot_type`     | str  | `"potpawPBE"`   |
| `sub_pot_type` | str  | `""`            |
| `gw`           | bool | `False`         |
| `gamma_kpoint` | bool | `False`         |
| `cell_type`    | str  | `""`            |
| `input_file`   | str  | `"olcao.skl"`   |

No list / tuple / dict defaults.  No axis-, atom-,
element-, or species-indexed data.  **OK** -- matches
the `pdb2sklrc.py` / `condenserc.py` /
`make_reactionsrc.py` shape documented in the
structure_control audit's rc-file summary.

## Non-findings (verified consistent)

- **`recordCLP` command-file writer** at
  `skl2vasp.py:328` iterates `for argument in sys.argv`
  -- a 0-indexed iteration over `sys.argv` strings,
  matching Perl's `foreach $argument (0..$#ARGV)` at
  `skl2vasp:188`.  No indexing concerns.

- **`ScriptSettings.reconcile`** at `skl2vasp.py:318`
  does plain scalar assignment from the argparse
  namespace onto the instance.  No list plumbing.

- **`get_cell_type_from_sg`** at `skl2vasp.py:378` uses
  `int(f.readline().split()[0])` -- a 0-indexed access
  on the split-words result of a single file line.  The
  Perl equivalent at `skl2vasp:279-281` calls `prepLine`
  on the file handle twice, then reads `$values[0]`.
  Both read the first whitespace-delimited token of the
  second line of the space-group file (the numeric
  group number).  **OK.**

- **`pre_sort_skl` header parse** at
  `skl2vasp.py:471-472` -- `values = line.split();
  num_atoms = int(values[1])`.  The Perl equivalent at
  `skl2vasp:227-228` uses
  `prepLine("", $line, '\s+')` followed by
  `$numAtoms = $values[1]`.  The `prepLine("", $line,
  ...)` form of `prepLine` splits the passed-in string
  directly without prepending a sentinel (the sentinel
  convention only applies to certain file-handle-mode
  callers; confirmed by reading the Perl `prepLine` in
  `StructureControl.pm`).  So both Perl `$values[1]` and
  Python `values[1]` read the second whitespace-delimited
  token of the `frac N` / `cart N` header line, which is
  the atom count.  **OK.**  (Flagged here because the
  slot-convention interaction with `prepLine` is not
  obvious on first reading.)

- **`assign_rc_defaults`** at `skl2vasp.py:151` assigns
  seven scalar values from the rc dict to instance
  attributes.  No lists.

- **`add_parser_arguments`** at `skl2vasp.py:207` defines
  argparse options for seven scalar settings.  No
  indexing concerns.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every subsystem (A through F).

## Fix-order recommendation

**None required.**  No findings.

## Notes on scope gap (not an audit finding)

The Python port of `print_vasp` only writes POSCAR; the
Perl original also wrote POTCAR, KPOINTS, and INCAR.  This
is a documented intentional gap in
`structure_control.py:2988` (`print_vasp` docstring)
because POTCAR generation depends on an external
pseudopotential database tree referenced via environment
variables (`$VASPPOT_DIR`, `$VASPPOT_PAW_PBE`, etc.) and
shell-level decompression.  The five settings parsed into
`ScriptSettings` that are only consumed by the unported
files (`pot_type`, `sub_pot_type` after
`resolve_sub_pot_type`, `gw`, `gamma_kpoint`, `job_type`)
are stored but not currently used.  This does not affect
the 1-indexing audit but is flagged here for context: when
the POTCAR / KPOINTS / INCAR writers are eventually
ported, a Pass 2 audit should re-check any per-element
iteration inside those writers.

## Fix status

No fixes needed.  No tests exercise `skl2vasp.py` directly
(confirmed by absence of `test_skl2vasp*` files under
`tests/`).  Regression verification post-Pass-1 is by
inspection only; no edits were made.
