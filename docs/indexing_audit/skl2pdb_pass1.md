# skl2pdb.py — Pass 1 Audit

**Status: zero findings.** Thin CLI wrapper with no indexing-
sensitive surface area of its own.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | five clusters (SP-A through SP-E) |

## Scope and method

- Perl original: `src/scripts/skl2pdb` (109 lines)
- Python port:   `src/scripts/skl2pdb.py` (~242 lines)
- rc companion:  `src/scripts/skl2pdbrc.py` (38 lines)

The job of the script is narrow: read an OLCAO skeleton (`.skl`)
input file, then write it back out as a Protein Data Bank (`.pdb`)
file with atoms reordered by element. All of the real work — parsing
the skeleton, sorting atoms by element, and producing the PDB text
— happens inside `StructureControl`. The script itself is a thin
CLI wrapper around two method calls.

## Summary table

| # | Topic                     | Symbols                         | Verdict |
|---|---------------------------|---------------------------------|---------|
| A | CLI / settings class      | `ScriptSettings.*`              | OK      |
| B | Main entry point          | `main()`                        | OK      |
| C | `sc.*` calls              | `read_input_file`, `print_pdb`  | OK      |
| D | Local bookkeeping / loops | (none)                          | OK      |
| E | rc companion defaults     | `parameters_and_defaults`       | OK      |

## Cluster SP-A — CLI / settings class (OK)

Perl `skl2pdb` lines 42-46 and the `parseCommandLine` sub (lines
76-108): declares four scalar globals (`$sklTypes`, `$pdbFile`,
`$sklFile`, `$mapFile`) and walks `@ARGV` to override them. Python
`skl2pdb.py:45-202` mirrors this with argparse plumbing. The four
instance variables `input_file`, `output_file`, `map_file`, and
`skl_types` are all scalars. **OK.**

## Cluster SP-B — Main entry point (OK)

Perl `skl2pdb:54-58`:
```perl
StructureControl::readInputFile($sklFile,$sklTypes);
StructureControl::printPDB($pdbFile,"molecule");
```
Python `skl2pdb.py:209`:
```python
sc = StructureControl()
sc.read_input_file(
    settings.input_file,
    use_file_species=settings.skl_types,
)
sc.print_pdb(filename=settings.output_file)
```
Both calls take only scalar arguments — filename strings and a
boolean flag. No per-atom data flows in either direction.

## Cluster SP-C — `sc.*` call boundaries (OK)

- `structure_control.py:1495`:
  `def read_input_file(self, filename=None, use_file_species=True)`
- `structure_control.py:3085`:
  `def print_pdb(self, filename=None)`

Both signatures take scalars only. The 1-indexed rule does not
apply at this boundary.

## Cluster SP-D — Local bookkeeping / loops (OK)

A complete read-through of `skl2pdb.py` finds:

- Zero `range(...)` loops over any atom / axis / element /
  species index space.
- Zero `.append([...])` sites building nested per-atom or
  per-vector data.
- Zero `[None]`-sentineled lists.
- The only `for` loop in the file walks `sys.argv` inside
  `recordCLP()` to echo the invocation into the `command`
  provenance file — legitimately 0-indexed argv.

## Cluster SP-E — rc companion defaults (OK)

`skl2pdbrc.py:16-33` returns four scalars:
- `"input_file": "olcao.skl"`
- `"output_file": "olcao.pdb"`
- `"map_file": "sklPDB.map"`
- `"skl_types": False`

All scalars. No list default crosses a 1-indexed `sc.*` call
boundary. Joins the scalar-only rc class.

## Non-indexing observations (out of scope)

- **Map file is accepted but never written.** Both the Perl
  original and the Python port parse the `-m` / `--map` option and
  store it in `$mapFile` / `self.map_file`, but neither actually
  produces the mapping file. The Perl source has a commented-out
  call to `printOLCAOMap` (line 61 of `skl2pdb`); the Python port
  drops the stub.
- **`print_pdb` mode argument dropped.** Perl calls
  `printPDB($pdbFile,"molecule")` with `"molecule"` as the mode /
  group name. The Python `print_pdb(filename=...)` no longer
  accepts that argument. Signature change, not indexing.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every cluster

## Fix order

None required.
