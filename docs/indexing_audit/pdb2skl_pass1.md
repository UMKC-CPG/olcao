# pdb2skl.py — Pass 1 audit

## Scope and method

`src/scripts/pdb2skl.py` (~274 lines) is a Python port of the Perl
`pdb2skl` script (~125 lines).  Both are thin wrappers that parse
three file-name arguments plus a buffer and a `--pdbtypes` flag,
load the named PDB into a `StructureControl` instance, and emit an
OLCAO skeleton (.skl) file plus a PDB↔SKL atom-index map file.

Because the script is so short and contains no per-atom, per-axis,
per-element, or per-species array accesses of its own, the audit
reduces to two questions:

1. Does the CLI parser preserve any Perl-side indexing conventions?
2. Are the three `StructureControl` calls (`read_input_file`,
   `print_olcao`, `print_olcao_map`) used with correct argument
   conventions after the Cluster A–I fixes in that module?

## Findings

**Zero indexing findings.**

- No `range(N)` loops, no `[0]` / `[1]` / `[2]` / `[3]` subscripts,
  no list-slicing of coordinate data, no lattice-vector math, no
  per-atom loop.  A `Grep` sweep for `range(`, `[1]`, `[2]`, `[3]`,
  `[0]`, and `[None]` across the entire file returns no matches.
- `ScriptSettings` stores scalar values only (file names, a float
  buffer, a bool flag).  No list or tuple layout choices to audit.
- The three `StructureControl` method calls in `main` each receive
  only scalar arguments:

  ```python
  sc.read_input_file(
      settings.input_file,
      use_file_species=settings.pdb_types,
  )
  sc.print_olcao(
      filename=settings.output_file,
      title="Generated from PDB file.",
      style="cartesian",
  )
  sc.print_olcao_map(filename=settings.map_file)
  ```

  All three methods were audited clean in Section 1.10 of the
  `structure_control.py` audit and require no coordinate vectors,
  slot-sentinel lists, or axis indices from the caller.

## Consistency / non-indexing observations

Not findings, but worth flagging for the separate quality-pass
backlog:

- **`--buffer` / `-b` is accepted but never applied.**  Both the
  Perl original and the Python port parse `-b <value>` into a
  script-local variable (`$buffer` / `settings.buffer`) and never
  pass it through to `StructureControl.set_buffer`.  Since
  `set_buffer` exists on the Python `StructureControl` and is
  required for non-crystalline PDBs to get a sensible bounding
  box (via `compute_crystal_parameters`), the value silently
  defaults to zero.  This is a pre-existing Perl bug faithfully
  preserved in the port; **not an indexing issue**.
- **Docstring drift** between the script's top module docstring
  and `parse_command_line`'s `description_text` — the top text
  mentions `--help` while the argparse parser relies on the
  default help handler.  Cosmetic only.

Both observations are outside the indexing-audit scope.  If you
want them tracked, add them to a separate "pdb2skl follow-up"
TODO item.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: everything

No fix phase is needed for `pdb2skl.py`.
