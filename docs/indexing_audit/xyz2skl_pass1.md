# xyz2skl.py — Pass 1 audit

## Scope and method

`src/scripts/xyz2skl.py` (~292 lines) is a Python port of the Perl
`xyz2skl` script (~133 lines).  It is a thin CLI driver: read an XYZ
structure file through `StructureControl`, prepend a provenance note
to the title, write an OLCAO skeleton file in Cartesian form, and
write an atom-index map between the XYZ ordering and the resulting
`.skl` ordering.

Pass 1 approach mirrors the earlier `pdb2skl_pass1.md` /
`bond_analysis_pass1.md` sweeps:

- Read the Perl original (`src/scripts/xyz2skl`) end-to-end.
- Read the Python port (`src/scripts/xyz2skl.py`) end-to-end.
- Read the companion rc file (`src/scripts/xyz2sklrc.py`) end-to-end.
- Cross-check every `sc.*` call against `structure_control.py` for
  signature / indexing convention compatibility.
- Inspect `read_xyz` (`structure_control.py:1807`) for how the XYZ
  file populates per-atom 1-indexed arrays and the title state.

The Perl original has **no per-atom, per-axis, per-element, or
per-species data of its own** — every indexed array lives inside the
`StructureControl` module and is manipulated entirely through calls
into it.  Accordingly the port surface is equally narrow: one
`read_input_file`, one title read, one `print_olcao`, one
`print_olcao_map`.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings | `ScriptSettings.*`, rc loader | OK |
| B | Read dispatch | `sc.read_input_file` | OK |
| C | Title plumbing | `sc.system_title[1]` read vs. `getTitle` | OK (see note) |
| D | Skeleton output | `sc.print_olcao` | OK |
| E | Map output | `sc.print_olcao_map` | OK |
| F | rc file | `xyz2sklrc.parameters_and_defaults` | OK |

**Zero indexing findings.**  Everything crossing an `sc.*` boundary
is either a scalar (filename, bool, title string, style token) or a
1-indexed read from `sc.system_title` that matches the sentinel
layout in `structure_control.py`.  The port is strictly narrower than
`bond_analysis.py` or `condense.py` — there is simply no per-atom
nested data being marshalled at this level.

## Per-cluster detail

### Cluster X2-A — `ScriptSettings` / rc file (Pass 1 items 1 and 5)

**Perl** `xyz2skl:73-82, 84-125`:

```perl
sub setDefaults {
   $xyzFile  = "model.xyz";
   $sklFile  = "olcao.skl";
   $mapFile  = "sklXYZ.map";
   $xyzTypes = 0;
}
```

Four scalar defaults.  No per-atom/axis/element/species arrays.
Command-line parsing walks `@ARGV` by numeric index; no indexing
crosses into `StructureControl`.

**Python** `xyz2skl.py:59-222` and `xyz2sklrc.py:16-34`:

```python
def parameters_and_defaults():
    param_dict = {
        "input_file":  "model.xyz",
        "output_file": "olcao.skl",
        "map_file":    "sklXYZ.map",
        "xyz_types":   False,
    }
    return param_dict
```

Same four scalars.  The rc dictionary contains only strings and a
bool; no list-valued defaults and therefore no 1-indexed container
boundaries to preserve.  Matches the rc-audit pattern seen for
`pdb2sklrc.py`, `condenserc.py`, and `make_reactionsrc.py`.

**Severity**: **OK**.  No action needed.

### Cluster X2-B — `sc.read_input_file` call

**Perl** `xyz2skl:59`:

```perl
StructureControl::readInputFile($xyzFile, $xyzTypes);
```

Two scalars (filename string, Boolean flag).  No 1-indexed container
crosses the call boundary.

**Python** `xyz2skl.py:253-256`:

```python
sc.read_input_file(
    settings.input_file,
    use_file_species=settings.xyz_types,
)
```

Same two scalars.  `read_input_file` at `structure_control.py:1495`
dispatches on filename substring to `read_xyz`
(`structure_control.py:1807`), which populates the 1-indexed per-atom
arrays (`direct_xyz`, `atom_element_name`, etc.) internally with the
correct `[None]` sentinel layout.  No indexing data crosses the call
boundary.

**Severity**: **OK**.

### Cluster X2-C — Title plumbing (indexing OK, data-flow drift noted)

**Perl** `xyz2skl:62-63`:

```perl
$title = StructureControl::getTitle;
$title = "Generated from XYZ file: $title";
```

`getTitle` at `StructureControl.pm:239-240` returns the module-level
scalar `$title`, which `readXYZ` populates from the first token of
line 2 of the XYZ file (`StructureControl.pm:1772`:
`$title = $values[0];`).  Scalar in, scalar out — no indexed data.

**Python** `xyz2skl.py:261-264`:

```python
title = ""
if sc.system_title and len(sc.system_title) > 1:
    title = sc.system_title[1].strip()
title = f"Generated from XYZ file: {title}"
```

The Python port reads **`sc.system_title[1]`** — slot 1 of the
1-indexed title array (slot 0 holds the `[None]` sentinel;
meaningful data lives in slots 1..N).  The 1-indexed read is
*correct in form* and matches the `system_title = [None]`
initialisation at `structure_control.py:458`.

However, `read_xyz` at `structure_control.py:1807-1896` writes its
title to **`self.title`** (a separate scalar attribute, set at line
1852) and does **not** append to `self.system_title`.  After an XYZ
read, `self.system_title` is still the bare `[None]` sentinel list
from `_init_state`.  The `len(...) > 1` guard therefore evaluates
False and `title` stays as the empty string, so the generated
output title is always literally `"Generated from XYZ file: "`
with nothing after the colon.

- `xyz2skl.py:263` reads the right slot of the right kind of
  container, so it is **not** an indexing finding.
- The real issue is that `read_xyz` stores into `self.title`
  (scalar) while `xyz2skl.py` reads from `self.system_title`
  (1-indexed list) — a *data-flow* drift, not an *indexing* drift.

**Severity (indexing)**: **OK**.  `system_title[1]` is the right
way to read slot 1 of a 1-indexed list.

**Severity (data flow)**: out of audit scope.  Flag for a separate
follow-up: either teach `read_xyz` to also
`system_title.append(...)` the comment-line content (matching how
`read_olcao_skl` at `structure_control.py:2451` and `read_pdb` at
`structure_control.py:2159` populate `self.title`, and how the
OLCAO reader populates `self.system_title` around
`structure_control.py:1097`), or have `xyz2skl.py` fall back to
`sc.title` when `system_title` has no meaningful slot-1 entry.
Non-blocking and does not affect indexing compliance.

### Cluster X2-D — `sc.print_olcao` call

**Perl** `xyz2skl:66`:

```perl
StructureControl::printOLCAO(\*SKL, $title, "cartesian");
```

Filehandle, scalar title, scalar style keyword.  No indexed data.

**Python** `xyz2skl.py:270-274`:

```python
sc.print_olcao(
    filename=settings.output_file,
    title=title,
    style="cartesian",
)
```

`print_olcao` signature at `structure_control.py:2819` is
`(self, filename=None, title='', style='frac')` — three scalars.
The method internally reads every 1-indexed array it needs
(`system_title[1:]`, `direct_xyz[atom][1..3]`, element/species
arrays) using correct sentinel-aware slicing and loop bounds.
No indexed argument crosses the call boundary.

**Severity**: **OK**.

### Cluster X2-E — `sc.print_olcao_map` call

**Perl** `xyz2skl:69`:

```perl
StructureControl::printOLCAOMap(\*MAP);
```

Filehandle only.

**Python** `xyz2skl.py:281`:

```python
sc.print_olcao_map(filename=settings.map_file)
```

`print_olcao_map` signature at `structure_control.py:2936` is
`(self, filename=None)` — a single scalar.  The method walks
`range(1, self.num_atoms + 1)` internally against the 1-indexed
per-atom arrays (`molecule_seq_num`, `residue_seq_num`,
`residue_name`, `atom_element_name`, `atom_tag`).  No indexed
argument at the boundary.

**Severity**: **OK**.

## Non-findings (verified consistent)

- **No local bookkeeping lists.**  Unlike `bond_analysis.py` (which
  builds `bond_info_temp`, `histogram`, `reduced_bond_info`, etc.)
  or `condense.py` (which builds `BondData.hooke_bond_coeffs`,
  `AngleData.hooke_angle_coeffs`, etc.), `xyz2skl.py` owns *no*
  per-atom, per-axis, per-element, or per-species arrays of its
  own.  Every list that could have an indexing convention lives
  inside `StructureControl`.

- **No range() loops over indexed data.**  The Python port has
  zero `range(...)` expressions touching atom/axis/element/species
  data.  The only loop is a pure `for argument in sys.argv`
  iteration in `recordCLP` for provenance, which is 0-indexed in
  both Perl and Python and does not cross into
  `StructureControl`.

- **No slot-0 sentinel hazards.**  Nothing is constructed with a
  `[None]` prefix because nothing 1-indexed is constructed.

- **rc list defaults.**  `xyz2sklrc.py:16-34` returns a dict
  containing three strings and one bool.  No list defaults of any
  kind — no indexing boundary to preserve.  Matches the clean
  profile of `pdb2sklrc.py`, `condenserc.py`, and
  `make_reactionsrc.py`.

- **`scalar(@ARGV)` / `foreach $argument (0..$#ARGV)` in Perl**
  (`xyz2skl:91, 121`).  These walk `@ARGV` which is 0-indexed in
  both Perl and Python; no 1-indexed convention applies.  The
  Python port's `sys.argv` loop at `xyz2skl.py:220` is the natural
  equivalent.

## Non-indexing observations (out of audit scope)

- **Title data-flow drift** (described under Cluster X2-C above):
  `read_xyz` writes to `self.title` while `xyz2skl.py` reads from
  `self.system_title[1]`.  The output is always
  `"Generated from XYZ file: "` with an empty trailing field.
  Not an indexing issue; flag for a separate fix alongside the
  `read_xyz` cleanup.

- **`print_olcao_map` default filename drift.**  Perl opens the
  map file at `$mapFile` (caller-supplied, default `sklXYZ.map`
  at `xyz2skl:78`) and writes into it directly.  The Python
  `print_olcao_map` default at `structure_control.py:2974` is
  `'olcao.map'` — different from `sklXYZ.map` — but `xyz2skl.py`
  always passes `filename=settings.map_file`, so the default
  never takes effect.  Harmless; noted only for completeness.

- **Documentation policy compliance.**  The Python port preserves
  the full Perl help text as a module docstring and replicates
  the option explanations in argparse help strings, matching the
  project's documentation-preservation rule.

## Running tally

- **BUG**: 0
- **DRIFT**: 0 (indexing)
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every subsystem.

One non-indexing data-flow drift flagged for a separate follow-up
(`read_xyz` → `self.title` vs. `xyz2skl.py` → `self.system_title
[1]`).

## Fix-order recommendation

No fix phase required for indexing compliance.  The script is a
thin CLI→`StructureControl` wrapper with no internal indexed data
structures; Pass 1 is complete with zero findings.

If the title data-flow drift is addressed later, the minimal fix
is one of:

1. Teach `read_xyz` to `self.system_title.append(...)` the
   comment-line content alongside the existing `self.title`
   assignment, matching `read_olcao_skl`'s pattern.  Downstream
   callers that already read `system_title[1]` then work
   uniformly.
2. Alternatively, in `xyz2skl.py`, fall back to `sc.title` when
   `sc.system_title` has no meaningful slot-1 entry.

Option 1 is preferred because it normalises the title state
across all input formats.

## Fix status

Not applicable — no indexing findings to fix.
