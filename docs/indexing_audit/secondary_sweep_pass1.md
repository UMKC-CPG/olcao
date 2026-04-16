# Secondary Sweep Pass 1 — structure_control intersections in
# non-ported Python files

## Scope

This sweep covers Python source files in `src/scripts/` that have
**no direct Perl ancestor** (i.e. they are new Python code, not a
port of a `StructureControl.pm`-era Perl script) but that might
still touch the `structure_control.py` API.  The goal is narrow:
for each file we only inspect whether and how it crosses the
1-indexed boundary enforced by `structure_control.py`, not the
full internal indexing bookkeeping of the file itself.

Reference rules: `docs/indexing_audit/README.md`.
Enforcement reference: `src/scripts/structure_control.py`.
Format reference: `docs/indexing_audit/bond_analysis_pass1.md`.

Audit approach:

1. `Grep` every target file for `structure_control`,
   `StructureControl`, and the common alias pattern `\bsc\.`.
2. `Grep` the import block of each file to confirm there is no
   indirect import of `structure_control` through a wrapper
   module.
3. Spot-check any file whose imports route through another local
   module (`olcaoObject`, `olcaoStructures`, `org_radwavefn`,
   `osrecurintglib`) to verify that wrapper does not itself
   import `structure_control`.

## Top-of-file summary

| File | Imports sc? | Uses `sc.*`? | Findings |
|---|---|---|---|
| `ab2pf.py` | No | No | none |
| `graspElems.py` | No | No | none |
| `hdf5Comp.py` | No | No | none |
| `interFit.py` | No | No | none |
| `lmp2skl.py` | No | No | none |
| `make_veusz_graph.py` | No | No | none |
| `nanoTube.py` | No | No | none |
| `olcaoObject.py` | No | No | none |
| `olcaoStructures.py` | No | No | none |
| `org_radwavefn.py` | No | No | none |
| `osrecurintg.py` | No | No | none |
| `osrecurintgana.py` | No | No | none |
| `osrecurintglib.py` | No | No | none |
| `osrecurintg_makenum.py` | No | No | none |
| `osrecurintgnum.py` | No | No | none |
| `pot2plot.py` | No | No | none |
| `viewCell.py` | No | No | none |
| `element_pct.py` | No | No | none |
| `dead-md.py` | No | No | none |
| `dpmts.py` | No | No | none |
| `XYZ.py` | No | No | none |
| `XYZrc.py` | No | No | none |
| `olcao.py` | No | No | none (Py2 leftover, see note) |

**No file in the secondary sweep imports `structure_control` or
instantiates `StructureControl`.**  Consequently there are no
Pass-1, Pass-2, Pass-3, or Pass-4 `sc.*` call-boundary findings
to record for any of these files.

### Finding tally by severity

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK (no sc intersection) | 23 |

No critical bugs surfaced in this sweep.  Every file in scope is
either a self-contained computational tool, a plotting / I/O
utility, or a wrapper around an unrelated subsystem
(`element_data`, `org_radwavefn`, `osrecurintglib`, HDF5, vedo,
matplotlib, etc.).

## Per-file notes

### `ab2pf.py`
Imports: `argparse`, `gzip`, `os`, `sys`, `subprocess`.  Zero
references to `structure_control` or `StructureControl`.  Zero
uses of the `sc.` alias pattern.  **No `sc.*` usage — nothing
to audit.**

### `graspElems.py`
Imports: `argparse`, `os`, `sys`, `subprocess`, `re`, `datetime`,
and `element_data.ElementData` only.  **No `sc.*` usage.**

### `hdf5Comp.py`
Imports: `argparse`, `os`, `sys`, `math`, `h5py`, `numpy`,
`matplotlib.pyplot`, `datetime`.  Pure HDF5 comparison tool.
**No `sc.*` usage.**

### `interFit.py`
Imports: `argparse`, `numpy`, `scipy`, `time`, `copy`, `os`,
`sys`, `math`, `pandas`, `datetime`, `org_radwavefn`,
`element_data.ElementData`, and `scipy.linalg.lu_solve`.  It
imports `org_radwavefn` for its symbol tables
(`lktol`, `sktol`, `ltosk`, `ltok`, `ltolsym`) — and
`org_radwavefn.py` itself does **not** import
`structure_control` (verified by direct search).  **No `sc.*`
usage.**

### `lmp2skl.py`
Imports: `os` only.  Stub / trivial file.  **No `sc.*` usage.**

### `make_veusz_graph.py`
Imports: `time`, `os`, `pandas`, `org_radwavefn`.  Plotting
helper.  `org_radwavefn` does not touch `structure_control`.
**No `sc.*` usage.**

### `nanoTube.py`
Imports: `os`, `sys`, `numpy`, `time`, `copy`.  Self-contained
nanotube geometry generator that emits coordinates directly
(does not route through `structure_control.py`).  **No `sc.*`
usage.**

### `olcaoObject.py`
Imports: `os`, `sys`, `time`, and `from olcaoStructures import
*`.  The star-import brings symbols from `olcaoStructures.py`,
which itself does not import `structure_control` (verified).
`olcaoObject.py` defines an `olcaoJob` class that drives job
I/O but does not cross the `sc.*` boundary.  **No `sc.*`
usage.**

### `olcaoStructures.py`
Imports: `os`, `sys` only.  Self-contained job / directory
bookkeeping layer used by `olcaoObject.py`.  **No `sc.*`
usage.**

### `org_radwavefn.py`
Imports: `time`, `numpy`, `math`, `copy`, `pandas`,
`matplotlib.pyplot`.  Radial-wavefunction organiser; operates
on orbital-indexed tables, not on atom/axis/element data from
`structure_control`.  **No `sc.*` usage.**

### `osrecurintg.py`
Imports: `argparse`, `os`, `sys`, `datetime`, `osrecurintgana`,
`osrecurintgnum`.  Thin top-level driver for the
analytic / numerical recursive-integral generators.  Both
submodules verified (below) to be clean.  **No `sc.*` usage.**

### `osrecurintgana.py`
Imports: `osrecurintglib` only.  Analytic branch; pure symbolic
math over the OS integral recursion relations.  **No `sc.*`
usage.**

### `osrecurintglib.py`
Imports: nothing shown in the header sweep (no top-level
`import` / `from` statements).  Pure utility library for
`osrecurintgana.py` / `osrecurintgnum.py`.  **No `sc.*`
usage.**

### `osrecurintg_makenum.py`
Imports: `sympy`, `argparse`, `os`, `sys`, `re`,
`osrecurintglib`.  Symbolic code generator for the numerical
branch.  **No `sc.*` usage.**

### `osrecurintgnum.py`
Imports: `osrecurintglib` only.  Numerical branch of the
recursive OS integral generator.  **No `sc.*` usage.**

### `pot2plot.py`
Imports: `argparse`, `os`, `sys`, `math`, `numpy`.  Potential
plotter.  **No `sc.*` usage.**

### `viewCell.py`
Imports: `argparse`, `os`, `sys`, `datetime`, `numpy`, `math`,
`vedo`.  3D cell visualiser.  Reads coordinates directly from
files (or its own parsers) and hands them to `vedo` for
rendering — does not instantiate `StructureControl`.  **No
`sc.*` usage.**

### `element_pct.py`
Imports: `sys` only.  Elemental-percentage counter; operates on
raw element-name tokens.  **No `sc.*` usage.**

### `dead-md.py`
Imports: `os`, `numpy`, `random`, `subprocess`.  MD-related
helper; does not touch `structure_control`.  **No `sc.*`
usage.**

### `dpmts.py`
Imports: `numpy`, `matplotlib.pyplot`.  Data plotter.  **No
`sc.*` usage.**

### `XYZ.py`
Imports: `argparse`, `os`, `sys`, `datetime`.  XYZ-format
helper / driver.  **No `sc.*` usage.**  Note: the file is
small and does *not* route through `structure_control` even
though it ships alongside a rc companion.

### `XYZrc.py`
Imports: `os` only.  Resource-control defaults for `XYZ.py`.
No list defaults that would cross a `structure_control` call
boundary (no `sc.*` consumer to cross into).  **No `sc.*`
usage.**

### `olcao.py`
Imports: `os`, `sys`, and `from olcaoObject import olcaoJob`.
**No `sc.*` usage.**

Status note (as the task description requested this be
explicitly captured): `olcao.py` is a **pre-Python-3 leftover**.
Line 9 reads `print "Help Stuff"` and line 25 reads
`print cmdParams` — both use the Python 2 print-statement
syntax that is a hard `SyntaxError` under Python 3.  The file
would not even parse, let alone run, on any modern interpreter.
It is therefore dead code relative to the current repository
and contributes nothing to the indexing audit.  Flag for a
separate cleanup decision (delete vs. port), but it is
**not** an indexing issue and has no `structure_control`
interaction whatsoever.

## Non-findings (verified)

- **No indirect `structure_control` imports.**  The only files
  in scope that route imports through a local module are
  `olcaoObject.py` → `olcaoStructures.py`, `olcao.py` →
  `olcaoObject.py`, `interFit.py` → `org_radwavefn.py`,
  `make_veusz_graph.py` → `org_radwavefn.py`, and the
  `osrecurintg*.py` family → `osrecurintglib.py`.  Each of
  these wrapper modules was searched directly for
  `structure_control` / `StructureControl` references and
  returned no matches.

- **No `*rc.py` companion list defaults touch sc.**
  `XYZrc.py` is the only rc companion in the secondary sweep,
  and it contains no list-valued defaults at all (its parent
  `XYZ.py` does not instantiate `StructureControl`), so there
  is no 1-indexed boundary for its defaults to cross.

- **The three "old / new" alternates** (`interFit.new.py`,
  `interFit.old.py`, `make_veusz_graph.old.py`,
  `org_radwavefn.old.py`) were skipped per task instructions
  (primary-file-only sweep).

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK (no sc intersection)**: 23 / 23

No fix phase is required for this sweep — there is nothing to
fix.  The secondary-sweep files form a clean partition: they
are either self-contained computational / plotting utilities or
wrappers around subsystems (`element_data`, `org_radwavefn`,
`osrecurintglib`, HDF5 / vedo / matplotlib / sympy) that are
disjoint from the `structure_control.py` 1-indexed API.

## Follow-up observations (out of audit scope)

- **`olcao.py` Python 2 leftover.**  As noted above, `olcao.py`
  uses Python 2 `print` statements and would not parse under
  Python 3.  Not an indexing concern, but worth flagging for a
  separate housekeeping pass (delete or port to Python 3).

- **`lmp2skl.py` is a one-line stub** (`import os` only).  Not
  an indexing concern, but if a real implementation is ever
  added it will need to be re-swept.

- **`osrecurintglib.py` has no top-level imports.**  Unusual but
  not suspicious — it is pure utility code called only from the
  sibling `osrecurintg*.py` generators.  Again, not an indexing
  concern.
