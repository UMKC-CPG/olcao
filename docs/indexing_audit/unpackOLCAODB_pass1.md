# unpackOLCAODB.py — Pass 1 Audit

**Status: zero findings.** Smallest converted script in the
priority chain; carries no scientific data structures.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | all three subsystems |

## Scope and method

- Perl original: `src/scripts/unpackOLCAODB` (42 lines)
- Python port:   `src/scripts/unpackOLCAODB.py` (~210 lines)
- rc companion:  none

The script is a pure filesystem / subprocess utility: `chdir`s
into `$OLCAO_DATA`, extracts five `.tgz` database archives with
`tar -xzf`, then runs the `remake` script inside `spaceDB/` to
regenerate the individual space-group files.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings | `ScriptSettings.__init__`, `parse_command_line`, `recordCLP` | OK |
| B | Database unpacking | `unpack_databases` | OK |
| C | Entry point | `main`, `__main__` guard | OK |

## Coverage verification

### 1. Public entry-point params / returns

- `ScriptSettings.__init__(self)`: consistency scaffold; no
  parameters. **OK.**
- `parse_command_line(self)`: no arguments beyond `-h/--help`.
- `recordCLP(self)`: writes `sys.argv` to the `command` file; no
  atom / axis / element / species data flows through.
- `unpack_databases()`: no parameters, no return value.
- `main()`: no parameters, no return value.

No public method accepts or returns an atom / axis / element /
species / vector index.

### 2. `sc.*` / `structure_control` calls

**None.** The full import list is `argparse`, `os`, `subprocess`,
`sys`, `datetime.datetime` (lines 32-36). **OK.**

### 3. Local bookkeeping lists

Exactly **one** list is constructed in the entire script
(`unpack_databases`, `unpackOLCAODB.py:149-155`):

```python
databases = [
    "atomicBDB.tgz",
    "atomicPDB.tgz",
    "precursorDB.tgz",
    "spaceDB.tgz",
    "sybdDB.tgz",
]
```

Flat 0-indexed list of five filename strings, consumed only by a
`for db in databases:` value iteration at `unpackOLCAODB.py:158`.
Not atom / axis / element / species / vector keyed. The Perl
original hard-codes the same five filenames as separate literal
`system()` calls — no 1-indexed `@databases` array exists in Perl
either. **OK.**

### 4. Loop bounds

Exactly **two** loops:

- `for db in databases:` at `unpackOLCAODB.py:158` — value
  iteration. No index arithmetic.
- `for argument in sys.argv:` at `unpackOLCAODB.py:112` — value
  iteration over the runtime argv. Not atom / axis / element /
  species / vector.

No `range(...)` call appears anywhere.

## Non-findings

- **No scientific data path.** The script's purpose is extracting
  tarballs — a pre-computation setup step.
- **`chdir` sequence matches Perl exactly.**
- **Error handling** is the Python enhancement (Perl bare `die` →
  Python descriptive `sys.exit` with captured stderr). Not
  indexing.
- **`ScriptSettings`** is a scaffold with no user-configurable
  parameters.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: all three topics (A-C)

## Fix order

None required.
