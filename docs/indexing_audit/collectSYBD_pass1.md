# collectSYBD.py — Pass 1 Audit

**Status: zero indexing findings.** `collectSYBD.py` is *not* a Perl
port — no `collectSYBD` ancestor exists in `src/scripts/`. It is a
179-line greenfield Python stub with no `structure_control`
intersection and no per-atom/axis/element/species arrays.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | all subsystems |

## Scope

- Perl original: **none** (greenfield Python, not a port)
- Python port:   `src/scripts/collectSYBD.py` (179 lines)
- rc companion:  `src/scripts/collectSYBDrc.py` (15 lines)

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings | `ScriptSettings.*` | OK |
| B | File-list builder | `make_SYBD_file_list` | OK (+ functional bug, out of scope) |
| C | hskp row extractor | `get_hskp_data` | OK (+ **syntax error**, out of scope) |
| D | Driver | `main` | OK |
| E | rc companion | `collectSYBDrc.parameters_and_defaults` | OK |

## Per-topic detail

### A. `ScriptSettings` (lines 11-103)

Holds exactly one user-facing parameter, `self.hskp`, a single-
letter high-symmetry k-point name (default `"G"` from rc). No
list attributes, no slot-0 sentinels, no indices. **OK.**

### B. `make_SYBD_file_list` (lines 108-122)

Walks the directory tree with `os.walk` and appends matching paths
to a flat list `sybd_file_list`. The list is a bag of filesystem
paths — there is no atom/axis/element/species dimension, so the
1-indexed convention does not apply. Matches the natural Perl
idiom. **OK.**

*Out-of-scope functional bug*: `return sybd_file_list.sort()`
returns `None` (Python `list.sort()` is in-place). `main()` then
passes `None` into `get_hskp_data`, which would raise `TypeError`.
Not an indexing finding.

### C. `get_hskp_data` (lines 125-152)

Parses each SYBD file line by line looking for a row whose last
whitespace token equals the target hskp name, then reads the
second-to-last token as the hskp position. The two list accesses
`line[len(line)-1]` and `line[len(line)-2]` are 0-indexed reads of a
whitespace-split token list — the universal idiom in both Perl and
Python. **OK.**

*Out-of-scope syntax error*: lines 140 and 144 contain an `if`
body-less statement and a bare `for` keyword. The file will not
`ast.parse`; the pandas code below is dead until this is fixed.

### D. `main` (lines 155-170)

Three-line orchestrator: build settings, build file list, call
extractor. No per-atom / per-axis / per-element / per-species
arguments. **OK.**

### E. `collectSYBDrc.py` (15 lines)

A single scalar string `"hskp": "G"` in a dict. No list defaults,
no indexed data, no 1-indexed boundary. **OK.**

## Pass 1 coverage checklist

1. Public entry-point params/returns: all scalar. **OK.**
2. `sc.*` calls: **zero.**
3. Local bookkeeping lists — slot 0 sentinelled: no atom/axis/
   element/species lists exist anywhere. **OK.**
4. Loop bounds: `os.walk` and `for sybd_file in sybd_file_list`.
   No `range()` calls. **OK.**
5. rc companion: single scalar `"hskp"`. **OK.**

## Non-indexing observations (out of scope)

- Stub is unrunnable (`get_hskp_data` incomplete statements at
  lines 140, 144; whole file fails `ast.parse`).
- `make_SYBD_file_list` returns `None` via `.sort()`.
- Unused / incomplete pandas consumption.
- **No Perl original**. The Perl → Python port framing does not
  apply to this file.

## Running tally

- **BUG**: 0 (indexing)
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every subsystem

## Fix order

No indexing fixes required. Two functional problems (syntax error,
`.sort()` return) are flagged for a separate follow-up, not part
of the 1-indexing sweep.
