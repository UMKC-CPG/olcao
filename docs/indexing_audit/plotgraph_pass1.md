# plotgraph.py — Pass 1 Audit

**Status: zero findings.** `plotgraph.py` is a pure plotting utility
with no `structure_control` imports or calls, and no per-atom, per-
axis, per-element, or per-species array access anywhere in the file.
It is the second priority-file (after `pdb2skl.py`) to pass clean.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK (verified) | all relevant constructs |

## Scope of audit

- Perl original: `src/scripts/plotgraph`
- Python port:   `src/scripts/plotgraph.py`
- rc companion:  `src/scripts/plotgraphrc.py`

## Verifications (all clean)

### PG-N1 — No structure_control intersection

A grep for `sc\.`, `StructureControl`, and `structure_control`
returns zero matches across `plotgraph.py` and `plotgraphrc.py`.
The script has no call-boundary into `structure_control`, so it
cannot drift against the 1-indexed convention at an sc.* edge.

### PG-N2 — `y_col` / `x_col` hold gnuplot/user-facing column numbers

These settings hold 1-indexed data-file column numbers (the values a
user types at the command line) inside a flat 0-indexed Python
container. This matches the Perl original exactly: `@ycol` in the
Perl script is likewise a 0-indexed container whose *values* are
1-indexed column numbers destined for gnuplot. The Python port
translates these to pandas indices at four well-defined crossing
points via `col - 1` / `col - 2` offsets:

- `read_data_headers` (~line 692, 698)
- `apply_x_limits` (~line 822)
- `apply_y_limits` (~line 867)

The user-facing 1-indexing is preserved; only the internal pandas
lookup is decremented. No DRIFT.

### PG-N3 — List-valued visual settings are flat 0-indexed

`curve_width`, `curve_style`, `curve_color`, `curve_mark`,
`curves_per_subplot`, `subplots_per_fig`, `curve_separation` are flat
0-indexed Python lists with no Perl ancestor (the Perl original used
a single `$lineStyle` string). They never cross a 1-indexed
boundary and do not represent per-atom / per-axis / per-element /
per-species data. No DRIFT.

### PG-N4 — rc companion is scalars + flat 0-indexed lists only

`plotgraphrc.py` contains only scalars and flat 0-indexed list
defaults. No nested list defaults, no sc consumer, no 1-indexed
boundary to violate.

### PG-N5 — Matplotlib-forced 1-indexed counters

`curr_subplot_curve` and `curr_subplot` are 1-indexed because
`plt.subplot(nrows, ncols, index)` takes a 1-indexed `index`
argument. These are matplotlib API requirements, not domain indices,
so they do not participate in the 1-indexing convention at all.

### PG-N6 — Local iterator counters

`wrap_subplot`, `wrap_fig`, `curve` are 0-indexed local iterators
over flat Python lists. They are not domain indices (atom / axis /
element / species) and do not cross a 1-indexed boundary.

### PG-N7 — Loop bounds

All `range()` loops in `plotgraph.py` walk 0-indexed containers that
are either Python-native (matplotlib figure/axes objects, pandas
columns) or flat style lists. None of them touch a 1-indexed
container.

## Non-indexing observation (out of scope for this audit)

`settings.x_col` is defined as a scalar in `plotgraphrc.py`, but
`apply_x_limits` (~line 822) reads it as `settings.x_col[curve + idx]`
when `multi_x_cols == 2`. This would raise a latent `TypeError` the
first time a multi-x-column invocation runs. Not an indexing-audit
finding — flagged here only for separate follow-up.

## Fix order

No fixes required — the file is clean.

## Running tally

- Methods / subsystems audited: full file
- Findings: 0
- Status: **CLEAN**
