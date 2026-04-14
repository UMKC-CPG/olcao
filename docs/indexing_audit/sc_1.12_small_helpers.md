# structure_control.py — Section 1.12: Small-helper utilities

## 1.12.a `_prep_line(f, full_line=False)` — OK (dead code)

**Perl** `prepLine` at `StructureControl.pm:8821`: accepts
`(fileHandle, line, splitter)`, reads a line if the file handle is
non-empty, splits with a custom regex, and shifts off a leading
empty token if the split produced one. Returns a 0-indexed
`@values`.

**Python** `_prep_line` at `structure_control.py:9205`: simplified
signature `(f, full_line=False)`; drops the `line` and `splitter`
parameters because callers either already have the string (and use
`str.split()` directly) or only need whitespace tokenization.
Returns a 0-indexed list.

**Severity**: Signature change, not indexing drift. Both return
0-indexed lists; the inner convention matches Perl.

**Dead-code note**: `_prep_line` has **no callers** anywhere in
`structure_control.py`. Every in-file caller that would have used
`prepLine` calls `line.strip().split()` inline or uses a regex
directly. The method is preserved only for API symmetry with the
Perl source. Fine to leave in place but not a live concern.

## 1.12.b `_find_number(values, target)` — signature change (dead code)

**Perl** `findNumber($itemToFind, $array, $startIndex, $endIndex)`
at `StructureControl.pm:8846`: bounded linear scan between
`startIndex` and `endIndex` (both 1-indexed), returns the loop
index (1-indexed position) on match, 0 on not found.

**Python** `_find_number(values, target)` at
`structure_control.py:9247`: drops `startIndex`/`endIndex` — callers
pass a slice if they want a bounded scan. Returns `i + 1`
(1-indexed position) on match, 0 on not found.

**Severity**: Signature change (start/end dropped). Indexing is
preserved: the return value is 1-indexed in both languages.

**Dead-code note**: No callers in `structure_control.py`.

## 1.12.c `_find_string(values, target)` — signature change (dead code)

Same pattern as `_find_number`. Signature change (start/end
dropped), 1-indexed return preserved. No callers.

## 1.12.d `_int_remainder(a, b)` — behavior fix (dead code)

**Perl** `intRemainder` at `StructureControl.pm:8896`: has a latent
bug — returns `$num2 - int($num2)` (fractional part of `num2`),
completely ignoring `num1`. The commented-out alternative computes
the fractional part of the *quotient* `num1/num2`, which is also
not a conventional integer remainder.

**Python** `_int_remainder` at `structure_control.py:9314`: returns
`a % b`, the conventional integer remainder, which matches the
function's name.

**Severity**: Behavior fix of a Perl bug, not indexing drift.

**Dead-code note**: No callers in either Perl or Python.

## 1.12.e `_min(a, b)` — OK

Trivial two-argument min. No indexing concerns.

## 1.12.f `stable_sort(values, keys)` — DRIFT + signature change (dead code in this file)

**Perl** `stableSort($value_ref, $index_ref, $numItems)` at
`StructureControl.pm:8788`: **modifies `$value_ref` and
`$index_ref` in place**, both 1-indexed, sort from `$i=2` to
`$numItems` (note: Perl also has a latent bug using `$numAtoms`
module global instead of the passed `$numItems` parameter).

**Python** `stable_sort(values, keys)` at `structure_control.py:9107`:
**returns a new list of 0-indexed integers**, computed as
`[i for i, _ in sorted(enumerate(values), key=lambda x: keys[x[0]])]`.
No in-place modification. No `numItems` parameter.

**Severity**: DRIFT (return values are 0-indexed while Perl
operated on 1-indexed arrays in place) AND major signature change
(in-place → functional, drops `numItems` parameter). The Python
version also fixes the latent `$numAtoms`-vs-`$numItems` bug.

**Dead-code note in this file**: The only callers inside
`structure_control.py` are the `sort_atoms` method, which does not
actually call `stable_sort` — it builds its own `sorted(...)` call
inline. However, `stable_sort` is part of the public class API and
may be called by external scripts (`makeinput.py`, `condense.py`,
etc.); those callers should be checked during the fix phase.

**Recommended fix**: Either delete `stable_sort` entirely (if no
external callers exist) or restore the 1-indexed in-place signature
to match Perl. Choose after Pass-2 cross-file call-site sweep.

## 1.12.g `gaussian_broaden(data, sigma)` — DRIFT + signature change

**Perl** `gaussianBroaden` at `StructureControl.pm:8908`: signature
`($sigma, $start, $end, $points_ref, $graph_ref)`. Accepts explicit
`$start`/`$end` loop bounds on 1-indexed arrays, reads
`$points_ref->[$start..$end]`, accumulates into `$graph_ref->[$start..$end]`
(the output is written back into the caller's ref). Uses
`$graph[$graphValue] +=` on an *uninitialised* local array, then
copies `@graph[$start..$end]` into `$graph_ref` at the end.

**Python** `gaussian_broaden(data, sigma)` at
`structure_control.py:9136`: drops `start`/`end`/`graph_ref`; loops
over `range(len(data))` (0-indexed). Returns a new list of the
same length. Documented: *"Callers should ensure no indices need
to be skipped before calling."*

**Severity**: DRIFT + signature change.

- DRIFT: the Python version operates on a 0-indexed `data` array
  (index 0 is live), while Perl's `points_ref` was 1-indexed (with
  start=1 in the typical caller).
- Signature change: start/end bounds and the output ref are
  removed; the Python version always processes the full array.

**Live caller**: `compute_rpdf` at `structure_control.py:7752`:

```python
broadened = self.gaussian_broaden(raw, sigma)
```

The caller builds `raw` explicitly before this line; the current
layout of `raw` (1-indexed or 0-indexed) decides whether the DRIFT
is a local bug. This needs to be confirmed during the fix pass.

**Recommended fix**: Restore the 1-indexed layout. Match Perl's
signature more closely, or at minimum document that `data` is
1-indexed with `data[0] = None`.

## Section 1.12 verdict

7 utility helpers audited. 2 DRIFT (`stable_sort`, `gaussian_broaden`),
2 signature changes with preserved indexing (`_find_number`,
`_find_string`), 1 signature change on tokens (`_prep_line`), 1
behavior fix (`_int_remainder`), 1 trivial OK (`_min`). Four of
the seven (`_prep_line`, `_find_number`, `_find_string`,
`_int_remainder`) are dead code within `structure_control.py` but
preserved for API symmetry — external callers in the other Python
scripts need verification during Pass 2.
