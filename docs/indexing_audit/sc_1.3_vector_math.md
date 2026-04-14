# structure_control.py — Section 1.3: Vector math helpers

Pure-math 3-vector routines (dot product, cross product, plane
normal, vector angle). In Perl these accept 1-indexed array refs
(`$vector_ref->[1..3]`) and return 1-indexed lists (slot 0 is an
explicit `0.0` sentinel). In Python they have been uniformly
converted to 0-indexed sequences `[x, y, z]`, with the docstrings
explicitly acknowledging the convention change.

## 1.3.a `dot_product(a, b)` — DRIFT

**Perl** `dotProduct` at `StructureControl.pm:8282`: reads
`$vector_ref->[1] * $vector_ref->[1]` etc. for indices `[1..3]`.

**Python** `dot_product` at `structure_control.py:8585`: reads
`a[0]*b[0] + a[1]*b[1] + a[2]*b[2]`. Docstring explicitly notes the
convention change ("Perl equivalent: dotProduct (line 8282), which
operated on 1-indexed array references").

## 1.3.b `cross_product(a, b)` — DRIFT

**Perl** `crossProduct` at `StructureControl.pm:8299`: writes
`$product[0] = 0.0;` as explicit sentinel, then computes
`$product[1..3]` from `$vector_ref->[1..3]`. Returns a 4-element
array.

**Python** `cross_product` at `structure_control.py:8603`: reads
`a[0..2]`, `b[0..2]`, returns `[cx, cy, cz]` (3-element, no
sentinel). Docstring explicitly notes the convention change.

## 1.3.c `cross_product_mag(a, b)` — DRIFT

Inherits drift from `cross_product` (which it calls internally).
Perl computed `sqrt($cpv[1]**2 + $cpv[2]**2 + $cpv[3]**2)`; Python
uses `c[0]**2 + c[1]**2 + c[2]**2`.

## 1.3.d `normalized_cross_product(a, b)` — DRIFT

Same pattern: Perl reads indices `[1..3]` and normalizes in place on
the 1-indexed result; Python normalizes the 0-indexed result.

## 1.3.e `get_vector_angle(a, b)` — DRIFT

**Perl** `getVectorAngle` at `StructureControl.pm:8363`: reads
`$vector_ref->[1..3]`.

**Python** `get_vector_angle` at `structure_control.py:8683`: reads
`a[0..2]`, `b[0..2]`. Docstring notes "the Perl original uses
1-indexed array references; here a and b are 0-indexed sequences
[x, y, z]."

Also includes a Python-only improvement: the cosine is clamped to
`[-1, 1]` before `acos` to guard floating-point rounding. Documented
as a Python addition. Behaviorally safer; not an indexing issue.

## 1.3.f `get_plane_normal(p1, p2, p3)` — DRIFT

**Perl** `getPlaneNormal` at `StructureControl.pm:8235`: three
1-indexed point refs. Returns a 4-element array with `$normal[1..3]`
holding the result.

**Python** `get_plane_normal` at `structure_control.py:8556`: three
0-indexed length-3 lists. Returns a length-3 list. Docstring:
"Perl used 1-indexed [None, x, y, z] refs; Python uses plain lists."

Also omits the dead-code centroid calculation that Perl computed
but never used — a good cleanup, not an indexing issue.

## Section 1.3 verdict

6 DRIFT findings. The vector math helpers are a *consistent* island
of 0-indexing: every method was converted uniformly, and every
method's docstring openly acknowledges the change. They are
therefore internally self-consistent but collectively violate the
rule. All six must be restored to 1-indexed together, with their
callers updated to pass `[None, x, y, z]` refs and read `[1..3]` on
returns.
