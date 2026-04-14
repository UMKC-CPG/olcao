# structure_control.py — Section 1.4: Rotation helpers

## 1.4.a `define_rot_matrix(axis, angle_deg)` — MIXED

**Perl** `defineRotMatrix` at `StructureControl.pm:2859`: takes
`$rotAngle, $rot_ref` (note: opposite order from Python); reads
`$rot_ref->[1..3]` and writes `$rotMatrix[1..3][1..3]`. Both input
and output are 1-indexed.

**Python** `define_rot_matrix` at `structure_control.py:5797`: takes
`(axis, angle_deg)` — **argument order is inverted** relative to
Perl. Reads `axis[0..2]` (0-indexed input) but writes
`self.rot_matrix[1..3][1..3]` (1-indexed output).

**Severity**: MIXED — the output matrix is correctly 1-indexed but
the input axis is 0-indexed. Within a single method, the convention
is inconsistent. Also: argument order is swapped (`(angle, axis)` in
Perl vs `(axis, angle)` in Python).

## 1.4.b `rotate_one_point(point)` — DRIFT + signature change

**Perl** `rotateOnePoint` at `StructureControl.pm:2907`: takes
`$point_ref, $orig_ref` (both 1-indexed), mutates `$point_ref` in
place to hold the rotated result.

**Python** `rotate_one_point(point)` at `structure_control.py:5848`:
takes `point` (0-indexed) only — the `orig_ref` parameter is
**dropped** because the Python rotation pipeline always rotates
about the origin. Returns a new list instead of mutating in place.

**Severity**: DRIFT on indexing (0-indexed input/output). The
missing `orig` parameter is a documented signature change (the Perl
caller always passed `(0,0,0,0)` anyway, so behavior is preserved).
The mutate-in-place → return-new change is a Python idiom
improvement.

Inside the method, the matrix is read as `self.rot_matrix[1..3][k+1]`
with `k in range(3)` — a deliberate `k+1` offset that correctly
bridges the 0-indexed loop variable into the 1-indexed matrix.

## 1.4.c `rotate_arb_axis(axis, angle_deg)` — signature + semantics change

**Perl** `rotateArbAxis` at `StructureControl.pm:2838`: low-level
helper taking `(point_ref, unit_vect_ref, rot_angle_RADIANS)`.
Rotates a single point; the caller handled looping over atoms.

**Python** `rotate_arb_axis(axis, angle_deg)` at
`structure_control.py:5769`: high-level entry point that normalizes
the axis, builds the rotation matrix, and rotates all atoms via
`rotate_all_atoms()`. Angle is now in **degrees**.

**Severity**: Not an indexing issue per se — this is a major
refactor. The Perl low-level helper role is filled by
`rotate_one_point` in Python. Both behavior changes are documented
in the docstring. Flag for awareness but outside the indexing
audit.

## Section 1.4 verdict

1 MIXED (1.4.a), 1 DRIFT (1.4.b), 1 refactor (1.4.c, not strictly
an indexing finding). The rotation machinery has been substantially
reorganized beyond the indexing question.
