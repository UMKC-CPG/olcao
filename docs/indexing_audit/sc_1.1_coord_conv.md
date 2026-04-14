# structure_control.py — Section 1.1: Coordinate conversion helpers

Four helper methods that convert between Cartesian (XYZ) and
fractional or direct-space (ABC) coordinate representations. In Perl
all four accept 1-indexed refs and return 1-indexed results; in
Python they have been converted to 0-indexed length-3 lists,
docstrings openly acknowledging the change.

## 1.1.a `direct_xyz2fract_abc(xyz)` — DRIFT

**Perl** `directXYZ2fractABC` at `StructureControl.pm:4736`

- Takes a ref `$posXYZ_ref` to a 1-indexed array. Reads
  `$posXYZ_ref->[$xyzAxis]` for `$xyzAxis in 1..3`.
- Returns `(@posABC[1..3])` — a 3-element list assembled from slots
  `[1..3]` of a 1-indexed local.
- Single caller at `StructureControl.pm:3515` passes
  `\@{$newDirectXYZ[$newNumAtoms]}` (a ref to a 1-indexed atom row)
  and immediately `unshift`s `""` onto the returned list to
  re-establish the 1-indexed convention before iterating `1..3`.

**Python** `direct_xyz2fract_abc` at `structure_control.py:3985`

- Takes `xyz`, a 3-element 0-indexed list. Accesses `xyz[0..2]`.
- Returns a 3-element 0-indexed list `abc[0..2]`.
- Caller at `structure_control.py:4826` manually constructs a
  0-indexed slice `[new_xyz[1], new_xyz[2], new_xyz[3]]` to bridge
  from the 1-indexed `new_xyz` into this helper's 0-indexed
  interface, then iterates the result with `range(3)`.

**Assessment**: End-to-end behavior is correct (the slice conversion
and loop bounds line up), but the helper's interface convention was
silently switched from 1-indexed to 0-indexed. Every other coordinate
method in the file (`get_direct_xyz`, `get_direct_abc`,
`get_fract_abc`, `get_min_max_xyz`, etc.) preserves the 1-indexed
convention, so this helper is an island of 0-indexing surrounded by
1-indexed code. That forces callers to perform manual slice
conversions at the boundary, which is exactly the kind of friction
the user's rule was written to prevent.

**Recommended fix**: Restore the 1-indexed interface. Accept a
1-indexed input (read `[1..3]`) and return a 1-indexed list (slot 0
holds `None` or `""` sentinel, data in `[1..3]`). Update the caller
at `structure_control.py:4826` to pass `new_xyz` directly (already
1-indexed) and read `temp_fract[1..3]` in the subsequent loop.

## 1.1.b `direct_xyz2direct_abc(xyz)` — DRIFT

**Perl** `directXYZ2directABC` at `StructureControl.pm:4767`

Same structure as `directXYZ2fractABC`: 1-indexed input ref, returns
`(@posABC[1..3])`. Identical argument convention.

**Python** `direct_xyz2direct_abc` at `structure_control.py:4016`

- Takes 0-indexed `xyz` (length 3).
- Internally delegates to `direct_xyz2fract_abc(xyz)` (line 4046),
  inheriting that helper's 0-indexed interface.
- Returns 0-indexed `abc[0..2]`.

**Assessment**: Same DRIFT as 1.1.a — inherited through the
delegation chain. Fix 1.1.a first; then this method's interface
should also be restored to 1-indexed (input and output).

## 1.1.c `fract_abc2direct_xyz(fract)` — DRIFT

**Perl** `fractABC2directXYZ` at `StructureControl.pm:4801`

- Takes a ref `$posABC_ref` to a 1-indexed array. Reads
  `$posABC_ref->[$abcAxis]` for `$abcAxis in 1..3`.
- Note: the Perl version writes into `@posXYZ[1..3]` but does **not**
  return the result — this is almost certainly a latent Perl bug
  (unused computation), but that is a separate concern from the
  indexing audit.

**Python** `fract_abc2direct_xyz` at `structure_control.py:4052`

- Takes 0-indexed `fract` (length 3).
- Returns 0-indexed `xyz[0..2]`.

**Assessment**: DRIFT. Same class of change as 1.1.a/1.1.b. The
Python port appears to have *also* fixed the Perl bug of not
returning the computed result, which is a valuable behavioral
improvement — but the interface convention should still be
1-indexed to match the ambient style.

## 1.1.d `get_direct_xyz_point(fract)` — DRIFT (mild)

**Perl** `getDirectXYZPoint` at `StructureControl.pm:4659`

- Takes 3 scalar arguments via `$_[0]`, `$_[1]`, `$_[2]` (the
  positional-args convention is effectively 0-indexed from Perl's
  `@_`), then immediately copies them into 1-indexed local slots
  `$Pabc[1..3]`.
- Returns `(@Pxyz)` — the full Perl array including the undefined
  slot `[0]`, with data in `[1..3]`. Caller reads `$values[1..3]`.

**Python** `get_direct_xyz_point(fract)` at
`structure_control.py:3918`

- Takes `fract`, a 0-indexed length-3 list. Reads `fract[0..2]`.
- Returns a 0-indexed length-3 list `xyz[0..2]`.
- Caller `get_direct_xyz` at line 3913 reads `xyz[0], xyz[1],
  xyz[2]` and stores them into the 1-indexed
  `self.direct_xyz[atom][1..3]`.

**Assessment**: *Milder* DRIFT than 1.1.a-c because the Perl
original already had a mixed convention at this boundary (positional
args in, 1-indexed return out). The Python port settled on
fully-0-indexed. The consistent alternative is 1-indexed at both
ends (with a sentinel at slot 0), matching the surrounding file
convention.

## Section 1.1 verdict

4 DRIFT findings. All four helpers will be restored to 1-indexed
input and output, and their call sites updated to remove the manual
slicing scaffolds.
