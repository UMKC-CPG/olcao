# structure_control.py — Section 1.7: Remaining Pass-1 methods

## 1.7.a `apply_perturbation(magnitude)` — OK

**Perl** `applyPerturbation` at `StructureControl.pm:5294`:
`@perturbRThetaPhi[1..3]`, `@perturbXYZ[1..3]`, atom loop
1..numAtoms (though note the latent `my $numAtoms` shadowing bug
that collapsed the loop — Python fixes this, documented).

**Python** at `structure_control.py:5501`: atom loop 1-indexed,
writes `direct_xyz[atom][1..3]`. Temporary `dx/dy/dz` scalars
rather than `[None, ...]` lists, but since they are never passed to
a helper and are consumed immediately by the `+=` writes into the
1-indexed `direct_xyz`, there is no drift.

## 1.7.b `apply_filter(min_dist_filter)` — OK (with mechanism change)

**Perl** `applyFilter` at `StructureControl.pm:5335`: 1-indexed
everywhere; uses `@rejectedAtoms[1..numRejected]` as a linear list
and a `foreach` search to test membership.

**Python** at `structure_control.py:5561`: 1-indexed everywhere;
uses `rejected = set()` for O(1) membership testing instead of
Perl's linear search. A mechanism improvement, not an indexing
change.

Also builds each new `fract_abc_local` row as `[None, a, b, c]` to
preserve the 1-indexed inner convention when appending. Correct.

## 1.7.c `compute_qn(bridges, ions)` — REVERSE DRIFT

**Perl** `computeQn` at `StructureControl.pm:7456`:

- `$atomQn[1..numAtoms]` is 1-indexed on atoms (matching the
  per-atom-array convention).
- `$absoluteSysQn[0..8]`, `$fractionalSysQn[0..8]` are **0-indexed**:
  the loop is `foreach $n (0..8)` and the storage indexes directly
  by the Q^n value `n`. This is a deliberate Perl choice — `n` is
  literally the array index, there is no sentinel slot, and the
  arrays genuinely use 0 as a live data index.
- `$netConnQn` accumulates `$n * $fractionalSysQn[$n]` for
  `n in 0..8`.

**Python** `compute_qn` at `structure_control.py:7902`:

- `atom_qn[1..num_atoms]` 1-indexed on atoms — OK.
- `absolute_sys_qn` and `fractional_sys_qn` were **forced to
  1-indexed** with an `[n+1]` offset: `absolute_sys_qn = [None] + [0]*9`,
  then writes at index `atom_qn[atom] + 1`, reads at `n+1` in a
  `range(9)` loop.

**Severity**: REVERSE DRIFT. This is the inverse of the Section 1.1
pattern — there, the Python port pulled Perl's 1-indexed helpers
down to 0-indexed; here, it pushed a Perl 0-indexed array up to
1-indexed. Both directions violate the rule. The docstring
explicitly acknowledges the change: *"Index 0 is the None
placeholder; Perl stored these 0-indexed directly by n, but Python
uses the standard 1-indexed convention."*

**Recommended fix**: Revert to Perl's 0-indexed storage. Allocate
`absolute_sys_qn = [0] * 9` and `fractional_sys_qn = [0.0] * 9`,
remove the `+1` offsets everywhere, and let `n` be the direct
index. The `range(9)` loops become `for n in range(9)`.

This finding is particularly important as a precedent: the user's
rule is **preserve Perl conventions in both directions**, not
"standardize on 1-indexing throughout Python."

## 1.7.d `save_ring(ring_len, curr_ring)` — OK

**Perl** `saveRing` at `StructureControl.pm:7931`: mixed convention
correctly — rings outer dimension is 1-indexed
(`foreach $ring (1..ringCounts-1)`), but the atom list *within* a
ring is 0-indexed (`foreach $atom (0..$targetRingLen-1)`, slice
`@{$currRing_ref}[0..$targetRingLen-1]`).

**Python** at `structure_control.py:8349`: matches the mixed
convention. Outer ring loop uses `range(1, ring_counts[len])`;
inner sort and comparison use the 0-indexed list slice
`curr_ring[:ring_len]` and direct list equality.

## 1.7.e `compute_ring_level_delta(...)` — OK

Pure integer arithmetic on `target_ring_len`, `ring_midpoint`,
`curr_ring_max_index`. No indexing concerns. Matches Perl.

## 1.7.f `obtain_atomic_interaction(...)` — OK

**Perl** `obtainAtomicInteraction` at `StructureControl.pm:6912`:
`foreach $item1 (1..numItems1)`, `foreach $item2 (1..numItems2)`,
reads `$item1XYZ_ref->[$item1][1..3]`, builds local `@diff[1..3]`,
passes `\@diff` to `recordInteractionData`.

**Python** at `structure_control.py:7268`: `range(1, num_items1+1)`,
`range(1, num_items2+1)`, reads `item1_xyz[item1][1..3]`. Calls
`record_interaction_data(..., [None, dx, dy, dz], ...)` —
**correctly builds a 1-indexed 4-element list** to pass the diff
vector. Good.

## 1.7.g `initialize_interaction_data(...)` — OK

Both Perl and Python preserve 1-indexing throughout. Python
`initialize_interaction_data` at line 7393 uses `[None] + [0] * N`
patterns for every per-item array (`num_bonds`, `bonded_ext`,
`bond_length_ext`, `rpdf`, `q_bar`, `ext_scan_dist`, etc.), and all
inner arrays that get `.append()`-ed are pre-initialized with a
`[None]` sentinel slot so that `append()` lands at index 1 and
matches Perl's `$bonded_ext[$item1][++$numBonds[$item1]]`
assignment. This is the correct pattern — `append` after `[None]`
initialization is equivalent to Perl's `[1..N]` assignment.

Note the BOO `q_bar` layout is 3-deep 1-indexed: `q_bar[item][m_idx]`
where `m_idx` is 1..2*l+1, and each entry is a `[real, imag]` pair
where the 0,1 sub-indices are genuinely 0-indexed (real and imag
are pair components, not atom/axis-style indices). Matches Perl
intent.

## 1.7.h `record_interaction_data(...)` — OK (with 1.7.i consequence)

**Perl** `recordInteractionData` at `StructureControl.pm:7178`:
receives `\@diff` (1-indexed), passes to `sphericalAngles` which
reads `$diff_ref->[1..3]`.

**Python** at `structure_control.py:7575`: receives `diff` as
1-indexed `[None, dx, dy, dz]`, reads `diff[1..3]` correctly. Calls
`.append()` on pre-initialized `[None]` lists — matches Perl
conventions.

*However*, at line 7670 there is a call-site consequence:

```python
theta, phi = self.spherical_angles([diff[1], diff[2], diff[3]])
```

A manual 0-indexed 3-element list is built here to bridge into the
drifted `spherical_angles` helper. This is a symptom of 1.7.i
below, not a local bug.

## 1.7.i `spherical_angles(vector)` — DRIFT

**Perl** `sphericalAngles` at `StructureControl.pm:8077`: takes
`$diff_ref` (1-indexed), reads `$diff_ref->[1..3]`, mutates in
place with epsilon clamping.

**Python** `spherical_angles` at `structure_control.py:8771`: takes
`vector` (0-indexed length-3), reads `vector[0..2]`, returns
without mutating (behavior improvement: Perl version mutated the
caller's array, Python does not).

**Severity**: DRIFT. Same pattern as Section 1.3 — throwaway vector
helper forced from 1-indexed to 0-indexed. Docstring explicitly
acknowledges the change: *"The Perl version accepts a 1-indexed
reference [$x,$y,$z] and clamps in-place; here we clamp via a list
comprehension so the caller's array is unchanged."*

**Recommended fix**: Restore 1-indexed input. The non-mutating
behavior can be kept as a Python-side improvement by building a
local 1-indexed copy `clamped = [None, ...]` rather than clamping
in place.

## Section 1.7 verdict

9 methods audited. 1 REVERSE DRIFT (`compute_qn` — a
precedent-setting finding), 1 DRIFT (`spherical_angles`), 7 OK.
All the interaction-data machinery (`obtain_*`, `initialize_*`,
`record_*`) cleanly preserves Perl's 1-indexing with correct use
of `[None]`-sentinel + append.
