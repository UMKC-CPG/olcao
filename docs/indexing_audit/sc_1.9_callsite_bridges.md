# structure_control.py — Section 1.9: Call-site bridge scaffolding

Several methods that are themselves *internally* correct
nevertheless contain call-site scaffolding that manually bridges
between the 1-indexed persistent arrays and the 0-indexed drifted
helpers. These bridge sites will disappear once the helpers in
Sections 1.1, 1.3, 1.4, and 1.7.i are restored to 1-indexed. Noting
them here so they are fixed in the same pass.

## 1.9.a `rotate_all_atoms()` — 4 bridge sites

**Python** at `structure_control.py:5894`. Contains four places
where a 0-indexed length-3 list is manually built from 1-indexed
`direct_xyz[atom][1..3]` values:

- Lines 5943–5945 (pass 1 input):
  `pt = [direct_xyz[atom][1], direct_xyz[atom][2], direct_xyz[atom][3]]`
- Lines 5947–5949 (pass 1 output):
  `direct_xyz[atom][1] = rpt[0]; [2] = rpt[1]; [3] = rpt[2]`
- Lines 5987–5989 (pass 2 input): same pattern
- Lines 5991–5993 (pass 2 output): same pattern

All four exist **only** because `rotate_one_point` (Section 1.4.b)
was drifted to 0-indexed I/O. Once 1.4.b is fixed, these become:

```python
rpt = self.rotate_one_point(direct_xyz[atom])
direct_xyz[atom] = rpt
```

or equivalent in-place assignment using 1-indexed slices.

## 1.9.b `record_interaction_data()` — 1 bridge site

**Python** at `structure_control.py:7670`:

```python
theta, phi = self.spherical_angles([diff[1], diff[2], diff[3]])
```

Manually builds a 0-indexed 3-element list to bridge into the
drifted `spherical_angles` (Section 1.7.i). Once 1.7.i is fixed,
this becomes `self.spherical_angles(diff)`.

## 1.9.c `create_min_dist_matrix` / surrounding callers — verify after fix

The `obtain_atomic_interaction` caller at line 7384 builds
`[None, dx, dy, dz]` — **already correctly 1-indexed**. This is
evidence that the author *can* pass 1-indexed lists when the helper
expects them; the drifts in 1.1/1.3/1.7 are convention choices at
the helper side, not an inability at the call sites.

## 1.9.d `direct_xyz2fract_abc` caller at line 4826 — 1 bridge site

Already noted in Section 1.1.a. Currently:

```python
temp_fract = self.direct_xyz2fract_abc(
        [new_xyz[1], new_xyz[2], new_xyz[3]])
```

After fixing 1.1.a, becomes:

```python
temp_fract = self.direct_xyz2fract_abc(new_xyz)
```

## 1.9.e `rotate_arb_axis` → `define_rot_matrix` — 1 bridge site

**Python** at `structure_control.py:5792–5794`:

```python
mag = math.sqrt(sum(v*v for v in axis))
unit = [v / mag for v in axis]
self.define_rot_matrix(unit, angle_deg)
```

`unit` is built as a 0-indexed length-3 list from the 0-indexed
`axis` parameter and passed to `define_rot_matrix`, which has the
0-indexed-input MIXED problem in Section 1.4.a. Added during Pass 2.

Once Section 1.4.a is fixed to accept a 1-indexed `axis` input,
this becomes:

```python
# axis is already 1-indexed [None, x, y, z]
mag = math.sqrt(sum(v*v for v in axis[1:]))
unit = [None] + [v / mag for v in axis[1:]]
self.define_rot_matrix(unit, angle_deg)
```

The `rotate_arb_axis` public signature (Section 1.4.c) still
needs revisiting separately — it currently accepts a 0-indexed
`axis` parameter.

## Section 1.9 totals

7 call-site bridge sites catalogued (6 in 1.9.a–d during Pass 1 +
1 in 1.9.e during Pass 2). All will be removed when the underlying
helpers (Sections 1.1, 1.3, 1.4, 1.7.i) are restored to 1-indexed
I/O.
