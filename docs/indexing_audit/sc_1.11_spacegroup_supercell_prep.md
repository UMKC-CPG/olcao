# structure_control.py — Section 1.11: Space-group, supercell, prep-surface

## 1.11.a `apply_space_group()` — OK

**Python** at `structure_control.py:4213`. Calls the external
`applySpaceGroup` Fortran binary via subprocess to expand the
asymmetric unit.

- `real_lattice[axis][1..3]` — 1-indexed on both dimensions.
- `rl = self.real_lattice[axis]; rl[1..3]` — local alias of a
  1-indexed row.
- Atom loop `range(1, self.num_atoms + 1)` — 1-indexed.
- Reset arrays initialized as `[None]` sentinel lists; `.append()`
  extends them at index 1..num_atoms. ✓
- `self.fract_abc.append([None, fa1, fa2, fa3])` — 1-indexed inner. ✓
- Token arrays from `readline().split()` (`vals[0..4]`) are
  legitimately 0-indexed parsing tokens matching Perl.

**Verdict**: OK.

## 1.11.b `apply_supercell()` — OK

**Python** at `structure_control.py:4362`. Replicates the cell along
each abc axis.

- `self.supercell[1..3]`, `self.sc_mirror[1..3]` — 1-indexed reads.
- `self.mag[abc_axis] *= n` and
  `self.real_lattice[abc_axis][xyz_axis] *= n` — 1-indexed writes.
- `new_fa = [None, 0.0, 0.0, 0.0]` — 1-indexed inner row built per
  new atom.
- `new_fract_abc.append(new_fa)` and analogous for
  `new_atom_element_name`, `new_atom_element_id`, etc. — correct
  1-indexed builders with `[None]` sentinel.
- Post-loop re-allocation of `direct_xyz` and `direct_abc` uses the
  `[None] + [[None, 0.0, 0.0, 0.0] for _ in range(self.num_atoms)]`
  pattern — 1-indexed outer, 1-indexed inner.

**Verdict**: OK.

## 1.11.c `prep_surface(hkl)` — MANY call-site bridges (no original drift)

**Python** at `structure_control.py:4488`. This is a large (~500
line) method that orients and cuts the cell to expose an (hkl)
surface plane. It is not itself introducing new drift, but it is
the single largest consumer of the drifted helpers from Sections
1.1, 1.3, 1.4, and 1.7.i. Every use of `dot_product`,
`cross_product_mag`, `get_plane_normal`, `get_vector_angle`,
`define_rot_matrix`, `rotate_one_point`, and `direct_xyz2fract_abc`
inside this method has to manually bridge between the 1-indexed
persistent `real_lattice`/`direct_xyz` arrays and the 0-indexed
temporary vectors expected by those helpers.

Bridge sites inside `prep_surface` (non-exhaustive — these all go
away automatically when the helpers are fixed):

- Line 4586: `diagonal = [0.0, 0.0, 0.0]` — 0-indexed local.
- Line 4589: `diagonal[xyz_axis - 1] += self.real_lattice[abc][xyz]`
  — `-1` bridge.
- Line 4640: `uvw_normal = [0.0, 0.0, 0.0]` — 0-indexed.
- Line 4643: `uvw_normal[xyz_axis - 1] += ...` — `-1` bridge.
- Lines 4666–4669: `lp = [0.0, 0.0, 0.0]; lp[m-1] += ...` —
  0-indexed point in `cell_points`.
- Line 4694–4712: `new_real_lattice[i][1..3] = tv[0..2]` — writing
  1-indexed slots from a 0-indexed temp.
- Lines 4740–4760: `new_cell_center`, `face_center[1..6]` — outer
  1-indexed, inner 0-indexed. Mixed.
- Lines 4770–4780: `rl1, rl2, rl3`, `origin_v`, `face_normal[1..6]`
  — inner 0-indexed.
- Line 4783: `f2c = [... for i in range(3)]` — bridge into
  `dot_product`.
- Line 4805–4808: `new_xyz = [None, ...]` 1-indexed — good.
- Line 4815: `f2a = [new_xyz[ax] - face_center[face][ax-1] for ax
  in range(1, 4)]` — mixed 1-indexed read + `-1` bridge.
- Line 4826–4827: `direct_xyz2fract_abc([new_xyz[1], new_xyz[2],
  new_xyz[3]])` — already noted as Section 1.9.d.
- Line 4828: `for i in range(3)` reading 0-indexed `temp_fract`.
- Lines 4869, 4870, 4873: `x_axis = [1.0, 0.0, 0.0]`, etc. —
  0-indexed axes.
- Line 4885: `define_rot_matrix(rot_axis, ...)` — 0-indexed
  `rot_axis` input (Section 1.4.a).
- Lines 4889–4892, 4900–4903, 4934–4937, 4945–4948: four copies of
  the `pt = [real_lattice/direct_xyz[...][j] for j in range(1, 4)]`
  / `rpt = rotate_one_point(pt)` / `[...][j+1] = rpt[j]` pattern.
- Lines 4919–4921: `yz_proj = [0.0, real_lattice[2][2],
  real_lattice[2][3]]` — 0-indexed 3-vector for use with
  `get_vector_angle`.

**Severity**: Not a finding in its own right — `prep_surface` is a
victim, not a perpetrator. All bridges disappear once Sections 1.1,
1.3, 1.4.a/b, and 1.7.i are fixed.

**Recommended fix sequence**: after restoring the helpers in
Sections 1.1/1.3/1.4/1.7.i to 1-indexed I/O, do a single sweep of
`prep_surface` to:

1. Replace `diagonal = [0.0, 0.0, 0.0]` and similar with
   `diagonal = [None, 0.0, 0.0, 0.0]` and drop the `-1` offsets.
2. Replace `face_center[face]`, `face_normal[face]`, `rl1/rl2/rl3`,
   `uvw_normal`, `x_axis`, `y_axis`, `origin_v`, `rot_axis`,
   `yz_proj` with 1-indexed `[None, x, y, z]` layout.
3. Pass `direct_xyz[atom]` directly to `rotate_one_point` and
   assign the 1-indexed return back to `direct_xyz[atom]`.
4. Pass `real_lattice[abc_axis]` directly to `rotate_one_point`
   and assign back.
5. Pass `new_xyz` directly to `direct_xyz2fract_abc`; read
   `temp_fract[1..3]` in the subsequent loop.

## Section 1.11 verdict

2 OK methods (`apply_space_group`, `apply_supercell`), 1 large
consumer of drifted helpers (`prep_surface` — ~25 bridge sites
inside one method). No new drift introduced by Section 1.11 methods
themselves.
