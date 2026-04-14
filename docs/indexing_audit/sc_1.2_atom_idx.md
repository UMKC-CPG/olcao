# structure_control.py — Section 1.2: Atom-index accepting methods

All methods in this section accept one or more atom indices (or an
inclusive atom range) as parameters. For each, the verdict is
whether atom indices, axis indices, and per-atom coordinate arrays
are 1-indexed to match Perl.

| # | Method | File:line | Verdict |
|---|--------|-----------|---------|
| 1.2.a | `get_min_max_xyz(atom1, atom2)` | `structure_control.py:3872` | OK — atom 1-indexed, axis 1-indexed, `max_pos`/`min_pos` use `[None, x, y, z]` sentinel layout. |
| 1.2.b | `get_direct_xyz(atom)` | `structure_control.py:3900` | OK — atom 1-indexed; reads `fract_abc[atom][1..3]` and writes `direct_xyz[atom][1..3]`. Only drift is at the helper `get_direct_xyz_point` (Section 1.1.d). |
| 1.2.c | `get_direct_abc(atom)` | `structure_control.py:3948` | OK — atom 1-indexed, axis 1-indexed throughout. |
| 1.2.d | `get_fract_abc(atom)` | `structure_control.py:4083` | OK — atom 1-indexed, axis 1-indexed, slot 0 set to `None` to match Perl's `""` placeholder. |
| 1.2.e | `shift_xyz_center(atom1, atom2)` | `structure_control.py:5110` | OK — atom range 1-indexed; axis 1-indexed; local `center_pos`, `center_cell`, `center_diff` all use `[None, ...]` 1-indexed layout. Note: `shift_xyz_center` now calls `get_min_max_xyz` internally (documented); Perl required the caller to do so beforehand. Not an indexing issue. |
| 1.2.f | `translate_atoms(atom1, atom2, displacement)` | `structure_control.py:5161` | OK — `displacement` parameter is explicitly documented as 1-indexed `[None, dx, dy, dz]` to match Perl's `$translation_ref->[1..3]`. |
| 1.2.g | `check_bounding_box(atom1=None, atom2=None)` | `structure_control.py:5650` | OK on indexing — atom range and axes both 1-indexed. **Signature change** (not indexing): the Perl `checkBoundingBox` took a third parameter `$shiftStyle` (0=translation, 1=rotation); the Python version only implements the translation case and documents that the rotation case is handled inside `rotate_all_atoms` differently. The caller `translate_atoms` correctly omits the 3rd arg. Flag for awareness but not in audit scope. |
| 1.2.h | `cut_block(zone, abcxyz_flag, block_border)` | `structure_control.py:5264` | OK — `block_border` is explicitly 1-indexed: `block_border[axis] = [None, low, high]` for axis in 1..3. Axis loops use `range(1, 4)`. Survivor arrays all built with `[None]` sentinel. |
| 1.2.i | `cut_sphere(zone, abcxyz_flag, sphere_rad, target_atom, sphere_loc)` | `structure_control.py:5347` | OK — `target_atom` is 1-indexed; `sphere_loc` is 1-indexed `[None, c1, c2, c3]`. Reads `direct_xyz[target_atom][1..3]` and `fract_abc[atom][1..3]`/`direct_xyz[atom][1..3]`. |
| 1.2.j | `insert_vacuum(vac_axis, vac_amt)` | `structure_control.py:5229` | OK — `vac_axis` is 1-indexed (1=a, 2=b, 3=c) per docstring. |
| 1.2.k | `make_ortho()` | `structure_control.py:5442` | OK — internal `range(1, 4)` loops over abc/xyz axes; `real_lattice[abc][xyz]` access matches Perl. |
| 1.2.l | `apply_sort(order)` | `structure_control.py:5721` | OK on indexing — `order[1..num_atoms]` is 1-indexed, reorders each array at indices 1..num_atoms. **Signature change** (not indexing): the Perl `applySort` was called once per list to sort; the Python version bundles all per-atom arrays internally into a single call. Documented. |
| 1.2.m | `apply_filter(min_dist_filter)` | `structure_control.py:5561` | See Section 1.7.b. |
| 1.2.n | `apply_perturbation(magnitude)` | `structure_control.py:5501` | See Section 1.7.a. |

## Section 1.2 verdict

11 methods OK on indexing, 0 drift. Two methods with documented
Python-side signature changes that are not indexing issues (1.2.g
`check_bounding_box` drops a parameter, 1.2.l `apply_sort` reshapes
to a bundled call). 1.2.m and 1.2.n covered in Section 1.7.
