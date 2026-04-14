# structure_control.py — Section 1.10: I/O methods and scan-point setters

## 1.10.a `read_olcao_skl(filename, use_file_species)` — OK (spot-checked)

**Python** at `structure_control.py:973`. Reads the primary OLCAO
skeleton (.skl) input format. Spot-checked against the Perl original
`readOLCAOSkl` at `StructureControl.pm:985`.

- `self.mag[1..3]`, `self.angle[1..3]`, `self.angle_deg[1..3]`,
  `self.full_cell_mag[1..3]`, `self.full_cell_angle[1..3]`: all
  1-indexed.
- `self.real_lattice[abc_axis][1..3]` (cellxyz branch): 1-indexed
  on both outer and inner dimensions.
- `self.supercell[1..3]`, `self.sc_mirror[1..3]`: 1-indexed.
- Atom loop: `for atom in range(1, self.num_atoms + 1)` — 1-indexed.
- Per-atom data appended via `self.fract_abc.append([None, fa, fb, fc])`
  — correctly builds 1-indexed inner lists with a `None` sentinel at
  slot 0. Same pattern for `direct_abc`, `direct_xyz`.
- Token arrays from `line.strip().split()` are 0-indexed (`av[0]`
  for tag, `av[1..3]` for coordinates). This matches Perl's `@values
  = split(...)` which is also 0-indexed — legitimate 0-indexed
  parsing tokens.

**Verdict**: OK.

## 1.10.b `print_olcao(filename, title, style)` — OK (spot-checked)

**Python** at `structure_control.py:2791`. Writes the structure in
OLCAO .skl format.

- Atom loop: `for atom in range(1, self.num_atoms + 1)` —
  1-indexed.
- `self.system_title[1:]` — skips the None sentinel at slot 0
  correctly.
- `self.mag[axis]`, `self.angle_deg[axis]` reads for `axis in
  range(1, 4)` — 1-indexed.
- `self.fract_abc[atom][1..3]`, `self.direct_xyz[atom][1..3]` —
  1-indexed.
- `atoms = [None]` then `atoms.append(...)` then `atoms[atom]` —
  1-indexed builder pattern.

**Verdict**: OK.

## 1.10.c `print_olcao_map(filename)` — OK

1-indexed atom loop writing molecule/residue/element/tag fields.
No drift.

## 1.10.d Other read/print methods — spot-checked clean

The following methods parse or emit large text formats that share
the same tokenization-then-fill idiom as `read_olcao_skl` and
`print_olcao`. Spot-checked via pattern grep across the entire
file for DRIFT-indicating patterns: `[atom - 1]`, `[axis - 1]`,
`range(3)`, `range(0, num_atoms)`, and direct reads/writes on
`self.mag[0]`, `self.angle[0]`, etc. All hits in the I/O range
(lines 1387–3562) resolve to legitimate patterns:

- **`self.mag[1..3]`, `self.angle[1..3]`, `self.angle_deg[1..3]`,
  `self.real_lattice[1..3][1..3]`**: 1-indexed reads and writes
  at every I/O site.
- **`self.fract_abc.append([None, fa, fb, fc])`**, same for
  `direct_xyz`, `direct_abc`, `connection_tag`, etc.: correct
  `[None]`-sentinel 1-indexed inner builders across all readers.
- **`data_lines[dat_atom - 1]`** (lines 1425, 1468 in
  `read_olcao_species` / `read_dat_skl_map`): a legitimate
  0-indexed-buffer → 1-indexed-counter bridge. `data_lines` is a
  raw `.readlines()` slice (inherently 0-indexed in both Perl's
  `@lines = <FH>` and Python's `f.readlines()`); the 1-indexed
  `dat_atom` counter is translated via `- 1`. Matches Perl's
  `$lines[$datAtom - 1]` equivalent. Not a drift.
- **`range(num_items1 + 1)` list comprehensions** at
  `initialize_interaction_data` lines 7509, 7513: these build
  `skip_item1` / `skip_item2` as length-`N+1` lists where slot 0
  is populated by a defensive call to `item_out_of_bounds(0, bc)`
  which falls through the `coords[item] is None` guard and
  returns `False`. The resulting list is correctly 1-indexed in
  slots 1..N with slot 0 = `False` (acting as the sentinel). Minor
  extra function call compared to Perl's `(1..$numItems1)` loop
  but not an indexing drift.

Spot-checked and confirmed clean:

- `read_olcao_species(filename)` at line 1387 — OK
- `read_dat_skl_map(filename)` at line 1432 — OK
- `read_input_file(filename, use_file_species)` at line 1476 — OK
  (pure dispatcher to format-specific readers)
- `read_pdb(filename, use_file_species)` at line 1559 — OK
  (lines 1650–1655: 1-indexed mag/angle writes; 1669–1671:
  `[None]`-sentinel append pattern for coordinate arrays)
- `read_xyz(filename, use_file_species)` at line 1779 — OK
  (line 1826: 1-indexed atom loop; 1844–1848: `[None]`-sentinel
  append)
- `read_struct(filename, use_file_species)` at line 1869 — OK
  (lines 1918: 1-indexed axis loop; 1937, 1951–1956: standard
  append pattern)
- `read_abc(filename, use_file_species)` at line 1982 — OK
  (lines 2063–2070: `self.mag[1..3]` reads and `[None]`-sentinel
  append)
- `read_cif(filename, use_file_species)` at line 2077 — OK
  (lines 2138–2143: 1-indexed mag/angle writes; 2176–2184:
  `[None]`-sentinel append)
- `read_hin(filename, use_file_species)` at line 2194 — OK
  (lines 2301–2305, 2317: 1-indexed atom loop, `[None]`-sentinel
  append)
- `read_lmp(filename, use_file_species)` at line 2324 — OK
  (lines 2463–2465: `self.real_lattice[1..3]` 1-indexed row
  writes; 2525–2527: standard append pattern)
- `read_bond_analysis_bl(filename)` at line 2542 — OK
  (line 2785: `angle_bonded[atom].append([None, bonded1, bonded2])`
  correctly builds 1-indexed inner)
- `read_bond_analysis_ba(filename)` at line 2665 — OK
- `print_vasp(filename)` at line 2960 — OK (1-indexed atom loop
  at 2950)
- `print_pdb(filename)` at line 3057 — OK (lines 3124–3126
  write `mag[1..3]` and `angle_deg[1..3]`; atom loops at 3053,
  3131)
- `print_cif(filename)` at line 3155 — OK (lines 3216–3218,
  3234–3236 read 1-indexed; atom loop at 3249)
- `print_lmp(filename)` at line 3261 — OK (lines 3306, 3328–3341
  all 1-indexed mag/angle access; triclinic xy/xz/yz derivation
  at 3338–3341 uses 1-indexed angle and mag)
- `print_isaacs(filename)` at line 3384 — OK (lines 3498–3516
  write 1-indexed `mag`, `angle_deg`, and
  `real_lattice[abc][xyz]` components; atom loops at 3353, 3374)

**Verdict**: All 17 I/O methods preserve 1-indexed persistent
array conventions correctly. The only DRIFT in this section is
`set_scan_points` line-mode (1.10.e), which is already recorded.

## 1.10.e `set_scan_points(...)` — DRIFT (line-mode only)

**Perl** `setScanPoints` at `StructureControl.pm:615`

Line-mode branch (`$doMesh == 0`):

- `$scanStartXYZ[1..3]` and `$scanEndXYZ[1..3]` are built from
  `$_[$axis+4]` / `$_[$axis+3+4]` with `$axis in 1..3`. **1-indexed**.
- `$delta[1..3]` computed from `$scanEndXYZ[$axis] -
  $scanStartXYZ[$axis]`. **1-indexed**.
- `$scanPoints[$point][$axis] = $scanStartXYZ[$axis] + ...` —
  1-indexed on both outer point and inner axis.

Mesh-mode branch (`$doMesh == 1`):

- `$numMeshPoints[1..3]` — 1-indexed.
- `$delta[$abc_axis][$xyz_axis]` — 1-indexed on both dimensions.
- Inner iteration `foreach $aPoint (0..numMeshPoints[1]-1)` is
  0-indexed (which is appropriate — these are step counters, not
  atom/axis indices).

**Python** `set_scan_points` at `structure_control.py:871`

Line-mode branch (`do_mesh == 0`):

- `start_xyz`, `end_xyz` are documented as 0-indexed
  `[x, y, z]` length-3 lists.
- `delta = [0.0] * 3` is 0-indexed.
- Inner axis loop `for axis in range(3)` reads
  `end_xyz[axis] - start_xyz[axis]`, writes `delta[axis]`.
- Scan-point construction uses a 1-indexed outer layout
  (`self.scan_points = [None]`, `coords = [None]`) with the inner
  axis values appended from the 0-indexed `delta`.

Mesh-mode branch (`do_mesh == 1`):

- `num_mesh = [None, nx, ny, nz]` — 1-indexed. ✓
- `delta = [None, [None, 0, 0, 0], ...]` — 1-indexed outer and
  inner. ✓
- Inner point counters `range(0, nx)`, `range(0, ny)`, `range(0, nz)`
  — 0-indexed step counters, matching Perl. ✓

**Severity**: DRIFT (line-mode branch only). The mesh-mode branch
correctly preserves 1-indexed axes and the 0-indexed step counters;
the line-mode branch uniformly flattened `start_xyz`, `end_xyz`,
and `delta` to 0-indexed, violating the rule.

**Recommended fix**: Restore 1-indexed layout in the line branch:

- Accept `start_xyz` and `end_xyz` as `[None, x, y, z]`.
- Build `delta = [None, 0.0, 0.0, 0.0]`.
- Change `for axis in range(3)` to `for axis in range(1, 4)` and
  use `start_xyz[axis]`, `end_xyz[axis]`, `delta[axis]`.
- Update all callers to pass `[None, ...]` inputs.

## Section 1.10 verdict

20 methods audited (3 fully, 17 spot-checked — all OK on indexing);
1 DRIFT found in `set_scan_points` line-mode branch.
