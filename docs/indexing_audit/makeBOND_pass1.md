# makeBOND.py -- Pass 1 audit

## Scope and method

`src/scripts/makeBOND.py` (~3086 lines) is a Python port of the
Perl `makeBOND` script (~2580 lines).  It parses OLCAO bond
order output, reads a structure file, optionally applies bond
control and group control filters, and emits one of several
visualisation formats: scatter-plot data, openDX model files,
3-D openDX bond order mesh, bond profile curve, or VTK legacy
model files.

Unlike most of the converted scripts audited so far,
`makeBOND.py` does **not** import or use `structure_control.py`.
The original Perl script only called `StructureControl::prepLine`
and `StructureControl::getBohrRad` (two scalar helpers).  The
Python port inlines these (`line.split()` and a hard-coded
`BOHR_RAD` constant).  Pass-1 section "2. sc.* calls" is
therefore N/A for this file.

The script is organised around two classes, `ScriptSettings`
and `BondData`, plus ~20 module-level functions that operate
on a `BondData` instance.  Most of these are output writers;
the structural data arrays are populated in
`get_system_structure`, `read_control_files`, `read_data`,
`implicit_data`, `get_bond_positions`, `get_extra_atom_info`,
and `get_extended_positions`.

Pass-1 approach:

- `Grep` sweep for drift signals: `range(3)`, `range(0, ...)`,
  `[atom - 1]`, `[i + 1]`, `.append([`, `[None]` sentinels.
- Full read of all ~3086 Python lines against the Perl original.
- Side-by-side verification of every per-atom, per-bond,
  per-axis, per-filter, per-limit loop bound and index access.
- Special attention to the rc defaults (`makeBONDrc.py`) and
  their crossing of the 1-indexed boundary.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings | `ScriptSettings.*` | **1 DRIFT** (MBND-A: `num_mesh_points` 0-indexed inner) |
| B | BondData init | `BondData.__init__` | OK (all persistent arrays `[None]` sentineled) |
| C | Structure read | `get_system_structure` | OK |
| D | Control files | `read_control_files` | OK |
| E | Raw data read | `read_data` | **1 minor DRIFT** (MBND-B: `element_z` append skip on "3c") |
| F | Implicit data | `implicit_data` | OK |
| G | Bond midpoints | `get_bond_positions` | OK |
| H | Per-atom stats | `get_extra_atom_info` | OK |
| I | Extended cell | `get_extended_positions` | OK |
| J | Check requirements | `check_requirements` | OK |
| K | Scatter output | `print_scatter_data` | OK |
| L | openDX model | `print_model_data`, `print_lattice_dx`, `create_dx_list`, `add_dx_atom`, `add_dx_bond` | OK |
| M | 3-D mesh | `print_3d_mesh_data`, `make_dx_mesh`, `eval_dx_mesh`, `print_dx_mesh` | **cluster consumer** (MBND-A `num_mesh_points` drift) |
| N | Bond profile | `print_bond_profile` | OK |
| O | VTK output | `print_model_vtk` | OK |
| P | rc defaults | `makeBONDrc.py` | **1 DRIFT** (MBND-A source) |

**Running tally: 0 BUGs, 2 DRIFT findings.**  The only drift
that crosses any consumer boundary is MBND-A (`num_mesh_points`).
MBND-B is a dead-storage drift (the skipped slots are never
read back).  Everything else preserves the 1-indexed Perl
convention faithfully, including the doubly-1-indexed 3x3
lattice tensors (`real_lattice`, `delta`), the
1-indexed-outer-plus-1-indexed-inner bond record arrays
(`bond_atom`, `bond_atom_ext`, `bond_atom_dx`, `bo_limits`,
`bl_limits`), the 1-indexed 3-vectors (`a`, `b`, `c`,
`x_factor`, `y_factor`, `z_factor`), and the per-atom /
per-bond / per-filter / per-group scalar lists.

---

## Finding MBND-A -- `num_mesh_points` inner axis 0-indexed

**Perl** `src/scripts/makeBOND`:

```perl
# default init (lines 464-466)
$numMeshPoints[1] = 10;
$numMeshPoints[2] = 10;
$numMeshPoints[3] = 10;

# CLI parse (lines 557-559)
$numMeshPoints[1] = $ARGV[$number+2];
$numMeshPoints[2] = $ARGV[$number+3];
$numMeshPoints[3] = $ARGV[$number+4];

# makeDXMesh (line 1444)
$numScanPoints = $numMeshPoints[1]*$numMeshPoints[2]*$numMeshPoints[3];
# (line 1451-1457)
$delta[1][$axis] = $realLattice[1][$axis] / ($numMeshPoints[1]-1);
$delta[2][$axis] = $realLattice[2][$axis] / ($numMeshPoints[2]-1);
$delta[3][$axis] = $realLattice[3][$axis] / ($numMeshPoints[3]-1);
# (line 1463-1467)
foreach $aPoint (0..$numMeshPoints[1]-1)
foreach $bPoint (0..$numMeshPoints[2]-1)
foreach $cPoint (0..$numMeshPoints[3]-1)

# printDXMesh (line 1588, 1595, 1601-1605)
print ODX "object 1 class gridpositions counts $numMeshPoints[1]"
      . " $numMeshPoints[2] $numMeshPoints[3]\n";
print ODX "object 2 class gridconnections counts $numMeshPoints[1]"
      . " $numMeshPoints[2] $numMeshPoints[3]\n";
foreach $aPoint (1..$numMeshPoints[1])
foreach $bPoint (1..$numMeshPoints[2])
foreach $cPoint (1..$numMeshPoints[3])
```

Perl `@numMeshPoints` is unambiguously 1-indexed: slots 1, 2,
3 carry the a-, b-, c-axis point counts.

**Python port** -- rc default at `makeBONDrc.py:135`:

```python
"num_mesh_points": [10, 10, 10],
```

Stored into `ScriptSettings.num_mesh_points` at
`makeBOND.py:265` and `:526` via
`list(rc["num_mesh_points"])` / `list(args.num_mesh_points)`.
The `argparse` entry at `makeBOND.py:466-473`:

```python
parser.add_argument(
    '-steps', dest='num_mesh_points',
    nargs=3, type=int,
    default=self.num_mesh_points,
    ...
)
```

Consumers -- all read slots `[0]`, `[1]`, `[2]` instead of the
Perl `[1]`, `[2]`, `[3]`:

- `make_dx_mesh` at `makeBOND.py:2428-2456`:
  ```python
  nmp = s.num_mesh_points
  bd.num_scan_points = nmp[0] * nmp[1] * nmp[2]
  # delta[axis=1] <- real_lattice[1]/(nmp[0]-1)
  bd.delta[1][axis] = bd.real_lattice[1][axis] / (nmp[0] - 1)
  bd.delta[2][axis] = bd.real_lattice[2][axis] / (nmp[1] - 1)
  bd.delta[3][axis] = bd.real_lattice[3][axis] / (nmp[2] - 1)
  for a_pt in range(nmp[0]):
      for b_pt in range(nmp[1]):
          for c_pt in range(nmp[2]):
  ```

- `print_dx_mesh` at `makeBOND.py:2550-2585`:
  ```python
  nmp = s.num_mesh_points
  f.write(f"object 1 class gridpositions counts "
          f"{nmp[0]} {nmp[1]} {nmp[2]}\n")
  f.write(f"object 2 class gridconnections counts "
          f"{nmp[0]} {nmp[1]} {nmp[2]}\n")
  for _ in range(nmp[0]):
      for _ in range(nmp[1]):
          for _ in range(nmp[2]):
  ```

**Severity**: DRIFT (inner axis only).  All consumers are
internally consistent with the 0-indexed storage, so the mesh
output is correct end-to-end.  But the inner layout violates
the Perl 1-indexed convention: Perl's `$numMeshPoints[1]` is
the a-axis, Python's `nmp[0]` is the a-axis.  This is directly
analogous to the already-catalogued `makeinput.py` Cluster
MI-B finding (`kp_mesh_scf`/`kp_mesh_pscf`), where the fix was
to store the array as `[None, na, nb, nc]` and update every
consumer site.

**Recommended fix** (single cluster MBND-A):

1. `makeBONDrc.py:135` -- change default to
   `"num_mesh_points": [None, 10, 10, 10]`.
2. `makeBOND.py:265` / `:526` -- replace
   `list(rc["num_mesh_points"])` and
   `list(args.num_mesh_points)` with
   `[None] + list(...)`.
3. `makeBOND.py:466-473` argparse entry -- either rebuild the
   defaulted display to show three values (not four) and
   still prepend `None` when reading `args.num_mesh_points`,
   or accept three values and prepend on reconcile.  The
   simplest solution is to prepend on `reconcile`.
4. `make_dx_mesh`, `print_dx_mesh` -- change every `nmp[0]` ->
   `nmp[1]`, `nmp[1]` -> `nmp[2]`, `nmp[2]` -> `nmp[3]`.
5. No external caller exists; `num_mesh_points` is never
   passed to another module.

After the fix, `bd.delta[axis]` / `nmp[axis]` are both
1-indexed and match the Perl convention exactly.

---

## Finding MBND-B -- `element_z` append-skip on "3c" atoms

**Perl** `src/scripts/makeBOND:701-718`:

```perl
if (lc($elementName[$i]) eq "3c")  # Fake atom for a 3-center bond
{
   $atomRadius[$i] = 0.25 * $radiusFactor;
   $atomChargeTransfer[$i] = $atomQStar[$i];
   $atomColor[$i] = 100;
}
else
{
   $elementZ[$i] = ElementData::getElementZ(lc($elementName[$i]));
   $atomRadius[$i] = $covalRadii_ref->[$elementZ[$i]] * $radiusFactor;
   ...
}
```

Perl uses direct 1-indexed assignment (`$elementZ[$i]`,
`$atomRadius[$i]`, ...).  For a "3c" atom, `$elementZ[$i]`
remains `undef`, leaving a "hole" in the 1-indexed array at
slot `i`.  The slot alignment is preserved: every per-atom
array still has one entry per atom at the correct index.

**Python** `makeBOND.py:1148-1172`:

```python
if bd.element_name[i].lower() == "3c":
    # Fake atom for a 3-centre bond.
    bd.atom_radius.append(0.25 * s.radius_factor)
    bd.atom_charge_transfer.append(bd.atom_q_star[i])
    bd.atom_color.append(100)
else:
    z = ed.get_element_z(bd.element_name[i].lower())
    bd.element_z.append(z)
    bd.atom_radius.append(ed.coval_radii[z] * s.radius_factor)
    ct = bd.atom_q_star[i]
    for qn_l in range(ed.max_qn_l + 1):
        ct -= ed.vale_charge[z][qn_l]
    bd.atom_charge_transfer.append(ct)
    if not s.grey_scale:
        bd.atom_color.append(ed.color_dx[z])
    else:
        bd.atom_color.append(ed.grey_dx[z])
```

The three "other" arrays (`atom_radius`, `atom_charge_transfer`,
`atom_color`) are appended in **both** branches, so their
indices stay in lock-step with atom number.  But
`bd.element_z` is appended **only** in the else branch.  When
atom `i` is "3c", no append happens, so atom `i+1`'s Z lands
at Python index `i` instead of index `i+1`.  Every subsequent
atom shifts by one slot (by two if another "3c" appears,
etc.).

**Severity**: minor DRIFT (dead storage).  `bd.element_z` is
never read back outside `read_data`:

- grep for `element_z` in `makeBOND.py` finds only three
  hits: the `self.element_z = [None]` init at line 608, the
  `z = ed.get_element_z(...)` computation at line 1158
  (which stores into a **local** `z` variable, not into
  `bd.element_z`), and the `bd.element_z.append(z)` at
  line 1161.
- The 1-indexed-per-atom uses inside the else branch all read
  from the local `z`, not from `bd.element_z[i]`.

So the drift does not produce a runtime error or a wrong
answer in any currently-exercised path.  It is a latent
hazard: any future reader of `bd.element_z[i]` after a "3c"
atom would read the wrong Z.

**Recommended fix** (single line):

```python
if bd.element_name[i].lower() == "3c":
    bd.element_z.append(None)   # <-- add this line
    bd.atom_radius.append(0.25 * s.radius_factor)
    ...
```

This mirrors the Perl `undef` hole at the `[i]` slot while
preserving lock-step 1-indexing with all other per-atom
arrays.

Alternative: delete `self.element_z` entirely since nothing
reads it.  Either fix is acceptable.  The first is safer for
future use.

---

## Non-findings (verified consistent)

### B -- `BondData.__init__`

Every persistent 1-indexed array is initialised with a
`[None]` sentinel:

```python
self.element_name = [None]
self.element_z = [None]
self.atom_q_star = [None]
self.num_filters = [None]
self.filter_name = [None]
self.filter_value = [None]
self.bond_atom = [None]
self.bond_length = [None]
self.bond_order = [None]
self.num_bonds_by_atom = [None]
self.bond_atom_by_atom = [None]
self.bond_length_by_atom = [None]
self.bond_order_by_atom = [None]
self.atom_tag = [None]
self.bond_tag = [None]
self.atom_radius = [None]
self.atom_color = [None]
self.atom_charge_transfer = [None]
self.x_pos = [None];  self.y_pos = [None];  self.z_pos = [None]
self.a_pos = [None];  self.b_pos = [None];  self.c_pos = [None]
...
```

Doubly-1-indexed 3x3 lattice tensors use the nested
`[None, [None, 0, 0, 0], ...]` pattern:

```python
self.real_lattice = [
    None,
    [None, 0.0, 0.0, 0.0],
    [None, 0.0, 0.0, 0.0],
    [None, 0.0, 0.0, 0.0],
]
self.delta = [
    None,
    [None, 0.0, 0.0, 0.0],
    [None, 0.0, 0.0, 0.0],
    [None, 0.0, 0.0, 0.0],
]
self.scan_points = [None, [None], [None], [None]]
```

The 1-indexed 3-vector lattice edges and conversion factors
use the single-`[None]` pattern:

```python
self.a = [None, 0.0, 0.0, 0.0]
self.b = [None, 0.0, 0.0, 0.0]
self.c = [None, 0.0, 0.0, 0.0]
self.x_factor = [None, 0.0, 0.0, 0.0]
self.y_factor = [None, 0.0, 0.0, 0.0]
self.z_factor = [None, 0.0, 0.0, 0.0]
```

All consistent with Perl.

### C -- `get_system_structure`

Perl reads `$a[0..2]` from the line, rotates into 1-indexed
by `$a[1..3] = $a[0..2]; $a[0] = ""`, then computes
`$magA = sqrt($a[1]*$a[1] + $a[2]*$a[2] + $a[3]*$a[3])`.
Python writes directly into slots 1-3:

```python
for i in range(3):
    bd.a[i + 1] = float(vals[i]) * BOHR_RAD
bd.mag_a = math.sqrt(bd.a[1]**2 + bd.a[2]**2 + bd.a[3]**2)
```

Then `bd.real_lattice[1][1] = bd.a[1]; ...` through
`bd.real_lattice[3][3] = bd.c[3]`.  Matches Perl exactly.

Atom positions: Perl writes `$xPos[$i] = $values[3]*$bohrRad`
directly.  Python, starting from `[None]`, appends in order:

```python
for i in range(1, bd.num_atoms + 1):
    vals = pos.readline().split()
    x = float(vals[3]) * BOHR_RAD
    ...
    bd.x_pos.append(x)   # lands at index i since list was [None]
```

Because the loop walks `1..num_atoms` and each iteration
appends exactly one value, atom `i`'s x-position lands at
index `i` in the 1-indexed list.  Matches Perl.

### D -- `read_control_files`

`bo_limits` / `bl_limits` outer 1-indexed via append from
`[None]`; inner entries built as
`[None, elem1, elem2, minv, maxv]` -- doubly 1-indexed.
Matches Perl `$BOLimits[$i][1..4]` / `$BLLimits[$i][1..4]`.

`group_name` 1-indexed via append from `[None]`; both the
default-path `group_name.append("ELEMENT_NAME")` and the
file-reading `for i, name in enumerate(vals): ... append(name)`
produce a 1-indexed list.  `bd.num_groups = len(vals)` (Python)
equals Perl's `$numGroups = $#groupName` (last index of
1-indexed = count of real entries).

### E -- `read_data` main body (minus MBND-B)

First pass over the data file populates a temporary
`bd.element_name` via append; the list is then deliberately
overwritten by a re-allocation `[None] * (num_atoms + 1)`
before the per-atom loop runs:

```python
bd.element_name = [None] * (bd.num_atoms + 1)
for i in range(1, bd.num_atoms + 1):
    ...
    if vals[0] == "ELEMENT_NAME":
        bd.element_name[i] = vals[1]
```

This re-allocation ensures slot `i` is reserved for atom `i`
regardless of element-name order or duplicates.  Matches
Perl.

Per-atom `num_filters`, `filter_name`, `filter_value`
initialised via `append(0)` / `append([None])` / `append([None])`
in the pre-loop, then filled in the per-atom loop with
`filter_name[i].append(name)` / `filter_value[i].append(val)`,
building a doubly-1-indexed layout matching Perl
`$filterName[$i][$k]` / `$filterValue[$i][$k]`.

`num_bonds_by_atom`, `bond_atom_by_atom`, `bond_length_by_atom`,
`bond_order_by_atom` are similarly 1-indexed via
`[None]`-append pattern per atom.

Bond record building:

```python
idx = j + bd.num_bonds_by_atom[i] - num_bonds_temp
while len(bd.bond_atom_by_atom[i]) <= idx:
    bd.bond_atom_by_atom[i].append(None)
    ...
bd.bond_atom_by_atom[i][idx] = bonded_atom
```

Direct port of Perl
`$bondAtomByAtom[$i][$j + $numBondsByAtom[$i] - $numBondsTemp]`
with explicit Python-side grow-to-size (Perl autovivifies).
1-indexed correct because `j` walks `1..num_bonds_temp` and
`idx` stays >= 1.

Per-bond storage:

```python
bd.bond_atom.append([None, i, bonded_atom])
bd.bond_length.append(b_length)
bd.bond_order.append(b_order)
```

`bond_atom` outer 1-indexed (list started with `[None]`),
inner 1-indexed (`[None, a1, a2]`).  Matches Perl
`$bondAtom[$numSystemBonds][1..2]`.

Reject path (bond fails `check_requirements`):

```python
bd.num_bonds_by_atom[i] -= 1
bd.num_system_bonds -= 1
bd.bond_atom.pop()
bd.bond_length.pop()
bd.bond_order.pop()
bd.num_bonds_by_atom[bonded_atom] -= 1
```

`pop()` removes the just-appended slot, mirroring Perl's
implicit shrink on decrement.  OK.

### F -- `implicit_data`

Outer loops walk `range(1, bd.num_atoms + 1)` and
`range(1, bd.num_system_bonds + 1)`.  Inner group/filter
loops walk `range(1, bd.num_groups + 1)` and
`range(1, bd.num_filters[i] + 1)`.  All 1-indexed, matching
Perl `for ($i=1;$i<=$numAtoms;$i++)` etc.

`unique_atom_tags` / `unique_bond_tags` initialised to
`[None, atom_tag[1]]` then extended with subsequent unique
tags starting from index 2.  Matches Perl's
`$uniqueAtomTag[1] = $atomTag[1]; ... $uniqueAtomTag[++$count]`.

### G -- `get_bond_positions`

Per-bond loop `range(1, bd.num_system_bonds + 1)` reads
`bd.bond_atom[i][1]` and `bd.bond_atom[i][2]` as 1-indexed
atom refs and writes `bd.a_bond_pos.append(a_mid)` (etc.)
from the `[None]`-initialised list so idx `i` holds bond
`i`'s midpoint.  Matches Perl.

### H -- `get_extra_atom_info`

Outer loop `range(1, bd.num_atoms + 1)`, inner bond loop
`range(1, bd.num_bonds_by_atom[i] + 1)`, all 1-indexed reads
and writes.  Init arrays via append from `[None]` so atom
`i` stats land at idx `i`.  Matches Perl.

### I -- `get_extended_positions`

Initialisation of extended arrays as `list(bd.x_pos)` etc.
copies the `[None]` sentinel and `num_atoms` entries; the
resulting `bd.x_pos_ext` is 1-indexed `[None, p1, p2, ...]`.
Matches Perl `@xPosExt = @xPos` (1-indexed outer).

`bd.bond_atom_ext` built via an explicit loop:

```python
bd.bond_atom_ext = [None]
for i in range(1, bd.num_system_bonds + 1):
    bd.bond_atom_ext.append(
        [None, bd.bond_atom[i][1], bd.bond_atom[i][2]]
    )
```

Doubly 1-indexed, matching Perl `$bondAtomExt[$i][1..2]`.

Local scratch vectors `min_dist`, `min_x`, `min_y`, `min_z`
are 1-indexed `[None, ...]` triples with only slots 1 and 2
used (one per bond end):

```python
min_dist = [None, math.sqrt(...), 0.0]
min_dist[2] = min_dist[1]
min_x = [None, 0.0, 0.0]
min_y = [None, 0.0, 0.0]
min_z = [None, 0.0, 0.0]
```

Matches Perl `$minDistance[1]`, `$minX[1]`, `$minY[1]` etc.

The periodic-image scan uses `for m in (1, 2)` explicitly,
not `range(2)`, preserving the 1-indexed `m` variable for
`bd.bond_atom[i][m]` access.  Matches Perl
`foreach $m (1..2)`.

The `for j_idx in (1, 2)` loop that writes new extended
atoms reads `min_x[j_idx]` etc. with the 1-indexed scratch
refs.  Matches Perl `for ($j=1;$j<=2;$j++)`.

All the subsequent bond-midpoint computations are 1-indexed
reads and writes that exactly mirror the Perl block at
`makeBOND:2122-2180`.

One oddity worth noting for cleanup (**not an audit
finding**): lines 1635-1640 first compute `tz` using
`bd.x_pos[ba]` (a typo that was ported forward verbatim),
then immediately overwrite it with the correct
`bd.z_pos[ba]` expression at lines 1643-1648.  The Python
port flags the fix in a comment:

```python
tz = (bd.x_pos[ba] + j * bd.a[3] + k * bd.b[3] + l * bd.c[3])
# Correctly compute tz using z_pos, not x_pos.
tz = (bd.z_pos[ba] + j * bd.a[3] + k * bd.b[3] + l * bd.c[3])
```

Perl at `makeBOND:2032` has no such typo
(`$testZ[$m] = $zPos[$bondAtom[$i][$m]] + ...`).  The dead
first assignment is a cosmetic dust-devil, not an indexing
issue.

### J -- `check_requirements`

`n = bd.num_system_bonds` is the current (just-appended) bond
index.  Reads `bd.bond_length[n]`, `bd.bond_order[n]`,
`bd.bond_atom[n][1..2]` as 1-indexed.  Sweeps
`bd.bo_limits[i][1..4]` and `bd.bl_limits[i][1..4]` via
`for i in range(1, bd.num_bo_limits + 1)` and
`for i in range(1, bd.num_bl_limits + 1)`.  Returns a Python
`bool` instead of writing Perl's `$requirementsPassed`; the
read_data caller tests the return value directly.  Matches
Perl.

### K -- `print_scatter_data`

Outer loops `range(1, num_unique_atom + 1)` and
`range(1, num_unique_bond + 1)` (where `num_unique_atom = len(
unique_atom_tags) - 1`, i.e. the count of real entries past
the `[None]` sentinel).  Inner loops walk
`range(1, bd.num_atoms + 1)` and
`range(1, bd.num_system_bonds + 1)` reading 1-indexed
scalars/pairs.  Matches Perl.

### L -- openDX model family

`print_lattice_dx` reads `bd.a[1..3]`, `bd.b[1..3]`,
`bd.c[1..3]` -- 1-indexed triples.  OK.

`print_model_data` outer loop walks
`range(1, num_unique_tags + 1)`.  Inner atom-write loops walk
`range(1, bd.num_atoms_dx + 1)` and
`range(1, bd.num_atoms_ext_dx + 1)`, reading
`bd.x_pos_ext_dx[j]`, `bd.atom_radius_dx[j]`, etc.  1-indexed
throughout.

The openDX connection output subtracts 1 at the file-format
boundary:

```python
for j in range(1, bd.num_bonds_dx + 1):
    # Subtract 1 because openDX is C-based
    #   and array refs start from 0.
    fa = bd.bond_atom_dx[j][1] - 1
    sa = bd.bond_atom_dx[j][2] - 1
    f.write(f"{fa}  {sa}\n")
```

This is a legitimate 1-indexed-to-0-indexed bridge at the
output boundary (analogous to `bond_analysis.py`'s VTK/OpenDX
index translations -- see `bond_analysis_pass1.md` section
"Non-findings").  OK.

`create_dx_list` allocates `bd.atom_map_2_dx = [None] * (
bd.num_atoms_ext + 1)` -- 1-indexed, slot 0 reserved.
`bd.x_pos_ext_dx`, `bd.y_pos_ext_dx`, `bd.z_pos_ext_dx`,
`bd.atom_color_dx`, `bd.atom_radius_dx`,
`bd.atom_charge_transfer_dx` all reset to `[None]` sentinel.
`bd.bond_atom_dx`, `bd.bond_order_dx`, `bd.abc_bond_pos`,
`bd.xyz_bond_pos` likewise reset to `[None]`.  Outer sweeps
use `range(1, bd.num_atoms_ext + 1)` and
`range(1, bd.num_bonds_ext + 1)`.  Matches Perl.

`add_dx_atom` uses the 1-indexed `atom_number` parameter to
read `bd.x_pos_ext[atom_number]` etc. and writes
`bd.atom_map_2_dx[atom_number] = bd.num_atoms_ext_dx`.  The
per-atom DX arrays grow by append, so the new atom's data
lands at index `num_atoms_ext_dx` (which was just
incremented).  Matches Perl exactly.

The "extend the map if needed" guard loop at lines 2306-2307
uses `while len(bd.atom_map_2_dx) <= atom_number:
bd.atom_map_2_dx.append(None)` as a defensive padding in case
the atom_number exceeds the initial allocation -- harmless.

`add_dx_bond` builds doubly-1-indexed entries for
`bd.bond_atom_dx`, `bd.abc_bond_pos`, `bd.xyz_bond_pos`:

```python
bd.bond_atom_dx.append([
    None,
    bd.atom_map_2_dx[bd.bond_atom_ext[bond_number][1]],
    bd.atom_map_2_dx[bd.bond_atom_ext[bond_number][2]],
])
bd.abc_bond_pos.append([
    None,
    bd.a_bond_pos[bond_number],
    bd.b_bond_pos[bond_number],
    bd.c_bond_pos[bond_number],
])
bd.xyz_bond_pos.append([
    None,
    bd.x_bond_pos[bond_number],
    bd.y_bond_pos[bond_number],
    bd.z_bond_pos[bond_number],
])
```

Inner pattern `[None, v1, v2, v3]`, outer append from `[None]`
initial.  Matches Perl
`$bondAtomDX[$numBondsDX][1..2]` / `$abcBondPos[$numBondsDX][1..3]`
/ `$xyzBondPos[$numBondsDX][1..3]`.

### M -- 3-D mesh family (minus MBND-A consumer sites)

`print_3d_mesh_data`, `make_dx_mesh`, `eval_dx_mesh`,
`print_dx_mesh` use `bd.delta[1..3][1..3]`,
`bd.real_lattice[1..3][1..3]`, `bd.scan_points[1..3][1..N]`,
`bd.mesh_bo_dx[1..N]`, `bd.xyz_bond_pos[bond][1..3]`.  Each
is doubly-1-indexed, matching Perl.

Per-mesh-point contribution lists are 1-indexed with `[None]`
sentinel:

```python
gaussian_dist_factor = [None]
contrib_bond_number = [None]
...
gaussian_dist_factor.append(...)   # slot 1, 2, ...
contrib_bond_number.append(bond)
```

Consumer loop walks
`for cb in range(1, num_contrib_bonds + 1)`.  Matches Perl
`foreach $contribBond (1..$numContribBonds)`.

`bd.mesh_bo_dx = [None] * (bd.num_scan_points + 1)` -- 1-indexed.

Only the `nmp[0]`/`nmp[1]`/`nmp[2]` accesses violate the
convention, and those are all part of Cluster MBND-A above.

### N -- `print_bond_profile`

Inner scratch arrays preserve Perl's 0-indexed start (slot 0
is the sentinel `-0.01` for numerical-accuracy reasons as
explained in the Perl comment):

```python
axis_points = [-0.01]
cumul_bond_order = [0.0]
cumul_num_bonds = [0]
for i in range(num_points):
    axis_points.append(i * 0.01)
    cumul_bond_order.append(0.0)
    cumul_num_bonds.append(0)
```

This builds three length-`(num_points + 1)` lists where slot
0 = -0.01 and slots 1..num_points = 0.00, 0.01, ...,
`(num_points-1)*0.01`.  Matches Perl exactly:

```perl
$axisPoints[0] = -0.01;
for ($i=0;$i<=$numPoints-1;$i++) {
   $axisPoints[$i+1] = $i * 0.01;
   $cumulBondOrder[$i+1] = 0;
   $cumulNumBonds[$i+1]  = 0;
}
```

`atom_index = [None]` then per-atom append via
`for i in range(1, bd.num_atoms_ext + 1)` -- 1-indexed.
Matches Perl `$atomIndex[$i]`.

Bond-cross loop:

```python
for i in range(1, bd.num_bonds_ext + 1):
    left = atom_index[bd.bond_atom_ext[i][1]]
    right = atom_index[bd.bond_atom_ext[i][2]]
```

1-indexed bond, 1-indexed inner.  Matches Perl.

Output loop `for i in range(1, num_points)` skips slot 0
(the -0.01 sentinel) and walks 1..num_points-1.  Matches Perl
`for ($i=1;$i<=$numPoints-1;$i++)`.

### O -- `print_model_vtk`

Same structure as `print_model_data` with its own file
format.  The VTK connection output subtracts 1 at the
boundary -- legitimate 0-indexed VTK bridge:

```python
for j in range(1, bd.num_atoms_ext_dx + 1):
    ...
    k = j - 1
    f.write(f"1 {k}\n")
...
for j in range(1, bd.num_bonds_dx + 1):
    fa = bd.bond_atom_dx[j][1] - 1
    sa = bd.bond_atom_dx[j][2] - 1
    f.write(f"2 {fa} {sa}\n")
```

Matches Perl.  All input reads are 1-indexed.  OK.

### P -- `makeBONDrc.py`

Only one list-valued setting: `num_mesh_points` (see
MBND-A).  Every other key is a scalar (str, float, int,
bool) that cannot carry an indexing convention.
`bond_control_file`, `group_control_file`, `bond_profile_file`
are file-path strings -- no indexing concern.

---

## Non-indexing observations (out of audit scope)

### Dead first `tz` assignment in `get_extended_positions`

At `makeBOND.py:1635-1648`, the first `tz = ...` uses
`bd.x_pos[ba]` (a typo), then the same `tz = ...` is
immediately recomputed with the correct `bd.z_pos[ba]`.  The
first assignment is dead code.  Perl has no corresponding
bug at `makeBOND:2032` -- the Python port evidently wrote the
wrong line, noticed the mistake, and left the old line in
place with a "correctly compute" comment rather than deleting
it.  Cosmetic cleanup task, not an indexing issue.

### Unused `dipole` flag

`ScriptSettings.dipole` is wired through argparse and
reconcile but never consulted in `main()`'s dispatch.  The
Perl original also had `$dipole` but I did not trace whether
it was used; the Python port is in any case a faithful port
of the parse-and-ignore behaviour.  Not an indexing issue.

### Verbose progress prints in `eval_dx_mesh`

`print(num_contrib_bonds)` and `print(f"NO BONDS {point}")`
are left in from the Perl original (which had them as
`print STDOUT` lines).  These will emit one line per mesh
point -- noisy for a 10x10x10 mesh (1000 lines).  Matches
Perl behaviour, not an indexing issue.

### `element_z` never read

See MBND-B above: `bd.element_z` is a latent dead-store that
could be deleted entirely.  The proposed MBND-B fix (append
`None` on "3c") is a safer choice if the field may be read
in a future addition, but deleting the field is equally
acceptable.

---

## Running tally

- **BUG**: 0
- **DRIFT**: 2
  - **MBND-A**: `num_mesh_points` inner axis (rc default and
    all consumers in `make_dx_mesh` / `print_dx_mesh`).
    Analogous to `mi_pass1.md` Cluster MI-B.  ~6 touch points
    in `makeBOND.py` plus 1 in `makeBONDrc.py`.
  - **MBND-B**: `element_z` append-skip on "3c" atoms.
    Latent (the field is never read).  Single-line fix or
    delete the field.
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every other subsystem, covering
  `ScriptSettings`, `BondData.__init__`, structure read,
  control files, raw data, implicit data, bond midpoints,
  per-atom stats, extended cell, bond profile, all three
  openDX variants, and VTK.

## Fix-order recommendation

1. **Cluster MBND-A** -- `num_mesh_points` inner axis.
   Single coordinated edit across `makeBONDrc.py`,
   `ScriptSettings.assign_rc_defaults`,
   `ScriptSettings.reconcile`, `make_dx_mesh`, and
   `print_dx_mesh`.  Six consumer sites total.  No external
   callers; no test coverage (regression verification by
   inspection + full `structure_control` test suite, which is
   unaffected since `makeBOND.py` does not call into
   `structure_control`).
2. **Cluster MBND-B** -- add `bd.element_z.append(None)` to
   the "3c" branch of `read_data`, or delete the field
   entirely.  Single-line change, zero risk.

No dependency between the two clusters; they can be fixed in
either order.  Apply MBND-A first since it is the only
finding that crosses a real consumer boundary.

## Fix status

**Not yet applied.**  This audit file is Pass 1 only.  Apply
the two clusters in a follow-up commit per the recommendation
above.
