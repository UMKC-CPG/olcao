# condense.py — Pass 1 audit

## Scope and method

`src/scripts/condense.py` (~2426 lines) is a Python port of the
Perl `condense` script (~2156 lines).  It builds LAMMPS input
files for condensation simulations by:

1. Reading a keyword-driven input file (composition, stages,
   reactions).
2. Copying reaction templates from the precursor database.
3. Running `packmol` to pack molecules into a cubic cell.
4. Calling `pdb2skl.py` and `bond_analysis.py` as subprocesses.
5. Building `lammps.dat`, `lammps.in`, and a slurm submission
   file, with per-atom, per-bond, per-angle type normalization
   across the main data file and every pre/post reaction
   template.

The file has three classes (`BondData`, `AngleData`,
`ScriptSettings`) plus one workhorse class `Condense` with
methods for each workflow step.  Total ~15 methods under
audit.

Pass 1 approach:

- Read every method end-to-end.  Most drift signals are in
  nested-list *inner* layouts (storing ≥2-level structures).
- Compare against the Perl `src/scripts/condense` source
  directly — especially the `hookeBondCoeffs`,
  `hookeAngleCoeffs`, `orderedBondedAtoms`, `angleBondedAtoms`,
  `uniqueBondTypes`, `uniqueAngleTypes`, `uniqueBondCoeffs`,
  and `uniqueAngleCoeffs` arrays which are all 1-indexed on
  every axis in Perl.
- `Grep` sweep for `range(3)`, `range(0, …)`, `[atom - 1]`,
  `[i - 1]`, `[0] +=` style patterns.  **No matches** — every
  drift in this file is purely in the inner layouts of nested
  lists, never as an axis-offset bridge or a bare `range(3)`.

## Summary table

| # | Topic | Perl layout | Python layout | Verdict |
|---|-------|-------------|---------------|---------|
| C-A | `BondData.hooke_bond_coeffs` inner | `[1..4]` (Z1,Z2,k,r0) | `[0..3]` | DRIFT (inner axis) |
| C-B | `AngleData.hooke_angle_coeffs` inner | `[1..6]` (Z1,Zv,Z2,k,angle,tol) | `[0..5]` | DRIFT (inner axis) |
| C-C | `ordered_bonded_atoms` inner | `[1..2]` | `[0..1]` | DRIFT (inner axis) |
| C-D | `angle_bonded_atoms` inner | `[1..3]` | `[0..2]` | DRIFT (inner axis) |
| C-E | `ordered_species_pair_coeffs` inner | `[1..2]` | `[0..1]` | DRIFT (inner axis, inherited from `element_data.lj_pair_coeffs`) |
| C-F | `unique_bond_coeffs` inner | `[1..2]` | `[0..1]` | DRIFT (inner axis) |
| C-G | `unique_angle_coeffs` inner | `[1..2]` (or `[1..3]` in normalize_types) | `[0..1]` (or `[0..2]`) | DRIFT (inner axis) |
| C-H | `unique_bond_types` inner | `[1..2]` | `[0..1]` | DRIFT (inner axis) |
| C-I | `unique_angle_types` inner | `[1..4]` | `[0..3]` | DRIFT (inner axis) |
| — | `BondData`/`AngleData` file parse | `@values` 0-indexed | `values[0..]` 0-indexed | OK (legitimate file-token buffer) |
| — | `normalize_types` token field indices | `@vals[3..16]` 0-indexed | `vals[3..16]` 0-indexed | OK (legitimate LAMMPS line parse) |
| — | Stage / composition / reaction outer arrays | 1-indexed `[None, …]` | same | OK |
| — | Atom / molecule / family / num-species outer arrays | 1-indexed | same | OK |
| — | Call boundaries into `StructureControl` | 1-indexed coord vectors, scalar atom indices | same | OK |

Nine DRIFT findings, every single one on an **inner-axis**
flattening.  No runtime bugs (producer/consumer pairs are
self-consistent so the script functions correctly).  No outer
array drift, no `range(3)` bridges, no off-by-one hazards.

## Finding details

Each finding lists the producer site, the consumer site(s), and
the Perl equivalent.  All file:line references below are for
`src/scripts/condense.py` unless stated otherwise.

### C-A: `BondData.hooke_bond_coeffs` inner list

**Perl** (`src/scripts/condense` lines 887–897 via the
`BondData.pm` module):

```perl
foreach $hookeBond (1..$numHookeBonds)
{
   if (($atom1Z == $hookeBondCoeffs_ref->[$hookeBond][1]) and
       ($atom2Z == $hookeBondCoeffs_ref->[$hookeBond][2]))
   {
      $uniqueBondCoeffs[$id][1] = $hookeBondCoeffs_ref->[$hookeBond][3];
      $uniqueBondCoeffs[$id][2] = $hookeBondCoeffs_ref->[$hookeBond][4];
   }
}
```

`$hookeBondCoeffs_ref->[$hookeBond][1..4]` — slot 1 = Z1,
slot 2 = Z2, slot 3 = k, slot 4 = r0.  Fully 1-indexed inner.

**Python** `BondData.init_bond_data` at `condense.py:230–237`:

```python
for bond in range(1, self.num_hooke_bonds + 1):
    values = self._prep_line(f)
    self.hooke_bond_coeffs.append([
        int(values[0]),    # Z of first element
        int(values[1]),    # Z of second element
        float(values[2]),  # k spring constant
        float(values[3]),  # Rest bond length
    ])
```

The outer `self.hooke_bond_coeffs = [None]` is 1-indexed on the
bond index (good).  The inner entry is a plain 4-element
0-indexed list.

**Consumers**:

- `Condense.create_lammps_files` at `condense.py:1191–1198`:
  ```python
  for hb in range(1, bd.num_hooke_bonds + 1):
      hbc = bd.hooke_bond_coeffs[hb]
      if atom1_z == hbc[0] and atom2_z == hbc[1]:
          ...
          unique_bond_coeffs[found] = [hbc[2], hbc[3]]
  ```
- `Condense.normalize_types` at `condense.py:1987–1993`:
  ```python
  for hb in range(1, bd.num_hooke_bonds + 1):
      hbc = bd.hooke_bond_coeffs[hb]
      if (((atom1_z == hbc[0] and atom2_z == hbc[1]) or
           (atom2_z == hbc[0] and atom1_z == hbc[1]))):
          unique_bond_coeffs[ub] = [hbc[2], hbc[3]]
  ```

**Severity**: DRIFT (inner axis only).  Self-consistent; the
only reason to fix it is rule compliance.

**Recommended fix**: change the builder to
`self.hooke_bond_coeffs.append([None, Z1, Z2, k, r0])` and
change all consumers to read `hbc[1..4]`.

### C-B: `AngleData.hooke_angle_coeffs` inner list

**Perl** `src/scripts/condense:944–951`:

```perl
foreach $hookeAngle (1..$numHookeAngles)
{
   if (($atom1Z == $hookeAngleCoeffs_ref->[$hookeAngle][1]) and
       ($atom2Z == $hookeAngleCoeffs_ref->[$hookeAngle][2]) and
       ($atom3Z == $hookeAngleCoeffs_ref->[$hookeAngle][3]))
   {
      if ((abs($bondAnglesExt_ref->[$atom][$angle]
               - $hookeAngleCoeffs_ref->[$hookeAngle][5]))
          <= $hookeAngleCoeffs_ref->[$hookeAngle][6])
         ...
```

`$hookeAngleCoeffs_ref->[$hookeAngle][1..6]` — slots 1..6 are
Z1, Zv, Z2, k, angle_deg, tolerance.  Fully 1-indexed inner.

**Python** `AngleData.init_angle_data` at `condense.py:321–330`:

```python
for angle in range(1, self.num_hooke_angles + 1):
    values = self._prep_line(f)
    self.hooke_angle_coeffs.append([
        int(values[0]),    # Z first element
        int(values[1]),    # Z vertex element
        int(values[2]),    # Z second element
        float(values[3]),  # k angle spring const
        float(values[4]),  # Rest angle (degrees)
        float(values[5]),  # Tolerance (degrees)
    ])
```

Inner is a 6-element 0-indexed list.

**Consumers**:

- `create_lammps_files` at `condense.py:1253–1271` and
  `1291–1301` — reads `hac[0..5]` throughout (element match,
  tolerance check, coefficient extraction).
- `normalize_types` at `condense.py:2016–2023` — same pattern.

**Severity**: DRIFT (inner axis only).

**Recommended fix**: `append([None, Z1, Zv, Z2, k, angle, tol])`
and change every `hac[n]` read to `hac[n + 1]`.

### C-C: `ordered_bonded_atoms` inner pair

**Perl** `src/scripts/condense:882–883`:

```perl
$orderedBondedAtoms[$bondCount][1] = $atom;
$orderedBondedAtoms[$bondCount][2] = $bonded_ref->[$atom][$bond];
```

**Python** `condense.py:1185`:

```python
ordered_bonded_atoms.append([atom, bonded_atom])
```

Inner is a 2-element 0-indexed pair.

**Consumer** `condense.py:1436–1441`:

```python
for bond in range(1, bond_count + 1):
    bt = ordered_bond_type[bond]
    ba = ordered_bonded_atoms[bond]
    lmp.write(
        f"{bond} {bt} "
        f"{ba[0]} {ba[1]} # "
        f"{unique_bond_tags_local[bt]}\n"
    )
```

**Severity**: DRIFT (inner axis only).  Local bookkeeping with
no external boundary; purely a convention issue.

**Recommended fix**: append `[None, atom, bonded_atom]` and
read `ba[1]`, `ba[2]`.

### C-D: `angle_bonded_atoms` inner triple

**Perl** `src/scripts/condense:999–1001`:

```perl
$angleBondedAtoms[$angleCount][1]=$angleBonded_ref->[$atom][$angle][1];
$angleBondedAtoms[$angleCount][2]=$atom;
$angleBondedAtoms[$angleCount][3]=$angleBonded_ref->[$atom][$angle][2];
```

Note the Perl convention: the vertex atom is stored at slot
*2* (middle), not slot 1 or 3 — the slot index literally gives
the role (slot 1 = first end, slot 2 = vertex, slot 3 = second
end).

**Python** `condense.py:1303`:

```python
angle_bonded_atoms.append([a1, atom, a2])
```

Inner is a 3-element 0-indexed triple with the vertex at slot
*1* (second).

**Consumer** `condense.py:1446–1452`:

```python
for angle in range(1, angle_count + 1):
    at = ordered_angle_type[angle]
    aa = angle_bonded_atoms[angle]
    lmp.write(
        f"{angle} {at} "
        f"{aa[0]} {aa[1]} {aa[2]} # "
        f"{unique_angle_tags_local[at]}\n"
    )
```

**Severity**: DRIFT (inner axis only).  Self-consistent (the
LAMMPS `Angles` line expects the two-end-atoms / vertex order
and the Python builder-reader pair preserves it), but the
vertex slot index differs from Perl.

**Recommended fix**: append `[None, a1, atom, a2]` and read
`aa[1..3]`.

### C-E: `ordered_species_pair_coeffs` inner pair (inherited)

**Python** `condense.py:1121–1124`:

```python
ordered_species_pair_coeffs[sp] = [
    ed.lj_pair_coeffs[z][0],
    ed.lj_pair_coeffs[z][1],
]
```

**Consumer** `condense.py:1341–1345`:

```python
pc = ordered_species_pair_coeffs[sp]
lmp.write(
    f"{sp} {pc[0]} {pc[1]} "
    f"# {ordered_species_tag[sp]}\n"
)
```

This finding is **inherited** from `element_data.py`, which
stores `lj_pair_coeffs[element]` as a 2-element tuple
`(epsilon, sigma)` — a 0-indexed inner pair.  The Perl
`ElementData.pm:173–174` source has
`$ljPairCoeffs[$element][1] = $values[0];` /
`$ljPairCoeffs[$element][2] = $values[1];` — 1-indexed inner.

**Severity**: DRIFT (inner axis, inherited).  The
`condense.py`-local copy `ordered_species_pair_coeffs` can be
fixed to 1-indexed independently of the upstream drift, but
the cleanest fix is to repair `element_data.py` first and
then change condense's builder and reader in lock-step.

**Scope note**: the root `lj_pair_coeffs` drift belongs in the
separate `element_data.py` audit (pending).  This finding
remains recorded here so the condense-side consumers are not
forgotten when the upstream fix is made.

### C-F: `unique_bond_coeffs` inner pair

**Perl** `src/scripts/condense:894–897` and `1047–1050`:

```perl
$uniqueBondCoeffs[$bond][1] = $hookeBondCoeffs_ref->[$hookeBond][3];
$uniqueBondCoeffs[$bond][2] = $hookeBondCoeffs_ref->[$hookeBond][4];
...
print LMPDAT "$bond $uniqueBondCoeffs[$bond][1] " .
   "$uniqueBondCoeffs[$bond][2] " .
   "# $uniqueBondTags[$bond]\n";
```

**Python** `condense.py:1198` (create_lammps_files) and
`1993` (normalize_types):

```python
unique_bond_coeffs[found] = [hbc[2], hbc[3]]
unique_bond_coeffs[ub] = [hbc[2], hbc[3]]
```

Readers at `condense.py:1351–1355` and `2092–2100`:

```python
bc = unique_bond_coeffs[b]
lmp.write(f"{b} {bc[0]} {bc[1]} ...")
```

**Severity**: DRIFT (inner axis).

**Recommended fix**: store `[None, k, r0]` at every producer
site; read `bc[1]`, `bc[2]` at every consumer.

### C-G: `unique_angle_coeffs` inner list

**Perl** `src/scripts/condense:992–994` and `1056–1059`:

```perl
$uniqueAngleCoeffs[$angle][1] =
   $hookeAngleCoeffs_ref->[$hookeAngle][4];
$uniqueAngleCoeffs[$angle][2] =
   $hookeAngleCoeffs_ref->[$hookeAngle][5];
...
print LMPDAT "$angle $uniqueAngleCoeffs[$angle][1] " .
   "$uniqueAngleCoeffs[$angle][2] ...";
```

**Python** `condense.py:1301` (create_lammps_files):

```python
unique_angle_coeffs[found] = [hac[3], hac[4]]
```

And `condense.py:2023` (normalize_types) — with a **slightly
different shape**:

```python
unique_angle_coeffs[ua] = [hac[3], hac[4], ha]
```

Note the 3-element form in `normalize_types` carries an extra
`ha` Hooke-angle ID appended to the tail.  The reader at
`condense.py:2155–2166` only uses `ac[0]` and `ac[1]`, so the
3rd element is written but never read — a minor latent mess
but not an audit finding by itself.

**Severity**: DRIFT (inner axis).  The 2-vs-3 element
inconsistency across the two methods is cosmetic only and
does not change behaviour.

**Recommended fix**: store `[None, k, angle_deg]` (2 elems in
create_lammps_files; optionally `[None, k, angle_deg, ha]` in
normalize_types if the tail is intended to be preserved).

### C-H: `unique_bond_types` inner pair

**Perl** `src/scripts/condense:1561–1562`:

```perl
$uniqueBondTypes[$numUniqueBonds][1] = $bondType1;
$uniqueBondTypes[$numUniqueBonds][2] = $bondType2;
```

**Python** `condense.py:1832` and `1939`:

```python
unique_bond_types.append([bt1, bt2])
```

Readers throughout `normalize_types`:

- `condense.py:1824–1827`
- `condense.py:1930–1934`
- `condense.py:2096–2098`
- `condense.py:2297–2307`

Every read is `unique_bond_types[ub][0]` / `[1]`.

**Severity**: DRIFT (inner axis).

**Recommended fix**: `append([None, bt1, bt2])` and read
`[1]`, `[2]`.

### C-I: `unique_angle_types` inner quadruple

**Perl** `src/scripts/condense:1592–1595`:

```perl
$uniqueAngleTypes[$numUniqueAngles][1] = $angleType1;
$uniqueAngleTypes[$numUniqueAngles][2] = $angleType2;
$uniqueAngleTypes[$numUniqueAngles][3] = $angleType3;
$uniqueAngleTypes[$numUniqueAngles][4] = $angleType4;
```

**Python** `condense.py:1860` and `1971`:

```python
unique_angle_types.append([at1, at2, at3, at4])
```

Readers:

- `condense.py:1848–1855` (lammps.dat Angle Coeffs parse)
- `condense.py:1956–1966` (mol-file Angles parse)
- `condense.py:2159–2165` (lammps.dat rewrite)
- `condense.py:2191–2198` (lammps.dat rewrite, angle reorder)
- `condense.py:2345–2354` (mol-file rewrite, angle reorder)

All reads are `unique_angle_types[ua][0..3]`.

**Severity**: DRIFT (inner axis).

**Recommended fix**: `append([None, at1, at2, at3, at4])` and
read `[1..4]`.

## Non-findings (verified consistent)

- **Outer stage / composition / reaction arrays.**
  `squish_factor`, `squish_step_size`, `ensemble_type`,
  `ensemble_temp_start`, `ensemble_temp_end`,
  `ensemble_temp_damp`, `run_steps`, `molecule_name`,
  `family_name`, `num_molecules`, `num_mol_atoms`,
  `atom_molecule_id`, `atom_molecule_name`,
  `molecule_type_of_mol`, `rxn_probability` — every one of
  these is initialised to `[None]` and appended to within a
  `range(1, N + 1)` loop.  Matches Perl 1-indexed conventions
  exactly.  OK.
- **Two-sided reaction data.** `rxn_mol_name` and
  `rxn_binding` are each built as `[None, [None], [None]]`
  giving a 1-indexed outer (side ∈ {1,2}) and a 1-indexed
  inner (reaction index).  Matches Perl `$rxnMolName[$side][$rxn]`.
  OK.
- **`trans_dist` passed to `sc.translate_atoms`.**  Built at
  `condense.py:1391`, `1393`, `1395` as
  `[None, dx, dy, dz]` — 1-indexed as required by the
  post-Cluster-A `StructureControl` interface.  OK.
- **`set_lattice_from_mag_angle` call.**  `condense.py:986–988`
  passes `[None, cell_size, cell_size, cell_size]` and
  `[None, 90.0, 90.0, 90.0]` — 1-indexed as required.  OK.
- **`direct_xyz[atom][1..3]` read in the Atoms section.**
  `condense.py:1409–1414` — 1-indexed axis access.  OK.
- **`ab = angle_bonded[atom]; a1 = ab[angle_idx][1]; a2 =
  ab[angle_idx][2]`.**  `condense.py:1218–1225` — reads
  `structure_control`'s 1-indexed inner pair from
  `angle_bonded[atom][angle_idx]` correctly.  OK.
- **`bond_tag_id_local[atom]` / `angle_tag_id_local[atom]`.**
  Each outer slot is a 1-indexed list with a `[None]` sentinel
  inner, filled via `.append(...)` per atom/bond.  OK.
- **File-token buffers.**  `values[0]`, `vals[3..16]`, etc.
  throughout `BondData`/`AngleData` parsers, `parse_input_file`,
  and `normalize_types` are legitimate 0-indexed reads from
  `.split()` output — the Perl analogue uses `@values` with
  0-indexed access too.  OK.

## Non-indexing observations

Flagged for a separate backlog, outside the audit scope:

- **`unique_angle_coeffs` 2-vs-3 element inconsistency.**
  `create_lammps_files` stores a 2-element entry `[k, angle]`
  at `condense.py:1301`; `normalize_types` stores a 3-element
  entry `[k, angle, hooke_id]` at `condense.py:2023`.  The
  extra `hooke_id` tail is never read.  Cosmetic; could be
  harmonised when Cluster C-G is fixed.
- **`bond_tag_id_local` outer storage is opaque.**  The
  builder at `condense.py:1171` appends both `None` entries
  (for bonds where `atom > bonded_atom`) and real tag IDs.
  The array is never subsequently *read* anywhere in the
  current code path — only the flat `unique_bond_tags_local`
  list is used downstream.  Not an audit finding but a dead
  data structure that could be removed during a follow-up
  pass.  Same applies to `angle_tag_id_local`.
- **`found_angle == 0` fallback prints a slightly wrong
  message.**  `condense.py:1264–1268` prints the
  `angle_idx` as "Angle number" but the user-facing label
  should be the running `angle_count`.  Cosmetic and
  unrelated to indexing.

## Running tally

- **BUG**: 0
- **DRIFT** (inner-axis flattening of nested-list inner
  dimensions; Perl 1-indexed → Python 0-indexed): 9
  - C-A through C-I; all self-consistent.
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK** (outer arrays, axis passes into `StructureControl`,
  file-parse token buffers, stage/composition/reaction 1-D
  arrays): all other layouts.
- **Inherited drift** (not in condense.py's scope, remediate
  upstream): 1 — `element_data.lj_pair_coeffs` (tuple form).

## Fix-order recommendation

All nine findings are inner-axis flattenings of nested lists.
The builder and every consumer of each nested structure are
local to `condense.py`, so the fix for each cluster is a
closed-file edit.  Suggested ordering once the user approves:

1. **Cluster CD-A** — `BondData.hooke_bond_coeffs` inner
   layout.  Touches the `BondData.init_bond_data` builder
   and two consumers in `create_lammps_files` and
   `normalize_types`.
2. **Cluster CD-B** — `AngleData.hooke_angle_coeffs` inner
   layout.  Same shape as CD-A but with 6-element entries
   and two consumer sites.
3. **Cluster CD-C** — `ordered_bonded_atoms` inner pair.
   Single builder + single consumer.
4. **Cluster CD-D** — `angle_bonded_atoms` inner triple.
   Single builder + single consumer (note the vertex-slot
   reordering).
5. **Cluster CD-E** — `ordered_species_pair_coeffs` inner
   pair.  *Blocked* on the upstream `element_data.py` fix
   for `lj_pair_coeffs`; both edits should land together.
6. **Cluster CD-F** — `unique_bond_coeffs` inner pair.
   Two builder sites (create_lammps_files + normalize_types)
   and two consumer sites.
7. **Cluster CD-G** — `unique_angle_coeffs` inner list.
   Two builder sites + two consumer sites; optionally
   harmonise the 2-vs-3 element inconsistency at the same
   time.
8. **Cluster CD-H** — `unique_bond_types` inner pair.
   Two builder sites + four consumer sites.
9. **Cluster CD-I** — `unique_angle_types` inner quadruple.
   Two builder sites + five consumer sites.

There is **no test coverage for `condense.py`** in
`tests/`.  After the fixes, regression verification will
have to be by inspection plus the full `structure_control`
test suite (which exercises the `direct_xyz` / `mag` /
`angle` / `translate_atoms` / `set_lattice_from_mag_angle`
call boundaries that condense depends on).  A smoke test
invoking `condense.py` on a tiny molecule set would be a
natural follow-up but is out of the indexing-audit scope.

## Fix status

Eight of the nine clusters have been applied to
`src/scripts/condense.py`.  Cluster CD-E is intentionally
deferred until the separate `element_data.py` audit lands.
Summary of the edits:

- **Cluster CD-A** — `BondData.init_bond_data` builds
  `[None, Z1, Z2, k, r0]`.  Consumers in
  `create_lammps_files` and `normalize_types` now read
  `hbc[1..4]`.
- **Cluster CD-B** — `AngleData.init_angle_data` builds
  `[None, Z1, Zv, Z2, k, angle_deg, tolerance]`.  Consumers
  in both methods now read `hac[1..6]`.
- **Cluster CD-C** — `ordered_bonded_atoms.append([None,
  atom, bonded_atom])`; reader prints `ba[1]`, `ba[2]`.
- **Cluster CD-D** — `angle_bonded_atoms.append([None, a1,
  atom, a2])` with the vertex at slot 2 matching Perl;
  reader prints `aa[1]`, `aa[2]`, `aa[3]`.
- **Cluster CD-E** — **deferred** (blocked on upstream
  `element_data.lj_pair_coeffs` tuple drift).
- **Cluster CD-F** — `unique_bond_coeffs[ub] = [None, k,
  r0]` at both producer sites; both readers print
  `bc[1]`, `bc[2]`.
- **Cluster CD-G** — `unique_angle_coeffs[ua] = [None, k,
  angle_deg]` at both producer sites; readers print
  `ac[1]`, `ac[2]`.  The historical 3-element form in
  `normalize_types` (which trailed an unused hooke-angle
  ID) has been **harmonised** to the same 2-element
  `[None, k, angle_deg]` shape as `create_lammps_files`,
  fixing the minor latent inconsistency noted in the
  non-findings section.
- **Cluster CD-H** — `unique_bond_types.append([None, bt1,
  bt2])` at both producer sites; every reader in
  `normalize_types` (six sites across lammps.dat parse,
  mol-file parse, coefficient computation, and both
  rewrite phases) now reads slots `[1]` / `[2]`.
- **Cluster CD-I** — `unique_angle_types.append([None, at1,
  at2, at3, at4])` at both producer sites; every reader
  (eight sites across lammps.dat parse, mol-file parse,
  coefficient computation, and both rewrite phases) now
  reads slots `[1]` / `[2]` / `[3]` / `[4]`.

Verification: `python -c "import ast; ast.parse(...)"`
passes on `condense.py`; a smoke import of `BondData`,
`AngleData`, and `Condense` classes succeeds.  Full
`tests/` suite — **295 passed, 23 skipped** — unchanged
from the pre-fix baseline (as expected; no
`condense.py`-specific tests exist, so the suite only
covers the shared `structure_control` call boundaries
that `condense.py` uses).
