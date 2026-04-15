# mod_struct.py — Pass 1 audit

## Scope and method

`src/scripts/mod_struct.py` (~1061 lines) is a Python port of
the Perl `modStruct` script (~562 lines).  It applies a
user-specified sequence of structural operations (supercell,
translate, rotate, filter, add vacuum, cut sphere/block,
orthorhombic conversion, perturbation, surface preparation)
to an OLCAO structure file.  Each operation is parsed from
the command line into a dict, queued, and then dispatched in
order to a small `do_*` helper that calls the corresponding
`StructureControl` method.

The file has one class (`ScriptSettings`) plus ~12 top-level
`do_*` operation helpers plus a `main()` driver.

Pass 1 approach:

- `Grep` sweep for drift signals: `range(3)`, `range(0, …)`,
  `[atom - 1]`, `[axis - 1]`, `+1]`, `-1]`, nested-list
  append patterns.
- Read every operation parser branch plus every `do_*` helper
  and the `main` driver.  Cross-reference with Perl
  `src/scripts/modStruct` for every data structure that
  crosses a `StructureControl` call boundary.
- Run a direct smoke test against
  `sc.set_abc_xyz_assignment_order` with the defaults
  `rc["abc_order"]` to confirm at least one of the findings is
  a live crash.

## Summary table

| # | Topic | Perl layout | Python layout | Verdict |
|---|-------|-------------|---------------|---------|
| MS-A | `hkl` → `sc.prep_surface` | `$hkl[$op][1..3]` 1-indexed | `[h, k, l]` 0-indexed | **BUG** (IndexError on first read) |
| MS-B | `block_border` → `sc.cut_block` | `$blockBorder[1..3][1..2]` doubly 1-indexed | `[[f,t], [f,t], [f,t]]` doubly 0-indexed | **BUG** (IndexError on slot 3) |
| MS-C | `sphere_loc` → `sc.cut_sphere` | `$sphereLoc[1..3]` 1-indexed | `[x, y, z]` 0-indexed | **BUG** (IndexError on slot 3) |
| MS-D | `abc_order` / `xyz_order` → `sc.set_abc_xyz_assignment_order` | `$abc[1..3]` 1-indexed | `[1, 2, 3]` 0-indexed | **BUG** (IndexError on slot 3; fires on every invocation of `main()` before any operation runs) |
| MS-E | `op["sc"]` / `op["mirror"]` inner layout | `$sc[$op][1..3]` / `$mirror[$op][1..3]` 1-indexed | `[s1, s2, s3]` / `[m1, m2, m3]` 0-indexed | **DRIFT** (bridge via `i + 1` in `do_supercell` works; layout diverges from Perl) |
| — | `op["trans"]` | `$trans[$op][1..3]` 1-indexed | `[None, tx, ty, tz]` 1-indexed | OK |
| — | `op["rot"]`, `op["origin"]` | `$rot[$op][1..3]`, `$rotOrig[$op][1..3]` 1-indexed | `[None, rx, ry, rz]`, `[None, ox, oy, oz]` 1-indexed | OK |
| — | `p1`, `p2`, `p3` in `-rotP` | `@p1[1..3]` 1-indexed in Perl, but passed to `getPlaneNormal` via reference and then consumed element-wise | `[x, y, z]` 0-indexed local scratch for `_get_plane_normal` | OK — local helper takes 0-indexed per its docstring contract |

## Why `set_abc_xyz_assignment_order` crashes every invocation

The audit turned up a **critical live bug**: every invocation
of `mod_struct.py main()` crashes in `sc.set_abc_xyz_assignment_order`
*before any operation runs*.

**Reproduction** (self-contained smoke test):

```python
from structure_control import StructureControl
sc = StructureControl()
abc_order = [1, 2, 3]   # default from mod_structrc.py:26
xyz_order = [1, 2, 3]   # default from mod_structrc.py:27
sc.set_abc_xyz_assignment_order(abc_order, xyz_order)
# -> IndexError: list index out of range
```

Confirmed with the actual `StructureControl` class from
`src/scripts/structure_control.py`.

`set_abc_xyz_assignment_order` (`structure_control.py:775`)
reads `abc_order[axis]` for `axis in range(1, 4)`:

```python
for axis in range(1, 4):
    self.abc_order[axis] = abc_order[axis]
    self.xyz_order[axis] = xyz_order[axis]
```

With a 3-element `[1, 2, 3]`, slot 3 does not exist → the
third iteration raises `IndexError`.

**Impact**: any tool that imports `mod_struct` and calls its
`main()` (which is what the `main.reactions.py` condense
pipeline does via `subprocess.run(["mod_struct.py", ...])`)
dies on line 1019 before reading the input file.  This bug
has been latent since the initial port; no test exercises it
and the failing call sits inside an external-subprocess
boundary that the existing test suite never crosses.

## Finding details

### MS-A — `hkl` → `sc.prep_surface` (BUG)

**Perl** `modStruct:408–414`:

```perl
elsif ($ARGV[$number] eq "-rothkl")
{
   $ops[++$numOps] = 8;
   $hkl[$numOps][1] = $ARGV[++$number];
   $hkl[$numOps][2] = $ARGV[++$number];
   $hkl[$numOps][3] = $ARGV[++$number];
}
```

Perl stores `$hkl[$op][1..3]`.  Consumer at `modStruct:228`:

```perl
elsif ($ops[$op] == 8)
   {StructureControl::prepSurface(\@{$hkl[$op]});}
```

Passes the 1-indexed row by reference.

**Python** `scan_operations` at `mod_struct.py:506–516`:

```python
elif flag == "-rothkl":
    # -rothkl H K L
    hkl = [
        float(argv[n + 1]),
        float(argv[n + 2]),
        float(argv[n + 3]),
    ]
    n += 3
    self.ops.append({
        "op": OP_ROTHKL,
        "hkl": hkl,
    })
```

Builds a 3-element **0-indexed** `[h, k, l]`.  Consumer
`do_rothkl` at `mod_struct.py:955`:

```python
sc.prep_surface(op["hkl"])
```

`sc.prep_surface` (`structure_control.py:4519`) reads
`hkl[hkl_idx]` for `hkl_idx in range(1, 4)`:

```python
for hkl_idx in range(1, 4):
    if hkl[hkl_idx] != 0:
        product *= hkl[hkl_idx]
```

With a 3-element input, slot 3 does not exist → IndexError on
the first iteration that actually reads it (after loading the
structure and getting past the first two reads, which would
silently return k and l where h and k were expected).

**Severity**: **BUG**.  Any invocation with `-rothkl` crashes.

**Recommended fix**: store
`hkl = [None, float(argv[n+1]), float(argv[n+2]),
float(argv[n+3])]`.

### MS-B — `block_border` → `sc.cut_block` (BUG)

**Perl** `modStruct:457–462, 470–475`:

```perl
$blockBorder[1][1] = $ARGV[$number+2];  # from-a
$blockBorder[1][2] = $ARGV[$number+3];  # to-a
$blockBorder[2][1] = $ARGV[$number+4];
$blockBorder[2][2] = $ARGV[$number+5];
$blockBorder[3][1] = $ARGV[$number+6];
$blockBorder[3][2] = $ARGV[$number+7];
```

Doubly 1-indexed: outer axis 1..3, inner from/to 1..2.
Consumer at `modStruct:230` passes it by reference.

**Python** `scan_operations` at `mod_struct.py:557–561,
578–585, 591–598`:

```python
block_border = [
    [0.0, 0.0],
    [0.0, 0.0],
    [0.0, 0.0],
]
...
block_border = [
    [float(argv[n + 2]), float(argv[n + 3])],
    [float(argv[n + 4]), float(argv[n + 5])],
    [float(argv[n + 6]), float(argv[n + 7])],
]
```

Doubly 0-indexed.  Consumer `do_cutblock` at
`mod_struct.py:971`:

```python
sc.cut_block(
    op["zone"],
    op["abcxyz_flag"],
    op["block_border"],
)
```

`sc.cut_block` (`structure_control.py:5321`) reads
`block_border[axis][1]` and `block_border[axis][2]` for
`axis in range(1, 4)`:

```python
for axis in range(1, 4):
    ...
    if block_border[axis][2] == sentinels[axis]:
        block_border[axis][2] = self.mag[axis]
```

With doubly-0-indexed input, `block_border[1]` returns the
second row (b-axis), `block_border[1][1]` returns the to-value
rather than from, and `block_border[3]` raises IndexError.

**Severity**: **BUG**.  Any invocation with `-cutblock`
crashes (after silently mis-shifting the first two axes).

**Recommended fix**: store
`[None, [None, f1, t1], [None, f2, t2], [None, f3, t3]]`.

### MS-C — `sphere_loc` → `sc.cut_sphere` (BUG)

**Perl** `modStruct:519–521, 529–531`:

```perl
$sphereLoc[1] = $ARGV[$number+2];
$sphereLoc[2] = $ARGV[$number+3];
$sphereLoc[3] = $ARGV[$number+4];
```

1-indexed.  Passed to `cutSphere` via `\@sphereLoc` reference.

**Python** `scan_operations` at `mod_struct.py:618, 645–649,
655–659`:

```python
sphere_loc = [0.0, 0.0, 0.0]
...
sphere_loc = [
    float(argv[n + 2]),
    float(argv[n + 3]),
    float(argv[n + 4]),
]
```

3-element 0-indexed.  Consumer `do_cutsphere` at
`mod_struct.py:992`:

```python
sc.cut_sphere(
    op["zone"],
    op["abcxyz_flag"],
    op["sphere_rad"],
    op["target_atom"],
    op["sphere_loc"],
)
```

`sc.cut_sphere` (`structure_control.py:5404`) reads
`sphere_loc[1..3]`:

```python
distance = ((coords[1] - sphere_loc[1])**2 +
            (coords[2] - sphere_loc[2])**2 +
            (coords[3] - sphere_loc[3])**2) ** 0.5
```

Same IndexError / silent slot-shift pattern as MS-B.

**Note**: `cut_sphere` internally rebuilds `sphere_loc` as
1-indexed when `target_atom != 0`, but the `target_atom == 0`
path (which is taken for `-atxyz` and `-atabc`) passes the
caller's `sphere_loc` through unchanged.

**Severity**: **BUG**.  Any invocation with `-cutsphere
-atxyz` or `-cutsphere -atabc` crashes.  `-cutsphere -atom`
accidentally works because `cut_sphere` short-circuits the
bad `sphere_loc`.

**Recommended fix**: store
`[None, float(argv[n+2]), float(argv[n+3]), float(argv[n+4])]`.

### MS-D — `abc_order` / `xyz_order` → `sc.set_abc_xyz_assignment_order` (BUG)

**Perl** `modStruct:260–261, 289–294`:

```perl
$abc[1] = 1; $abc[2] = 2; $abc[3] = 3;
$xyz[1] = 1; $xyz[2] = 2; $xyz[3] = 3;
...
$abc[1] = $ARGV[++$number];
$xyz[1] = $ARGV[++$number];
$abc[2] = $ARGV[++$number];
$xyz[2] = $ARGV[++$number];
$abc[3] = abs($abc[1]+$abc[2]-6);
$xyz[3] = abs($xyz[1]+$xyz[2]-6);
```

1-indexed.  Consumer at `modStruct:200`:

```perl
StructureControl::setABCXYZAssignmentOrder(\@abc,\@xyz);
```

**Python** `mod_structrc.py:26–27`:

```python
"abc_order": [1, 2, 3],  # a=1, b=2, c=3
"xyz_order": [1, 2, 3],  # x=1, y=2, z=3
```

3-element 0-indexed.  `ScriptSettings.assign_rc_defaults` at
`mod_struct.py:227–228` copies them verbatim.  `reconcile`
at `mod_struct.py:718–729` writes slots `[0], [1], [2]` (also
0-indexed).  Consumer `main` at `mod_struct.py:1019`:

```python
sc.set_abc_xyz_assignment_order(
    settings.abc_order, settings.xyz_order
)
```

`sc.set_abc_xyz_assignment_order`
(`structure_control.py:775`) reads `abc_order[axis]` for
`axis in range(1, 4)`.  With a 3-element input, slot 3 does
not exist.

**Severity**: **BUG** — critical, fires on **every**
invocation of `mod_struct.py main()` before any operation
parser is consulted.  **Confirmed with a live smoke test**
(see the "Why … crashes every invocation" section above).

**Recommended fix**: update `mod_structrc.py` defaults to
`[None, 1, 2, 3]`, update `reconcile` to write slots
`[1], [2], [3]`.

### MS-E — `op["sc"]` / `op["mirror"]` inner layout (DRIFT)

**Perl** `modStruct:299–304`:

```perl
$sc[$numOps][1] = $ARGV[++$number];
$sc[$numOps][2] = $ARGV[++$number];
$sc[$numOps][3] = $ARGV[++$number];
$mirror[$numOps][1] = $ARGV[++$number];
$mirror[$numOps][2] = $ARGV[++$number];
$mirror[$numOps][3] = $ARGV[++$number];
```

1-indexed inner.  Consumer at `modStruct:212`:

```perl
StructureControl::applySupercell(
    @{$sc[$op]}[1..3],
    @{$mirror[$op]}[1..3]
);
```

Passes the slice `[1..3]` of each row.

**Python** `scan_operations` at `mod_struct.py:344–353`:

```python
sc_vals = [
    int(argv[n + 1]),
    int(argv[n + 2]),
    int(argv[n + 3]),
]
mirror_vals = [
    int(argv[n + 4]),
    int(argv[n + 5]),
    int(argv[n + 6]),
]
```

3-element 0-indexed.  Consumer `do_supercell` at
`mod_struct.py:804–806` bridges via `i + 1`:

```python
for i in range(3):
    sc.supercell[i + 1] = op["sc"][i]
    sc.sc_mirror[i + 1] = op["mirror"][i]
```

**Severity**: **DRIFT** (not a bug — the bridge works).
Inner layout diverges from Perl's 1-indexed convention.

**Recommended fix**: store
`[None, int(...), int(...), int(...)]` at the producer and
read `op["sc"][axis]` / `op["mirror"][axis]` for `axis in
range(1, 4)` in `do_supercell`.

## Non-findings (verified consistent)

- **`op["trans"]`** — built 1-indexed
  `[None, tx, ty, tz]` at `mod_struct.py:364–368`.  Consumer
  `do_translate` passes it straight to
  `sc.translate_atoms(1, sc.num_atoms, op["trans"])`, which
  per the `structure_control` Section 1.2.f audit expects a
  1-indexed displacement.  OK.
- **`op["rot"]`, `op["origin"]`** (both `-rotA` and `-rotP`
  parsers) — built 1-indexed `[None, …]` at
  `mod_struct.py:381–386, 400–413, 454, 474, 487–492`.
  Consumer `do_rotate` reads `origin[1..3]` and passes
  `axis = op["rot"]` straight to `sc.define_rot_matrix`
  which was fixed to 1-indexed during the
  `structure_control` Cluster B fix sweep.  OK.
- **`p1`, `p2`, `p3` in `-rotP`** — built 3-element 0-indexed
  locals (`mod_struct.py:426–440`), consumed only by
  `_get_plane_normal(p1, p2, p3)` at line 448.  The local
  helper's docstring explicitly documents the 0-indexed
  input contract and returns a 1-indexed
  `[None, nx, ny, nz]` vector.  The returned vector is then
  normalised into a 1-indexed rot axis and stored in
  `op["rot"]`.  OK.
- **`settings.ops` queue** — plain Python list built via
  `.append({...})`; each dict carries its per-operation
  parameters.  The outer queue has no index convention to
  choose.  OK.
- **`_get_plane_normal` 0→1 hybrid** — input 0-indexed,
  output 1-indexed.  This is the intentional interface that
  was established during the `structure_control` Cluster A
  fix sweep for `mod_struct.py`; the Perl
  `getPlaneNormal` in `StructureControl.pm` is fully
  1-indexed on both ends, but the Python local helper in
  `mod_struct.py` retains the 0-indexed input because the
  `-rotP` parser receives its points as flat 0-indexed
  tuples from argparse and the helper is purely local to
  `mod_struct.py`.  Documented as a deliberate convention.
  OK.

## Running tally

- **BUG**: 4
  - MS-A: `hkl` via `do_rothkl` → `sc.prep_surface`
  - MS-B: `block_border` via `do_cutblock` → `sc.cut_block`
  - MS-C: `sphere_loc` via `do_cutsphere` → `sc.cut_sphere`
  - MS-D: `abc_order` / `xyz_order` via `main` →
    `sc.set_abc_xyz_assignment_order`  (critical; fires on
    every invocation of `main()`)
- **DRIFT**: 1
  - MS-E: `op["sc"]` / `op["mirror"]` inner layout
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: everything else

## Fix-order recommendation

Four BUGs plus one DRIFT, all in `mod_struct.py`.  Suggested
order:

1. **Cluster MS-D** first — the `abc_order` / `xyz_order`
   crash is the most critical because it fires on every
   invocation, blocking all other operations.  Fix
   `mod_structrc.py` defaults to `[None, 1, 2, 3]`, update
   `reconcile` slot writes, update the internal recovery
   line `self.abc_order[2] = abs(self.abc_order[0] + self.abc_order[1] - 6)`
   to use slots `[1], [2], [3]`.
2. **Clusters MS-A / MS-B / MS-C** — fix the producer side
   of each to emit `[None, …]` / doubly `[None, …]` lists.
   No consumer-side changes needed because the SC methods
   already expect 1-indexed inputs.
3. **Cluster MS-E** — rewrite the `-sc` parser to emit
   `[None, sc1, sc2, sc3]` and `[None, m1, m2, m3]`, then
   update `do_supercell` to index via `op["sc"][axis]` for
   `axis in range(1, 4)`.

Regression testing: no `mod_struct.py`-specific tests exist
in `tests/`, so post-fix verification is by inspection plus:

- A standalone smoke test that imports mod_struct, builds a
  trivial `ScriptSettings`, and verifies
  `sc.set_abc_xyz_assignment_order` no longer raises.
- The full `structure_control` suite (unchanged expectations;
  the fixes do not touch `structure_control.py`).

A proper end-to-end smoke test that runs `mod_struct.py` as a
subprocess on a tiny .skl fixture with each operation in turn
would be a valuable follow-up but is out of the indexing-audit
scope.  Given the severity of MS-D and that make_reactions.py
invokes `mod_struct.py` as a subprocess in its condense
pipeline, this follow-up should be prioritized.

## Fix status

**All five clusters applied.**  Summary of the edits:

- **Cluster MS-D** (`mod_structrc.py` + `mod_struct.py
  reconcile`).  Defaults updated to
  `abc_order = [None, 1, 2, 3]` and
  `xyz_order = [None, 1, 2, 3]`.  `reconcile` now writes
  slots `[1]`, `[2]`, `[3]` (previously `[0]`, `[1]`, `[2]`),
  and the `|a+b−6|` recovery formula reads/writes the same
  slots.  Confirmed with a live smoke test that
  `sc.set_abc_xyz_assignment_order(rc["abc_order"],
  rc["xyz_order"])` no longer raises.
- **Cluster MS-A** (`scan_operations`, `-rothkl` branch).
  `hkl` is now built as `[None, float(...), float(...),
  float(...)]`.  Consumer `do_rothkl` is unchanged because
  `sc.prep_surface` already expects the 1-indexed layout.
- **Cluster MS-B** (`scan_operations`, `-cutblock` branch).
  `block_border` is now built as
  `[None, [None, f1, t1], [None, f2, t2], [None, f3, t3]]`
  at both the default-initialisation site and the two
  per-axis populators (`-xyz` and `-abc`).  Consumer
  `do_cutblock` is unchanged because `sc.cut_block` already
  expects doubly 1-indexed.
- **Cluster MS-C** (`scan_operations`, `-cutsphere` branch).
  `sphere_loc` is now built as `[None, x, y, z]` at the
  default and the two populators (`-atxyz` and `-atabc`).
  Consumer `do_cutsphere` is unchanged because
  `sc.cut_sphere` already expects the 1-indexed layout
  (except on the `target_atom != 0` short-circuit path,
  which rebuilds internally).
- **Cluster MS-E** (`scan_operations`, `-sc` branch, and
  `do_supercell`).  `sc_vals` and `mirror_vals` are now
  built as `[None, v1, v2, v3]`.  `do_supercell`'s axis loop
  is rewritten to `for axis in range(1, 4)` and reads
  `op["sc"][axis]` / `op["mirror"][axis]` directly, removing
  the former `i + 1` bridge.

Verification steps:

- `python -c "import ast; ast.parse(open('mod_struct.py').read())"`
  passes for both `mod_struct.py` and `mod_structrc.py`.
- The live smoke test for MS-D now succeeds (no IndexError);
  the full smoke script also verifies MS-A/B/C/E slot reads
  by constructing the new layouts and reading slots 1..3 on
  each.
- Full `tests/` suite — **295 passed, 23 skipped** —
  unchanged from the pre-fix baseline (as expected: no
  `mod_struct.py`-specific tests exist, and the
  `structure_control` methods that the fixed `mod_struct.py`
  calls into are unaffected by the caller-side convention
  flip).

A follow-up to add an actual end-to-end smoke test for
`mod_struct.py` (invoking it as a subprocess on a tiny
fixture with each operation in turn) remains recommended
for a later quality pass but is out of the indexing-audit
scope.
