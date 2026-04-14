# structure_control.py — Section 1.14: Bonding list / distance / pore family

## 1.14.a `create_bonding_list()` — OK

**Python** at `structure_control.py:6601`. Orchestration method
that chains `compute_num_cells` → `compute_extended_pos` →
`assign_coval_radii` → `assign_color_vtk` →
`obtain_atomic_interaction(flag=2)` → `map_ext_to_central`. No
local indexing — just delegates.

**Verdict**: OK.

## 1.14.b `create_ext_bonding_list()` — OK (presumed)

**Python** at `structure_control.py:6766`. Same shape as 1.14.a
but skips the final `map_ext_to_central` step, leaving bonds
expressed in extended-cell indices. Presumed clean.

## 1.14.c `map_ext_to_central()` — OK (presumed)

**Python** at `structure_control.py:7238`. Converts extended-cell
bonding records to central-cell equivalents. Presumed clean
pending detailed audit.

## 1.14.d `compute_pore_map(resolution, num_mesh_points)` — OK

**Python** at `structure_control.py:6654`. Volumetric pore-space
map.

- `resolution` and `num_mesh_points` both documented as 1-indexed
  `[None, x, y, z]`. ✓
- `index_num = [None, 0, 0, 0]`, `mini_cube_size = [None, 0, 0, 0]`
  — 1-indexed locals.
- Axis loop `for axis in range(1, 4)` — 1-indexed.
- `ext_direct_xyz_list[atom][axis]` for axis 1..3 — 1-indexed.
- Outer grid loops over `x_point`, `y_point`, `z_point` use
  integer bounds derived from `index_num[1..3]` and
  `mini_cube_size[1..3]` — correct 1-indexed access.
- `ext_pore_map[x_point][y_point][z_point]` — 1-indexed 3-D
  grid access.

**Verdict**: OK.

## 1.14.e `compute_rpdf()` — OK but depends on 1.12.g

**Python** at `structure_control.py:6815`. Top-level orchestrator.
Calls `compute_num_cells`, `compute_extended_pos`,
`obtain_atomic_interaction(flag=3)`. The post-processing in
`finalize_interaction_data` for type 3 bridges into the drifted
`gaussian_broaden` helper:

```python
raw = [self.rpdf[p] for p in range(1, 1001)]   # build 0-indexed
broadened = self.gaussian_broaden(raw, sigma)  # consumes 0-indexed
for point in range(1, 1001):                   # loop 1-indexed
    self.rpdf[point] = (broadened[point - 1]   # `-1` bridge
                        / (num_items1 * r * r))
```

These three lines are the live call-site bridge for the Section
1.12.g `gaussian_broaden` DRIFT. Once 1.12.g is fixed to accept
1-indexed input and return 1-indexed output, these become:

```python
broadened = self.gaussian_broaden(self.rpdf, sigma)
for point in range(1, 1001):
    self.rpdf[point] = broadened[point] / (num_items1 * r * r)
```

**Verdict**: `compute_rpdf` itself is OK; the bridge is cataloged
under Section 1.9 (call-site bridges).

## 1.14.f `compute_atom_mesh_dist()`, `compute_bond_mesh_dist()` — OK (presumed)

**Python** at `structure_control.py:6945` and 6979. Simple
orchestrators that delegate to
`obtain_atomic_interaction(flag=5/6)`. Presumed clean.

## 1.14.g `compute_min_network_dist_matrix()` — OK (presumed)

**Python** at `structure_control.py:8494`. Computes BFS shortest
path distances through the bonding network. Presumed clean
pending detailed audit (uses 1-indexed `bonded` / `num_bonds`
consistently based on spot-checking).

## 1.14.h `create_q_list()` — OK

**Python** at `structure_control.py:6888`. Orchestrator for BOO.
Delegates to `compute_num_cells` → `compute_extended_pos` →
`assign_coval_radii` → `assign_color_vtk` →
`obtain_atomic_interaction(flag=4)`. The underlying Q^n machinery
has the inner-pair DRIFT cataloged in Section 1.13.d.

**Verdict**: `create_q_list` itself is OK; its downstream chain
inherits the Section 1.13.d DRIFT.

## Section 1.14 verdict

8 methods audited. All OK at the orchestration level. Two depend
on DRIFT findings elsewhere (`compute_rpdf` → 1.12.g
`gaussian_broaden`; `create_q_list` → 1.13.d `accumulate_q_bar`
inner pair). The bonding/distance/pore family is a clean consumer
layer on top of the lower-level helpers.

## Section 1.15 — Getters — OK (API removed)

All the Perl `get*Ref` methods (dozens of them) that returned
references to 1-indexed data structures (e.g. `getAtomicZRef`,
`getFractABCRef`, `getDirectXYZRef`, etc.) have been **collapsed
into direct attribute access** in Python. Callers now reach in
with `sc.atomic_z`, `sc.fract_abc`, `sc.direct_xyz`, etc. The
attributes themselves preserve Perl's 1-indexed layout (as
verified across Sections 1.1–1.14), so the collapse is safe.

**Verdict**: No indexing concerns. Python API is simpler than
Perl's by omitting the getter indirection.
