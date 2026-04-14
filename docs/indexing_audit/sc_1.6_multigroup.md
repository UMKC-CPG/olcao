# structure_control.py — Section 1.6: Multi-group distance/interaction builders

## 1.6.a `create_min_dist_matrix(...)` — OK (with parameter reorder)

**Perl** `createMinDistMatrix` at `StructureControl.pm:6364`:
signature is `(numItems1, numItems2, items1ABC, items2ABC,
items1XYZ, items2XYZ)` — counts first, then all ABCs, then all XYZs.
All list refs 1-indexed; internal item loops `(1..numItems1)` and
`(1..numItems2)`.

**Python** `create_min_dist_matrix` at `structure_control.py:6384`:
signature is `(num_items1=None, items1_abc=None, items1_xyz=None,
num_items2=0, items2_abc=None, items2_xyz=None)` — per-group
interleaving, all optional. Docstring explicitly documents the
parameter reorder. All indexing remains 1-indexed.

**Severity**: Signature change, not indexing. The reorder is a
Python idiom improvement and is documented.

## 1.6.b `compute_num_cells(num_items=None, fract_abc=None)` — OK (with output-arg change)

**Perl** `computeNumCells` at `StructureControl.pm:6658`: takes
`($numItems, $fractABC_ref, $negBit_ref, $posBit_ref)` — note the
last two are output ref args. All 1-indexed internally.

**Python** `compute_num_cells` at `structure_control.py:7017`:
takes only `(num_items, fract_abc)`; the output arrays `neg_bit`,
`pos_bit` are stored as `self.neg_bit`, `self.pos_bit` instead of
being passed in. Documented.

All loops are `range(1, 4)` and `range(1, num_items + 1)`; all
accesses use 1-indexed subscripts; sentinel `[None]` is used when
initializing per-item rows. No indexing drift.

## 1.6.c `compute_extended_pos(num_items=None, fract_abc=None)` — OK (with output-arg change)

Same pattern as 1.6.b. The Perl version had seven output ref args
(`\@negBit, \@posBit, ..., \@itemsTotalExtXYZList,
\@itemsTotalExtABCList`); the Python version stores results as
`self.*` attributes or constructs local 1-indexed lists as needed.
Indexing is preserved; only the mechanism for returning results
differs.

## Section 1.6 verdict

3 OK findings. Signature/mechanism changes without indexing drift.
