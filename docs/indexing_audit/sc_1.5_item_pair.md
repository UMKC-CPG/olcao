# structure_control.py — Section 1.5: Item/pair bounds helpers

## 1.5.a `item_out_of_bounds(item, coords=None)` — OK (with expanded signature)

**Perl** `itemOutOfBounds` at `StructureControl.pm:8190`: takes
`$item, $itemCoords_ref`. `$borderZone` and `$borderLow`/`$borderHigh`
are package globals.

**Python** `item_out_of_bounds` at `structure_control.py:6472`:
takes `item` (1-indexed) and an optional `coords` array (1-indexed,
with the Python method auto-selecting `direct_xyz` / `direct_abc` /
`fract_abc` via `_border_coords()` if not provided). Axis loop
`range(1, 4)` matches Perl `(1..3)`. Accesses `coords[item][axis]`
just like the Perl version.

**Python addition** (not indexing): top-of-method guard
`if self.border_zone == 0: return False` — the Perl version relied
on callers to check `$borderZone` beforehand. Documented.

**Python addition** (not indexing): length-guard
`if item >= len(coords) or coords[item] is None: return False` —
again safer boundary handling absent in Perl.

Neither addition affects indexing.

## 1.5.b `item_not_requested(item, tag_list, request_list)` — OK

**Perl** `itemNotRequested` at `StructureControl.pm:8109`: iterates
`foreach $tag (1..$#$requestList)` and accesses `$tagList->[$item]`
vs `$requestList->[$tag]`. All 1-indexed.

**Python** `item_not_requested` at `structure_control.py:6525`:
iterates `for idx in range(1, len(request_list))` and accesses
`tag_list[item]` vs `request_list[idx]`. All 1-indexed. Empty-list
test uses `len(request_list) <= 1` to match Perl's
`$#$requestList < 0` (because the list has a `None` sentinel at
slot 0). Docstring explicitly calls out the correspondence.

## 1.5.c `pair_not_of_element(item1, item2)` — OK

**Perl** `pairNotOfElement` at `StructureControl.pm:8160`: reads
`$atomElementName[$item1]`, `$atomElementName[$item2]`,
`$rpdfSelectAtoms[1]`, `$rpdfSelectAtoms[2]`.

**Python** `pair_not_of_element` at `structure_control.py:6568`:
reads `self.atom_element_name[item1]`,
`self.atom_element_name[item2]`, `self.rpdf_select_atoms[1]`,
`self.rpdf_select_atoms[2]`. All 1-indexed. Full logic preserved.

## Section 1.5 verdict

3 OK findings. Pair/bounds helpers correctly preserved the
1-indexed convention.
