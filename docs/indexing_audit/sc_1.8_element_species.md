# structure_control.py — Section 1.8: Element/species-index methods

## 1.8.a `create_element_list()` — OK

**Perl** `createElementList` at `StructureControl.pm:4016`:
`$elementList[0] = ""`, `foreach $atom (1..numAtoms)`,
`foreach $element (1..numElements)`, assigns
`$atomElementName[$atom] = $values[0]` (where `$values` is the
0-indexed return from `prepLine` — which is a separate helper,
audited in Pass 4 below).

**Python** at `structure_control.py:6001`: `element_list = [None]`,
`atom_element_name = [None] * (num_atoms + 1)`, all loops
1-indexed, uses `re.split(r'[0-9]+', ...)[0]` to match Perl's
0-indexed-first behavior. Preserved correctly.

## 1.8.b `create_species_data(use_file_species=True)` — OK

**Perl** `createSpeciesData` at `StructureControl.pm:3949`:
`$speciesList[$elementID][$species]` — doubly 1-indexed (both outer
element and inner species). `$numSpecies[1..numElements]`,
`$atomSpeciesID[1..numAtoms]`.

**Python** at `structure_control.py:6050`:
`species_list = [None] + [[None] for _ in range(num_elements)]` —
doubly 1-indexed, matching Perl. Inner lists are extended via
`.append()` on a pre-`[None]` list so new species land at index 1,
2, ... as expected.

## 1.8.c `count_element_atoms()` — OK

`element_count = [None] + [0] * num_elements` — 1-indexed. Direct
increment via `element_count[atom_element_id[atom]] += 1`. Matches
Perl.

## 1.8.d `map_element_number()` — OK

Writes `atomic_z[1..num_atoms]` from `get_element_number`, which
delegates to `element_data.get_element_z()`. 1-indexed throughout.

## 1.8.e `get_element_number(element_name)` — OK (with mechanism change)

**Perl** `getElementNumber` at `StructureControl.pm:6324`: linear
scan through `$elementNames_ref->[1..N]` (1-indexed). Return value
is the 1-based loop index, which equals the atomic number Z.

**Python** at `structure_control.py:6295`: delegates to
`element_data.get_element_z()`, which performs the same scan. The
returned Z value is the 1-based index, matching Perl. Documented.

## 1.8.f `assign_coval_radii()` — OK

Both Perl and Python: `coval_radii[1..num_atoms]` 1-indexed,
populated from `element_data.coval_radii[z]` where `z` is also
1-indexed.

## 1.8.g `assign_color_vtk()` — OK (legitimate RGBA exception)

The `color_vtk` outer list is 1-indexed on atoms with `[None]`
sentinel, but each entry is an RGBA quadruple that is **genuinely
0-indexed** in both Perl and Python (`[0]=r, [1]=g, [2]=b, [3]=a`).
The Python docstring explicitly notes: *"The Perl implementation
copies the four components explicitly as indices [0]..[3]; Python
uses list() to copy all components at once."*

This is the one legitimate case where an inner list is 0-indexed
in both languages because the Perl original is 0-indexed.
Correctly preserved.

## Section 1.8 verdict

7 methods audited, all OK. The element/species infrastructure
matches Perl exactly, including the RGBA-specific 0-indexed
exception.
