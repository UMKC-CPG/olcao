# element_data.py — Pass 1 Audit

**Status: 1 DRIFT (ED-A).** **Confirmed: ED-A is the root of the
deferred CD-E cluster in `condense_pass1.md`.**

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 1 (ED-A) |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | all other tables (ED-B through ED-N) |

## Scope

- Perl original: `src/scripts/ElementData.pm`
- Python port:   `src/scripts/element_data.py`

This file is the per-element data provider for every converted
script. Any drift here propagates to every consumer — audit is
high priority.

## Finding ED-A — `lj_pair_coeffs` inner pair — DRIFT (root of CD-E)

**Perl** `ElementData.pm:173-174`:
```perl
$ljPairCoeffs[$element][1] = $epsilon;
$ljPairCoeffs[$element][2] = $sigma;
```
Inner axis 1-indexed with slot 0 unused.

**Python** `element_data.py:322-326`: builds a 0-indexed 2-tuple
`(epsilon, sigma)`. Outer element axis is correctly `[None]`-
padded 1-indexed, but inner pair is flattened to 0-indexed.

Docstring drift corroborates: `element_data.py:38, 81-83`
explicitly documents the tuple shape.

**Consumers**:
- `condense.py:1153-1156` (`ordered_species_pair_coeffs`
  producer, reads `[0]` / `[1]`).
- `condense.py:2161-2162` (`unique_lj_pair_coeffs` formatter,
  reads `[0]` / `[1]`).

No other scripts touch this table.

**Severity**: DRIFT only; self-consistent, no runtime hazard.

## CD-E resolution target — CONFIRMED

Cluster ED-A **is** the concrete upstream fix target that
`condense_pass1.md` Cluster CD-E was deferred for. The coordinated
fix is a five-touchpoint sweep:

1. `element_data.py:322-326` — builder →
   `[None, float(vals[0]), float(vals[1])]` (list, not tuple).
2. `element_data.py:38, 81-83, 147` — docstring and init comment
   updated to `[None, epsilon, sigma]`.
3. `condense.py:1153-1156` — producer reads
   `ed.lj_pair_coeffs[z][1]` / `[2]`.
4. `condense.py:2161-2162` — formatter reads `lj[1]` / `lj[2]`.
5. Simultaneously apply CD-E: store
   `ordered_species_pair_coeffs[sp] = [None, eps, sig]` and read
   `pc[1]` / `pc[2]` at `condense.py:1341-1345`.

**All five touchpoints must land in one commit** — no intermediate
state works, because the tuple → list change breaks both existing
`condense.py` reads.

## Non-findings (all verified OK)

- **spdf inner axes** (`core_orbitals`, `core_charge`,
  `vale_charge`, `num_terms_wf`, `vale_orbitals` inner) — all
  legitimately 0-indexed in both Perl (`0..$maxQN_l`) and Python.
  Match, not drift.
- **`color_vtk`** RGBA inner is legitimately 0-indexed per
  README.md's stated exception and SC §1.8.g.
- **`orbital_terms[qn_l][element]`** axis ordering matches Perl
  exactly.
- **`vale_orbitals`** outer basis axis is 1-indexed with
  `[None, None, None, None]` pre-pad.
- **All 1-D element-keyed tables** (`masses`, `radii`, `names`,
  `neut_scatt`, `num_uj_electrons`, `min/max_term_wf`,
  `num_terms_pot`, `min/max_term_pot`, `color_dx`, `grey_dx`,
  `element_names`, `element_full_names`) — 1-indexed with `[None]`
  sentinel.
- **`apply_bond_factor`** and **`get_element_z`** loops walk
  `1..numElements` correctly.

## Running tally

- **BUG**: 0
- **DRIFT**: 1 (ED-A)
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every other table

## Fix-order recommendation

**Highest priority**: ED-A + CD-E must be fixed together as a
single coordinated commit across `element_data.py` and
`condense.py`. This unblocks the deferred cluster that was noted
as blocking at the start of the audit sweep.

No tests for `element_data.py`; regression verification is
`ast.parse` + a smoke-load of `lj_pair_coeffs[1]` / `[14]` + the
full `structure_control` suite baseline + an optional `condense.py`
smoke run.
