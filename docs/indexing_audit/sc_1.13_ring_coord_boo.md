# structure_control.py — Section 1.13: Ring / coordination / BOO family

## 1.13.a `compute_ring_distribution(min_ring_len, max_ring_len)` — OK (presumed)

**Python** at `structure_control.py:8022`. Top-level driver for
ring identification. Not deeply audited yet but the supporting
methods (`compute_ring_level_delta`, `save_ring`) were audited
clean in Section 1.7.d/e and the method consumes the 1-indexed
`num_bonds` / `bonded` arrays.

**Verdict**: Presumed clean pending detailed audit.

## 1.13.b `compute_bond_angles_ext()` — OK

**Python** at `structure_control.py:8405`. Computes bond angles
via law of cosines.

- `self.num_bond_angles = [None] + [0 for _ in range(num_atoms)]`
  — 1-indexed atom-count list.
- `self.bond_angles_ext = [None] + [[None] for _ in range(num_atoms)]`
  — 1-indexed outer on atom, 1-indexed inner with `[None]`
  sentinel for the per-atom angle list (appended at index 1..).
- Atom loop `range(1, num_atoms + 1)` — 1-indexed.
- Bond loops `range(1, num_bonds[atom])` and `range(bond1+1,
  num_bonds[atom]+1)` — 1-indexed.
- `bond_length_ext[atom][bond1]`, `bonded_ext[atom][bond1]` —
  1-indexed.
- `ext_direct_xyz_list[ext1][ax]` for `ax in range(1, 4)` —
  1-indexed.

**Verdict**: OK.

## 1.13.c `init_ylm_coeff()` — OK

**Python** at `structure_control.py:8817`. Initialises the
spherical harmonic normalisation coefficients.

- `self.ylm_coeff = [None, ...]` — 1-indexed with sentinel.
  For `l=0`: 1 coefficient, list length 2 (slots 0 and 1).
  For `l=1`: 3 coefficients, list length 4 (slots 0 through 3).
  For `l=6`: 13 coefficients, list length 14 (slots 0 through 13).
- The coefficient at m-index `k` is accessed as `ylm_coeff[k]` for
  `k in 1..2l+1` — matching Perl `$Ylm_Coeff[1..2l+1]`.

**Verdict**: OK.

## 1.13.d `accumulate_q_bar(atom, theta, phi)` — DRIFT (inner [real, imag] pair)

**Perl** `accumulateQBar` at `StructureControl.pm:8526`: accesses
`$qBar[$atom][$k][1]` for real part and `$qBar[$atom][$k][2]` for
imaginary part. The innermost `[real, imag]` pair is **1-indexed**
with slots `[1]` and `[2]` (slot `[0]` is unused).

**Python** `accumulate_q_bar` at `structure_control.py:8874`: uses
`qb[k][0]` for real part and `qb[k][1]` for imaginary part. The
innermost pair is **0-indexed** `[0]=real, [1]=imag`.

**The docstring explicitly acknowledges the change**: *"Each entry
q_bar[atom][k] is a two-element list [real, imag] (Python indices
0/1; Perl used 1/2)."*

**Severity**: DRIFT. The `[real, imag]` inner pair was 1-indexed
in Perl (with slot 0 unused) and was flattened to 0-indexed
`[real, imag]` in Python. Every read and write of `qb[k][0]`,
`qb[k][1]` is affected.

**Affected sites** (all in `structure_control.py`):

- `initialize_interaction_data` at line 7544:
  `[[0.0, 0.0] for _ in range(m)]` — 0-indexed pair initialization.
- `accumulate_q_bar` at lines 8924–8925 (l=0 branch),
  8929–8935 (l=1 branch), 8962–8969 (l=6 branch) — all reads/writes
  of `qb[k][0]` and `qb[k][1]`.
- `q_bar_per_bond` at line 8996 onwards — divides `q_bar[atom][m][0]`
  and `q_bar[atom][m][1]`.
- `q_bar_normalize` — normalises real/imag at `[0]`/`[1]`.
- `q_bar_correlation` — accumulates dot product over `[0]` and `[1]`.

**Recommended fix**: Restore the 1-indexed inner layout. Initialise
as `[None, 0.0, 0.0]` in `initialize_interaction_data`. Change
every `qb[k][0]` → `qb[k][1]` (real) and every `qb[k][1]` →
`qb[k][2]` (imag). Update `q_bar_per_bond`, `q_bar_normalize`,
`q_bar_correlation` loops that iterate the pair to use
`for part in range(1, 3)`.

Consequence: the 5-method chain
(`accumulate_q_bar`, `q_bar_per_bond`, `q_bar_normalize`,
`q_bar_correlation`, `initialize_interaction_data` type-4 branch)
must all be edited together as a single DRIFT repair.

## 1.13.e `q_bar_per_bond()`, `q_bar_normalize()`, `q_bar_correlation()` — DRIFT (inner pair, inherited)

All three methods inherit the 0-indexed `[real, imag]` inner pair
from Section 1.13.d. Same recommended fix — update alongside
`accumulate_q_bar`.

## 1.13.f `create_coordination_list()`, `create_coordination_summary()` — OK (presumed)

**Python** at `structure_control.py:7782` and 7846. Not deeply
audited but consume the 1-indexed `num_bonds` / `bonded` arrays
established by `create_bonding_list`.

**Verdict**: Presumed clean pending detailed audit.

## Section 1.13 verdict

4 OK (2 confirmed, 2 presumed), 1 DRIFT (`accumulate_q_bar` with
inner `[real, imag]` pair, plus three inherited sites in
`q_bar_per_bond`/`q_bar_normalize`/`q_bar_correlation` and the
`initialize_interaction_data` type-4 initialization).

The DRIFT is particularly tricky because the `[real, imag]` pair
is not an atom or axis index — it is a component pair that looks
naturally 0-indexed to a Python programmer. But Perl's `qBar`
explicitly uses slots `[1]` and `[2]`, leaving `[0]` unused, so
per the user's rule this must be preserved.
