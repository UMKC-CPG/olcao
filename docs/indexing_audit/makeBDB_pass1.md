# makeBDB.py — Pass 1 Audit

**Status: 3 DRIFT, 0 BUG.** All three clusters are locally self-
consistent (producer/consumer agree); no runtime hazard.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 3 (MBDB-A, MBDB-B, MBDB-C) |
| REVERSE DRIFT | 0 |
| MIXED | 0 |

## Scope

- Perl original: `src/scripts/makeBDB`
- Python port:   `src/scripts/makeBDB.py`
- rc companion:  `src/scripts/makeBDBrc.py`

No `structure_control` dependency in `makeBDB.py` — the single
Perl `StructureControl::getAtomicBDB()` call was replaced with a
direct `$OLCAO_DATA/atomicBDB` lookup in `get_database_path()`.

## Finding MBDB-A — `orbital_qn_n` / `up_charge` / `dn_charge` inner per-orbital axis — DRIFT

**Perl** `makeBDB:582-607`: `$orbitalQN_n[$qn_l][$orbital]` walks
`$orbital (1..numValeOrbitals)` with slot 0 unused.

**Python** `makeBDB.py:512-556` (`get_element_data`): uses
`.append()` so the inner per-orbital axis is 0-indexed.

**Consumer**: `make_atom_scf` (`makeBDB.py:760-780`) walks
`range(n_orb)` — self-consistent with producer, no runtime hazard.

**Severity**: DRIFT.

**Recommended fix**: seed each inner orbital list with `[None]`
sentinel in `get_element_data`; update consumer to walk
`range(1, n_orb + 1)`.

## Finding MBDB-B — `contract_core` / `contract_vale` outer file-variant axis — DRIFT

**Perl** `makeBDB:508-524`: `$coreOrbitals[1]` / `[2]` and
`$valeOrbitals[1..2][...]` are 1-indexed on the file-variant axis
(`[1]` = with-core, `[2]` = nocore).

**Python** `makeBDB.py:565-570` (returned from
`get_element_data`): uses 0-indexed outer file-variant axis
(`[0]` = with-core, `[1]` = nocore).

**Severity**: DRIFT. Must be fixed together with MBDB-C (same
axis).

## Finding MBDB-C — `file_names` outer file-variant axis — DRIFT

**Perl** `makeBDB:409-423`: `@fileName[1..2]` + `foreach $file
(1..2)` — 1-indexed.

**Python** `makeBDB.py:955-984` (`make_contracts`): `file_names`
plus `for file_idx in range(2)` — 0-indexed.

**Severity**: DRIFT. Must be fixed together with MBDB-B.

**Recommended fix for B + C**: seed `file_names = [None, ..., ...]`
and the `contract_core` / `contract_vale` outer lists with `[None]`;
update the `make_contracts` loop to `range(1, 3)`; update all
consumers that read `[0]` / `[1]` to read `[1]` / `[2]`.

## Verified OK subsystems

- Element axis (outer-most loop over elements): 1-indexed.
- `element_data` reads: 1-indexed as shipped.
- `qn_l (0..3)` axes: legitimately 0-indexed in both Perl and
  Python (quantum number is the literal index).
- `basis (1..3)` dict keys: 1-indexed on both sides.
- `ORBITAL_TAGS`: legitimately 0-indexed constant table.
- `makeBDBrc.py`: scalar defaults only, OK.

## Fix order

1. **MBDB-A** (standalone): seed inner orbital lists with `[None]`,
   update one consumer loop.
2. **MBDB-B + MBDB-C** (must be fixed together): seed file-variant
   outer axes, rewrite the `range(2)` loop to `range(1, 3)`,
   update all `[0]` / `[1]` reads to `[1]` / `[2]`.

Validation: inspection + smoke run of `makeBDB.py` to regenerate
the basis database.

## Running tally

- **BUG**: 0
- **DRIFT**: 3 (MBDB-A, MBDB-B, MBDB-C)
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every other subsystem
