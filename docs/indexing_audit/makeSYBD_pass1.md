# makeSYBD.py — Pass 1 Audit

**Status: zero findings.** Joins `pdb2skl.py`, `make_reactions.py`
as a clean port.

## Summary

| Severity | Count |
|---|---|
| BUG | 0 |
| DRIFT | 0 |
| REVERSE DRIFT | 0 |
| MIXED | 0 |
| OK | seven clusters (A through G) |

## Scope

`makeSYBD.py` (~581 lines) is a Python port of Perl `makeSYBD`
(~210 lines). Post-processes three OLCAO SYBD files (`.dat`,
`.out`, `.raw`) into a single annotated `.plot`. The script has
**no calls into `structure_control.py`** at all — no atoms / axes
/ elements / species touchpoints anywhere. All bookkeeping is pure
file-parse scratch.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings | `ScriptSettings.*` | OK |
| B | K-point name reader | `read_kpoint_names` | OK |
| C | K-point position reader | `read_kpoint_positions` | OK |
| D | Path-break merger | `merge_path_breaks` | OK |
| E | Raw→plot writer | `write_plot_file` | OK |
| F | `main()` orchestrator | `main` | OK |
| G | rc companion | `makeSYBDrc.parameters_and_defaults` | OK |

## Detailed verification

### Section A — `ScriptSettings` (Python 73–280)

Six scalars (four file-name strings plus `basis`, `edge` tags). No
lists, no indices, no cross-script calls. Perl `makeSYBD:24–78`.
**OK.**

### Section B — `read_kpoint_names` (Python 287–338; Perl 86–114)

`@KNames` in Perl is genuinely 0-indexed: it starts empty, is grown
by `push`, and the later merge loop reads `$KValues[$kp]` vs
`$KValues[$kp-1]` starting at `$kp=1` (slot 1 vs slot 0), which
only makes sense if slot 0 holds real data. This is one of the
rare Perl cases of native 0-indexing (no atom/axis key → no
sentinel convention). Python uses a native 0-indexed list via
`.append()`. Matches Perl exactly. **OK.**

### Section C — `read_kpoint_positions` (Python 345–408; Perl 122–159)

Same story as `@KNames`. `k_values` is populated via `.append()`
and then overwritten in place using `current = 0..len(k_values)-1`.
`tokens[1]` and `tokens[2]` are legitimate 0-indexed `split()`
reads, matching Perl `$tempArray[1]` / `$tempArray[2]`. **OK.**

### Section D — `merge_path_breaks` (Python 415–458; Perl 167–181)

Perl: `foreach my $kp (1..$numHighSymKP-1)` comparing `$KValues[$kp]`
with `$KValues[$kp-1]`, guarded by `if ($kp > $#KValues) last`.
Python mirrors with `kp = 1; while kp < len(k_values)`. Slot-0
semantics preserved on both sides. **OK.**

**Non-indexing observation (out of audit scope)**: Perl's `foreach`
over a fixed integer range unconditionally increments `$kp` each
pass, so after a `splice` the element that slides into the just-
vacated slot is **skipped**. Python's `while` deliberately does
*not* increment after `pop()`, so cascaded merges are possible.
With three or more consecutive identical positions the two scripts
produce different labels. Not an indexing drift — a logic
divergence worth flagging for the developer.

### Section E — `write_plot_file` (Python 465–536; Perl 183–209)

Outer loop: `range(num_kpoints)` ≡ Perl `$i=0..$NumKPoints-1`
(0-indexed). Inner continuation-line loop:
`range(num_cont_lines)` where `num_cont_lines = num_states // 11`
≡ Perl `$j=1..$NumStates/11` — identical iteration count (the body
never references `$j`, it's a pure counter). `tokens[1]`/`tokens[2]`
header reads match Perl's `$tempArray[1]`/`$tempArray[2]`. **OK.**

### Section F — `main()` (Python 543–571)

Three-call orchestrator. Passes `k_names` into
`read_kpoint_positions` only for its length (via `len(k_names)`),
then parallel `k_names`/`k_values` through `merge_path_breaks` and
`write_plot_file`. No `sc.*` calls. **OK.**

### Section G — `makeSYBDrc.py` (rc companion)

Six scalar string settings (`dat_file`, `out_file`, `raw_file`,
`plot_file`, `basis`, `edge`). No list-valued defaults. Matches
the `make_reactionsrc.py` / `condenserc.py` / `pdb2sklrc.py` pure-
scalar pattern. **OK.**

## Non-findings verified

- No `sc.*` calls anywhere — `structure_control` is never imported.
- `tokens[-1]`, `tokens[1]`, `tokens[2]` reads are legitimate
  0-indexed `split()` results matching Perl `@tempArray`.
- All four loops walk genuinely 0-indexed structures; no
  `range(1, N+1)` is needed because no 1-indexed container exists.
- rc has no list defaults; only scalar strings.

## Non-indexing observations (out of audit scope)

1. `merge_path_breaks` iteration-semantics divergence
   (Section D): Perl skips, Python cascades. Different output
   only for 3+ consecutive equal positions.
2. `GAMMA` detection: Python's `if/else` avoids a latent Perl
   double-push bug where a line contains "GAMMA" as a substring
   but isn't the gamma label. Improvement, not drift.
3. `num_states_orig` vs `num_states` naming: Python cleanly
   separates raw header value from rounded-up line-count.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every subsystem (A–G)

## Fix order

None required.
