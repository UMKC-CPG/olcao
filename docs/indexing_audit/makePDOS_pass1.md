# makePDOS.py — Pass 1 audit

## Scope and method

`src/scripts/makePDOS.py` (~1141 lines) is a Python port of the
Perl `makePDOS` script (~714 lines).  It reads a raw PDOS data
file produced by the OLCAO DOS program, applies a user-supplied
(or internally-generated default) set of filter/column
definitions, accumulates the selected PDOS contributions into
output curves, and writes a whitespace-delimited plot file.

Unlike most other scripts in this audit series, `makePDOS.py`
makes **zero calls** into `structure_control.py`.  The script
is entirely self-contained: no atom, axis, element, or species
indices are passed across an external API boundary.  The
indexing audit therefore reduces to three local questions:

1. Are the internal bookkeeping lists (`energy_points`,
   `format_names`, `col_defs`, `pdos_col_data`, `unit_data`,
   `input_col_labels`) 1-indexed with a slot-0 sentinel, as
   their Perl originals are?
2. Do the loop bounds walk `1..N+1` for 1-indexed containers?
3. Does the rc companion (`makePDOSrc.py`) expose any
   list-valued defaults that would cross a 1-indexed boundary?

Pass 1 approach:

- Full read of both the Perl original (`src/scripts/makePDOS`)
  and the Python port side-by-side.
- Identify every nested-list data structure and compare its
  Perl producer/consumer convention to the Python port.
- Verify all loop ranges and `len()` expressions against the
  Perl `$#array` idioms.
- Cross-check the rc file for list-valued defaults.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / rc loading | `ScriptSettings.*` | OK |
| B | Header reader | `PDOSData.init_pdos` | OK |
| C | Control-file reader | `PDOSData.read_control` | **1 DRIFT** (`filter_terms` inner axis) |
| D | Default-control builder | `PDOSData.default_control` | **1 DRIFT** (`col_def_entry[2]` / `[3]` inner axis) |
| E | Data reader / accumulator | `PDOSData.read_data` | OK (consumes drifts C/D self-consistently) |
| F | Name assigner | `PDOSData._assign_default_name` | **1 DRIFT** (`col_defs[col_idx][1]` inner label list is 0-indexed) |
| G | Output writer | `PDOSData.print_results` | OK (consumes drift F self-consistently) |
| H | rc companion | `makePDOSrc.py` | OK (scalars only) |

**Totals**: 0 BUGs, 3 DRIFT, 0 REVERSE DRIFT, 0 MIXED.  Every
drift is confined to the *innermost* axis of the `col_defs`
nested-list structure and is self-consistent between its
producer (`read_control` / `default_control` /
`_assign_default_name`) and its single consumer in `read_data`
/ `print_results`.  The script functions correctly as-is; the
fixes are pure rule compliance.

---

## Finding MPD-A — `filter_terms` inner axis in `read_control`

**Perl** `src/scripts/makePDOS` lines 286–332 (`sub
readControl`):

```perl
# Read the column definition on this line and store its filter parts in
#   the definition array.
@values = StructureControl::prepLine("",$line,':');
for ($i=0;$i<=$#values;$i++)
{
   $numFilterTerms = 0;
   @filterTermTemp = split(/\s+/,$values[$i]);
   ...
   for ($j=0;$j<=$#filterTermTemp;$j++)
   {
      ...
      if ($filterTermTemp[$j] =~ /\-/)
      {
         @range = split(/\-/,$filterTermTemp[$j]);
         for ($k=$range[0];$k<=$range[1];$k++)
         {
            $numFilterTerms++;
            $colDef[$numColDefs][$i+1][$numFilterTerms] = $k;
         }
      }
      else
      {
         $numFilterTerms++;
         $colDef[$numColDefs][$i+1][$numFilterTerms] =
               $filterTermTemp[$j];
      }
   }
}
```

`$numFilterTerms` starts at 0 and is pre-incremented before
every store, so the first term lands at slot **1** of
`$colDef[$numColDefs][$i+1]`.  The innermost axis is
1-indexed with slot 0 unused.  The outer two axes are also
1-indexed (`$numColDefs` starts at 0 and is pre-incremented
before use at line 284; `$i+1` explicitly shifts the format
position to 1-based at lines 318 and 328).

**Python** `makePDOS.py:611–663` (`read_control`):

```python
parts = line.split(':')
parts = [p.strip() for p in parts]

col_def_entry = [None]  # index 0 unused

for part in parts:
    filter_terms = []               # <-- 0-indexed inner
    term_tokens = part.split()

    for token in term_tokens:
        if not token:
            continue
        if '-' in token:
            range_parts = token.split('-')
            try:
                range_start = int(range_parts[0])
                range_end = int(range_parts[1])
                for k in range(range_start, range_end + 1):
                    filter_terms.append(str(k))
                continue
            except (ValueError, IndexError):
                pass
        filter_terms.append(token)

    col_def_entry.append(filter_terms)

self.col_defs.append(col_def_entry)
```

- `self.col_defs = [None]` at line 597 makes the **outer**
  axis 1-indexed.  OK.
- `col_def_entry = [None]` at line 618 makes the **middle**
  axis (format position) 1-indexed.  OK.
- `filter_terms = []` at line 623 makes the **innermost**
  axis 0-indexed.  **Drifted** from Perl.

**Consumer** (`read_data` at lines 910–921):

```python
acceptable = self.col_defs[col_idx][fmt_idx]
matched_term = False
for acc_val in acceptable:
    if (current_term_value == acc_val
            or acc_val == 'all'):
        matched_term = True
        break
```

The Python consumer iterates with `for acc_val in acceptable`
which walks from index 0.  The Perl consumer at line 502 uses
`for ($l=1;$l<=$#{$colDef[$j][$k]};$l++)` which walks 1..N.
Python is self-consistent with its own 0-indexed producer and
works correctly — but the convention drifted.

The same pattern is used for the requested column labels at
line 942:

```python
last_fmt = len(self.format_names) - 1
requested = self.col_defs[col_idx][last_fmt]

for req_label in requested:
    ...
```

Perl at line 523 walks
`for ($k=1;$k<=$#{$colDef[$j][$#format]};$k++)`.

**Severity**: DRIFT (inner axis only, self-consistent).

**Recommendation**: switch `filter_terms = []` to
`filter_terms = [None]` and adjust the two consumer loops in
`read_data` to walk `range(1, len(acceptable))` and
`range(1, len(requested))` respectively.  Affected touch
points:

- `makePDOS.py:623` — `filter_terms = []` → `[None]`.
- `makePDOS.py:917` — `for acc_val in acceptable:` →
  `for acc_idx in range(1, len(acceptable)):` +
  `acc_val = acceptable[acc_idx]`.
- `makePDOS.py:942` — `for req_label in requested:` →
  `for req_idx in range(1, len(requested)):` +
  `req_label = requested[req_idx]`.

---

## Finding MPD-B — `col_def_entry[2]` / `[3]` inner axis in `default_control`

**Perl** `src/scripts/makePDOS` lines 364–392 (`sub
defaultControl`):

```perl
for ($i=1;$i<=$numColDefs;$i++)
{
   if ($xanes==0)
   {
      $colDef[$i][2][1] = $i;
      $colDef[$i][3][1] = "TOTAL";
   }
   else
   {
      if ($i<=$numUnits)
      {
         $colDef[$i][2][1] = $i;
         $colDef[$i][3][1] = "s";
         $colDef[$i][3][2] = "d";
      }
      else
      {
         $colDef[$i][2][1] = $i-$numUnits;
         $colDef[$i][3][1] = "p";
      }
   }
}
```

Every store uses an explicit `[1]` or `[2]` innermost
subscript.  The innermost axis is 1-indexed with slot 0
unused.

**Python** `makePDOS.py:720–751` (`default_control`):

```python
for i in range(1, self.num_col_defs + 1):
    col_def_entry = [
        None,   # index 0 unused
        [],     # [1] NAME — filled in during read
        [],     # [2] SEQUENCE_NUM
        [],     # [3] COL_LABELS
    ]

    if not self.settings.xanes:
        col_def_entry[2] = [str(i)]
        col_def_entry[3] = ['TOTAL']
    else:
        if i <= self.num_units:
            col_def_entry[2] = [str(i)]
            col_def_entry[3] = ['s', 'd']
        else:
            col_def_entry[2] = [
                str(i - self.num_units)
            ]
            col_def_entry[3] = ['p']

    self.col_defs.append(col_def_entry)
```

- `self.col_defs = [None]` at line 718 — outer 1-indexed.  OK.
- `col_def_entry` prelude with `None` at slot 0 — middle
  1-indexed.  OK.
- The list literals `[str(i)]`, `['TOTAL']`, `['s', 'd']`,
  `['p']` at lines 736–749 all lack the slot-0 sentinel.
  These become the **innermost** axis consumed by
  `read_data`.  **Drifted** from Perl's `[1]` / `[1..2]`
  innermost layout.

Consumer is the same `for acc_val in acceptable` /
`for req_label in requested` loops at lines 917 and 942,
already flagged under MPD-A.  The drift is self-consistent
because `default_control` and `read_control` feed the same
consumers and both producers use the same 0-indexed inner
convention.

**Severity**: DRIFT (inner axis only, self-consistent).

**Recommendation**: prepend `None` to every innermost literal
in `default_control`, i.e. `[None, str(i)]`,
`[None, 'TOTAL']`, `[None, 's', 'd']`, `[None, 'p']`.  This
change must be applied in lockstep with MPD-A so the
consumers can be fixed once.  Affected touch points:

- `makePDOS.py:736` — `[str(i)]` → `[None, str(i)]`.
- `makePDOS.py:737` — `['TOTAL']` → `[None, 'TOTAL']`.
- `makePDOS.py:743` — `[str(i)]` → `[None, str(i)]`.
- `makePDOS.py:744` — `['s', 'd']` → `[None, 's', 'd']`.
- `makePDOS.py:746–748` — `[str(i - self.num_units)]` →
  `[None, str(i - self.num_units)]`.
- `makePDOS.py:749` — `['p']` → `[None, 'p']`.

---

## Finding MPD-C — `col_defs[col_idx][1]` label list inner axis

**Perl** `src/scripts/makePDOS` lines 553–563 (inside
`readData`, the default-name branch):

```perl
if ($style == 1)
{
   $colDef[$j][1][1] = $unitHash{"ELEMENT_NAME"} . "_" .
                       $unitHash{"SPECIES_ID"} . "_" .
                       $unitHash{"TYPE_ID"} . $tailID;
}
elsif ($style == 2)
{
   $colDef[$j][1][1] = $unitHash{"ELEMENT_1_NAME"} .
         "_" . $unitHash{"ELEMENT_2_NAME"} . $tailID;
}
```

And `printResults` at line 594:

```perl
for ($i=1;$i<=$numColDefs;$i++)
   {printf OUTPUT ("%12s",$colDef[$i][1][1]);}
```

The innermost label slot is explicitly `[1]` — 1-indexed.

**Python** `makePDOS.py:1026` (`_assign_default_name`):

```python
self.col_defs[col_idx][1] = [name]
```

and `makePDOS.py:1058` (`print_results`):

```python
for col_idx in range(1, self.num_col_defs + 1):
    label = self.col_defs[col_idx][1][0]
    out_fh.write(f"{label:>12s}")
```

The label is stored as a 1-element 0-indexed list `[name]`
and read back as `[1][0]`, whereas Perl stores `[1][1]` and
reads `[1][1]`.

**Severity**: DRIFT (inner axis only, self-consistent).

**Note**: the *external* control-file path does not set this
slot at all in the current Python port.  In the Perl original,
`readControl` would read the user-supplied label as the
`$filterTermTemp[0]` of the first colon-separated field and
store it at `$colDef[$numColDefs][1][1]`.  In Python,
`read_control` appends the first filter-terms list
(`filter_terms`) exactly the same way as every other format
position, so the label ends up as
`self.col_defs[col_idx][1] = ['label_value']` — i.e. in slot
`[1][0]` after the MPD-A drift.  This means the
user-control path happens to be consistent with
`_assign_default_name` and with `print_results`, but the
consistency is coincidental: it's the same drift applied
twice.

**Recommendation**: together with MPD-A / MPD-B, change the
store at line 1026 to

```python
self.col_defs[col_idx][1] = [None, name]
```

and the read at line 1058 to

```python
label = self.col_defs[col_idx][1][1]
```

After the MPD-A inner-axis fix, the user-control path will
automatically produce `[None, 'label_value']` at this slot
because the same `filter_terms` builder is used.

---

## Non-findings (verified consistent)

- **`energy_points`** — built at `init_pdos:530–533` as
  `[0.0]` sentinel + appended values, and read at
  `print_results:1068` as `self.energy_points[energy_idx]`
  with `energy_idx in 1..num_energy_points`.  Matches Perl
  `$energyPoint[$i]` (1-indexed).  **OK**.

- **`format_names`** — built at `read_control:583–589` as
  `['']` sentinel + appended, and at `default_control:705`
  directly as `['', 'NAME', 'SEQUENCE_NUM', 'COL_LABELS']`.
  Both branches give slot 0 = '' (matching Perl's unused
  slot 0, since Perl explicitly stores `$format[$i+1]`).
  Consumers: `read_data:890` uses
  `len(self.format_names) - 1` (equivalent to Perl's
  `$#format`) and the for-loop `range(2, num_formats)`
  matches Perl's `for ($k=2;$k<=$#format-1;$k++)`.  The
  last-format lookup at `read_data:939`
  `len(self.format_names) - 1` equals Perl's `$#format`.
  **OK**.

- **`col_defs`** **outer axis** — `read_control:597` and
  `default_control:718` both initialise with `[None]`.
  Consumers at `read_data:877` and `print_results:1057`
  walk `range(1, self.num_col_defs + 1)`.  Matches Perl's
  1-indexed `@colDef`.  **OK** (outer only; inner axes are
  findings MPD-A, MPD-B, MPD-C above).

- **`pdos_col_data`** — `read_data:795–800` builds outer as
  `[None]` + appended per-column lists, each of which is
  `[0.0]` + `num_energy_points` zeros.  Both axes are
  1-indexed.  Read at `read_data:960–964` as
  `self.pdos_col_data[col_idx][m]` with
  `col_idx in 1..num_col_defs` and
  `m in 1..num_energy_points`.  Read at `print_results:1074`
  as `self.pdos_col_data[col_idx][energy_idx]`.  Matches
  Perl `$pdosColData[$i][$j]` (both 1-indexed).  **OK**.

- **`unit_data`** — `read_data:850–859` builds outer as
  `[None]` + appended per-energy row lists, each of which is
  `['']` + appended floats.  Both axes 1-indexed.  Read at
  `read_data:962–964` as `unit_data[m][label_idx]` with
  `m in 1..num_energy_points` and
  `label_idx in 1..len(input_col_labels)-1`.  Matches Perl
  `$unitData[$j][$l]` (both 1-indexed, after Perl's
  per-line `unshift ""`).  **OK**.

  *(Side note: the Perl original has a latent quirk where
  `unshift ""` is inside the
  `foreach my $dataLine (1..$numLines)` loop, so for
  multi-line column blocks it prepends one empty string per
  line.  The Python port correctly pre-allocates `['']` once
  per energy point, which produces the intended 1-indexed
  row regardless of how many lines the block spans.  Not an
  audit finding — the Python behaviour is strictly more
  correct.)*

- **`input_col_labels`** — `read_data:812` pre-allocates
  `['']` as the slot-0 sentinel; `read_data:836–839`
  appends the labels read from the file.  Read at
  `read_data:945–947` as
  `range(1, len(input_col_labels))`.  Matches Perl's
  1-indexed `@inputColLabels`.  **OK**.

- **`unit_hash`** — plain dictionary keyed by filter name.
  No indexing convention applies.  **OK**.

- **Per-unit reset** — Python relies on `unit_hash = {}`,
  `input_col_labels = ['']`, and `unit_data = [None]` being
  reinitialised at the top of each unit iteration
  (`read_data:811–812`, `850`).  Perl uses `undef` at the
  bottom of the loop (lines 577–579).  Both achieve the
  same reset.  **OK**.

- **`num_units`, `num_energy_points`, `num_col_defs`,
  `style`** — plain integer scalars; no indexing convention
  applies.  **OK**.

- **Loop bound sweep** — every `range()` expression in the
  port:

  | Line | Expression | Container | Perl equivalent | OK? |
  |---|---|---|---|---|
  | 531 | `range(self.num_energy_points)` | counter for appending to `energy_points` | `for ($i=1;$i<=$numEnergyPoints;$i++)` | OK (same net count) |
  | 645 | `range(range_start, range_end + 1)` | filter-term range expansion | `for ($k=$range[0];$k<=$range[1];$k++)` | OK |
  | 720 | `range(1, self.num_col_defs + 1)` | iterate `col_defs` | `for ($i=1;$i<=$numColDefs;$i++)` | OK |
  | 796 | `range(1, self.num_col_defs + 1)` | allocate `pdos_col_data` | `for ($i=1;$i<=$numColDefs;$i++)` | OK |
  | 798 | `range(self.num_energy_points)` | counter for appending to `col_data` | `for ($j=1;$j<=$numEnergyPoints;$j++)` | OK (same net count) |
  | 803 | `range(1, self.num_units + 1)` | iterate units | `for ($i=1;$i<=$numUnits;$i++)` | OK |
  | 836 | `range(num_lines)` | counter for file lines | `foreach my $dataLine (1..$numLines)` | OK (same net count) |
  | 851 | `range(1, self.num_energy_points + 1)` | iterate `unit_data` | `for ($j=1;$j<=$numEnergyPoints;$j++)` | OK |
  | 855 | `range(num_lines)` | counter for file lines | `foreach my $dataLine (1..$numLines)` | OK (same net count) |
  | 877 | `range(1, self.num_col_defs + 1)` | iterate col definitions | `for ($j=1;$j<=$numColDefs;$j++)` | OK |
  | 892 | `range(2, num_formats)` with `num_formats = len - 1` | walks format 2..N-1 | `for ($k=2;$k<=$#format-1;$k++)` | OK |
  | 945 | `range(1, len(input_col_labels))` | iterate input labels | `for ($l=1;$l<=$#inputColLabels;$l++)` | OK |
  | 956 | `range(1, self.num_energy_points + 1)` | iterate energy bins | `for ($m=1;$m<=$numEnergyPoints;$m++)` | OK |
  | 1057 | `range(1, self.num_col_defs + 1)` | iterate col labels for header | `for ($i=1;$i<=$numColDefs;$i++)` | OK |
  | 1063 | `range(1, self.num_energy_points + 1)` | iterate energy rows | `for ($i=1;$i<=$numEnergyPoints;$i++)` | OK |
  | 1071 | `range(1, self.num_col_defs + 1)` | iterate output columns | `for ($j=1;$j<=$numColDefs;$j++)` | OK |

  No loop-bound drift.

- **rc companion (`makePDOSrc.py`)** — returns a dict with
  only scalar values (three strings, two booleans).  No
  list-valued defaults cross a 1-indexed boundary.  **OK**.

- **No `structure_control` calls** — verified by a full
  read of the port and a grep for `sc.` / `StructureControl`
  / `structure_control`.  `makePDOS.py` imports only
  `argparse`, `math`, `os`, `sys`, `datetime`, and the rc
  module.  No external indexing contract exists.

## Non-indexing observations (out of audit scope)

- **`math` import unused** — `import math` at line 135 is
  never referenced.  Harmless but dead.

- **`col_def_entry` prelude pre-fill vs. overwrite** — in
  `default_control:726–731` the slots `[1]`, `[2]`, `[3]`
  are first allocated as empty `[]` lists and then
  overwritten by the filter-term assignments on the next
  few lines.  The empty list at `[1]` (NAME) is never
  populated by `default_control` itself — it is filled by
  `_assign_default_name` only when a matching unit is
  found in `read_data`.  If no unit matches (pathological
  input, e.g. `num_units = 0`), `print_results` will crash
  on `label = self.col_defs[col_idx][1][0]` with an
  IndexError.  Not reachable from normal usage but worth a
  comment or a guard.

- **Perl `$#format - 1` vs. Python `len - 1`** — both
  correctly exclude the last format slot (COL_LABELS) from
  the unit-matching loop.  The equivalence is exact:
  Perl's `$#format` is the highest 1-indexed slot, Python's
  `len(self.format_names) - 1` is the same quantity when
  slot 0 is a sentinel.  Verified.

- **`prep_line` helper** — local reimplementation of
  `StructureControl::prepLine`.  Returns a standard
  0-indexed Python list (matching what Perl's `prepLine`
  returns in list context).  This 0-indexed return is a
  legitimate parse buffer; every caller treats the result
  as an ordinary list.  **OK**, matches Perl.

## Running tally

- **BUG**: 0
- **DRIFT**: 3
  - MPD-A (`filter_terms` inner axis in `read_control`)
  - MPD-B (`col_def_entry[2]` / `[3]` inner axis in
    `default_control`)
  - MPD-C (`col_defs[col_idx][1]` label list inner axis in
    `_assign_default_name` and `print_results`)
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every other subsystem (rc companion, loop bounds,
  outer/middle axes of all nested lists, all of the
  scalar / dict / flat-list buffers).

## Fix-order recommendation

Single cluster, single file, minimal coordinated edit:

1. **Cluster MPD-A+B+C** — innermost axis of `col_defs`.
   The three drifts share a common consumer in `read_data`
   (the two `for … in acceptable` / `for … in requested`
   loops) and a common print-site in `print_results` (the
   `label = ...[1][0]` lookup), so they must be fixed
   together.  Touch points:

   - `read_control:623` — `filter_terms = []` → `[None]`.
   - `default_control:736,737,743,744,746–748,749` — prepend
     `None` to each innermost literal.
   - `_assign_default_name:1026` — `[name]` → `[None, name]`.
   - `read_data:917` — switch from `for acc_val in
     acceptable:` to an index-based walk starting at 1.
   - `read_data:942` — switch from `for req_label in
     requested:` to an index-based walk starting at 1.
   - `print_results:1058` — read `[1][1]` instead of
     `[1][0]`.

   No external callers (the script is a standalone
   post-processor).  No test coverage exists for
   `makePDOS.py` in `tests/`; regression verification after
   the fix will need to be by inspection plus a smoke run
   against a real `gs_dos-fb.p.raw` fixture if one is
   available.

## Fix status

**Not yet applied.**  This file is Pass 1 audit output only;
no code has been modified.
