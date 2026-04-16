# processPOPTC.py — Pass 1 audit

## Scope and method

`src/scripts/processPOPTC.py` (~892 lines) is a Python port of the
Perl `processPOPTC` script (~544 lines).  The program reads the raw
partial Epsilon-2 data file produced by an OLCAO poptc calculation
(`fort.250` or `fort.251`), extracts each partial pair's eps2 through
`makePDOS`, feeds it to `OLCAOkkc` for the Kramers-Kronig conversion,
and appends the resulting derived spectra (eps1, ELF, refractive
index, extinction coefficient, eps1i, reflectivity, absorption
coefficient) into eight accumulated "raw" output files.

The script is a thin orchestrator around two external binaries
(`makePDOS` and `OLCAOkkc`).  It does **not** call into
`structure_control.py` at all — the only `StructureControl` reference
in the original Perl is `StructureControl::prepLine`, which the
Python port has reimplemented locally as the module-level `prep_line`
helper.  Consequently there is no `sc.*` call boundary to audit, and
the audit reduces to internal bookkeeping lists, loop bounds, the
public-API parameters of the `POPTCData` methods, and the `rc`
companion's defaults.

Pass 1 approach:

- Read the entire Perl original (`src/scripts/processPOPTC`, 544
  lines) end-to-end.
- Read the entire Python port (`src/scripts/processPOPTC.py`, 892
  lines) end-to-end.
- Read the `rc` companion (`src/scripts/processPOPTCrc.py`) —
  it exposes a single scalar (`spin`) and has no list-valued
  defaults.
- Grep-verify there are no hidden `structure_control` / `sc.*`
  call sites (confirmed: zero).
- Compare every 1-indexed Perl array against its Python
  counterpart for slot-0 sentinel, loop bounds, and consumer
  reads.

## Summary table

| # | Topic | Functions | Verdict |
|---|-------|-----------|---------|
| A | CLI / settings | `ScriptSettings.*` | OK |
| B | File-name setup | `_init_file_names` | OK |
| C | Metadata reader | `get_eps2_raw_metadata` | OK |
| D | KKC control reader | `get_kkc_control_data` | OK |
| E | Header writer | `write_header` | OK |
| F | Per-sequence control / KKC drivers | `create_pdos_control`, `call_make_pdos`, `call_olcao_kkc` | OK |
| G | Per-sequence data copier | `copy_data` | OK |
| H | Main loop | `process_all_sequences`, `main` | OK |
| I | `prep_line` helper | `prep_line` | OK (file-handle bridge, values returned 0-indexed — matches Perl's `@values` return) |
| J | rc companion | `processPOPTCrc.py` | OK (single scalar default) |

**Zero findings.**  No BUG, no DRIFT, no REVERSE DRIFT, no MIXED.
The port is clean with respect to the 1-indexing convention.

## Verification details

### PP-A — `energy_points`, `pair_kkc_factor` (1-indexed, OK)

**Perl** `processPOPTC:67, 72, 225–229, 278–282`:

```perl
my @energyPoint;
my @pairKKCFactor;
...
for ($i=1;$i<=$numEnergyPoints;$i++)
{
   @values = StructureControl::prepLine(\*EPS2RAW,$line,'\s+');
   $energyPoint[$i] = $values[0];
}
...
for ($i=1;$i<=$numSequences;$i++)
{
   @values = StructureControl::prepLine(\*KKCCON,$line,'\s+');
   $pairKKCFactor[$i] = $values[$#values];
}
```

Both arrays are 1-indexed with slot 0 left undef.

**Python** `processPOPTC.py:362, 367, 476–478, 522–524`:

```python
self.energy_points = [None]       # Index 0 unused.
self.pair_kkc_factor = [None]     # Index 0 unused.
...
for i in range(1, self.num_energy_points + 1):
    values = prep_line(f)
    self.energy_points.append(values[0])
...
for i in range(1, self.num_sequences + 1):
    values = prep_line(f)
    self.pair_kkc_factor.append(values[-1])
```

Slot 0 is sentinelled with `None`.  Loops walk `1..N+1`.  Appends
produce entries at slots `1..N`.  Matches Perl exactly.  **OK.**

### PP-B — `element1`, `element2`, `type1`, `type2` (1-indexed, OK)

**Perl** `processPOPTC:68–71, 231–265`:

```perl
my @element1;
my @element2;
my @type1;
my @type2;
...
$element1[0] = "";
$element2[0] = "";
$type1[0] = "";
$type2[0] = "";
...
$j = 1;
while ($line=<EPS2RAW>)
{
   if ($line =~ /ELEMENT_1_NAME/)
   {
      ...
      $element1[$j]=$line;
      ...
      $element2[$j]=$line;
      ...
      $type1[$j]=$line;
      ...
      $type2[$j]=$line;
      $j = $j + 1;
   }
}
```

Four parallel 1-indexed arrays with explicit empty-string
sentinels at slot 0, populated by an incrementing `$j` that
starts at 1.

**Python** `processPOPTC.py:363–366, 486–505`:

```python
self.element1 = [""]         # Index 0 unused.
self.element2 = [""]         # Index 0 unused.
self.type1 = [""]            # Index 0 unused.
self.type2 = [""]            # Index 0 unused.
...
j = 1
for line in f:
    if "ELEMENT_1_NAME" in line:
        self.element1.append(line.strip())
        ...
        self.element2.append(line2)
        ...
        self.type1.append(line3)
        ...
        self.type2.append(line4)
        j += 1
```

Slot 0 is sentinelled with `""` (matching Perl's `""`
assignment), the `j` counter is preserved for fidelity with the
Perl, and `append` correctly lands each new record at slot `j`
(because the list has `j` elements before the append, the first
append lands at slot 1).  **OK.**

Note: `j` is no longer read after the loop (its value is never
used to size a downstream allocation), but keeping it mirrors the
Perl and costs nothing.

### PP-C — Consumer reads in `copy_data` (OK)

**Perl** `processPOPTC:408–412`:

```perl
print $fileHandleRaw ($element1[$i],"\n");
print $fileHandleRaw ($element2[$i],"\n");
print $fileHandleRaw ($type1[$i],"\n");
print $fileHandleRaw ($type2[$i],"\n");
```

where `$i` is the main loop variable running `1..$numSequences`.

**Python** `processPOPTC.py:700–714`:

```python
file_handle.write(f"SEQUENCE_NUM {seq_num}\n")
file_handle.write(f"{self.element1[seq_num]}\n")
file_handle.write(f"{self.element2[seq_num]}\n")
file_handle.write(f"{self.type1[seq_num]}\n")
file_handle.write(f"{self.type2[seq_num]}\n")
```

`seq_num` is the loop variable from `process_all_sequences` which
walks `range(1, self.num_sequences + 1)` — directly 1-indexed,
matches Perl's `1..$numSequences`.  **OK.**

### PP-D — Consumer read of `pair_kkc_factor[seq_num]` (OK)

**Perl** `processPOPTC:392`:

```perl
system("$OLCAO_BIN/OLCAOkkc $optcLines 1 1 $pairKKCFactor[$i]");
```

**Python** `processPOPTC.py:656–664`:

```python
subprocess.run(
    [
        kkc_cmd,
        str(optc_lines),
        "1",
        "1",
        str(self.pair_kkc_factor[seq_num]),
    ],
    check=True,
)
```

1-indexed index drawn from the outer `process_all_sequences` loop.
**OK.**

### PP-E — Consumer read of `energy_points[i]` in `write_header` (OK)

**Perl** `processPOPTC:298–299`:

```perl
for ($i=1;$i<=$numEnergyPoints;$i++)
   {printf $fileHandle ("%10.4f\n",$energyPoint[$i]);}
```

**Python** `processPOPTC.py:553–556`:

```python
for i in range(1, self.num_energy_points + 1):
    file_handle.write(
        f"{float(self.energy_points[i]):10.4f}\n"
    )
```

Loop bounds match exactly (`1..N+1`).  Read uses the 1-indexed
slot.  **OK.**

### PP-F — Main sequence loop (OK)

**Perl** `processPOPTC:161`:

```perl
for ($i=1;$i<=$numSequences;$i++)
```

**Python** `processPOPTC.py:787`:

```python
for i in range(1, self.num_sequences + 1):
```

Identical 1-indexed walk.  **OK.**

### PP-G — `prep_line` return convention (OK, matches Perl)

**Perl** `StructureControl::prepLine`: returns an ordinary Perl
array (`@values`) that is 0-indexed (the first token sits at
`$values[0]`).  Every caller in the Perl original reads either
`$values[0]` (first token) or `$values[$#values]` (last token) —
both 0-indexed accesses on the return value.

**Python** `processPOPTC.py:73–120`: `prep_line` returns a plain
Python list that is 0-indexed.  Every caller in the file reads
`values[0]` or `values[-1]`.

The return value being 0-indexed is **correct** — Perl's
`@values` is genuinely 0-indexed in this context (it is a
tokenized split buffer, not an atom/axis/element array), so a
0-indexed Python list is the faithful translation.  This matches
the legitimate-0-indexed-parse-buffer category documented in the
SC audit (e.g. `bondBOList` in `bond_analysis.py`).  **OK.**

### PP-H — `_init_file_names` (OK)

The file-name attribute setup at `processPOPTC.py:376–435` uses
a simple `if spin == 1 / else` branching and assigns scalar string
attributes directly.  No list indexing of any kind.  Matches the
Perl `if ($ARGV[0] == 1)` branches at lines 88–128.  **OK.**

### PP-I — `copy_data` internal `values` slice (OK)

**Perl** `processPOPTC:422–429`:

```perl
while ($line=<TEMPIN>)
{
   chomp ($line);
   @values = split (/\s+/,"$line");
   shift (@values);
   shift (@values);
   printf $fileHandleRaw ("@values\n");
}
```

The two `shift` calls drop the first two 0-indexed tokens of the
split result — a 0-indexed parse buffer exactly like `prep_line`.

**Python** `processPOPTC.py:731–739`:

```python
for line in temp_in:
    values = line.split()
    if len(values) > 2:
        values = values[2:]
    file_handle.write(
        " ".join(values) + "\n"
    )
```

Drops the first two 0-indexed tokens via slice.  Legitimate
0-indexed parse buffer.  **OK.**

## Non-findings (verified consistent)

- **No `sc.*` call sites.**  `Grep` over `processPOPTC.py` for
  `sc\.|StructureControl|structure_control` turns up only a
  single docstring mention (line 77) inside the `prep_line`
  helper comment.  There is nothing to bridge, nothing to pass
  1-indexed through a call boundary.

- **`num_sequences`, `num_energy_points`, `style`** are scalars;
  no indexing concerns.

- **Loop variable `j` in `get_eps2_raw_metadata`** is tracked
  only for Perl fidelity; it is never read after the loop, and
  incrementing it 1..N+1 on append is correct.

- **`raw_files` dict in `process_all_sequences`** keyed by
  string labels ("e1e2elf", "eps1", "elf", …), not by any
  atom/axis/element/species index.  Matches the Perl file
  handles.  Not subject to the 1-indexing rule.

- **`temp_files` cleanup list at lines 838–848** is a plain
  0-indexed iteration list used only as a `for x in list:` loop
  — standard Python flat list, matches the Perl shell command
  string of the same names.  No 1-indexing semantics.

- **`prep_line` return value** is legitimately 0-indexed (parse
  buffer), matching Perl.  See PP-G.

## rc companion audit

**`src/scripts/processPOPTCrc.py`** (45 lines):

```python
param_dict = {
    "spin": 1,  # Integer
}
```

Single scalar entry.  No list-valued defaults cross a 1-indexed
boundary.  No `structure_control` call sites anywhere in the
companion.  **OK.**  Falls into the same category as
`pdb2sklrc.py`, `condenserc.py`, and `make_reactionsrc.py` in
the README table.

## Running tally

- **BUG**: 0
- **DRIFT**: 0
- **REVERSE DRIFT**: 0
- **MIXED**: 0
- **OK**: every subsystem.

## Fix-order recommendation

None required.  The port is clean.  No test coverage exists for
`processPOPTC.py` itself (the script exercises external binaries
`makePDOS` and `OLCAOkkc` and would need a full OLCAO run to
validate end-to-end), but the 1-indexed bookkeeping conventions
are preserved faithfully throughout and no `structure_control`
call boundary exists to drift against.

## Notes for later passes

- **Pass 2 (cross-method calls)**: There are no internal
  cross-method calls that pass atom/axis/element/species
  indices — only a bare sequence number `seq_num` (or the loop
  variable `i`) passed to `create_pdos_control`,
  `call_olcao_kkc`, and `copy_data`.  All three consumers read
  the 1-indexed data arrays at `[seq_num]` correctly.  Nothing
  further to audit.
- **Pass 3 (slot-0 sentinels)**: Every 1-indexed list is
  initialised with an explicit slot-0 sentinel (`[None]` or
  `[""]`).  Already verified in PP-A and PP-B above.
- **Pass 4 (loop bounds)**: All three `for i in range(...)`
  loops that index 1-indexed data walk
  `range(1, N + 1)`.  Already verified in PP-A, PP-E, PP-F
  above.

The file can be marked complete after Pass 1 — no follow-on
passes are needed.
