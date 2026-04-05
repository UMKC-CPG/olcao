# OLCAO Project

> **Session setup:** Run `/color blue` and `/rename OLCAO` at the
> start of each session.

## Document Hierarchy

This project follows a five-level document chain. All design documents
live in `dev/`. Read them in order when starting work:

1. `dev/VISION.md` -- goals, principles, non-negotiables
2. `dev/ARCHITECTURE.md` -- layout, modules, dependencies
3. `dev/DESIGN.md` -- algorithms, data structures, math
4. `dev/PSEUDOCODE.md` -- language-agnostic algorithm specs
5. Source code in `src/`

`dev/TODO.md` tracks tasks organized by level. Use `/focus` to start a
session and `/refine` to check consistency across the chain.

## Level Awareness

During development conversations, the programmer may shift between levels
of the design chain without explicitly noticing. For example, a discussion
about a code fix may drift into questioning an algorithm's design, or a
design discussion may surface a conflict with a core principle.

When you notice the conversation has moved to a different level than where
it started, say so briefly. For example: "This sounds like it's becoming
an ARCHITECTURE question -- should we capture it there before continuing
with the code?" The goal is awareness, not interruption. Let the programmer
decide whether to switch context, propagate the change to the appropriate
document, or stay focused and defer.

Do not enforce rigid boundaries. The levels exist to organize thinking, not
to prevent it. A developer who is on a productive train of thought should
not be stopped -- but when the thought resolves, help them recognize which
documents it touches so nothing is left inconsistent.

## Coding Style

### Line Length

Lines MUST NOT exceed 80 characters. This is a hard limit. At the same
time, every line MUST be filled to at least 70 characters whenever the
content allows it. A 45-character line when 75 characters were available
is a formatting error just as serious as exceeding 80. Target the
70-80 character band on every line; treat both the floor and the
ceiling as mandatory.

Short statements whose content is naturally brief (`implicit none`,
`endif`, `integer :: i`) are fine at their natural length -- there is
nothing to fill them with. The rule applies when content is available:
long expository comments, argument lists, complex expressions, etc.

**Common failure mode -- do not do this:**
```fortran
   ! Guard each deallocation because
   !   this subroutine is called from
   !   multiple program paths.
```
Those three ~35-character lines have plenty of content to fill longer
lines. Write them as:
```fortran
   ! Guard each deallocation because this subroutine is called from multiple
   !   program paths.
```
The idea: let each line run toward 80 before wrapping. Do not break
at natural phrase boundaries when the line is only at 40-50 characters.
The result will be fewer, fuller lines -- not lines padded with filler.

### Documentation and Naming

CRITICAL: All program code must include rich, expressive documentation
so that students reading the source can easily follow what is happening.
Every function, class, and non-trivial block should carry clear
explanatory comments or docstrings that describe purpose, inputs,
outputs, and any relevant physics or math.

Variable names must be readable and self-documenting. Avoid cryptic
one- or two-letter abbreviations. Prefer concise but meaningful names
that a reader can understand without cross-referencing a legend.
Slightly-too-long names are far better than opaque short ones.

Good naming examples:
- `elec_mom` instead of `em` or `electron_momentum`
- `nuc_pot` instead of `np` or `nuclear_potential`
- `grid_spacing` instead of `gs` or `the_spacing_between_grid_points`

The goal is a middle ground: short enough to keep expressions tidy,
long enough that any student can read the code cold and follow the
logic without guessing what a variable holds.

## What This Is

OLCAO (Orthogonalized Linear Combination of Atomic Orbitals) is an academic
electronic structure code applicable to a broad range of material systems:
crystals, amorphous solids, nanoparticles, molecules, interfaces, grain
boundaries, and more. Key characteristics:
- **All-electron**: no pseudopotentials; core electrons are treated
  explicitly
- **Periodic boundary conditions**: used throughout, even for non-periodic
  systems
- **LCAO basis**: wavefunctions expressed as linear combinations of atomic
  orbitals
- **Orthogonalization**: valence orbitals are orthogonalized to the core,
  reducing the size of the secular (eigenvalue) equation
It has been in active development for many years.

## Working Convention

Use `/focus <name>` to load context for a specific subprogram or script
before working on it. This avoids reading the whole codebase unnecessarily.

## Repository Layout

```
src/                  Source code
  uolcao/             Primary production program (unified/HDF5-based)
  atomSCF/            Atomic SCF (generates basis/potential)
  makeKPoints/        k-point mesh generation
  gaussFit/           Gaussian fitting of charge density
  contract/           Basis set contraction
  applySpaceGroup/    Space group symmetry operations
  scripts/            Utility scripts (Perl, Python, bash)
  kinds.f90           Shared: precision kinds
  constants.f90       Shared: physical/mathematical constants
  elementData.f90     Shared: periodic table data
  readData.f90        Shared: input file parsing
  writeData.f90       Shared: output file writing
  radialGrid.f90      Shared: radial grid definitions
bin/                  Installed executables and scripts
build/                CMake out-of-source build directory
docs/                 Documentation
dev/                  Design document chain
share/                Basis set and potential databases
```

Note: `polcao/` is under active development and is a near-term priority to
enable. `olcao/` and `upolcao/` also exist in src but are currently disabled
in the build.

## Language and Build

- Fortran 90 for all programs; Python and Perl for scripts
- CMake build system; out-of-source build in `build/`
- Compiler set via `$FC` env var (must match the HDF5-compiled Fortran
  compiler)
- Install prefix set via `$OLCAO_DIR` env var
- HDF5 is used for primary I/O in `uolcao`
- Supported compilers: gfortran, ifort

## Perl-to-Python Conversion Status

- **StructureControl.pm -> structure_control.py**: Essentially complete
  (9,335 lines, 171 tests passing). One minor gap remains:
  `add_connection_atoms` for 2nd-4th tetrahedral connection atoms (rarely
  used). All other methods fully implemented.
- **makeinput -> makeinput.py**: Next target for conversion.

## Key Conventions

- All Fortran files use `implicit none` and the shared `kinds.f90`
  precision types
- Shared modules in `src/` are compiled once and used across subprograms
- Scripts are standalone tools; many have associated Perl `.pm` library
  files

## Documentation Policy

This is an academic codebase used by students who frequently need to read
and understand the source code.  When converting, refactoring, or writing
any code (scripts, Fortran programs, Python modules):
- **Preserve all existing documentation.**  Every comment block, usage
  note, option explanation, and conceptual description in the original
  source must be carried over to the new version.  Do not summarize or
  abbreviate.
- **Use appropriate formats.**  In Python, use module-level docstrings,
  class/method docstrings, argparse help text, and inline comments.  In
  Fortran, use header comment blocks and inline comments.
- **Explain the "why", not just the "what".**  Students benefit from
  knowing the physics/chemistry motivation, not just the code mechanics.
- This policy applies to all future conversions and new code, not just
  the current Perl-to-Python effort.
