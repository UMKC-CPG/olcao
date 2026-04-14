#!/usr/bin/env python3

"""bond_analysis.py - Geometric bond analysis of atomic structures.

This is the Python port of the Perl ``bondAnalysis`` script.  It
takes an atomic structure data file in one of many possible formats
and produces a new file that contains information about the nearest
neighbor relationships between atoms.  This can typically be
considered as a bond analysis, but the data is produced simply based
on geometrical considerations, not on any quantum mechanical
calculations such as the OLCAO bond order.

Options
=======

-o OUTFILE
    Designate a different name for the resultant output file rather
    than the default "bondAnalysis.xx" file.  The "xx" in the output
    depends on the data requested: for -bs, "xx"="bs"; for -ba,
    "xx"="ba"; and so on.  A special case occurs for the -frag
    option where the name specified here serves as an identifier
    for a set of output files, each with a different suffix.

-i INFILE
    Designate a different name for the input file rather than the
    default expected input file of "olcao.skl".  Other options
    include files such as the "structure.dat" OLCAO input file and
    the VASP input file "POSCAR".  (At present only the olcao.skl
    option works.)

-s
    Allow the use of OLCAO bond order results to support the bond
    analysis.  The default file name expected to contain BO data is
    "S1__S1.dat".  This file is made using makeBOND with a group
    control file containing only the filter "SYSTEM_NUM".  Often,
    one may add an additional option to the makeBOND command of
    "-minbo 0.0".

-dist LIMITDIST
    A cutoff distance in angstroms beyond which atom-atom distances
    should not be considered.  If this option is *not* given, then
    the program will assume a default cutoff of 4 angstroms.

Mutually exclusive analysis modes
---------------------------------

The -bs, -ba, -bl, -boo, -dx, -st, -co, -qn, -vtk, and -rn
options are all mutually exclusive.  The one given last (relatively)
on the command line will be the one that is performed.

-bs
    Print a file that can be used to produce a ball-and-stick figure
    through some third party software.  The format is the lattice
    parameters given in a,b,c x,y,z format followed by a line with
    the number of atoms, followed by a line for each atom.  These
    lines have the element, the species number, the atomic
    coordinates in either x,y,z direct, a,b,c direct, or a,b,c
    fractional, and the list of the atom numbers to which each atom
    bonds.

-vtk
    Print a file in legacy VTK format that shows a ball-and-stick
    model.  This file can be easily read by Paraview.

-ba
    A list of the bond angles for each pair of atoms bonded to each
    atom.  The output is divided into groups, one group for each
    atom.  The first line of a group is the number of bond angles in
    that group.  The remaining lines are the various bond angles.
    Each bond angle line has the form of two triplets and the bond
    angle in degrees.  The first triplet is the element name and
    species number of each atom in the triplet with the vertex atom
    in the middle.  The second triplet contains the atom number in
    the same order.  Then, the bond angle is given in degrees.

-bl
    A list of the bond lengths for each atom that is "bonded" to
    each other atom.  The output is divided into groups, one group
    for each atom.  The first line of a group contains the number of
    bonds in that group and the element name and species number of
    the atom that is common to all the bonds in that group.  The
    next lines contain the element name and species number of each
    atom bonded to the common atom along with their bond lengths and
    atom number.

-boo
    Compute the so-called bond orientational order number of each
    atom in the model.  The output is simply a list of the element
    name, species number, index number, and BOO number of each atom.
    The -Ylm_l suboption sets the l order of the Ylm spherical
    harmonic used to distinguish the bond orientational order of
    different atoms.  Currently only 0, 1, and 6 are allowed values.

-dx
    Produce a file that can be used to view the structure in openDX.
    The bonds are based simply on the covalent radii of each atom.
    The atom colors and sizes are pulled from the ElementData module.

-st
    Compute a set of statistics about each element, species, and
    atom and print them.  The statistics include average bond
    lengths, the bond length standard deviation, the average bond
    angle, and the bond angle standard deviation.

-co
    Compute (and print) the coordination for each atom and also
    produce and print a summary of the coordinations for each
    element.

-qn
    Compute the Q^n distribution for a crystal.  By default, the
    assumption is that the "bridging" atom is oxygen.  However, the
    -b BRIDGES parameter can be used to define a different bridging
    element.  The element(s) must be specified as lower case names
    from the periodic table.  If there is more than one name, then
    enclose in quotes and separate by spaces.  Similarly, by
    default, silicon is assumed to be the coordinated cation but the
    -c IONS option can be used to change that.  Again, if more than
    one cation should be considered, then enclose a space-separated
    list in quotes.

-3dfull
    Produce a single OpenSCAD model file that can be used to produce
    data for a 3D printer.  This model will contain all atoms and
    bonds of the complete system unioned together so that the whole
    model prints as one piece.  Suboptions:
      -a ASCALE  scale the radius of the atoms (default: covalent
                 radius)
      -b BSCALE  scale the diameter of the bond sticks (default: 10%
                 of the covalent radius of the largest atom)
      -c CSCALE  scale the overall size of the crystal lattice
      -r BONDRADIUSSCALE  ratio between bond rod radius and covalent
                          radius
    NOTE: bad choices for these values will result in bad looking
    models.

-3dparts
    Produce a collection of OpenSCAD model files that contain
    descriptions of the individual parts associated with an atomic
    scale model.  When the parts models are printed on a 3D printer
    they can be used to assemble a large complex 3D structural
    model.  The three sub-options work the same way as for -3dfull.
    Additional sub-options:
      -vol X Y Z   specify the printable volume (mm)
      -f FSCALE    font size factor (default 0.9 * bond cylinder
                   radius)

-rn MINRINGLEN MAXRINGLEN
    Compute the ring distribution between MINRINGLEN and MAXRINGLEN.

-bf BONDINGFACTOR
    A multiplicative factor for the determination of which atoms
    are bonded.  The default value is 1.1.  Atoms are considered to
    be bonded if their interatomic distance is <= (the sum of their
    covalent radii) * bondingFactor.  Be careful that the largest
    expected bond length does not exceed the limit distance set with
    -dist.

-xyz
    Print the atomic positions in x,y,z direct space (cartesian)
    instead of the default a,b,c fractional coordinates.  This
    option will affect the interpretation of the coordinates
    specified in the -box option and it is mutually exclusive with
    -fract.

-abc
    Print the atomic positions in a,b,c direct space coordinates.

-fract
    Print the atomic positions in a,b,c fractional coordinates
    instead of the default x,y,z direct space (cartesian).  This
    option will affect the interpretation of the coordinates
    specified in the -box option and it is mutually exclusive with
    -xyz.

-box MIN1 MAX1 MIN2 MAX2 MIN3 MAX3 ZONE
    Use only a limited portion of the system to compute bond length
    and bond angle statistics.  The min and max values are used to
    define the inclusive range over each axis within which atoms are
    to be considered as "inside" the box.  You may optionally use
    the word "max" instead of a specific number for MAX1, MAX2,
    MAX3.  This will automatically use the maximum value so you
    don't have to look it up or worry about rounding issues.  The
    minimum value of 0 is easy to remember and works for all cells
    so there is no "min" option.  The ZONE can be either 1 or 2.
    If it is 1, then only atoms inside the box are considered.
    If it is 2, then only atoms outside the box are considered.
    The default behavior if this option is not given is to use abc
    coordinates and include all atoms inside the entire cell.

-hist HISTDELTA
    Instead of printing the computed statistics in a single
    structured data file with exact numerical values, the output
    will also include a histogram with spacing HISTDELTA.

USAGE
=====

::

    bond_analysis.py [-o OUTFILE] [-i INFILE] [-s]
        [-dist LIMITDIST]
        [-bs | -ba | -bl | -boo [-Ylm_l L] | -dx | -st
         | -co | -qn [-b BRIDGES] [-c IONS] | -vtk
         | -rn MINRINGLEN MAXRINGLEN]
        [-3dfull  [-a ASCALE -b BSCALE -c CSCALE
                   -r BONDRADIUSSCALE]]
        [-3dparts [-a ASCALE -b BSCALE -c CSCALE
                   -f FSCALE -r BONDRADIUSSCALE
                   -vol X Y Z]]
        [-bf BONDINGFACTOR] [-xyz | -abc | -fract]
        [-box MIN1 MAX1 MIN2 MAX2 MIN3 MAX3 ZONE]
        [-hist HISTDELTA]
        [-help]

Conversion Notes (Perl bondAnalysis -> Python bond_analysis.py)
===============================================================

This script was converted from the Perl ``bondAnalysis`` script.
The following documents known differences between the two versions.

Intentional behavioral changes / bug fixes
-------------------------------------------

1. **-boo -Ylm_l off-by-one (Perl bug fixed).**
   The Perl version has an off-by-one error when parsing the
   ``-Ylm_l`` sub-option of ``-boo``.  After detecting the
   ``-Ylm_l`` flag, it does ``$Ylm_l = $ARGV[++$number]``
   which reads the flag name itself (``"-Ylm_l"``) rather
   than the numeric value following it.  The Python version
   correctly skips the flag and reads the value.

2. **-3dfull sub-option loop count (Perl bug fixed).**
   The Perl iterates ``foreach $index (1..3)`` for sub-options
   of ``-3dfull``, but there are 4 possible sub-options
   (``-a``, ``-b``, ``-c``, ``-r``).  The Perl can only parse
   3 of 4 in a single invocation.  The Python loops 4 times,
   allowing all sub-options to be specified.

3. **-3dparts -f sub-option (Perl omission fixed).**
   The Perl help text documents ``-f $fScale`` as a sub-option
   of ``-3dparts``, but the implementation omits it.  The
   Python version implements what the help text promises.

4. **-qn ions sub-option flag: ``-i`` renamed to ``-c``.**
   In the Perl version, the ions sub-option for ``-qn`` uses
   ``-i``.  This conflicts with the top-level ``-i`` (input
   file) flag.  The Python version uses ``-c`` (for "cation")
   instead.  This is a **backward-incompatible change**.
   Users who had ``-qn -i "si al"`` must now use
   ``-qn -c "si al"``.

5. **getBondsToShow off-by-one (Perl bug fixed).**
   The Perl deduplication loop iterates
   ``0..$numBondsExt-1``, which skips the last bond in the
   sorted list.  The Python processes all bonds.  This may
   result in one additional bond appearing in DX and OpenSCAD
   output compared to the Perl version.

6. **Bond sheet position in -3dparts (Perl bug fixed).**
   When laying out bonds on sheets for 3D printing, the Perl
   always uses the first bond's length for the y-position
   initialization on every sheet.  The Python correctly uses
   the starting bond's length for each new sheet.

7. **Debug print removed.**
   The Perl ``printBondLengths`` contains a diagnostic
   ``print STDOUT "min max numPts: ..."`` that appears to be
   a leftover debug statement.  It is not ported.

Minor output format differences
-------------------------------

8. **Histogram bin counts.**
   The Perl histogram initializes bins to ``0.0`` but prints
   them via Perl string interpolation, which renders
   integer-valued floats without a decimal point (e.g. ``5``).
   The Python uses integer bins and prints them as integers
   to match this behavior.

9. **VTK float precision.**
   Default float-to-string conversion differs between Perl
   and Python for coordinate, radii, and color values in VTK
   output.  The values are numerically identical but may have
   different decimal representations.

10. **OpenSCAD module declarations.**
    Some long OpenSCAD module declarations are wrapped to
    multiple lines in the Python output where the Perl output
    has them on a single line.  OpenSCAD is
    whitespace-insensitive, so this has no functional effect.

Structural differences (no behavioral effect)
----------------------------------------------

11. The Perl fetches ~15 ``get*Ref`` references into
    module-scoped variables after ``readInputFile``.  The
    Python accesses ``StructureControl`` attributes directly
    at point of use.

12. Several small Perl subroutines (``printFullModel``,
    ``scaleAtomicCoords``, ``printKey``, ``printSubCube``)
    are inlined into their callers in the Python version.

13. The Python accepts ``--help`` and ``-h`` in addition to
    the Perl's ``-help``.

14. The Python adds harmless defaults for ``min_ring_len``
    (3) and ``max_ring_len`` (8); the Perl leaves these
    uninitialized until the ``-rn`` flag sets them.

15. The command file records ``bond_analysis.py`` instead of
    ``bondAnalysis`` as the program name.
"""

import argparse as ap
import math
import os
import sys
from datetime import datetime

from structure_control import StructureControl
from element_data import ElementData


# -------------------------------------------------------------------
# Operation codes -- map analysis mode names to integers that mirror
# the original Perl script's numbering.
# -------------------------------------------------------------------
OP_BALL_AND_STICK = 1
OP_BOND_ANGLES = 2
OP_BOND_LENGTHS = 3
OP_BOND_OO = 4
OP_OPEN_DX = 5
OP_STATISTICS = 6
OP_COORDINATION = 7
OP_3D_FULL = 8
OP_3D_PARTS = 9
# OP_FRAG = 10  # disabled in Perl original
OP_VTK = 11
OP_QN = 12
OP_RING = 13

# Map operation codes to default output file suffixes.
_OP_SUFFIX = {
    OP_BALL_AND_STICK: ".bs",
    OP_BOND_ANGLES: ".ba",
    OP_BOND_LENGTHS: ".bl",
    OP_BOND_OO: ".boo",
    OP_OPEN_DX: ".bx",
    OP_STATISTICS: ".st",
    OP_COORDINATION: ".co",
    OP_3D_FULL: ".scad",
    OP_VTK: ".vtk",
    OP_QN: ".qn",
    OP_RING: ".rn",
}

# Coordinate type codes.
COORD_XYZ = 1    # Direct space x, y, z (cartesian).
COORD_ABC = 2    # Direct space a, b, c.
COORD_FRACT = 3  # Fractional a, b, c.


# -------------------------------------------------------------------
# ScriptSettings -- command-line parameters
# -------------------------------------------------------------------

class ScriptSettings:
    """Holds all user-controllable parameters for bond_analysis.

    On construction the object:
      1. Reads defaults from bond_analysisrc.py.
      2. Parses the command line (argparse).
      3. Reconciles CLI values with rc defaults.
      4. Records the command line to the ``command`` file.
    """

    def __init__(self):
        # Read default variables from the resource control file.
        rc = self._load_rc()
        self.assign_rc_defaults(rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc
        # defaults.
        self.reconcile(args)

        # At this point, the command line parameters are set
        # and accepted.  We use this as a good spot to record
        # the command line parameters that were used.
        self.record_clp()

    # ---------------------------------------------------------------
    # RC loading
    # ---------------------------------------------------------------

    @staticmethod
    def _load_rc():
        """Import bond_analysisrc and return its parameter dict.

        Search order: current working directory first (local
        override), then $OLCAO_RC.
        """
        cwd_rc = os.path.join(
            os.getcwd(), "bond_analysisrc.py"
        )
        if os.path.isfile(cwd_rc):
            sys.path.insert(0, os.getcwd())
        else:
            rc_dir = os.getenv("OLCAO_RC")
            if not rc_dir:
                sys.exit(
                    "Error: $OLCAO_RC is not set and no "
                    "local bond_analysisrc.py found. "
                    "See instructions."
                )
            sys.path.insert(0, rc_dir)
        from bond_analysisrc import parameters_and_defaults
        return parameters_and_defaults()

    # ---------------------------------------------------------------
    # Assign rc defaults
    # ---------------------------------------------------------------

    def assign_rc_defaults(self, rc):
        """Populate instance attributes from the rc dictionary."""

        # Input / output files.
        self.in_file = rc["in_file"]
        self.out_file = rc["out_file"]
        self.bond_order_file = rc["bond_order_file"]

        # Analysis operation.
        self.operation = rc["operation"]

        # Bonding criteria.
        self.limit_dist = rc["limit_dist"]
        self.bonding_factor = rc["bonding_factor"]
        self.use_olcao_bo = rc["use_olcao_bo"]

        # Coordinate type.
        self.coord_type = rc["coord_type"]

        # Bond orientational order.
        self.ylm_l = rc["ylm_l"]

        # 3D printing scale factors.
        self.a_scale = rc["a_scale"]
        self.c_scale = rc["c_scale"]
        self.b_scale = rc["b_scale"]
        self.f_scale = rc["f_scale"]
        self.bond_radius_scale = rc["bond_radius_scale"]

        # 3D printing volume.
        self.x_vol_3d = rc["x_vol_3d"]
        self.y_vol_3d = rc["y_vol_3d"]
        self.z_vol_3d = rc["z_vol_3d"]

        # Bounding box.
        # box_borders[axis][side]: axis in {1,2,3},
        # side 1=min, 2=max.
        self.box_borders = [
            None,
            [None, rc["box_min1"], rc["box_max1"]],
            [None, rc["box_min2"], rc["box_max2"]],
            [None, rc["box_min3"], rc["box_max3"]],
        ]
        self.zone = rc["zone"]

        # Histogram.
        self.do_hist = rc["do_hist"]
        self.hist_delta = rc["hist_delta"]

        # Q^n distribution.
        self.bridges = list(rc["bridges"])
        self.ions = list(rc["ions"])

        # Ring distribution.
        self.min_ring_len = rc["min_ring_len"]
        self.max_ring_len = rc["max_ring_len"]

    # ---------------------------------------------------------------
    # Command-line parsing
    # ---------------------------------------------------------------

    def parse_command_line(self):
        """Build the argparse parser and parse the command line.

        Top-level flags (-o, -i, -s, -dist, -bf, -xyz, -abc,
        -fract, -box, -hist) are handled by argparse directly.

        Analysis mode flags (-bs, -ba, -bl, -boo, -dx, -st,
        -co, -qn, -vtk, -rn, -3dfull, -3dparts) and their
        context-dependent sub-options are collected from the
        unrecognized arguments and parsed separately, because
        short flags like -a, -b, -c have different meanings
        under different modes.
        """

        prog_name = "bond_analysis"

        description_text = """\
Geometric bond analysis of atomic structures.

Takes an atomic structure file and produces information
about the nearest-neighbor relationships between atoms
based on geometric (covalent radii) criteria.
"""

        epilog_text = """\
ANALYSIS MODES (mutually exclusive)
====================================

-bs          Ball-and-stick output.
-ba          Bond angle list.
-bl          Bond length list (default).
-boo         Bond orientational order.
               Sub-option: -Ylm_l L  (default 6).
-dx          OpenDX visualization.
-st          Bond statistics (lengths + angles).
-co          Coordination data.
-qn          Q^n distribution.
               Sub-options: -b BRIDGES  -c IONS
-vtk         VTK ball-and-stick (for Paraview).
-rn MIN MAX  Ring distribution.
-3dfull      OpenSCAD 3D full model.
               Sub-options: -a ASCALE  -b BSCALE
                            -c CSCALE  -r RSCALE
-3dparts     OpenSCAD 3D parts model.
               Sub-options: -a ASCALE  -b BSCALE
                            -c CSCALE  -r RSCALE
                            -f FSCALE
                            -vol X Y Z

Note: analysis mode flags and their sub-options
should appear after all top-level flags on the
command line.

Defaults are given in ./bond_analysisrc.py or
$OLCAO_RC/bond_analysisrc.py.
"""

        parser = ap.ArgumentParser(
            prog=prog_name,
            formatter_class=ap.RawDescriptionHelpFormatter,
            description=description_text,
            epilog=epilog_text,
        )
        self._add_parser_arguments(parser)
        args, remaining = parser.parse_known_args()

        # Parse mode flags and their sub-options from the
        # unrecognized arguments.
        args.mode_op = None
        args.mode_extra = {}
        self._parse_mode(args, remaining)

        return args

    def _add_parser_arguments(self, parser):
        """Register top-level CLI flags with argparse.

        Analysis mode flags (-bs, -ba, -bl, -boo, -dx, -st,
        -co, -qn, -vtk, -rn, -3dfull, -3dparts) and their
        sub-options are parsed separately from the remaining
        argv because they use context-dependent short flags
        (-a, -b, -c) that would collide between modes.
        """

        # ---- Input / output ----
        # The -o option designates a different name for the
        # resultant output file rather than the default
        # "bondAnalysis.xx" file.
        parser.add_argument(
            "-o", dest="out_file", type=str,
            default=None,
            help="Output file base name. Default: "
                 f"'{self.out_file}' + mode suffix."
        )
        # The -i option designates a different name for the
        # input file rather than the default "olcao.skl".
        parser.add_argument(
            "-i", dest="in_file", type=str,
            default=None,
            help="Input structure file. "
                 f"Default: '{self.in_file}'."
        )

        # ---- OLCAO bond order support ----
        # The -s flag enables the use of OLCAO bond order
        # results to support the bond analysis.
        parser.add_argument(
            "-s", dest="use_olcao_bo",
            action="store_true", default=False,
            help="Use OLCAO bond order data from "
                 f"'{self.bond_order_file}'."
        )

        # ---- Cutoff distance ----
        # A cutoff distance in angstroms beyond which
        # atom-atom distances should not be considered.
        parser.add_argument(
            "-dist", dest="limit_dist", type=float,
            default=None,
            help="Max atom-atom distance (angstroms). "
                 f"Default: {self.limit_dist}."
        )

        # ---- Bonding factor ----
        # A multiplicative factor for the determination of
        # which atoms are bonded.  Atoms are considered bonded
        # if distance <= (sum of covalent radii) * factor.
        parser.add_argument(
            "-bf", dest="bonding_factor", type=float,
            default=None,
            help="Multiplicative bonding factor. "
                 f"Default: {self.bonding_factor}."
        )

        # ---- Coordinate type (mutually exclusive) ----
        coord = parser.add_mutually_exclusive_group()
        # Print the atomic positions in x,y,z direct space
        # (cartesian) instead of the default.
        coord.add_argument(
            "-xyz", dest="coord_type",
            action="store_const", const=COORD_XYZ,
            help="Use x,y,z cartesian coordinates."
        )
        # Print the atomic positions in a,b,c direct space
        # coordinates.
        coord.add_argument(
            "-abc", dest="coord_type",
            action="store_const", const=COORD_ABC,
            help="Use a,b,c direct space coordinates."
        )
        # Print the atomic positions in a,b,c fractional
        # coordinates (default).
        coord.add_argument(
            "-fract", dest="coord_type",
            action="store_const", const=COORD_FRACT,
            help="Use a,b,c fractional coordinates "
                 "(default)."
        )

        # ---- Bounding box ----
        # Use only a limited portion of the system.  The
        # min/max values define the inclusive range over each
        # axis.  ZONE is 1 (inside) or 2 (outside).  Use
        # 'max' for automatic maximum.
        parser.add_argument(
            "-box", dest="box", nargs=7,
            metavar=(
                "MIN1", "MAX1", "MIN2", "MAX2",
                "MIN3", "MAX3", "ZONE"
            ),
            default=None,
            help="Bounding box: min/max for each axis "
                 "+ zone (1=inside, 2=outside). "
                 "Use 'max' for automatic maximum."
        )

        # ---- Histogram ----
        # Produce a histogram alongside the output with the
        # given bin width.
        parser.add_argument(
            "-hist", dest="hist_delta", type=float,
            default=None,
            help="Produce a histogram with the given "
                 "bin width."
        )

    @staticmethod
    def _parse_mode(args, remaining):
        """Parse analysis mode flags from unrecognized args.

        This handles mode flags and their context-dependent
        sub-options that cannot be registered as top-level
        argparse arguments because the short flags (-a, -b,
        -c, etc.) have different meanings under different
        modes.

        The parsing mirrors the original Perl sequential
        parser: each mode flag greedily consumes its
        recognized sub-options and stops at the first
        unrecognized flag.
        """
        n = 0
        while n < len(remaining):
            arg = remaining[n]

            if arg == "-bs":
                args.mode_op = OP_BALL_AND_STICK

            elif arg == "-vtk":
                args.mode_op = OP_VTK

            elif arg == "-ba":
                args.mode_op = OP_BOND_ANGLES

            elif arg == "-bl":
                args.mode_op = OP_BOND_LENGTHS

            elif arg == "-boo":
                args.mode_op = OP_BOND_OO
                # Check for optional -Ylm_l sub-option.
                if (n + 1 < len(remaining)
                        and remaining[n + 1] == "-Ylm_l"):
                    n += 2
                    args.mode_extra["ylm_l"] = int(
                        remaining[n]
                    )

            elif arg == "-dx":
                args.mode_op = OP_OPEN_DX

            elif arg == "-st":
                args.mode_op = OP_STATISTICS

            elif arg == "-co":
                args.mode_op = OP_COORDINATION

            elif arg == "-qn":
                args.mode_op = OP_QN
                # Parse up to two sub-options: -b and -c.
                for _ in range(2):
                    if n + 1 >= len(remaining):
                        break
                    nxt = remaining[n + 1]
                    if nxt == "-b":
                        n += 2
                        args.mode_extra["bridges"] = (
                            remaining[n].lower().split()
                        )
                    elif nxt == "-c":
                        n += 2
                        args.mode_extra["ions"] = (
                            remaining[n].lower().split()
                        )
                    else:
                        break

            elif arg == "-rn":
                args.mode_op = OP_RING
                n += 1
                min_r = int(remaining[n])
                n += 1
                max_r = int(remaining[n])
                if min_r < 0:
                    print(
                        f"{min_r} must be a "
                        "positive integer."
                    )
                    sys.exit(1)
                if max_r < 0:
                    print(
                        f"{max_r} must be a "
                        "positive integer."
                    )
                    sys.exit(1)
                args.mode_extra["min_ring_len"] = min_r
                args.mode_extra["max_ring_len"] = max_r

            elif arg == "-3dfull":
                args.mode_op = OP_3D_FULL
                # Parse up to 4 sub-options.
                for _ in range(4):
                    if n + 1 >= len(remaining):
                        break
                    nxt = remaining[n + 1]
                    if nxt == "-a":
                        n += 2
                        args.mode_extra["a_scale"] = (
                            float(remaining[n])
                        )
                    elif nxt == "-b":
                        n += 2
                        args.mode_extra["b_scale"] = (
                            float(remaining[n])
                        )
                    elif nxt == "-c":
                        n += 2
                        args.mode_extra["c_scale"] = (
                            float(remaining[n])
                        )
                    elif nxt == "-r":
                        n += 2
                        args.mode_extra[
                            "bond_radius_scale"
                        ] = float(remaining[n])
                    else:
                        break

            elif arg == "-3dparts":
                args.mode_op = OP_3D_PARTS
                # Parse up to 6 sub-options.
                for _ in range(6):
                    if n + 1 >= len(remaining):
                        break
                    nxt = remaining[n + 1]
                    if nxt == "-a":
                        n += 2
                        args.mode_extra["a_scale"] = (
                            float(remaining[n])
                        )
                    elif nxt == "-b":
                        n += 2
                        args.mode_extra["b_scale"] = (
                            float(remaining[n])
                        )
                    elif nxt == "-c":
                        n += 2
                        args.mode_extra["c_scale"] = (
                            float(remaining[n])
                        )
                    elif nxt == "-r":
                        n += 2
                        args.mode_extra[
                            "bond_radius_scale"
                        ] = float(remaining[n])
                    elif nxt == "-vol":
                        n += 2
                        args.mode_extra["x_vol_3d"] = (
                            float(remaining[n])
                        )
                        n += 1
                        args.mode_extra["y_vol_3d"] = (
                            float(remaining[n])
                        )
                        n += 1
                        args.mode_extra["z_vol_3d"] = (
                            float(remaining[n])
                        )
                    elif nxt == "-f":
                        n += 2
                        args.mode_extra["f_scale"] = (
                            float(remaining[n])
                        )
                    else:
                        break

            else:
                print(
                    f"UNKNOWN COMMAND LINE PARAMETER "
                    f"{arg}. ABORTING."
                )
                sys.exit(1)

            n += 1

    # ---------------------------------------------------------------
    # Reconcile CLI with rc defaults
    # ---------------------------------------------------------------

    def reconcile(self, args):
        """Merge parsed CLI arguments into the settings.

        Top-level argparse arguments override rc defaults only
        when explicitly provided (i.e. not None).  Mode flags
        override the operation code and any mode-specific
        parameters.
        """

        # Top-level flags.
        if args.out_file is not None:
            self.out_file = args.out_file
        if args.in_file is not None:
            self.in_file = args.in_file
        if args.use_olcao_bo:
            self.use_olcao_bo = True
        if args.limit_dist is not None:
            self.limit_dist = args.limit_dist
        if args.bonding_factor is not None:
            self.bonding_factor = args.bonding_factor
        if args.coord_type is not None:
            self.coord_type = args.coord_type
        if args.hist_delta is not None:
            self.do_hist = True
            self.hist_delta = args.hist_delta

        # Bounding box.
        if args.box is not None:
            self.box_borders[1][1] = args.box[0]
            self.box_borders[1][2] = args.box[1]
            self.box_borders[2][1] = args.box[2]
            self.box_borders[2][2] = args.box[3]
            self.box_borders[3][1] = args.box[4]
            self.box_borders[3][2] = args.box[5]
            self.zone = int(args.box[6])

        # Analysis mode.
        if args.mode_op is not None:
            self.operation = args.mode_op

        # Mode-specific extras.
        extra = args.mode_extra
        if "ylm_l" in extra:
            self.ylm_l = extra["ylm_l"]
        if "bridges" in extra:
            self.bridges = extra["bridges"]
        if "ions" in extra:
            self.ions = extra["ions"]
        if "min_ring_len" in extra:
            self.min_ring_len = extra["min_ring_len"]
        if "max_ring_len" in extra:
            self.max_ring_len = extra["max_ring_len"]
        if "a_scale" in extra:
            self.a_scale = extra["a_scale"]
        if "b_scale" in extra:
            self.b_scale = extra["b_scale"]
        if "c_scale" in extra:
            self.c_scale = extra["c_scale"]
        if "bond_radius_scale" in extra:
            self.bond_radius_scale = (
                extra["bond_radius_scale"]
            )
        if "f_scale" in extra:
            self.f_scale = extra["f_scale"]
        if "x_vol_3d" in extra:
            self.x_vol_3d = extra["x_vol_3d"]
        if "y_vol_3d" in extra:
            self.y_vol_3d = extra["y_vol_3d"]
        if "z_vol_3d" in extra:
            self.z_vol_3d = extra["z_vol_3d"]

    # ---------------------------------------------------------------
    # Record command line
    # ---------------------------------------------------------------

    def record_clp(self):
        """Append the command line with timestamp to 'command'.

        This follows the XYZ.py convention of recording the
        date and the full command that was used so that users
        can track what parameters produced each output.
        """
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime(
                "%b. %d, %Y: %H:%M:%S"
            )
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write("Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")

    # ---------------------------------------------------------------
    # Output file name resolution
    # ---------------------------------------------------------------

    def resolve_out_file(self):
        """Determine the output file name.

        If the user did not specify -o, append the appropriate
        suffix to the default base name.
        """
        if self.out_file == "bondAnalysis":
            suffix = _OP_SUFFIX.get(self.operation, "")
            self.out_file = self.out_file + suffix


# -------------------------------------------------------------------
# BondAnalysis -- main analysis engine
# -------------------------------------------------------------------

class BondAnalysis:
    """Performs geometric bond analysis on an atomic structure.

    This class orchestrates reading the input structure, computing
    bonding information, and writing the requested output.  It
    mirrors the sequential flow of the original Perl script.
    """

    def __init__(self, settings):
        self.s = settings

        # The StructureControl object holds the atomic structure
        # and all bonding/angle data.
        self.sc = StructureControl()

        # Computed statistics arrays (populated by
        # compute_statistics when operation == OP_STATISTICS).
        self.avg_elem_bl = []
        self.avg_elem_ba = []
        self.avg_spec_bl = []
        self.avg_spec_ba = []
        self.avg_atom_bl = []
        self.avg_atom_ba = []
        self.std_elem_bl = []
        self.std_elem_ba = []
        self.std_spec_bl = []
        self.std_spec_ba = []
        self.std_atom_bl = []
        self.std_atom_ba = []

    # ---------------------------------------------------------------
    # Main execution flow
    # ---------------------------------------------------------------

    def run(self):
        """Execute the full bond analysis workflow."""
        s = self.s
        sc = self.sc

        # Apply the bonding factor to the element database.
        sc.element_data.apply_bond_factor(s.bonding_factor)

        # Read the input file.
        sc.read_input_file(s.in_file, 1)

        # Set the limit for atom interaction (maximum bond limit).
        sc.set_limit_dist(s.limit_dist)

        # Establish the bounding box that constrains which atoms
        # to consider.
        self._set_box_border()

        # Create the list of which atoms are bonded to which atoms.
        sc.create_bonding_list()

        # Compute the bond angles.  For each atom in the central
        # cell, obtain a list of bond angles between all pairs of
        # bonded atoms.
        sc.compute_bond_angles_ext()

        # Get coordinate data in the requested format.
        self._setup_coords()

        # Perform operation-specific setup.
        if s.operation == OP_BOND_OO:
            # Set the order of the spherical harmonics to use.
            sc.set_ylm_l(s.ylm_l)
            # Create the bond orientational order list.
            sc.create_q_list()

        elif s.operation == OP_STATISTICS:
            self._compute_statistics()

        elif s.operation == OP_COORDINATION:
            # Create the coordination list and the coordination
            # summary.
            sc.create_coordination_list()
            sc.create_coordination_summary()

        elif s.operation == OP_QN:
            # Create the Q^n lists (atomic, absolute, and
            # fractional).  Also compute the Q^n network
            # connectivity.
            sc.compute_qn(s.bridges, s.ions)

        elif s.operation == OP_RING:
            # Create the ring membership lists.
            sc.compute_ring_distribution(
                s.min_ring_len, s.max_ring_len
            )

        # Resolve the output file name.
        s.resolve_out_file()

        # Dispatch to the appropriate output function.
        self._dispatch_output()

    # ---------------------------------------------------------------
    # Coordinate setup
    # ---------------------------------------------------------------

    def _setup_coords(self):
        """Obtain coordinates in the user-requested format.

        Also sets self.coord_label for output headers.
        """
        sc = self.sc
        s = self.s

        if s.coord_type == COORD_XYZ:
            self.coords = sc.direct_xyz
            self.coord_label = "(x,y,z) direct"
        elif s.coord_type == COORD_ABC:
            self.coords = sc.direct_abc
            self.coord_label = "(a,b,c) direct"
        else:  # COORD_FRACT
            self.coords = sc.fract_abc
            self.coord_label = "(a,b,c) fractional"

    # ---------------------------------------------------------------
    # Output dispatch
    # ---------------------------------------------------------------

    def _dispatch_output(self):
        """Call the appropriate print function for the operation."""
        s = self.s
        op = s.operation

        if op == OP_BALL_AND_STICK:
            self._print_ball_and_stick()
        elif op == OP_BOND_ANGLES:
            self._print_bond_angles()
        elif op == OP_BOND_LENGTHS:
            self._print_bond_lengths()
        elif op == OP_BOND_OO:
            self._print_bond_oo()
        elif op == OP_OPEN_DX:
            self._print_bond_dx()
        elif op == OP_STATISTICS:
            self._print_bond_stats()
        elif op == OP_COORDINATION:
            self._print_coordination_data()
        elif op == OP_3D_FULL:
            self._print_openscad(do_parts=False)
        elif op == OP_3D_PARTS:
            self._print_openscad(do_parts=True)
        elif op == OP_VTK:
            self._print_vtk_ball_and_stick()
        elif op == OP_QN:
            self._print_qn_distribution()
        elif op == OP_RING:
            self._print_ring_distribution()

    # ---------------------------------------------------------------
    # Box border setup
    # ---------------------------------------------------------------

    def _set_box_border(self):
        """Resolve 'max' tokens in box borders and call setBorder.

        For any end points given as "max" values, we fill them here
        depending on which type of coordinates we are using.
        """
        sc = self.sc
        s = self.s

        for axis in range(1, 4):
            if s.box_borders[axis][2] == "max":
                if s.coord_type in (COORD_XYZ, COORD_ABC):
                    # Direct space: max is the lattice vector
                    # magnitude.
                    s.box_borders[axis][2] = sc.mag[axis]
                else:
                    # Fractional: max is 1.0.
                    s.box_borders[axis][2] = 1.0
            else:
                s.box_borders[axis][2] = float(
                    s.box_borders[axis][2]
                )
            s.box_borders[axis][1] = float(
                s.box_borders[axis][1]
            )

        # Set the border so that we include only the requested
        # atoms.  Both low- and high-bound triples are passed as
        # 1-indexed [None, v1, v2, v3] lists to match the set_border
        # parameter convention.
        sc.set_border(
            s.zone,
            [None,
             s.box_borders[1][1], s.box_borders[2][1],
             s.box_borders[3][1]],
            [None,
             s.box_borders[1][2], s.box_borders[2][2],
             s.box_borders[3][2]],
            s.coord_type,
        )

    # ---------------------------------------------------------------
    # Ball and stick output
    # ---------------------------------------------------------------

    def _print_ball_and_stick(self):
        """Print ball-and-stick data to the output file."""
        sc = self.sc
        s = self.s
        coords = self.coords

        with open(s.out_file, "w") as out:
            # Print the cell parameters.
            out.write("Cell parameters:\n")
            rl = sc.real_lattice
            out.write(
                f"{rl[1][1]} {rl[1][2]} {rl[1][3]}"
                " ax ay az\n"
            )
            out.write(
                f"{rl[2][1]} {rl[2][2]} {rl[2][3]}"
                " bx by bz\n"
            )
            out.write(
                f"{rl[3][1]} {rl[3][2]} {rl[3][3]}"
                " cx cy cz\n"
            )

            # Print the number of atoms in the model.
            out.write(
                f"Number of atoms:  {sc.num_atoms}\n"
            )

            # Print header.
            out.write(
                f" Atom # "
                f" Elem. "
                f"            {self.coord_label}"
                f"             Bonded atom #s\n"
            )

            for atom in range(1, sc.num_atoms + 1):
                # Print the atom name and coordinate position.
                out.write(f"{atom:8d}")
                name = sc.atom_element_name[atom]
                spec = sc.atom_species_id[atom]
                c = coords[atom]
                out.write(
                    f"{name:>3s}{spec:<5d}"
                    f"{c[1]:12.5e} "
                    f"{c[2]:12.5e} "
                    f"{c[3]:12.5e}"
                )

                # Print the atom numbers that this atom is bonded
                # to.
                for bond in range(
                    1, sc.num_bonds[atom] + 1
                ):
                    out.write(
                        f" {sc.bonded[atom][bond]:5d}"
                    )

                out.write("\n")

    # ---------------------------------------------------------------
    # VTK ball and stick output
    # ---------------------------------------------------------------

    def _print_vtk_ball_and_stick(self):
        """Print legacy VTK format ball-and-stick model."""
        sc = self.sc
        s = self.s
        coords = self.coords
        bohr_rad = sc.get_bohr_rad()

        with open(s.out_file, "w") as out:
            # Print the header.
            out.write(
                "# vtk DataFile Version 4.2\n"
                "This is a model produced by the OLCAO "
                "program suite: bondAnalysis program\n"
                "ASCII\n"
                "DATASET POLYDATA\n"
                "\n"
            )

            # Print the positions of all atoms in the cell.
            out.write(f"POINTS {sc.num_atoms} float\n")
            for atom in range(1, sc.num_atoms + 1):
                ax = coords[atom][1] / bohr_rad
                ay = coords[atom][2] / bohr_rad
                az = coords[atom][3] / bohr_rad
                out.write(f"{ax} {ay} {az}\n")
            out.write("\n")

            # Now, list the connections between the data points
            # (POINTS) as LINES in the vtk file.  The header of
            # this section must include the number of lines,
            # followed by the total number of data points in the
            # section used to express those lines.  So, if we are
            # only using the number of points per line, and the
            # start and end points, then we have 3x the number of
            # bonds for the total number.
            vtk_bond_list = []
            for atom in range(1, sc.num_atoms + 1):
                for bond in range(
                    1, sc.num_bonds[atom] + 1
                ):
                    # Avoid double listing bonds by only printing
                    # those where the index number for atom one is
                    # less than the index number for atom 2.
                    bonded_atom = sc.bonded[atom][bond]
                    if atom < bonded_atom:
                        # -1 because VTK arrays start at 0.
                        a1 = atom - 1
                        a2 = bonded_atom - 1
                        vtk_bond_list.append(
                            f"2 {a1} {a2}\n"
                        )

            total_num_bonds = len(vtk_bond_list)
            total_bond_points = total_num_bonds * 3
            out.write(
                f"LINES {total_num_bonds} "
                f"{total_bond_points}\n"
            )
            for line in vtk_bond_list:
                out.write(line)
            out.write("\n")

            # Now here is where we need to include the additional
            # information to display the model in legacy vtk
            # format.  This next section marks the atom positions
            # from the POINTS section as vertices in vtk.  Marking
            # them as vertices allows programs like Paraview or
            # just VTK to apply glyphs to them, and for the
            # purposes of this program, we will need to apply
            # spherical glyphs using a filter made for Paraview
            # once this program has run.  The VERTICES header must
            # include the number of vertices and the number of data
            # points used to express them, which is only the number
            # of vertices followed by the index number of each atom
            # on every line, so 2x the number of vertices.
            total_vert_points = 2 * sc.num_atoms
            out.write(
                f"VERTICES {sc.num_atoms} "
                f"{total_vert_points}\n"
            )
            for atom in range(1, sc.num_atoms + 1):
                # For our purposes, we treat "cells" and "points"
                # as the same thing.  For now, we do not need to
                # specify "CELL_DATA" in this calculation.
                # We need to attach a vertex to each point/atom.
                # Since C arrays start at zero, we must start from
                # there.  Attempts to use in Paraview can sometimes
                # be "hit or miss" when applying glyphs to
                # vertices.  If spherical glyphs are not applying
                # to one or more of your vertices, change the "1"
                # to a higher number as a quick fix.
                atom_c_index = atom - 1
                out.write(f"1 {atom_c_index}\n")
            out.write("\n")

            # This final section adds the specified color table to
            # the bonds and atoms.  The color tables for VTK may
            # be left as default, or you may construct a color
            # table and specify how it is used.  Each entry in the
            # lookup table is an rgba array (red green blue alpha)
            # where alpha is transparency.  To save time, we use
            # the default color table.  But in the future a color
            # table should be generated in each file, specifying
            # what color is associated with what type of atom or
            # what strength of bond.
            out.write(f"POINT_DATA {sc.num_atoms}\n")

            # This section will attempt to apply a scaling factor
            # to spherical glyphs applied to the vertices, based
            # on the atomic radii.
            out.write("SCALARS atomic_radius float\n")
            out.write("LOOKUP_TABLE default\n")
            for atom in range(1, sc.num_atoms + 1):
                out.write(
                    f"{sc.coval_radii[atom]}\n"
                )
            out.write("\n")

            # This section stores the color data as point data.
            out.write("COLOR_SCALARS atom_colors 3\n")
            for atom in range(1, sc.num_atoms + 1):
                r = sc.color_vtk[atom][0] / 255.0
                g = sc.color_vtk[atom][1] / 255.0
                b = sc.color_vtk[atom][2] / 255.0
                out.write(f"{r} {g} {b}\n")

    # ---------------------------------------------------------------
    # Bond angles output
    # ---------------------------------------------------------------

    def _print_bond_angles(self):
        """Print bond angle data and optional histogram."""
        sc = self.sc
        s = self.s
        coords = self.coords

        max_ba = 0.0
        min_ba = 360.0

        with open(s.out_file, "w") as out:
            # For each atom print the bond angles between each
            # possible pair of atoms, and also print the atom ID
            # numbers for each pair.  Also, gather statistics that
            # will be useful for printing a histogram if asked.
            for atom in range(1, sc.num_atoms + 1):
                # Compare the position of this atom with the box
                # the user defined (or the default).  This also
                # considers if the user wanted only atoms inside
                # the box or only outside.  The return value is 1
                # if the atom is out of the desired region and 0
                # if the atom is inside the desired region.
                if sc.item_out_of_bounds(
                    atom, coords
                ):
                    continue

                # Print the number of bond angles for this atom.
                out.write(
                    f"{atom} Num bond angles:  "
                    f"{sc.num_bond_angles[atom]}\n"
                )

                # Initialize a count of the current bond angle.
                current_ba_count = 0

                nb = sc.num_bonds[atom]
                for bond1 in range(1, nb):
                    # Get the central cell atom number of this
                    # bonded atom.
                    bonded_atom1 = (
                        sc.ext_to_central_item_map[
                            sc.bonded_ext[atom][bond1]
                        ]
                    )

                    for bond2 in range(bond1 + 1, nb + 1):
                        # Increment the count of the number of
                        # bonds.
                        current_ba_count += 1

                        # Get the central cell atom number of
                        # this bonded atom.
                        bonded_atom2 = (
                            sc.ext_to_central_item_map[
                                sc.bonded_ext[atom][bond2]
                            ]
                        )

                        # Print the element name and species
                        # number of each atom in this bond set
                        # with the middle listed atom being the
                        # one at the vertex.
                        n1 = sc.atom_element_name[
                            bonded_atom1
                        ]
                        s1 = sc.atom_species_id[
                            bonded_atom1
                        ]
                        na = sc.atom_element_name[atom]
                        sa = sc.atom_species_id[atom]
                        n2 = sc.atom_element_name[
                            bonded_atom2
                        ]
                        s2 = sc.atom_species_id[
                            bonded_atom2
                        ]
                        out.write(
                            f"{n1}{s1} {na}{sa} "
                            f"{n2}{s2} "
                        )

                        # Print the index number of the atoms in
                        # this bond with the middle listed atom
                        # being the one at the vertex.  Note that
                        # the number is in reference to the atoms
                        # in the central cell.
                        out.write(
                            f"{bonded_atom1:7d}"
                            f"{atom:7d}"
                            f"{bonded_atom2:7d} "
                        )

                        # Finally, print the bond angle.
                        ba = sc.bond_angles_ext[atom][
                            current_ba_count
                        ]
                        out.write(f"{ba:16.8f}\n")

                        # Compare this bond angle to the min and
                        # max observed so far.
                        if ba < min_ba:
                            min_ba = ba
                        if ba > max_ba:
                            max_ba = ba

        # Once the traditional bond angle list has been created,
        # now we compute the histogram if asked.
        if s.do_hist:
            self._write_angle_histogram(min_ba, max_ba)

    def _write_angle_histogram(self, min_ba, max_ba):
        """Write a histogram of bond angles."""
        sc = self.sc
        s = self.s
        coords = self.coords

        # Compute the number of points in the histogram.
        num_points = (
            int((max_ba - min_ba) / s.hist_delta) + 1
        )

        # Initialize the histogram.
        histogram = [0] * (num_points + 1)

        # Revisit all the atoms to produce the histogram of every
        # bond angle.
        for atom in range(1, sc.num_atoms + 1):
            if sc.item_out_of_bounds(atom, coords):
                continue

            current_ba_count = 0
            nb = sc.num_bonds[atom]

            # Consider one bond.
            for bond1 in range(1, nb):
                # Consider each other bond.
                for bond2 in range(bond1 + 1, nb + 1):
                    # Get the bond angle of the current bond.
                    current_ba_count += 1
                    ba = sc.bond_angles_ext[atom][
                        current_ba_count
                    ]

                    # Compute which bucket in the histogram
                    # should be incremented.
                    bucket = (
                        int((ba - min_ba) / s.hist_delta)
                        + 1
                    )

                    # Increment the histogram.
                    histogram[bucket] += 1

        # Now that the histogram has been computed, print it.
        with open(s.out_file + ".hist", "w") as out:
            start_point = int(min_ba)
            for bucket in range(1, num_points + 1):
                current_point = (
                    start_point
                    + s.hist_delta * (bucket - 1)
                )
                out.write(
                    f"{current_point} "
                    f"{histogram[bucket]}\n"
                )

    # ---------------------------------------------------------------
    # Bond lengths output
    # ---------------------------------------------------------------

    def _print_bond_lengths(self):
        """Print bond length data and optional histogram."""
        sc = self.sc
        s = self.s
        coords = self.coords

        # Initialize the min/max BLs for potential use in making
        # a histogram.
        min_bl = 100.0
        max_bl = 0.0

        with open(s.out_file, "w") as out:
            for atom in range(1, sc.num_atoms + 1):
                # Compare the position of this atom with the box
                # the user defined (or the default).
                if sc.item_out_of_bounds(
                    atom, coords
                ):
                    continue

                # Print the ID data, the coordination, and the
                # number of bonds for this atom.
                name = sc.atom_element_name[atom]
                spec = sc.atom_species_id[atom]
                nb = sc.num_bonds[atom]
                out.write(
                    f"{name}{spec}_{atom} "
                    f"Num_bonds:  {nb}\n"
                )

                # Print each bond for this atom.
                for bond in range(1, nb + 1):
                    # Obtain the central cell atom number for the
                    # atom at the other end of the current bond.
                    bonded_atom = (
                        sc.ext_to_central_item_map[
                            sc.bonded_ext[atom][bond]
                        ]
                    )

                    bl = sc.bond_length_ext[atom][bond]

                    bn = sc.atom_element_name[bonded_atom]
                    bs = sc.atom_species_id[bonded_atom]
                    out.write(
                        f"    {bn}{bs}_{bonded_atom}"
                        f" {bl:8.5f}"
                    )

                    if bond % 4 == 0:
                        out.write("\n")

                    # Now that the bond has been printed, compute
                    # the bounds needed in case a histogram was
                    # asked for.
                    if bl > max_bl:
                        max_bl = bl
                    if bl < min_bl:
                        min_bl = bl

                if nb % 4 != 0:
                    out.write("\n")

        # Once the traditional bond length list has been created,
        # now we compute the histogram if asked.
        if s.do_hist:
            self._write_length_histogram(min_bl, max_bl)

    def _write_length_histogram(self, min_bl, max_bl):
        """Write a histogram of bond lengths."""
        sc = self.sc
        s = self.s
        coords = self.coords

        # Compute the number of points in the histogram.
        num_points = (
            int((max_bl - min_bl) / s.hist_delta) + 1
        )

        # Initialize the histogram.
        histogram = [0] * (num_points + 1)

        # Revisit all the atoms to produce the histogram of every
        # bond length.
        for atom in range(1, sc.num_atoms + 1):
            if sc.item_out_of_bounds(atom, coords):
                continue

            # Consider one bond.
            for bond in range(
                1, sc.num_bonds[atom] + 1
            ):
                # Get the bond length of the current bond.
                bl = sc.bond_length_ext[atom][bond]

                # Compute which bucket in the histogram should
                # be incremented.
                bucket = (
                    int((bl - min_bl) / s.hist_delta) + 1
                )

                # Increment the histogram.
                histogram[bucket] += 1

        # Now that the histogram has been computed, print it.
        with open(s.out_file + ".hist", "w") as out:
            start_point = int(min_bl)
            for bucket in range(1, num_points + 1):
                current_point = (
                    start_point
                    + s.hist_delta * (bucket - 1)
                )
                out.write(
                    f"{current_point} "
                    f"{histogram[bucket]}\n"
                )

    # ---------------------------------------------------------------
    # Coordination output
    # ---------------------------------------------------------------

    def _print_coordination_data(self):
        """Print coordination data for each atom and element."""
        sc = self.sc
        s = self.s
        coords = self.coords

        with open(s.out_file, "w") as out:
            for atom in range(1, sc.num_atoms + 1):
                # Compare the position of this atom with the box
                # the user defined (or the default).
                if sc.item_out_of_bounds(
                    atom, coords
                ):
                    continue

                # Print the ID data and the coordination for
                # this atom.
                name = sc.atom_element_name[atom]
                spec = sc.atom_species_id[atom]
                coord = sc.coordination[atom]
                out.write(
                    f"{name}{spec}_{atom} {coord}\n"
                )

            out.write(
                f"\nELEMENT DATA for "
                f"{sc.num_elements} elements\n"
            )

            for elem in range(1, sc.num_elements + 1):
                ename = sc.element_list[elem]
                summary = sc.coordination_summary[elem]
                out.write(f"{ename} {summary}\n")

    # ---------------------------------------------------------------
    # Q^n distribution output
    # ---------------------------------------------------------------

    def _print_qn_distribution(self):
        """Print the Q^n distribution."""
        sc = self.sc
        s = self.s

        with open(s.out_file, "w") as out:
            out.write(
                f"ATOM Q^n {sc.num_qn_atoms}\n"
            )
            for atom in range(1, sc.num_atoms + 1):
                if sc.atom_qn[atom] != -1:
                    name = sc.atom_element_name[atom]
                    qn = sc.atom_qn[atom]
                    out.write(f"{atom} {name} {qn}\n")

            out.write("\nAbsolute Q^n\n")
            for n in range(0, 9):
                out.write(
                    f"Q^{n} "
                    f"{sc.absolute_sys_qn[n]}\n"
                )

            out.write("\nFractional Q^n\n")
            for n in range(0, 9):
                out.write(
                    f"Q^{n} "
                    f"{sc.fractional_sys_qn[n]}\n"
                )

            out.write("\nNetwork Connectivity\n")
            out.write(f"{sc.net_conn_qn}\n")

    # ---------------------------------------------------------------
    # Ring distribution output
    # ---------------------------------------------------------------

    def _print_ring_distribution(self):
        """Print the ring distribution."""
        sc = self.sc
        s = self.s

        with open(s.out_file, "w") as out:
            out.write("LABEL : SEQUENCE_NUM : COL_LABELS\n")
            for ring_len in range(
                s.min_ring_len, s.max_ring_len + 1
            ):
                for ring in range(
                    1, sc.ring_counts[ring_len] + 1
                ):
                    atoms = " ".join(
                        str(a) for a in
                        sc.rings[ring_len][ring]
                    )
                    out.write(
                        f"{ring_len}_{ring} : "
                        f"{atoms} : TOTAL\n"
                    )

    # ---------------------------------------------------------------
    # Bond orientational order output
    # ---------------------------------------------------------------

    def _print_bond_oo(self):
        """Print bond orientational order data."""
        sc = self.sc
        s = self.s

        with open(s.out_file, "w") as out:
            # Print the number of atoms in the model.
            out.write(f"{sc.num_atoms}\n")

            # Print the data for each atom.
            for atom in range(1, sc.num_atoms + 1):
                # Print the atom number and bond orientational
                # order parameter.
                out.write(
                    f"{atom:10d} "
                    f"{sc.q_order[atom]:12.8f}\n"
                )

    # ---------------------------------------------------------------
    # OpenDX output
    # ---------------------------------------------------------------

    def _print_bond_dx(self):
        """Print an openDX format visualization file.

        The openDX file contains all atom positions and bond
        connections for the periodic system, along with the data
        values needed for visualization.
        """
        sc = self.sc
        s = self.s

        # Get the extended cell coordinates.
        coords_ext = sc.ext_direct_xyz_list

        # Get the list of bonds to show.
        reduced_bond_info, reduced_num_bonds = (
            self._get_bonds_to_show()
        )

        with open(s.out_file, "w") as out:
            # Print the atomic positions of all relevant atoms in
            # the periodic system.
            out.write(
                "object 1 class array type float rank 1 "
                f"shape 3 items {sc.num_atoms_ext} "
                "data follows\n"
            )
            for atom in range(1, sc.num_atoms_ext + 1):
                out.write(
                    f"{coords_ext[atom][1]:15.11f} "
                    f"{coords_ext[atom][2]:15.11f} "
                    f"{coords_ext[atom][3]:15.11f} \n"
                )
            out.write("\n")

            # Reduce all the atom numbers associated with each
            # bond by 1 because OpenDX starts counting everything
            # from 0.
            for bond in range(1, reduced_num_bonds + 1):
                reduced_bond_info[bond][1] -= 1
                reduced_bond_info[bond][2] -= 1

            # Print the bond pairs for all the relevant bonds in
            # the periodic system.
            out.write(
                "object 2 class array type int rank 1 "
                f"shape 2 items {reduced_num_bonds} "
                "data follows\n"
            )
            for bond in range(1, reduced_num_bonds + 1):
                out.write(
                    f"{reduced_bond_info[bond][1]} "
                    f"{reduced_bond_info[bond][2]}\n"
                )
            out.write(
                'attribute "ref" string "positions"\n'
                'attribute "element type" '
                'string "lines"\n'
                "\n"
            )

            # Now we list the values for those connections.  The
            # values are all constant so that the bond lines are
            # the same color.
            out.write(
                "object 3 class array type float rank 0"
                f" items {reduced_num_bonds} "
                "data follows\n"
            )
            for _ in range(1, reduced_num_bonds + 1):
                out.write("5\n")  # Arbitrary value.
            out.write(
                'attribute "dep" string "connections"\n'
                "\n"
            )

            # Finally we combine the components into a field with
            # the name "bo" so that the bond order part of the
            # openDX program can read it.
            out.write(
                'object "bo" class field\n'
                'component "positions" value 1\n'
                'component "connections" value 2\n'
                'component "data" value 3\n'
                "\n"
            )

            # The DX file is complete now so we can end it.  It is
            # important to note that the DX file format requires a
            # blank line at the end.  So the \\n\\n *must* be
            # there.
            out.write("end\n\n")

    # ---------------------------------------------------------------
    # Bond statistics output
    # ---------------------------------------------------------------

    def _print_bond_stats(self):
        """Print bond length and bond angle statistics."""
        sc = self.sc
        s = self.s

        with open(s.out_file, "w") as out:
            # Print out element statistics.
            out.write("Element Statistics:\n")
            out.write(
                "   Element    "
                "   Avg. BL    "
                "   Avg. BA    "
                "  StdDV. BL   "
                "  StdDV. BA\n"
            )
            for elem in range(1, sc.num_elements + 1):
                ename = sc.element_list[elem]
                out.write(f"   {ename:<10s}")
                out.write(
                    f"{self.avg_elem_bl[elem]:13.4e} "
                )
                out.write(
                    f"{self.avg_elem_ba[elem]:13.4e}"
                )
                out.write(
                    f"{self.std_elem_bl[elem]:13.4e} "
                )
                out.write(
                    f"{self.std_elem_ba[elem]:13.4e}"
                )
                out.write("\n")

            # Make space for the next data segment.
            out.write("\n\n")

            out.write("Species Statistics:\n")
            out.write(
                "   Species    "
                "   Avg. BL    "
                "   Avg. BA    "
                "  StdDV. BL   "
                "  StdDV. BA\n"
            )
            for elem in range(1, sc.num_elements + 1):
                ns = sc.num_species[elem]
                for spec in range(1, ns + 1):
                    sname = sc.species_list[elem][spec]
                    abl = self.avg_spec_bl[elem][spec]
                    aba = self.avg_spec_ba[elem][spec]
                    sbl = self.std_spec_bl[elem][spec]
                    sba = self.std_spec_ba[elem][spec]
                    out.write(f"   {sname:<10s}")
                    out.write(f"{abl:13.4e} ")
                    out.write(f"{aba:13.4e}")
                    out.write(f"{sbl:13.4e} ")
                    out.write(f"{sba:13.4e}")
                    out.write("\n")

            # Make space for the next data segment.
            out.write("\n\n")

            out.write("Atom Statistics:\n")
            out.write(
                "   Atom       "
                "   Avg. BL    "
                "   Avg. BA    "
                "  StdDV. BL   "
                "  StdDV. BA\n"
            )
            for atom in range(1, sc.num_atoms + 1):
                if self.avg_atom_bl[atom] != -1:
                    name = sc.atom_element_name[atom]
                    abl = self.avg_atom_bl[atom]
                    aba = self.avg_atom_ba[atom]
                    sbl = self.std_atom_bl[atom]
                    sba = self.std_atom_ba[atom]
                    out.write(f"   {name:<10s}{atom}")
                    out.write(f"{abl:13.4e} ")
                    out.write(f"{aba:13.4e}")
                    out.write(f"{sbl:13.4e} ")
                    out.write(f"{sba:13.4e}")
                    out.write("\n")

    # ---------------------------------------------------------------
    # OpenSCAD output (full and parts)
    # ---------------------------------------------------------------

    def _print_openscad(self, do_parts):
        """Generate OpenSCAD model file(s).

        The basic approach that is used in this subroutine is to
        simultaneously allow for two different model printing
        methods.

        The first approach creates a complete model (cell, atoms,
        bonds) that is all conjoined as one giant union in one
        single file.  The first approach may not be easy to actually
        print although it looks nice on a computer screen.

        The second approach creates the model in a set of files
        where the components of the model are disassembled and
        printed as a set of sheets for manual assembly.  The second
        approach should be easier to print but it will require
        manual assembly.

        Algorithm for the first approach:
          Use elementary modules for depositing a sphere at any
          point in space and for depositing a rod (cylinder) at any
          point in space.  Then create a module that will deposit
          every bond and each of the eight cell vectors in space and
          join them as a giant union (even though they may not be
          touching).  The bonds and vectors are arranged as they
          would be found in the crystal model.  Do a similar thing
          for the atoms (spheres).  Finally, make a union of the
          bonds, vectors, and atoms.  The output is all in one
          single file.

        Algorithm for the second approach:
          Here, the atoms and bonds will be built and laid out in a
          set of files.  There will be a bunch of files for the
          atoms and a bunch for the bonds.  After the 3D printer
          prints the model it will be assembled by gathering all the
          pieces together and putting the rods in the holes.

          The tricks are that (1) when the atoms are printed in
          sheets, they will all have to have holes in them in
          exactly the right places (angles); (2) after printing
          there will be a bunch of pieces that all look almost
          exactly the same.  Thus, there will need to be labels on
          the pieces to distinguish them.  There will need to be
          three types of labels.  Each atom will need a label for
          the atom as a whole, then each hole on each atom will need
          a label.  Finally, each bond (rod) will need a label on
          each end to identify which atom and which hole on which
          atom that end of the bond is supposed to go in to.

          The approach that is used to make the holes in the atoms
          is to create the same structure of rods used in the "first
          approach" described above.  That is, describe a structure
          with just the bonds and the cell vectors in their
          appropriate positions in the crystalline cell in space.
          Then, position each atom (one at a time) into that model
          in its appropriate place and take a difference so that
          holes are established in the sphere at the right angles.
          Then, move the atoms to their place in a sheet.

          For the labeling of the hole, an index number is used
          (1, 2, 3, ...).  The label is placed near the hole as an
          indentation in the sphere.

          For labeling the atoms, an index symbol triplet is used
          (AAA, AAB, AAC, ..., AAZ, AAa, AAb, AAc, ..., AAz, ABA,
          ABB, ABC, ..., zzz).  This will accommodate up to
          52^3 = 140608 number of atoms in a given model.

          For labeling the bonds, the combination of atom label and
          atom-hole label are used.

        Args:
            do_parts: If False, produce one complete model file.
                      If True, produce sheets of parts for assembly.
        """
        sc = self.sc
        s = self.s
        coords = self.coords

        # Scale the atomic coordinates.
        scaled_coords = [None] * (sc.num_atoms + 1)
        scaled_coords_ext = [None] * (
            sc.num_atoms_ext + 1
        )
        coords_ext = sc.ext_direct_xyz_list

        for atom in range(1, sc.num_atoms + 1):
            scaled_coords[atom] = [None, 0.0, 0.0, 0.0]
            for xyz in range(1, 4):
                scaled_coords[atom][xyz] = (
                    s.c_scale * coords[atom][xyz]
                )
        for atom in range(1, sc.num_atoms_ext + 1):
            scaled_coords_ext[atom] = [
                None, 0.0, 0.0, 0.0
            ]
            for xyz in range(1, 4):
                scaled_coords_ext[atom][xyz] = (
                    s.c_scale * coords_ext[atom][xyz]
                )

        # Compute the maximum and minimum atom radius.
        min_atom_radius = 10000000.0
        max_atom_radius = 0.0
        scaled_atom_radius = [None] * (sc.num_atoms + 1)
        for atom in range(1, sc.num_atoms + 1):
            # Compute the scaled radius of the atom.
            scaled_atom_radius[atom] = (
                sc.coval_radii[atom] * s.a_scale
            )
            if scaled_atom_radius[atom] < min_atom_radius:
                min_atom_radius = scaled_atom_radius[atom]
            if scaled_atom_radius[atom] > max_atom_radius:
                max_atom_radius = scaled_atom_radius[atom]

        # Compute the radius of all of the bond cylinders.
        bond_radius = 0.0
        for atom in range(1, sc.num_atoms + 1):
            r = (
                s.bond_radius_scale
                * sc.coval_radii[atom]
                * s.a_scale
            )
            if r > bond_radius:
                bond_radius = r

        # Compute the final scaled radius.
        bond_radius = bond_radius * s.b_scale

        # Get the list of bonds to show.
        reduced_bond_info, reduced_num_bonds = (
            self._get_bonds_to_show()
        )

        # Scale the bond lengths.
        for bond in range(1, reduced_num_bonds + 1):
            reduced_bond_info[bond][3] *= s.c_scale

        # Create an ordered sequence of the extended atom numbers
        # associated with the set of bonds being shown.
        atom_num_index = [0] * (sc.num_atoms_ext + 1)
        atom_count = 0
        prev1 = 0
        for bond in range(1, reduced_num_bonds + 1):
            if reduced_bond_info[bond][1] != prev1:
                atom_count += 1
                atom_num_index[
                    reduced_bond_info[bond][1]
                ] = atom_count
                prev1 = reduced_bond_info[bond][1]
        prev2 = 0
        for bond in range(1, reduced_num_bonds + 1):
            if reduced_bond_info[bond][2] != prev2:
                atom_count += 1
                atom_num_index[
                    reduced_bond_info[bond][2]
                ] = atom_count
                prev2 = reduced_bond_info[bond][2]

        # Compute the scaled maximum bond length.
        max_bond_length = 0.0
        for bond in range(1, reduced_num_bonds + 1):
            if reduced_bond_info[bond][3] > max_bond_length:
                max_bond_length = (
                    reduced_bond_info[bond][3]
                )

        # Define fixed model variables.
        angle = 53.7  # Fix the angle for the key.

        # Print out the models in different ways according to
        # whether we want the full model or the parts.
        if not do_parts:
            # Everything will be printed in one file.
            with open(s.out_file, "w") as out:
                # Print the module for making the rods.
                self._write_rod_module(out)

                # Print the module for making a sphere that is
                # translated to a given point.
                self._write_trans_sphere_module(out)

                # Print the module for making the bonds and cell
                # vectors.  Note that the do_corner_spheres=True
                # flag requests that corner spheres be included
                # in the print.  These should not be included
                # when printing the bonds and vectors (BAV) as a
                # skeleton for subtracting holes (as is done in
                # the other "parts" section).
                self._write_bav_module(
                    out, min_atom_radius,
                    max_atom_radius, scaled_coords,
                    scaled_coords_ext, bond_radius,
                    reduced_bond_info,
                    reduced_num_bonds,
                    atom_num_index,
                    do_corner_spheres=True,
                    do_parts=False,
                )

                # Print the module for printing all of the atoms
                # as a big union.
                self._write_union_atoms_module(
                    scaled_coords,
                    scaled_atom_radius, out,
                )

                # Print the model.
                out.write("union()\n{\n")
                out.write("   bondsAndVectors();\n")
                out.write("   allAtoms();\n")
                out.write("}\n")
        else:
            # Compute the font size for all labels.
            font_size = bond_radius * s.f_scale

            # Print the sheets of atom hemispheres first.
            x_step = (
                2.0 * max_atom_radius
                + 0.5 * max_atom_radius
            )
            y_step = (
                2.0 * max_atom_radius
                + 0.5 * max_atom_radius
            )

            # Determine the maximum number of atom hemispheres
            # in the x and y directions.
            max_x = int(math.floor(s.x_vol_3d / x_step))
            max_y = int(math.floor(s.y_vol_3d / y_step))

            # Print the sheets of atom hemispheres.  We do the
            # top first and then the bottom.
            for hemisphere in (1, 2):
                starting_item = 1

                sheet = 0
                while starting_item <= sc.num_atoms:
                    sheet += 1

                    # Open the file for this sheet.
                    if hemisphere == 1:
                        fname = (
                            f"atomTopSheet{sheet}.scad"
                        )
                    else:
                        fname = (
                            f"atomBotSheet{sheet}.scad"
                        )

                    with open(fname, "w") as out:
                        # Set default facet values.
                        out.write(
                            "$fa=0.5; // default "
                            "minimum facet angle\n"
                        )
                        out.write(
                            "$fs=0.5; // default "
                            "minimum facet size\n"
                        )

                        # Initialize the position for the
                        # hemispheres on the sheet.
                        sheet_pos = [
                            None,
                            max_atom_radius,
                            max_atom_radius,
                            0.0,
                        ]

                        x_grid = 1
                        y_grid = 1

                        self._write_rod_module(out)
                        self._write_sphere_label_module(
                            out, font_size,
                        )
                        self._write_trans_sphere_module(
                            out,
                        )
                        self._write_trans_key_module(
                            out, angle,
                        )
                        self._write_trans_sub_cube_module(
                            out,
                        )
                        self._write_bav_module(
                            out, min_atom_radius,
                            max_atom_radius,
                            scaled_coords,
                            scaled_coords_ext,
                            bond_radius,
                            reduced_bond_info,
                            reduced_num_bonds,
                            atom_num_index,
                            do_corner_spheres=False,
                            do_parts=True,
                        )

                        last_printed = starting_item
                        for atom in range(
                            starting_item,
                            sc.num_atoms + 1,
                        ):
                            self._add_atom_to_sheet(
                                out, atom,
                                scaled_coords[atom],
                                sheet_pos,
                                scaled_atom_radius[
                                    atom
                                ],
                                angle, hemisphere,
                            )

                            x_grid += 1
                            sheet_pos[1] += x_step

                            if x_grid > max_x:
                                x_grid = 1
                                sheet_pos[1] = (
                                    max_atom_radius
                                )
                                y_grid += 1
                                sheet_pos[2] += y_step

                            if y_grid > max_y:
                                y_grid = 1
                                sheet_pos[2] = (
                                    max_atom_radius
                                )
                                last_printed = atom
                                break

                            last_printed = atom

                    starting_item = last_printed + 1

            # Now we do the bonds.
            sheet_pos1 = [0, 0.0, 0.0, bond_radius]
            sheet_pos2 = [0, 0.0, 0.0, bond_radius]

            # Select the grid step size.
            x_step = 2.0 * bond_radius + 0.5 * bond_radius
            y_step = (
                max_bond_length + 0.05 * max_bond_length
            )

            max_x = int(
                math.floor(s.x_vol_3d / x_step)
            )
            max_y = int(
                math.floor(s.y_vol_3d / y_step)
            )

            # Initialize the starting item number for the first
            # sheet.
            starting_item = 1

            # Initialize the number of bonds that each atom has
            # had printed so far.
            atom_bond_count = [0] * (sc.num_atoms + 1)

            sheet = 0
            while starting_item <= reduced_num_bonds:
                sheet += 1

                fname = f"bondSheet{sheet}.scad"
                with open(fname, "w") as out:
                    # Set default facet values.
                    out.write(
                        "$fa=0.5; // default minimum "
                        "facet angle\n"
                    )
                    out.write(
                        "$fs=0.5; // default minimum "
                        "facet size\n"
                    )

                    # Initialize the end points for the first
                    # bond on the sheet in the x-y axis.
                    sheet_pos1[1] = bond_radius
                    sheet_pos1[2] = bond_radius
                    sheet_pos2[1] = bond_radius
                    if starting_item <= reduced_num_bonds:
                        sheet_pos2[2] = (
                            bond_radius
                            + reduced_bond_info[
                                starting_item
                            ][3]
                        )

                    x_grid = 1
                    y_grid = 1

                    self._write_rod_module(out)
                    self._write_rod_label_module(
                        out, font_size,
                    )

                    last_printed = starting_item
                    for bond in range(
                        starting_item,
                        reduced_num_bonds + 1,
                    ):
                        a1 = reduced_bond_info[bond][1]
                        a2 = reduced_bond_info[bond][2]
                        atom_bond_count[a1] += 1
                        atom_bond_count[a2] += 1

                        # Print a scaled rod between these
                        # two cartesian points.
                        self._write_rod(
                            sheet_pos1, sheet_pos2,
                            bond_radius, 0, 0, out,
                            atom_num_index[a1],
                            atom_num_index[a2],
                            atom_bond_count[a1],
                            atom_bond_count[a2],
                            label_code=2,
                        )

                        x_grid += 1
                        sheet_pos1[1] += x_step
                        sheet_pos2[1] += x_step

                        if x_grid > max_x:
                            x_grid = 1
                            sheet_pos1[1] = bond_radius
                            sheet_pos2[1] = bond_radius
                            y_grid += 1
                            sheet_pos1[2] += y_step
                            sheet_pos2[2] += y_step

                        if y_grid > max_y:
                            y_grid = 1
                            sheet_pos1[2] = bond_radius
                            sheet_pos2[2] = (
                                bond_radius
                                + reduced_bond_info[
                                    bond
                                ][3]
                            )
                            last_printed = bond
                            break

                        last_printed = bond

                starting_item = last_printed + 1

    # ---------------------------------------------------------------
    # OpenSCAD helper: rod module
    # ---------------------------------------------------------------

    @staticmethod
    def _write_rod_module(out):
        """Write the OpenSCAD module for making individual rods."""
        out.write(
            "module rod(xPos,yPos,zPos,angle,"
            "xNormal,yNormal,zNormal,shift,"
            "length,radius)\n"
            "{\n"
            "   translate(v=[xPos,yPos,zPos]) {\n"
            "      rotate(a=angle, "
            "v=[xNormal,yNormal,zNormal]) {\n"
            "         translate(v=[0,0,shift]) {\n"
            "            cylinder(h=length,"
            "r=radius);}}}\n"
            "}\n"
        )

    @staticmethod
    def _write_rod_label_module(out, font_size):
        """Write the OpenSCAD module for labeled rods."""
        fs = font_size
        out.write(
            "module labeledRod(xPos,yPos,zPos,angle,"
            "xNormal,yNormal,"
            "zNormal,shift,length,radius,"
            "atomLabel1,atomLabel2,"
            "label1,label2)\n"
            "{\n"
            "   difference() {\n"
            "      translate(v=[xPos,yPos,zPos]) {\n"
            "         rotate(a=angle,"
            "v=[xNormal,yNormal,zNormal]) {\n"
            "            translate(v=[0,0,shift]) {\n"
            "               cylinder(h=length,"
            "r=radius);}}}\n"
        )
        # Label 1 near top.
        out.write(
            "      translate(v=[xPos,"
            "yPos+length-length*1/6,"
            "zPos+radius/2]) {\n"
            "         rotate(a=-90,v=[0,0,1]) {\n"
            "            linear_extrude(radius) {\n"
            '               text(atomLabel1,'
            f'font="Hack",size={fs},'
            'halign="center", '
            'valign="center");}}}\n'
        )
        # Label 2 near bottom.
        out.write(
            "      translate(v=[xPos,"
            "yPos+length*1/6,"
            "zPos+radius/2]) {\n"
            "         rotate(a=-90,v=[0,0,1]) {\n"
            "            linear_extrude(radius) {\n"
            '               text(atomLabel2,'
            f'font="Hack",size={fs},'
            'halign="center", '
            'valign="center");}}}\n'
        )
        # Number label 1.
        out.write(
            "      translate(v=[xPos,"
            "yPos+length-length*1/3,"
            "zPos+radius/2]) {\n"
            "         rotate(a=-90,v=[0,0,1]) {\n"
            "            linear_extrude(radius) {\n"
            '               text(label1,'
            f'font="Hack",size={fs},'
            'halign="center", '
            'valign="center");}}}\n'
        )
        # Number label 2.
        out.write(
            "      translate(v=[xPos,"
            "yPos+length*1/3,"
            "zPos+radius/2]) {\n"
            "         rotate(a=-90,v=[0,0,1]) {\n"
            "            linear_extrude(radius) {\n"
            '               text(label2,'
            f'font="Hack",size={fs},'
            'halign="center", '
            'valign="center");}}}\n'
        )
        out.write("   }\n}\n")

    @staticmethod
    def _write_sphere_label_module(out, font_size):
        """Write the OpenSCAD module for sphere labels.

        This module is essentially the same as the rod module
        except that instead of printing a cylinder, this will
        create a requested number at the radius of the targeted
        atom.  Thus, the number will appear on the surface of an
        atom that has this number differenced out of it.
        """
        fs = font_size
        out.write(
            "module sphereLabel(xPos,yPos,zPos,angle,"
            "xNormal,yNormal,"
            "zNormal,shift,distance,thickness,"
            "radius,label)\n"
            "{\n"
            "   translate(v=[xPos,yPos,zPos]) {\n"
            "      rotate(a=angle, "
            "v=[xNormal,yNormal,zNormal]) {\n"
            "         translate(v=[radius,radius,"
            "distance-thickness/2]) {\n"
            "            linear_extrude(thickness) {\n"
            '               text(label,'
            f'font="Hack",size={fs},'
            'halign="center",'
            'valign="center");}}}}\n'
            "}\n"
        )

    @staticmethod
    def _write_trans_sphere_module(out):
        """Write the OpenSCAD module for translated spheres."""
        out.write(
            "module transSphere(xPos,yPos,zPos,"
            "radius)\n"
            "{\n"
            "   translate(v=[xPos,yPos,zPos]) "
            "{sphere(r=radius);}\n"
            "}\n"
        )

    @staticmethod
    def _write_trans_key_module(out, angle):
        """Write the OpenSCAD module for translated keys."""
        out.write(
            "module transKey(xPos,yPos,zPos,length)\n"
            "{\n"
            "   translate(v=[xPos,yPos,zPos]) {\n"
            f"      rotate([0,0,{angle}]) "
            "{cube(size=length,center=true);}\n"
            "   }\n"
            "}\n"
        )

    @staticmethod
    def _write_trans_sub_cube_module(out):
        """Write the OpenSCAD module for translated cubes."""
        out.write(
            "module transSubCube(xPos,yPos,zPos,"
            "length)\n"
            "{\n"
            "   translate(v=[xPos,yPos,zPos]) "
            "{cube(size=length,center=true);}\n"
            "}\n"
        )

    def _write_bav_module(
        self, out, min_atom_radius, max_atom_radius,
        scaled_coords, scaled_coords_ext,
        bond_radius, reduced_bond_info,
        reduced_num_bonds, atom_num_index,
        do_corner_spheres, do_parts,
    ):
        """Write the bonds-and-vectors (BAV) OpenSCAD module.

        This module contains all the bonds and lattice vectors
        combined as a single union.
        """
        sc = self.sc
        s = self.s

        out.write("module bondsAndVectors()\n{\n")
        out.write("   union()\n   {\n")

        # Create a scaled set of rods that define the cell.  The two
        # endpoint vectors are built in the 1-indexed [None, a, b, c]
        # layout that get_direct_xyz_point expects so that the modular
        # axis-rotation trick below can be expressed with 1-indexed
        # axis values directly.
        for abc in range(1, 4):
            for i in range(2):
                for j in range(2):
                    # This will make ABC fractional points for all cell
                    # edges.  The (abc + k) % 3 formula cycles through
                    # the three axes; we translate the 0-based cycle back
                    # into 1-indexed slots by adding one.
                    p1_abc = [None, 0.0, 0.0, 0.0]
                    p2_abc = [None, 0.0, 0.0, 0.0]
                    p1_abc[((abc - 1 + 0) % 3) + 1] = 0
                    p1_abc[((abc - 1 + 1) % 3) + 1] = i
                    p1_abc[((abc - 1 + 2) % 3) + 1] = j
                    p2_abc[((abc - 1 + 0) % 3) + 1] = 1
                    p2_abc[((abc - 1 + 1) % 3) + 1] = i
                    p2_abc[((abc - 1 + 2) % 3) + 1] = j

                    # Convert from ABC fractional to XYZ
                    # cartesian.  get_direct_xyz_point now takes and
                    # returns 1-indexed [None, v1, v2, v3] vectors.
                    p1_xyz = sc.get_direct_xyz_point(p1_abc)
                    p2_xyz = sc.get_direct_xyz_point(p2_abc)

                    # Scale the locations (slots 1..3 hold the xyz data).
                    for xyz in range(1, 4):
                        p1_xyz[xyz] *= s.c_scale
                        p2_xyz[xyz] *= s.c_scale

                    # Print a scaled rod.
                    self._write_rod(
                        p1_xyz, p2_xyz, bond_radius,
                        min_atom_radius,
                        max_atom_radius, out,
                        0, 0, 0, 0,
                        label_code=0,
                    )

                    if do_corner_spheres:
                        self._write_sphere(
                            p1_xyz,
                            min_atom_radius, out,
                        )
                        self._write_sphere(
                            p2_xyz,
                            min_atom_radius, out,
                        )

        # In the case that the model is not printed by parts,
        # simply create a bond in space for each bond.
        if not do_parts:
            for bond in range(
                1, reduced_num_bonds + 1
            ):
                p1 = [None, 0.0, 0.0, 0.0]
                p2 = [None, 0.0, 0.0, 0.0]
                for xyz in range(1, 4):
                    p1[xyz] = scaled_coords_ext[
                        reduced_bond_info[bond][1]
                    ][xyz]
                    p2[xyz] = scaled_coords_ext[
                        reduced_bond_info[bond][2]
                    ][xyz]

                self._write_rod(
                    p1, p2, bond_radius,
                    min_atom_radius,
                    max_atom_radius, out,
                    0, 0, 0, 0,
                    label_code=0,
                )
        else:
            # If the model is being printed in parts, labels
            # need to be printed in addition to simply printing
            # the bond in space.  These labels are the labels
            # that will appear on the spheres.
            atom_bond_count = [0] * (sc.num_atoms + 1)

            for bond in range(
                1, reduced_num_bonds + 1
            ):
                p1 = [None, 0.0, 0.0, 0.0]
                p2 = [None, 0.0, 0.0, 0.0]
                for xyz in range(1, 4):
                    p1[xyz] = scaled_coords_ext[
                        reduced_bond_info[bond][1]
                    ][xyz]
                    p2[xyz] = scaled_coords_ext[
                        reduced_bond_info[bond][2]
                    ][xyz]

                a1 = reduced_bond_info[bond][1]
                a2 = reduced_bond_info[bond][2]
                atom_bond_count[a1] += 1
                atom_bond_count[a2] += 1

                label1 = atom_bond_count[a1]
                label2 = atom_bond_count[a2]

                self._write_rod(
                    p1, p2, bond_radius,
                    min_atom_radius,
                    max_atom_radius, out,
                    atom_num_index[a1],
                    atom_num_index[a2],
                    label1, label2,
                    label_code=1,
                )

        # Close the union of bond and cell vector parts.
        out.write("   }\n")
        # Close the BAV module.
        out.write("}\n")

    def _write_union_atoms_module(
        self, scaled_coords, scaled_atom_radius, out
    ):
        """Write the allAtoms OpenSCAD module."""
        sc = self.sc

        out.write("module allAtoms()\n{\n")
        out.write("   union()\n   {\n")

        # Print all of the atoms.
        for atom in range(1, sc.num_atoms + 1):
            self._write_sphere(
                scaled_coords[atom],
                scaled_atom_radius[atom], out,
            )

        # Close the union of atoms.
        out.write("   }\n")
        # Close the atom module.
        out.write("}\n")

    def _write_rod(
        self, point1, point2, rod_radius,
        min_sphere_radius, max_sphere_radius,
        out, atom_label1, atom_label2,
        label1, label2, label_code,
    ):
        """Write a single rod (bond/vector) to the OpenSCAD file.

        Args:
            point1, point2: End points [None, x, y, z].
            rod_radius: Radius of the cylinder.
            min_sphere_radius: Min atomic sphere radius.
            max_sphere_radius: Max atomic sphere radius.
            out: File handle.
            atom_label1, atom_label2: Atom index labels.
            label1, label2: Bond count labels.
            label_code: 0=no labels, 1=labels on spheres,
                        2=labels on rods.
        """
        sc = self.sc

        # Compute the rod length between defined points.
        rod_length = 0.0
        for xyz in range(1, 4):
            rod_length += (
                (point2[xyz] - point1[xyz]) ** 2
            )
        rod_length = math.sqrt(rod_length)

        # Adjust the rod length so that it penetrates only so
        # far into the sphere.  Presently this is hard coded to
        # go 1/3 of the way into the two atoms at either end of
        # the bond.
        rod_length -= min_sphere_radius / 3.0 * 2.0

        # Compute the shift of the rod position so that the full
        # assembly does not have touching rods.
        rod_shift = min_sphere_radius / 3.0

        # Compute the positions of the rod's endpoints when one
        # endpoint is centered at the origin.
        shifted_p2 = [None, 0.0, 0.0, 0.0]
        for xyz in range(1, 4):
            shifted_p2[xyz] = (
                point2[xyz] - point1[xyz]
            )

        # Compute the angle of rotation for the rod (noting that
        # rotations are positive in the counterclockwise
        # direction).
        if (abs(shifted_p2[1]) < 0.00000001
                and abs(shifted_p2[2]) < 0.00000001):
            angle = 0
        else:
            angle = (
                -90.0
                + math.atan2(
                    shifted_p2[3],
                    math.sqrt(
                        shifted_p2[1] ** 2
                        + shifted_p2[2] ** 2
                    ),
                ) * 180.0 / math.pi
            )

        # Obtain the vector that is normal to the plane formed
        # by three points:  (1) The origin; (2) A point on the
        # z-axis at a distance equal to the length of the rod;
        # (3) The shifted end point of the rod.
        normal = sc.get_plane_normal(
            [0, 0, 0, 0],
            shifted_p2,
            [0, 0, 0, rod_length],
        )

        # Use the rod module to print this rod.  This is used
        # for labelCode 0 and 1.
        if label_code in (0, 1):
            out.write(
                f"      rod({point1[1]},{point1[2]},"
                f"{point1[3]},{angle},{normal[1]},"
                f"{normal[2]},{normal[3]},"
                f"{rod_shift},{rod_length},"
                f"{rod_radius});\n"
            )

        # If the labelCode is 1, print the label number next to
        # the rod (so that it will appear on the sphere when
        # subtracted away).
        if label_code == 1:
            # Compute the label thickness.
            label_thickness = (
                (max_sphere_radius - min_sphere_radius)
                * 1.3
            )
            if label_thickness == 0:
                label_thickness = rod_radius

            # Compute the distance from the center of the sphere
            # that the label should be positioned.
            label_distance = (
                (max_sphere_radius + min_sphere_radius)
                / 2.0
            )

            # Print the labels for the holes for this bond rod.
            out.write(
                f"sphereLabel({point1[1]},"
                f"{point1[2]},{point1[3]},"
                f"{angle + 10},{normal[1]},"
                f"{normal[2]},{normal[3]},0,"
                f"{label_distance},"
                f"{label_thickness},"
                f'{rod_radius},"{label1}");\n'
            )
            out.write(
                f"sphereLabel({point2[1]},"
                f"{point2[2]},{point2[3]},"
                f"{angle + 190},{normal[1]},"
                f"{normal[2]},{normal[3]},0,"
                f"{label_distance},"
                f"{label_thickness},"
                f'{rod_radius},"{label2}");\n'
            )

            # In the special case that the first bond for this
            # atom is being printed, we must also print the atom
            # label.
            if label1 == 1:
                out.write(
                    f"sphereLabel({point1[1]},"
                    f"{point1[2]},{point1[3]},"
                    f"{angle - 20},{normal[1]},"
                    f"{normal[2]},{normal[3]},0,"
                    f"{label_distance},"
                    f"{label_thickness},"
                    f'{rod_radius},"{atom_label1}");\n'
                )
            if label2 == 1:
                out.write(
                    f"sphereLabel({point2[1]},"
                    f"{point2[2]},{point2[3]},"
                    f"{angle + 160},{normal[1]},"
                    f"{normal[2]},{normal[3]},0,"
                    f"{label_distance},"
                    f"{label_thickness},"
                    f'{rod_radius},"{atom_label2}");\n'
                )

        # Print labels onto the rods if requested.
        if label_code == 2:
            out.write(
                f"labeledRod({point1[1]},"
                f"{point1[2]},{point1[3]},"
                f"{angle},{normal[1]},{normal[2]},"
                f"{normal[3]},{rod_shift},"
                f"{rod_length},{rod_radius},"
                f'"{atom_label1}","{atom_label2}",'
                f'"{label1}","{label2}");\n'
            )

    @staticmethod
    def _write_sphere(point, radius, out):
        """Write a single sphere to the OpenSCAD file."""
        out.write(
            f"      transSphere({point[1]},"
            f"{point[2]},{point[3]},{radius});\n"
        )

    @staticmethod
    def _add_atom_to_sheet(
        out, atom, atom_pos, sheet_pos,
        radius, angle, hemisphere,
    ):
        """Add an atom hemisphere to the sheet."""
        if hemisphere == 1:
            rotation = 0
        else:
            rotation = 180

        # Compute the translation vector for the atom from its
        # position within the cell to its position on the sheet.
        trans = [
            None,
            sheet_pos[1] - atom_pos[1],
            sheet_pos[2] - atom_pos[2],
            sheet_pos[3] - atom_pos[3],
        ]

        # Now we need to translate the difference object that is
        # taken between the atom and the skeleton.  The
        # difference object will have holes in all the right
        # places and it will be translated to the correct spot
        # on the sheet.
        out.write(f"rotate([0,{rotation},0]) {{\n")
        out.write(
            f"   translate(v=[{trans[1]},"
            f"{trans[2]},{trans[3]}]) {{\n"
        )
        out.write("      difference() {\n   ")
        out.write("         union() {\n   ")
        out.write(
            f"      transSphere({atom_pos[1]},"
            f"{atom_pos[2]},{atom_pos[3]},"
            f"{radius});\n"
        )

        # Print key.
        key_size = 1.00
        key_x = (
            atom_pos[1]
            + radius * math.cos(angle * math.pi / 180.0)
        )
        key_y = (
            atom_pos[2]
            + radius * math.sin(angle * math.pi / 180.0)
        )
        key_z = atom_pos[3]
        out.write(
            f"      transKey({key_x},{key_y},"
            f"{key_z},{key_size});\n"
        )

        out.write("         }\n")
        out.write("         bondsAndVectors();\n")

        # Print the subtraction cube.
        cube_x = atom_pos[1]
        cube_y = atom_pos[2]
        if hemisphere == 1:
            cube_z = atom_pos[3] - radius
        else:
            cube_z = atom_pos[3] + radius
        cube_size = radius * 2.0 + 0.01 * radius
        out.write(
            f"         transSubCube({cube_x},"
            f"{cube_y},{cube_z},{cube_size});\n"
        )

        out.write("      }\n")
        out.write("   }\n")
        out.write("}\n")

    # ---------------------------------------------------------------
    # Get bonds to show (with de-duplication and optional BO filter)
    # ---------------------------------------------------------------

    def _get_bonds_to_show(self):
        """Build the reduced (de-duplicated) list of bonds.

        Prior to assembling the ordered list of bonded atoms with
        additional meta data, we will read in the bond order data
        if requested.  This data will be used to further reduce the
        list of bonds to show.

        Returns:
            (reduced_bond_info, reduced_num_bonds):
                reduced_bond_info is a 1-indexed list where each
                entry is [None, ext_atom1, ext_atom2, bond_length].
                reduced_num_bonds is the count.
        """
        sc = self.sc
        s = self.s

        # Read in bond order data if requested.
        bond_bo_list = []
        if s.use_olcao_bo:
            with open(s.bond_order_file) as f:
                # Read past the BO header.
                f.readline()
                for line in f:
                    values = line.split()
                    bond_bo_list.append(
                        (int(values[2]), int(values[3]))
                    )

        # Extract all the bonds to be shown and count them up
        # along the way.
        bond_info_temp = []
        for atom in range(1, sc.num_atoms + 1):
            # Get the atom number of the current atom in the
            # extended cell list.
            bonded_atom1 = (
                sc.central2ext_item_map[atom]
            )

            # Consider each bonded atom in turn.
            for bond in range(
                1, sc.num_bonds[atom] + 1
            ):
                # Note that initially this double counts bonds
                # between atoms that are both in the central
                # cell.
                bonded_atom2 = sc.bonded_ext[atom][bond]

                # Demand that the small atom number is listed
                # first for comparison later to remove the
                # double counts.
                if bonded_atom1 < bonded_atom2:
                    a1, a2 = bonded_atom1, bonded_atom2
                else:
                    a1, a2 = bonded_atom2, bonded_atom1

                # Record bond length info for this bond pair.
                bond_info_temp.append(
                    [a1, a2,
                     sc.bond_length_ext[atom][bond]]
                )

        # Sort the bonds according to the second bonded atom
        # followed by a stable sort according to the first
        # bonded atom so that bonds that are between identical
        # atom pairs are next door.
        bond_info_temp.sort(key=lambda x: x[1])
        bond_info_temp.sort(key=lambda x: x[0])

        # Now we need to remove duplicate bonds by simply
        # recording any non-dups.  At the same time, we will
        # only keep bonds that correspond to bonds listed in the
        # bond order file if asked.
        reduced_bond_info = [None]  # 1-indexed.
        reduced_num_bonds = 0

        for entry in bond_info_temp:
            a1, a2, bl = entry

            # If we are limited only to bonds that are present
            # in the bond order list, make sure this bond is one
            # of them.
            if s.use_olcao_bo:
                c1 = sc.ext_to_central_item_map[a1]
                c2 = sc.ext_to_central_item_map[a2]
                found = False
                for bo_a1, bo_a2 in bond_bo_list:
                    if c1 == bo_a1 and c2 == bo_a2:
                        found = True
                        break
                if not found:
                    continue

            # Search the reduced bond list to make sure that the
            # current bond isn't a duplicate.
            found = False
            for rb in range(1, reduced_num_bonds + 1):
                if (a1 == reduced_bond_info[rb][1]
                        and a2 == reduced_bond_info[rb][2]):
                    found = True
                    break

            # If the bond was found, it's a duplicate and we
            # should not use it.  Otherwise, add it to the
            # reduced list.
            if not found:
                reduced_num_bonds += 1
                reduced_bond_info.append(
                    [None, a1, a2, bl]
                )

        return reduced_bond_info, reduced_num_bonds

    # ---------------------------------------------------------------
    # Statistics computation
    # ---------------------------------------------------------------

    def _compute_statistics(self):
        """Compute bond length and bond angle statistics.

        The goal here is to compute the mean and 1st standard
        deviation of the bond lengths (BLs) and bond angles (BAs)
        for each element, species, and atom.
        """
        sc = self.sc
        s = self.s
        coords = self.coords
        ne = sc.num_elements
        na = sc.num_atoms

        # Initialize counts and accumulators.  Initialize all
        # average values to -1 (sentinel for "not computed").
        num_elem_bonds = [0] * (ne + 1)
        num_elem_ba = [0] * (ne + 1)
        acc_elem_bl = [0.0] * (ne + 1)
        acc_elem_ba = [0.0] * (ne + 1)
        self.avg_elem_bl = [-1] * (ne + 1)
        self.avg_elem_ba = [-1] * (ne + 1)

        # Species arrays are 2D: [element][species].
        num_spec_bonds = [None] * (ne + 1)
        num_spec_ba = [None] * (ne + 1)
        acc_spec_bl = [None] * (ne + 1)
        acc_spec_ba = [None] * (ne + 1)
        self.avg_spec_bl = [None] * (ne + 1)
        self.avg_spec_ba = [None] * (ne + 1)

        for elem in range(1, ne + 1):
            ns = sc.num_species[elem]
            num_spec_bonds[elem] = [0] * (ns + 1)
            num_spec_ba[elem] = [0] * (ns + 1)
            acc_spec_bl[elem] = [0.0] * (ns + 1)
            acc_spec_ba[elem] = [0.0] * (ns + 1)
            self.avg_spec_bl[elem] = [-1] * (ns + 1)
            self.avg_spec_ba[elem] = [-1] * (ns + 1)

        # Atom arrays.
        num_atom_bonds = [0] * (na + 1)
        num_atom_ba = [0] * (na + 1)
        acc_atom_bl = [0.0] * (na + 1)
        acc_atom_ba = [0.0] * (na + 1)
        self.avg_atom_bl = [-1] * (na + 1)
        self.avg_atom_ba = [-1] * (na + 1)

        # Pass 1: Compute averages.
        for atom in range(1, na + 1):
            if sc.item_out_of_bounds(atom, coords):
                continue

            eid = sc.atom_element_id[atom]
            sid = sc.atom_species_id[atom]

            # Accumulate bond lengths.
            for bond in range(
                1, sc.num_bonds[atom] + 1
            ):
                num_atom_bonds[atom] += 1
                num_elem_bonds[eid] += 1
                num_spec_bonds[eid][sid] += 1
                bl = sc.bond_length_ext[atom][bond]
                acc_atom_bl[atom] += bl
                acc_elem_bl[eid] += bl
                acc_spec_bl[eid][sid] += bl

            # Accumulate bond angles.
            num_curr_ba = 0
            nb = sc.num_bonds[atom]
            for bond1 in range(1, nb):
                for bond2 in range(bond1 + 1, nb + 1):
                    num_curr_ba += 1
                    num_atom_ba[atom] += 1
                    num_elem_ba[eid] += 1
                    num_spec_ba[eid][sid] += 1
                    ba = sc.bond_angles_ext[atom][
                        num_curr_ba
                    ]
                    acc_atom_ba[atom] += ba
                    acc_elem_ba[eid] += ba
                    acc_spec_ba[eid][sid] += ba

        # Compute element averages.
        for elem in range(1, ne + 1):
            if num_elem_bonds[elem] == 0:
                self._element_warning(elem, 1)
            else:
                self.avg_elem_bl[elem] = (
                    acc_elem_bl[elem]
                    / num_elem_bonds[elem]
                )
            if num_elem_ba[elem] == 0:
                self._element_warning(elem, 2)
            else:
                self.avg_elem_ba[elem] = (
                    acc_elem_ba[elem]
                    / num_elem_ba[elem]
                )

            for spec in range(
                1, sc.num_species[elem] + 1
            ):
                if num_spec_bonds[elem][spec] == 0:
                    self._species_warning(elem, spec, 1)
                else:
                    self.avg_spec_bl[elem][spec] = (
                        acc_spec_bl[elem][spec]
                        / num_spec_bonds[elem][spec]
                    )
                if num_spec_ba[elem][spec] == 0:
                    self._species_warning(elem, spec, 2)
                else:
                    self.avg_spec_ba[elem][spec] = (
                        acc_spec_ba[elem][spec]
                        / num_spec_ba[elem][spec]
                    )

        # Compute atom averages.
        for atom in range(1, na + 1):
            if num_atom_bonds[atom] == 0:
                self._atom_warning(atom, 1)
            else:
                self.avg_atom_bl[atom] = (
                    acc_atom_bl[atom]
                    / num_atom_bonds[atom]
                )
            if num_atom_ba[atom] == 0:
                self._atom_warning(atom, 2)
            else:
                self.avg_atom_ba[atom] = (
                    acc_atom_ba[atom]
                    / num_atom_ba[atom]
                )

        # Pass 2: Compute standard deviations.
        # Reinitialize accumulators.
        num_elem_bonds = [0] * (ne + 1)
        num_elem_ba = [0] * (ne + 1)
        acc_elem_bl = [0.0] * (ne + 1)
        acc_elem_ba = [0.0] * (ne + 1)
        self.std_elem_bl = [-1] * (ne + 1)
        self.std_elem_ba = [-1] * (ne + 1)

        num_spec_bonds = [None] * (ne + 1)
        num_spec_ba = [None] * (ne + 1)
        acc_spec_bl = [None] * (ne + 1)
        acc_spec_ba = [None] * (ne + 1)
        self.std_spec_bl = [None] * (ne + 1)
        self.std_spec_ba = [None] * (ne + 1)

        for elem in range(1, ne + 1):
            ns = sc.num_species[elem]
            num_spec_bonds[elem] = [0] * (ns + 1)
            num_spec_ba[elem] = [0] * (ns + 1)
            acc_spec_bl[elem] = [0.0] * (ns + 1)
            acc_spec_ba[elem] = [0.0] * (ns + 1)
            self.std_spec_bl[elem] = [-1] * (ns + 1)
            self.std_spec_ba[elem] = [-1] * (ns + 1)

        num_atom_bonds = [0] * (na + 1)
        num_atom_ba = [0] * (na + 1)
        acc_atom_bl = [0.0] * (na + 1)
        acc_atom_ba = [0.0] * (na + 1)
        self.std_atom_bl = [-1] * (na + 1)
        self.std_atom_ba = [-1] * (na + 1)

        for atom in range(1, na + 1):
            if sc.item_out_of_bounds(atom, coords):
                continue

            eid = sc.atom_element_id[atom]
            sid = sc.atom_species_id[atom]

            # Accumulate squared differences for bond lengths.
            for bond in range(
                1, sc.num_bonds[atom] + 1
            ):
                num_atom_bonds[atom] += 1
                num_elem_bonds[eid] += 1
                num_spec_bonds[eid][sid] += 1
                bl = sc.bond_length_ext[atom][bond]
                acc_atom_bl[atom] += (
                    (bl - self.avg_atom_bl[atom]) ** 2
                )
                acc_elem_bl[eid] += (
                    (bl - self.avg_elem_bl[eid]) ** 2
                )
                acc_spec_bl[eid][sid] += (
                    (bl - self.avg_spec_bl[eid][sid])
                    ** 2
                )

            # Accumulate squared differences for bond angles.
            num_curr_ba = 0
            nb = sc.num_bonds[atom]
            for bond1 in range(1, nb):
                for bond2 in range(bond1 + 1, nb + 1):
                    num_curr_ba += 1
                    num_atom_ba[atom] += 1
                    num_elem_ba[eid] += 1
                    num_spec_ba[eid][sid] += 1
                    ba = sc.bond_angles_ext[atom][
                        num_curr_ba
                    ]
                    acc_atom_ba[atom] += (
                        (ba - self.avg_atom_ba[atom])
                        ** 2
                    )
                    acc_elem_ba[eid] += (
                        (ba - self.avg_elem_ba[eid])
                        ** 2
                    )
                    acc_spec_ba[eid][sid] += (
                        (ba - self.avg_spec_ba[eid][sid])
                        ** 2
                    )

        # Compute element standard deviations.
        for elem in range(1, ne + 1):
            if num_elem_bonds[elem] == 0:
                self._element_warning(elem, 1)
            else:
                self.std_elem_bl[elem] = math.sqrt(
                    acc_elem_bl[elem]
                    / num_elem_bonds[elem]
                )
            if num_elem_ba[elem] == 0:
                self._element_warning(elem, 2)
            else:
                self.std_elem_ba[elem] = math.sqrt(
                    acc_elem_ba[elem]
                    / num_elem_ba[elem]
                )

            for spec in range(
                1, sc.num_species[elem] + 1
            ):
                if num_spec_bonds[elem][spec] == 0:
                    self._species_warning(elem, spec, 1)
                else:
                    self.std_spec_bl[elem][spec] = (
                        math.sqrt(
                            acc_spec_bl[elem][spec]
                            / num_spec_bonds[elem][spec]
                        )
                    )
                if num_spec_ba[elem][spec] == 0:
                    self._species_warning(elem, spec, 2)
                else:
                    self.std_spec_ba[elem][spec] = (
                        math.sqrt(
                            acc_spec_ba[elem][spec]
                            / num_spec_ba[elem][spec]
                        )
                    )

        # Compute atom standard deviations.
        for atom in range(1, na + 1):
            if num_atom_bonds[atom] == 0:
                self._atom_warning(atom, 1)
            else:
                self.std_atom_bl[atom] = math.sqrt(
                    acc_atom_bl[atom]
                    / num_atom_bonds[atom]
                )
            if num_atom_ba[atom] == 0:
                self._atom_warning(atom, 2)
            else:
                self.std_atom_ba[atom] = math.sqrt(
                    acc_atom_ba[atom]
                    / num_atom_ba[atom]
                )

    # ---------------------------------------------------------------
    # Warning helpers
    # ---------------------------------------------------------------

    def _element_warning(self, element, warn_type):
        """Print a warning about an element with no bonds/angles."""
        sc = self.sc
        s = self.s
        ename = sc.element_list[element]
        bf = s.bonding_factor
        if warn_type == 1:
            print(
                f"Element {ename} has no bonds for "
                f"bondingFactor = {bf}"
            )
        else:
            print(
                f"Element {ename} has no bond angles "
                f"for bondFactor = {bf}"
            )

    def _species_warning(
        self, element, species, warn_type
    ):
        """Print a warning about a species with no bonds/angles."""
        sc = self.sc
        s = self.s
        sname = sc.species_list[element][species]
        bf = s.bonding_factor
        if warn_type == 1:
            print(
                f"Species {sname} has no bonds for "
                f"bondFactor = {bf}"
            )
        else:
            print(
                f"Species {sname} has no bond angles "
                f"for bondFactor = {bf}"
            )

    def _atom_warning(self, atom, warn_type):
        """Print a warning about an atom with no bonds/angles."""
        s = self.s
        bf = s.bonding_factor
        if warn_type == 1:
            print(
                f"Atom {atom} has no bonds for "
                f"bondFactor = {bf}"
            )
        elif warn_type == 2:
            print(
                f"Atom {atom} has no bond angles for "
                f"bondFactor = {bf}"
            )


# -------------------------------------------------------------------
# Entry point
# -------------------------------------------------------------------

def main():
    """Parse settings and run the bond analysis."""
    settings = ScriptSettings()
    analysis = BondAnalysis(settings)
    analysis.run()


if __name__ == "__main__":
    main()
