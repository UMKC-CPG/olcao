#!/usr/bin/env python3

"""mod_struct.py -- Modify an atomic structure input file.

This is the Python port of the Perl ``modStruct`` script.  It reads
an OLCAO structure file (olcao.skl by default) and applies a
user-specified sequence of modifications such as supercell generation,
translation, rotation, bond-valence analysis, filtering, vacuum
insertion, sphere/block cutting, orthorhombic conversion, perturbation,
and surface preparation.

The script follows the XYZ.py / XYZrc.py pattern:
  1. Load defaults from mod_structrc.py ($OLCAO_RC or cwd).
  2. Parse the command line with argparse.
  3. Reconcile CLI args with rc defaults.
  4. Execute the main workflow.

Operations are *modular and sequential*: the -sc, -trans, -rotA,
-rotP, -rothkl, -filter, -addvac, -cutsphere, -cutblock, -ortho,
and -perturb options can each be repeated any number of times and
are applied in the order given on the command line.


KEY OPTIONS
===========

-i INFILE
    Alternate input file (default: olcao.skl).

-o OUTFILE
    Alternate output file (default: olcao.skl.new).

-t TITLE
    Title string for the output file.  Use quotes if it contains
    whitespace.

-abcxyz ABC1 XYZ1 ABC2 XYZ2
    Define the orientation of the a,b,c lattice vectors with respect
    to the x,y,z Cartesian axes.  ABC1 is aligned with XYZ1 and ABC2
    is in the XYZ1,XYZ2 plane so that ABC3 can be in an arbitrary
    direction.  The default (1 1 2 2) causes a to align with x and b
    to be in the xy plane.  This option can only be applied once.

-sc SC1 SC2 SC3 MIRROR1 MIRROR2 MIRROR3
    Create a supercell by replicating the cell SC1, SC2, SC3 times
    along a, b, c respectively.  MIRROR values (0 or 1) control
    whether every other replica is a mirror image along that axis.
    Overrides any supercell line in the input file.

-trans TX TY TZ
    Translate all atomic positions by (TX, TY, TZ) in direct XYZ
    coordinates.  The special value "0 0 0" requests centering the
    model within the cell (applicable to molecular systems).

-rotA RX RY RZ ANGLE [-orig OX OY OZ]
    Rotate all atoms about an axis defined by the vector (RX,RY,RZ)
    in direct XYZ coordinates by ANGLE degrees.  The default origin
    for the rotation axis is (0,0,0); use -orig to change it.
    Periodic boundary conditions are retained: atoms rotated outside
    the cell are wrapped back in as if from a neighboring replica.

-rotP P1X P1Y P1Z P2X P2Y P2Z P3X P3Y P3Z [-angle ANGLE]
                                              [-orig OX OY OZ]
    Rotate all atoms about the normal of the plane defined by three
    points P1, P2, P3.  The rotation angle is either given explicitly
    with -angle or computed as the angle between the 2nd and 3rd
    points.  The default origin is (0,0,0).  PBCs are retained.

-rothkl H K L
    Reorient and resize the crystal lattice to prepare for surface
    creation.  The Miller-index direction (H,K,L) (in units of
    lattice-vector magnitudes) is rotated to coincide with the
    Cartesian z-axis.  The simulation cell is then enlarged to
    retain periodicity.
    PRESENTLY ONLY WORKS FOR SIMPLE CASES.

-filter MINDIST
    Remove atoms that are within MINDIST of another atom.
    THIS DOES NOT REALLY WORK YET.

-addvac AXIS AMOUNT
    Increase the cell size along AXIS (a, b, or c) by AMOUNT
    angstroms to create a vacuum region.  Atomic direct-space
    coordinates are retained.  A translation should typically be
    performed beforehand to choose the surface plane.

-cutsphere [-atom ATOM | -atxyz X Y Z | -atabc A B C]
           [-radius R] [-zone in|out]
    Remove atoms inside (default) or outside a sphere of radius R
    (default 4 A) centered on an atom number, an xyz coordinate,
    or an abc coordinate.

-cutblock [-zone in|out]
          <-xyz FROMX TOX FROMY TOY FROMZ TOZ |
           -abc FROMA TOA FROMB TOB FROMC TOC>
    Remove atoms inside (default) or outside a rectangular block
    defined by coordinate ranges in either xyz or abc space.

-ortho
    Convert the lattice to orthorhombic type.  Presently only
    hexagonal cells (after appropriate supercell doubling) are
    supported.

-perturb MAGNITUDE
    Apply a random shift of up to MAGNITUDE angstroms to every
    atom in the model.
"""

import argparse as ap
import math
import os
import socket
import sys
from datetime import datetime


# OP codes — one per queued operation type.
OP_SUPERCELL = 1
OP_TRANSLATE = 2
OP_ROTATE = 3
OP_FILTER = 4
OP_ADDVAC = 5
OP_ORTHO = 6
OP_PERTURB = 7
OP_ROTHKL = 8
OP_CUTBLOCK = 9
OP_CUTSPHERE = 10


def _get_plane_normal(p1, p2, p3):
    """Return the normal vector to the plane defined by three points.

    Computes the cross product of (p2-p1) x (p3-p1).  The result
    is NOT normalised to unit length; the caller must normalise.

    Args:
        p1, p2, p3: Lists of three floats [x, y, z].

    Returns:
        List of three floats [nx, ny, nz].
    """
    d1 = [p2[i] - p1[i] for i in range(3)]
    d2 = [p3[i] - p1[i] for i in range(3)]
    return [
        d1[1] * d2[2] - d1[2] * d2[1],
        d1[2] * d2[0] - d1[0] * d2[2],
        d1[0] * d2[1] - d1[1] * d2[0],
    ]


# ------------------------------------------------------------------
# ScriptSettings
# ------------------------------------------------------------------

class ScriptSettings:
    """Holds all user-controllable parameters for modStruct.

    On construction the object:
      1. Reads defaults from mod_structrc.py.
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

        # Reconcile the command line arguments with the rc defaults.
        self.reconcile(args)

        # Record the command line that was used.
        self.record_clp()

    # --------------------------------------------------------------
    # RC loading
    # --------------------------------------------------------------

    @staticmethod
    def _load_rc():
        """Import mod_structrc and return its parameter dictionary.

        Search order: cwd first (local override), then $OLCAO_RC.
        """
        cwd_rc = os.path.join(os.getcwd(), "mod_structrc.py")
        if os.path.isfile(cwd_rc):
            sys.path.insert(0, os.getcwd())
        else:
            rc_dir = os.getenv("OLCAO_RC")
            if not rc_dir:
                sys.exit(
                    "Error: $OLCAO_RC is not set and no local "
                    "mod_structrc.py found. See instructions."
                )
            sys.path.insert(0, rc_dir)

        from mod_structrc import parameters_and_defaults
        return parameters_and_defaults()

    # --------------------------------------------------------------
    # Assign rc defaults
    # --------------------------------------------------------------

    def assign_rc_defaults(self, rc):
        """Populate instance attributes from the rc dictionary."""

        # Input / output files.
        self.in_file = rc["in_file"]
        self.out_file = rc["out_file"]

        # Title: default is hostname + cwd (matches Perl behavior).
        hostname = socket.gethostname()
        self.title = hostname + os.getcwd()

        # ABC-XYZ axis assignment order.
        self.abc_order = list(rc["abc_order"])
        self.xyz_order = list(rc["xyz_order"])

        # Operations queue.  Each entry is a dict whose contents
        # depend on the operation type.  The list preserves the
        # command-line ordering so operations execute in sequence.
        self.ops = []

    # --------------------------------------------------------------
    # Command-line parsing
    # --------------------------------------------------------------

    def parse_command_line(self):
        """Build the argparse parser and parse the command line.

        Because modStruct operations are positional (order matters)
        and repeatable, we cannot use standard argparse actions for
        the operation flags.  Instead we do a manual pass over
        sys.argv to build the ordered operations list, while still
        using argparse for the simple scalar flags (-i, -o, -t,
        -abcxyz).
        """

        # --- simple flags via argparse ---
        prog_name = "mod_struct.py"

        description_text = """\
Modify an atomic structure input file.

Applies a user-specified sequence of operations (supercell, translate,
rotate, filter, add vacuum, cut sphere/block, orthorhombic conversion,
perturbation, surface preparation) to an OLCAO structure file.

Operations are modular and sequential: they can be repeated any number
of times and are applied in the order given on the command line.
"""

        epilog_text = """\
Defaults are given in ./mod_structrc.py or $OLCAO_RC/mod_structrc.py.
"""

        parser = ap.ArgumentParser(
            prog=prog_name,
            formatter_class=ap.RawDescriptionHelpFormatter,
            description=description_text,
            epilog=epilog_text,
        )
        self.add_parser_arguments(parser)

        # We parse *known* args only because the repeatable
        # operation flags are handled manually below.
        args, _ = parser.parse_known_args()

        # --- repeatable operations via manual argv scan ---
        self.scan_operations()

        return args

    def add_parser_arguments(self, parser):
        """Add the simple (non-repeatable) arguments to the
        argparse parser."""

        parser.add_argument(
            "-i", dest="in_file", type=str,
            default=self.in_file,
            help=(
                "Alternate input file. "
                f"Default: {self.in_file}"
            ),
        )

        parser.add_argument(
            "-o", dest="out_file", type=str,
            default=self.out_file,
            help=(
                "Alternate output file. "
                f"Default: {self.out_file}"
            ),
        )

        parser.add_argument(
            "-t", dest="title", type=str,
            default=self.title,
            help=(
                "Title for the output file. Use quotes if "
                "it contains whitespace. "
                f"Default: hostname+cwd"
            ),
        )

        parser.add_argument(
            "-abcxyz", nargs=4, dest="abcxyz",
            type=int, default=None,
            help=(
                "Define orientation of a,b,c lattice vectors "
                "with respect to x,y,z axes. Four integers: "
                "ABC1 XYZ1 ABC2 XYZ2. Default: 1 1 2 2"
            ),
        )

    def scan_operations(self):
        """Scan sys.argv manually to build the ordered operations
        list.

        This is necessary because operations are repeatable and
        order-dependent — standard argparse cannot handle this
        pattern naturally.  Each recognised operation flag
        consumes its arguments and appends a dict to self.ops.
        """
        argv = sys.argv[1:]
        n = 0

        while n < len(argv):
            flag = argv[n]

            if flag == "-sc":
                # -sc SC1 SC2 SC3 MIRROR1 MIRROR2 MIRROR3
                sc_vals = [
                    int(argv[n + 1]),
                    int(argv[n + 2]),
                    int(argv[n + 3]),
                ]
                mirror_vals = [
                    int(argv[n + 4]),
                    int(argv[n + 5]),
                    int(argv[n + 6]),
                ]
                n += 6
                self.ops.append({
                    "op": OP_SUPERCELL,
                    "sc": sc_vals,
                    "mirror": mirror_vals,
                })

            elif flag == "-trans":
                # -trans TX TY TZ (1-indexed for
                #   translate_atoms)
                trans = [
                    None,
                    float(argv[n + 1]),
                    float(argv[n + 2]),
                    float(argv[n + 3]),
                ]
                n += 3
                self.ops.append({
                    "op": OP_TRANSLATE,
                    "trans": trans,
                })

            elif flag == "-rotA":
                # -rotA RX RY RZ ANGLE [-orig OX OY OZ]
                rot = [
                    float(argv[n + 1]),
                    float(argv[n + 2]),
                    float(argv[n + 3]),
                ]
                angle_deg = float(argv[n + 4])
                angle_rad = angle_deg * math.pi / 180.0
                n += 4

                # Normalize the rotation axis to a unit vector.
                mag = math.sqrt(
                    rot[0] ** 2 + rot[1] ** 2 + rot[2] ** 2
                )
                rot = [r / mag for r in rot]

                # Default origin is (0,0,0) (1-indexed
                #   for translate_atoms).
                origin = [None, 0.0, 0.0, 0.0]

                # Check for the -orig sub-option (up to 2
                # passes, matching original Perl logic).
                for _ in range(2):
                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-orig"):
                        n += 1
                        origin = [
                            None,
                            float(argv[n + 1]),
                            float(argv[n + 2]),
                            float(argv[n + 3]),
                        ]
                        n += 3

                self.ops.append({
                    "op": OP_ROTATE,
                    "rot": rot,
                    "angle": angle_rad,
                    "origin": origin,
                })

            elif flag == "-rotP":
                # -rotP P1X P1Y P1Z P2X P2Y P2Z P3X P3Y P3Z
                #       [-angle ANGLE] [-orig OX OY OZ]
                p1 = [
                    float(argv[n + 1]),
                    float(argv[n + 2]),
                    float(argv[n + 3]),
                ]
                p2 = [
                    float(argv[n + 4]),
                    float(argv[n + 5]),
                    float(argv[n + 6]),
                ]
                p3 = [
                    float(argv[n + 7]),
                    float(argv[n + 8]),
                    float(argv[n + 9]),
                ]
                n += 9

                # Compute the plane normal (axis of rotation)
                # using the cross product of (p2-p1) x (p3-p1).
                # This will be normalized below.
                rot = _get_plane_normal(p1, p2, p3)

                # Normalize the rotation axis to a unit vector.
                mag = math.sqrt(
                    rot[0] ** 2 + rot[1] ** 2 + rot[2] ** 2
                )
                rot = [r / mag for r in rot]

                # Default angle: computed as the angle between
                # the 2nd and 3rd points.  Rotations are positive
                # in the counterclockwise direction.
                angle_deg = (
                    -90.0
                    + math.atan2(
                        p3[2] - p2[2],
                        math.sqrt(
                            (p3[0] - p2[0]) ** 2
                            + (p3[1] - p2[1]) ** 2
                        ),
                    )
                    * 180.0
                    / math.pi
                )

                # Default origin is (0,0,0) (1-indexed
                #   for translate_atoms).
                origin = [None, 0.0, 0.0, 0.0]

                # Check for sub-options (up to 3 passes,
                # matching original Perl logic).
                for _ in range(3):
                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-angle"):
                        n += 1
                        angle_deg = float(argv[n + 1])
                        n += 1
                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-orig"):
                        n += 1
                        origin = [
                            None,
                            float(argv[n + 1]),
                            float(argv[n + 2]),
                            float(argv[n + 3]),
                        ]
                        n += 3

                # Convert angle to radians.
                angle_rad = angle_deg * math.pi / 180.0

                self.ops.append({
                    "op": OP_ROTATE,
                    "rot": rot,
                    "angle": angle_rad,
                    "origin": origin,
                })

            elif flag == "-rothkl":
                # -rothkl H K L
                hkl = [
                    float(argv[n + 1]),
                    float(argv[n + 2]),
                    float(argv[n + 3]),
                ]
                n += 3
                self.ops.append({
                    "op": OP_ROTHKL,
                    "hkl": hkl,
                })

            elif flag == "-filter":
                # -filter MINDIST
                min_dist = float(argv[n + 1])
                n += 1
                self.ops.append({
                    "op": OP_FILTER,
                    "min_dist": min_dist,
                })

            elif flag == "-addvac":
                # -addvac AXIS AMOUNT
                axis_str = argv[n + 1]
                amount = float(argv[n + 2])
                n += 2
                if axis_str == "a":
                    axis = 1
                elif axis_str == "b":
                    axis = 2
                elif axis_str == "c":
                    axis = 3
                else:
                    sys.exit(
                        "Unacceptable axis specification. "
                        "Need a, b, or c."
                    )
                self.ops.append({
                    "op": OP_ADDVAC,
                    "axis": axis,
                    "amount": amount,
                })

            elif flag == "-cutblock":
                # -cutblock [-zone in|out]
                #   <-xyz FROMX TOX FROMY TOY FROMZ TOZ |
                #    -abc FROMA TOA FROMB TOB FROMC TOC>

                # Default zone is "in" (0).
                zone = 0
                abcxyz_flag = 1  # Default to xyz.
                block_border = [
                    [0.0, 0.0],
                    [0.0, 0.0],
                    [0.0, 0.0],
                ]

                # Scan sub-options (up to 3 passes,
                # matching original Perl logic).
                for _ in range(3):
                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-zone"):
                        zone_str = argv[n + 2]
                        if zone_str == "in":
                            zone = 0
                        elif zone_str == "out":
                            zone = 1
                        n += 2

                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-xyz"):
                        abcxyz_flag = 1
                        block_border = [
                            [float(argv[n + 2]),
                             float(argv[n + 3])],
                            [float(argv[n + 4]),
                             float(argv[n + 5])],
                            [float(argv[n + 6]),
                             float(argv[n + 7])],
                        ]
                        n += 7

                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-abc"):
                        abcxyz_flag = 0
                        block_border = [
                            [float(argv[n + 2]),
                             float(argv[n + 3])],
                            [float(argv[n + 4]),
                             float(argv[n + 5])],
                            [float(argv[n + 6]),
                             float(argv[n + 7])],
                        ]
                        n += 7

                self.ops.append({
                    "op": OP_CUTBLOCK,
                    "zone": zone,
                    "abcxyz_flag": abcxyz_flag,
                    "block_border": block_border,
                })

            elif flag == "-cutsphere":
                # -cutsphere [-atom ATOM | -atxyz X Y Z |
                #             -atabc A B C]
                #            [-radius R] [-zone in|out]

                # Defaults.
                zone = 0          # "in"
                target_atom = 0   # No specific atomic site.
                sphere_rad = 4.0  # Default radius.
                abcxyz_flag = 1   # Default to xyz.
                sphere_loc = [0.0, 0.0, 0.0]

                # Scan sub-options (up to 3 passes,
                # matching original Perl logic).
                for _ in range(3):
                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-zone"):
                        zone_str = argv[n + 2]
                        if zone_str == "in":
                            zone = 0
                        elif zone_str == "out":
                            zone = 1
                        n += 2

                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-radius"):
                        sphere_rad = float(argv[n + 2])
                        n += 2

                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-atom"):
                        target_atom = int(argv[n + 2])
                        n += 2

                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-atxyz"):
                        abcxyz_flag = 1
                        sphere_loc = [
                            float(argv[n + 2]),
                            float(argv[n + 3]),
                            float(argv[n + 4]),
                        ]
                        n += 4

                    if (n + 1 < len(argv)
                            and argv[n + 1] == "-atabc"):
                        abcxyz_flag = 0
                        sphere_loc = [
                            float(argv[n + 2]),
                            float(argv[n + 3]),
                            float(argv[n + 4]),
                        ]
                        n += 4

                self.ops.append({
                    "op": OP_CUTSPHERE,
                    "zone": zone,
                    "abcxyz_flag": abcxyz_flag,
                    "sphere_rad": sphere_rad,
                    "target_atom": target_atom,
                    "sphere_loc": sphere_loc,
                })

            elif flag == "-ortho":
                self.ops.append({"op": OP_ORTHO})

            elif flag == "-perturb":
                # -perturb MAGNITUDE
                mag = float(argv[n + 1])
                n += 1
                self.ops.append({
                    "op": OP_PERTURB,
                    "magnitude": mag,
                })

            else:
                # Skip flags handled by argparse (-i, -o, -t,
                # -abcxyz, -h, --help) and their arguments.
                # argparse already consumed them; we just need
                # to advance past them here.
                if flag in ("-i", "-o", "-t"):
                    n += 1  # Skip the value argument.
                elif flag == "-abcxyz":
                    n += 4  # Skip 4 integer arguments.
                elif flag in ("-h", "--help"):
                    pass  # argparse handles this.
                else:
                    print(
                        f"UNKNOWN COMMAND LINE PARAMETER "
                        f"{flag}. ABORTING."
                    )
                    sys.exit(1)

            n += 1

    # --------------------------------------------------------------
    # Reconcile
    # --------------------------------------------------------------

    def reconcile(self, args):
        """Merge parsed CLI arguments into the settings object."""

        # Input / output files.
        self.in_file = args.in_file
        self.out_file = args.out_file

        # Title.
        self.title = args.title

        # ABC-XYZ axis assignment order.
        if args.abcxyz is not None:
            self.abc_order[0] = args.abcxyz[0]
            self.xyz_order[0] = args.abcxyz[1]
            self.abc_order[1] = args.abcxyz[2]
            self.xyz_order[1] = args.abcxyz[3]
            # The third axis is determined by the first two:
            # |1+2-6|=3; |1+3-6|=2; |2+3-6|=1
            self.abc_order[2] = abs(
                self.abc_order[0] + self.abc_order[1] - 6
            )
            self.xyz_order[2] = abs(
                self.xyz_order[0] + self.xyz_order[1] - 6
            )

    # --------------------------------------------------------------
    # Record command line
    # --------------------------------------------------------------

    def record_clp(self):
        """Append the command line used to the ``command`` file."""
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write("Cmnd:")
            for arg in sys.argv:
                cmd.write(f" {arg}")
            cmd.write("\n\n")


# ------------------------------------------------------------------
# Operation execution
# ------------------------------------------------------------------

def execute_operations(settings, sc):
    """Loop through all queued operations and apply them in order.

    Each operation dict in settings.ops contains an 'op' key
    identifying the operation type plus whatever parameters that
    operation requires.

    Args:
        settings: ScriptSettings instance with the operations queue.
        sc: StructureControl instance holding the model.
    """
    dispatch = {
        OP_SUPERCELL: do_supercell,
        OP_TRANSLATE: do_translate,
        OP_ROTATE: do_rotate,
        OP_FILTER: do_filter,
        OP_ADDVAC: do_addvac,
        OP_ORTHO: do_ortho,
        OP_PERTURB: do_perturb,
        OP_ROTHKL: do_rothkl,
        OP_CUTBLOCK: do_cutblock,
        OP_CUTSPHERE: do_cutsphere,
    }

    for i, op in enumerate(settings.ops, 1):
        op_code = op["op"]
        handler = dispatch.get(op_code)
        if handler is None:
            print(f"Unknown operation code {op_code}. "
                  f"Skipping.")
            continue
        print(f"  Operation {i}: "
              f"{handler.__doc__.split(chr(10))[0]}")
        handler(op, sc)


# ------------------------------------------------------------------
# Individual operation helpers
# ------------------------------------------------------------------

def do_supercell(op, sc):
    """Apply a supercell expansion.

    Sets the supercell and mirror parameters on the StructureControl
    instance and then calls apply_supercell().  The supercell and
    sc_mirror arrays are 1-indexed (index 0 is a placeholder).

    Args:
        op: Operation dict with keys 'sc' ([sc1,sc2,sc3]) and
            'mirror' ([m1,m2,m3]).
        sc: StructureControl instance.
    """
    for i in range(3):
        sc.supercell[i + 1] = op["sc"][i]
        sc.sc_mirror[i + 1] = op["mirror"][i]
    sc.apply_supercell()


def do_translate(op, sc):
    """Translate all atoms by a displacement vector.

    Calls StructureControl.translate_atoms for every atom in the
    model.  The displacement is in direct XYZ (Cartesian)
    coordinates.  A displacement of [0,0,0] is treated specially
    by translate_atoms: it centers the model within the cell.

    Args:
        op: Operation dict with key 'trans' ([tx,ty,tz]).
        sc: StructureControl instance.
    """
    sc.translate_atoms(1, sc.num_atoms, op["trans"])


def do_rotate(op, sc):
    """Rotate all atoms about an axis by a given angle.

    The rotation is performed about an arbitrary origin by first
    translating all atoms so the origin becomes (0,0,0), applying
    the rotation, and then translating back.  This matches the
    Perl behavior where rotateOnePoint subtracts and re-adds the
    origin for each point.

    Periodic boundary conditions are retained by the two-pass
    algorithm inside rotate_all_atoms().

    Args:
        op: Operation dict with keys 'rot' ([rx,ry,rz]),
            'angle' (radians), and 'origin'
            ([None,ox,oy,oz], 1-indexed).
        sc: StructureControl instance.
    """
    origin = op["origin"]
    axis = op["rot"]
    angle_rad = op["angle"]

    # Convert radians to degrees for define_rot_matrix
    # which expects degrees.
    angle_deg = angle_rad * 180.0 / math.pi

    # If the origin is not (0,0,0), translate all atoms
    # so the rotation origin becomes the coordinate
    # origin.  Check axes 1..3 (origin is 1-indexed).
    has_origin = any(
        abs(origin[ax]) > 1e-12 for ax in range(1, 4)
    )
    if has_origin:
        neg_origin = [
            None,
            -origin[1], -origin[2], -origin[3],
        ]
        sc.translate_atoms(
            1, sc.num_atoms, neg_origin
        )

    # Build the rotation matrix and apply it to all
    # atoms.
    sc.define_rot_matrix(axis, angle_deg)
    sc.rotate_all_atoms()

    # Translate back if we shifted the origin.
    if has_origin:
        sc.translate_atoms(
            1, sc.num_atoms, origin
        )


def do_filter(op, sc):
    """Remove atoms that are too close to other atoms.

    Calls StructureControl.apply_filter with the minimum distance
    threshold.  Note: as noted in the original Perl script, this
    feature does not really work yet.

    Args:
        op: Operation dict with key 'min_dist'.
        sc: StructureControl instance.
    """
    sc.apply_filter(op["min_dist"])


def do_addvac(op, sc):
    """Insert vacuum along a lattice axis.

    Calls StructureControl.insert_vacuum.  The axis is specified
    as 1, 2, or 3 for a, b, or c respectively.  The amount is in
    angstroms.  Note that the ABC vectors in XYZ format will be
    reobtained after extending the ABC magnitudes only.  At
    present, the directions of ABC may change.  The assumption is
    that a is in the x direction, b is in the xy plane, etc.

    Args:
        op: Operation dict with keys 'axis' (1,2,3 for a,b,c)
            and 'amount' (angstroms).
        sc: StructureControl instance.
    """
    sc.insert_vacuum(op["axis"], op["amount"])


def do_ortho(op, sc):
    """Convert the lattice to orthorhombic type.

    This should only be used for cells that have been carefully
    analyzed.  At present only hexagonal cells can be converted
    to orthorhombic type and then only after they have been made
    into a supercell that doubles the cell size along the
    direction of the axis to be changed.  Normally, that is a
    1 2 1 supercell.

    Args:
        op: Operation dict (no extra keys needed).
        sc: StructureControl instance.
    """
    sc.make_ortho()


def do_perturb(op, sc):
    """Apply random perturbations to all atomic positions.

    Calls StructureControl.apply_perturbation with the maximum
    magnitude of the random shift (in angstroms).

    Args:
        op: Operation dict with key 'magnitude' (angstroms).
        sc: StructureControl instance.
    """
    sc.apply_perturbation(op["magnitude"])


def do_rothkl(op, sc):
    """Prepare a surface by rotating an (h,k,l) direction to z.

    Calls StructureControl.prep_surface.  The Miller-index
    direction (h,k,l) is specified in units of lattice vector
    magnitudes.  This direction is rotated so as to be coincident
    with the Cartesian z-axis, and the simulation cell is then
    enlarged to retain periodicity.

    PRESENTLY ONLY WORKS FOR SIMPLE CASES.

    Args:
        op: Operation dict with key 'hkl' ([h,k,l]).
        sc: StructureControl instance.
    """
    sc.prep_surface(op["hkl"])


def do_cutblock(op, sc):
    """Remove atoms inside or outside a rectangular block.

    Calls StructureControl.cut_block.  The block is defined by
    coordinate ranges along each of the three axes in either xyz
    or abc space.

    Args:
        op: Operation dict with keys 'zone' (0=in,1=out),
            'abcxyz_flag' (0=abc,1=xyz), and
            'block_border' ([[from1,to1],[from2,to2],[from3,to3]]).
        sc: StructureControl instance.
    """
    sc.cut_block(
        op["zone"],
        op["abcxyz_flag"],
        op["block_border"],
    )


def do_cutsphere(op, sc):
    """Remove atoms inside or outside a sphere.

    Calls StructureControl.cut_sphere.  The sphere is centered
    on either an atom, an xyz coordinate, or an abc coordinate
    with a specific radius.  The zone controls whether atoms
    inside or outside the sphere are removed.

    Args:
        op: Operation dict with keys 'zone' (0=in,1=out),
            'abcxyz_flag' (0=abc,1=xyz), 'sphere_rad',
            'target_atom', and 'sphere_loc' ([x,y,z]).
        sc: StructureControl instance.
    """
    sc.cut_sphere(
        op["zone"],
        op["abcxyz_flag"],
        op["sphere_rad"],
        op["target_atom"],
        op["sphere_loc"],
    )


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    print("\n\nmod_struct.py script executing.\n")

    # Parse command line and load defaults.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(f"\nStarting environment setup at........{ts}.")
    settings = ScriptSettings()

    # Import and create the StructureControl instance.
    from structure_control import StructureControl
    sc = StructureControl()

    # Set the lattice vector orientations according to the
    # command line.
    sc.set_abc_xyz_assignment_order(
        settings.abc_order, settings.xyz_order
    )

    # Read the input file.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(
        f"\nReading input file at................{ts}."
    )
    sc.read_input_file(settings.in_file, use_file_species=True)

    # Execute all queued operations in order.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(
        f"\nApplying operations at...............{ts}."
    )
    execute_operations(settings, sc)

    # Write the output file.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(
        f"\nWriting output file at...............{ts}."
    )
    sc.print_olcao(
        filename=settings.out_file,
        title=settings.title,
        style="cart",
    )

    # Done.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(
        f"\nmod_struct.py script complete at......{ts}.\n"
    )


if __name__ == "__main__":
    # Everything before this point was a subroutine definition or
    # a request to import information from external modules.  Only
    # now do we actually start running the program.  The purpose
    # of this is to allow another python program to import *this*
    # script and call its functions internally.
    main()
