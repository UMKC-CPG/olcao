#!/usr/bin/env python3

"""porosity.py -- Calculate porosity of a solid-state material.

PROGRAM: porosity
Last Modified: (original Perl: undated; Python port: 2026)

PURPOSE:
    Calculate the porosity of a solid-state material by
    sampling the unit cell on a regular mesh and determining
    which mesh points fall inside (solid) or outside (void /
    pore space) the atomic radii.  The porosity is defined as
    the ratio of pore volume to total cell volume.

USAGE:
    porosity.py  -m  numA  numB  numC
                [-l  limitDist]
                [-i  inFile]
                [-sf scaleFactor]
                [-r  resolution]
                [-help]

OPTIONS:
    -m numA numB numC
        Define the number of sampling points along each a, b,
        c lattice direction.  Mutually exclusive with -r.

    -r resolution
        Define the target resolution (Angstroms) for the
        sampling mesh.  The actual number of mesh points along
        each axis is computed from this resolution and the cell
        dimensions, then the resolution is refined so that an
        integer number of points fits evenly.  Mutually
        exclusive with -m.

    -l limitDist
        Radius (Angstroms) of the sphere within which atoms
        can contribute to the pore map.  If not given, the
        default value of 4 A is used.

    -i inFile
        Name of the input file.  If not specified, the default
        name "olcao.skl" will be used.

    -sf scaleFactor
        Multiplicative factor applied to the covalent radii
        from the element database.  The scaled radius defines
        the pore boundary: mesh points within this radius of
        any atom centre are considered solid (not void).

    -help
        Print this friendly message and exit.

HOW IT WORKS:
    1. Read the structure (lattice vectors, atom positions)
       from an OLCAO skeleton file via StructureControl.
    2. Either compute the per-axis resolution from the given
       mesh counts (-m), or compute the mesh counts from the
       given target resolution (-r).
    3. Scale the covalent radii by the given factor and set the
       interaction distance limit.
    4. Set up the 3-D mesh of sampling points.
    5. For every atom, mark all mesh points that fall within
       its scaled covalent radius as "solid" (not pore).
    6. Count the remaining unmarked (void) mesh points.
    7. Compute porosity = pore_volume / cell_volume.
"""

import argparse as ap
import math
import os
import sys
from datetime import datetime

from element_data import ElementData
from structure_control import StructureControl


# ============================================================
#  ScriptSettings -- command-line / rc-file parameters
# ============================================================
class ScriptSettings():
    """The instance variables of this object are the user
    settings that control the program.  The variable values are
    pulled from a list that is created within a resource control
    file and that are then reconciled with command-line
    parameters."""

    def __init__(self):
        """Define default values for the program parameters by
        pulling them from the resource control file in the
        default location: $OLCAO_RC/porosityrc.py or from the
        current working directory if a local copy of
        porosityrc.py is present."""

        # Read default variables from the resource control
        #   file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit(
                "Error: $OLCAO_RC is not set. "
                "See instructions."
            )
        sys.path.insert(1, rc_dir)
        from porosityrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command-line arguments with the rc
        #   file defaults.
        self.reconcile(args)

        # Record the command line that was used.
        self.record_clp()

    def assign_rc_defaults(self, rc):
        """Assign all default values from the rc dictionary."""

        # Input file name.
        self.in_file = rc["in_file"]

        # Mesh / resolution control.
        self.num_mesh_points = list(rc["num_mesh_points"])
        self.resolution = rc["resolution"]

        # Physical parameters.
        self.limit_dist = rc["limit_dist"]
        self.scale_factor = rc["scale_factor"]

        # Internal flags for which mode was selected.
        #   do_mesh: user gave explicit mesh counts (-m).
        #   do_res:  user gave a target resolution (-r).
        self.do_mesh = False
        self.do_res = False

    def parse_command_line(self):
        """Build the argparse parser and return parsed args."""

        prog_name = "porosity.py"

        description_text = """
Calculate the porosity of a solid-state material by sampling
the unit cell on a regular mesh.  Mesh points that fall within
a scaled covalent radius of any atom are considered solid; the
remaining points are pore (void) space.

Either -m (explicit mesh) or -r (target resolution) must be
given.  They are mutually exclusive.
"""

        epilog_text = """
Defaults are given in ./porosityrc.py or $OLCAO_RC/porosityrc.py.
"""

        parser = ap.ArgumentParser(
            prog=prog_name,
            formatter_class=ap.RawDescriptionHelpFormatter,
            description=description_text,
            epilog=epilog_text,
        )

        self.add_parser_arguments(parser)

        return parser.parse_args()

    def add_parser_arguments(self, parser):
        """Add all command-line arguments to the parser."""

        # ---- Mesh definition (mutually exclusive) ----
        mesh_group = parser.add_mutually_exclusive_group()

        mesh_group.add_argument(
            '-m', dest='num_mesh_points',
            nargs=3, type=int,
            default=None,
            help=(
                'Number of sampling points along the '
                'a, b, c directions.'
            ),
        )
        mesh_group.add_argument(
            '-r', dest='resolution', type=float,
            default=None,
            help=(
                'Target resolution (Angstroms) for '
                'the sampling mesh.'
            ),
        )

        # ---- Physical parameters ----
        parser.add_argument(
            '-l', dest='limit_dist', type=float,
            default=self.limit_dist,
            help=(
                'Radius of the inclusion sphere '
                '(Angstroms). '
                f'Default: {self.limit_dist}'
            ),
        )
        parser.add_argument(
            '-sf', dest='scale_factor', type=float,
            default=self.scale_factor,
            help=(
                'Scale factor for covalent radii. '
                f'Default: {self.scale_factor}'
            ),
        )

        # ---- Input file ----
        parser.add_argument(
            '-i', dest='in_file', type=str,
            default=self.in_file,
            help=(
                'Input structure file. '
                f'Default: {self.in_file}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command-line arguments with rc defaults.

        Determines whether the user selected mesh mode (-m) or
        resolution mode (-r) and sets the appropriate flags.
        If neither was given, the program will abort later with
        a clear message.
        """

        self.in_file = args.in_file
        self.limit_dist = args.limit_dist
        self.scale_factor = args.scale_factor

        if args.num_mesh_points is not None:
            self.do_mesh = True
            self.do_res = False
            self.num_mesh_points = (
                [None] + list(args.num_mesh_points)
            )
        elif args.resolution is not None:
            self.do_res = True
            self.do_mesh = False
            self.resolution = args.resolution
        # else: neither given -- will be caught in main().

        # Compute the total number of space points if mesh
        #   mode was selected.  (For resolution mode this is
        #   computed after refinement.)
        if self.do_mesh:
            self.num_space_points = (
                self.num_mesh_points[1]
                * self.num_mesh_points[2]
                * self.num_mesh_points[3]
            )
        else:
            self.num_space_points = 0

    def record_clp(self):
        """Append the command line used to the 'command'
        file."""
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


# ============================================================
#  Core processing functions
# ============================================================

def compute_resolution(settings, sc):
    """Compute the per-axis resolution from the mesh counts.

    When the user specifies the number of mesh points along
    each axis with -m, the per-axis resolution (Angstroms per
    mesh step) is derived by dividing the cell magnitude along
    each axis by the number of mesh points on that axis.

    Parameters
    ----------
    settings : ScriptSettings
        The program settings containing num_mesh_points.
    sc : StructureControl
        The structure data (provides sc.mag[1..3]).

    Returns
    -------
    axis_resolution : list of float
        1-indexed list [None, res_a, res_b, res_c] giving
        the resolution along each lattice direction.
    """
    nmp = settings.num_mesh_points

    # axis_resolution is 1-indexed: [None, res_a, res_b, res_c]
    #   as expected by StructureControl.compute_pore_map().
    axis_resolution = [None, 0.0, 0.0, 0.0]
    for axis in range(1, 4):
        axis_resolution[axis] = (
            sc.mag[axis] / nmp[axis]
        )
    return axis_resolution


def compute_mesh_and_refine_res(settings, sc):
    """Compute mesh counts from a target resolution and
    refine.

    When the user specifies a target resolution with -r,
    this function:
    1. Divides each cell magnitude by the target resolution
       to get a tentative number of points per axis.
    2. Rounds that number to the nearest integer.
    3. Recomputes the actual per-axis resolution so that an
       integer number of mesh steps fits exactly into the
       cell.

    Parameters
    ----------
    settings : ScriptSettings
        The program settings (resolution is read; mesh
        counts and num_space_points are written back).
    sc : StructureControl
        The structure data (provides sc.mag[1..3]).

    Returns
    -------
    axis_resolution : list of float
        1-indexed list [None, res_a, res_b, res_c] giving
        the refined resolution along each lattice direction.
    """
    axis_resolution = [None, 0.0, 0.0, 0.0]

    # Refine the resolution for each axis.
    for axis in range(1, 4):
        points_per_cell = (
            sc.mag[axis] / settings.resolution
        )

        # Round to the nearest integer.
        if (points_per_cell - int(points_per_cell)) > 0.5:
            points_per_cell = int(points_per_cell) + 1
        else:
            points_per_cell = int(points_per_cell)

        axis_resolution[axis] = (
            sc.mag[axis] / points_per_cell
        )

    # Compute the mesh with the new resolution.
    for axis in range(1, 4):
        settings.num_mesh_points[axis] = int(
            sc.mag[axis] / axis_resolution[axis]
        )

    settings.num_space_points = (
        settings.num_mesh_points[1]
        * settings.num_mesh_points[2]
        * settings.num_mesh_points[3]
    )

    return axis_resolution


def calc_void_points(settings, sc):
    """Count the number of mesh points that are in pore
    (void) space.

    After compute_pore_map() has been called on the
    StructureControl object, ext_pore_map[a][b][c] holds
    a value for every mesh point: zero means the point is
    in a pore (no atom is close enough), non-zero means
    the point is within an atom's scaled covalent radius
    (solid).

    This function loops over every mesh point and counts
    the number of pore points (value == 0).

    Parameters
    ----------
    settings : ScriptSettings
        Program settings (provides num_mesh_points).
    sc : StructureControl
        Structure data (provides ext_pore_map).

    Returns
    -------
    pore_point_count : int
        The total number of mesh points that are in pore
        space.
    """
    nmp = settings.num_mesh_points

    # Print a progress note.
    print("Counting number of pore points.")

    # Begin computing if each mesh point is in or out of a
    #   pore.  ext_pore_map[a][b][c] == 0 means pore (void),
    #   non-zero means solid (within an atom's radius).
    pore_point_count = 0
    for a_pt in range(1, nmp[1] + 1):
        for b_pt in range(1, nmp[2] + 1):
            for c_pt in range(1, nmp[3] + 1):
                if sc.ext_pore_map[a_pt][b_pt][c_pt] \
                        == 0:
                    pore_point_count += 1

    # Print porosity information.
    print(
        f"number of pore points = "
        f"{pore_point_count}"
    )

    return pore_point_count


def compute_porosity(settings, sc, pore_point_count):
    """Compute and print the porosity of the system.

    Porosity is defined as:

        porosity = pore_volume / cell_volume

    where:
    - cell_volume  = mag_a * mag_b * mag_c  (in A^3,
      using the magnitudes of the lattice vectors).
    - vol_per_point = cell_volume / num_space_points
      (the volume represented by each mesh point).
    - pore_volume  = vol_per_point * pore_point_count
      (the total volume of all void mesh points).

    Parameters
    ----------
    settings : ScriptSettings
        Program settings (provides num_space_points).
    sc : StructureControl
        Structure data (provides sc.mag[1..3]).
    pore_point_count : int
        Number of void mesh points from calc_void_points().
    """
    # Calculate cell volume.
    cell_volume = sc.mag[1] * sc.mag[2] * sc.mag[3]

    # Calculate volume per point.
    vol_per_point = cell_volume / settings.num_space_points

    # Calculate pore volume.
    pore_volume = vol_per_point * pore_point_count

    # Calculate porosity.
    porosity = pore_volume / cell_volume

    print(f"total pore volume = {pore_volume} A^3")
    print(f"porosity = {porosity}")


# ============================================================
#  Main program
# ============================================================

def main():
    """Main execution flow for porosity.py.

    The overall flow is:
    1. Initialise environment and parse command line.
    2. Set up the element database.
    3. Read the input structure file.
    4. Compute the mesh resolution (from -m) or the mesh
       counts (from -r).
    5. Apply the scale factor to covalent radii.
    6. Set the interaction distance limit.
    7. Set up the 3-D mesh of sampling points.
    8. Compute the pore map (mark solid vs. void).
    9. Count the void points.
    10. Compute and print the porosity.
    """

    print("\n\nScript Executing.\n")

    # ---- Step 1: Initialise and parse command line ----
    settings = ScriptSettings()

    # ---- Step 2: Set up the element database ----
    #   Initialise the element data from the database.
    #   This provides covalent radii, element names, etc.
    ed = ElementData()
    ed.init_element_data()

    # ---- Step 3: Read the input file ----
    sc = StructureControl()
    sc.read_input_file(settings.in_file, 0)

    # Get the cell magnitudes.  sc.mag is a 1-indexed list:
    #   sc.mag[1] = |a|, sc.mag[2] = |b|, sc.mag[3] = |c|.

    # ---- Step 4: Compute resolution or mesh ----
    if settings.do_mesh:
        axis_resolution = compute_resolution(settings, sc)
    elif settings.do_res:
        axis_resolution = compute_mesh_and_refine_res(
            settings, sc
        )
    else:
        sys.exit("Use the -r or -m option!")

    # ---- Step 5: Apply scale factor to database ----
    #   Scale the covalent radii by the user-specified
    #   factor.  The scaled radii define the pore boundary.
    ed.apply_bond_factor(settings.scale_factor)

    # ---- Step 6: Set the limit for atom interaction ----
    #   (Radius of the inclusion sphere.)
    sc.set_limit_dist(settings.limit_dist)

    # ---- Step 7: Set up the sampling points ----
    nmp = settings.num_mesh_points
    sc.set_xyz_mesh_points(nmp[1], nmp[2], nmp[3])

    # ---- Step 8: Compute the pore map ----
    #   Visit every atom and mark all mesh points near each
    #   atom as not void space.
    sc.compute_pore_map(axis_resolution, nmp)

    # ---- Step 9: Count void points ----
    pore_point_count = calc_void_points(settings, sc)

    # ---- Step 10: Compute and print porosity ----
    compute_porosity(settings, sc, pore_point_count)

    print("\nScript complete.\n")


if __name__ == '__main__':
    # Everything before this point was a class/function
    #   definition or an import.  Only now do we actually
    #   start running the program.  This structure allows
    #   another Python program to import this script and
    #   call its functions internally.
    main()
