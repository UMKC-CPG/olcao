#!/usr/bin/env python3

"""makeBOND.py -- Parse OLCAO bond order output and produce
visualisation-ready files.

Program: makeBOND
Last Modified:  (original Perl: March 30, 2012; Python port: 2026)
Purpose:  To parse the bond order output file and create a number
          of other files that can be easily plotted by either
          openDX, Paraview/VisIT, Origin, or similar tools.

It is assumed that the calculated bond order and Q* data will be
   in the gs_bond-mb.raw file.  This default can be overridden
   with the -data option.

It is also assumed that the atomic positions will be in the
   structure.dat file.  This default can be overridden with the
   -pos option.  In any case the format of the position file must
   be that of the current OLCAO output file.  This has the
   positions sorted by types.

The -bc option is used to define a file that will filter bonds
   from the raw output file based on bond length and bond order
   for each atomic element.  The format of the file is a series of
   consecutive lines of five space-separated columns of data.  The
   first column is a tag of either "BO" or "BL" which indicates
   that the line represents a restriction for the bond order or
   the bond length.  The second and third values are the element
   names from the periodic table of the elements.  The fourth and
   fifth values are the min and max ranges for the bond order or
   the bond length for that element pair.  You can provide as many
   lines as you wish and they are applied IN ADDITION to the
   command line options.

The -gc option is used to define a file that will group the
   output in a specific way that often corresponds to the
   definitions in a HIN file.  The format of the file is a series
   of space-separated names that correspond to names that are
   available in the bond order/Q* raw data file.  Each name that
   is given is used to construct the tag for all atoms and bonds.
   For example if the filter names are:
   "ELEMENT_NAME SPECIES_ID" each atom will then have a tag
   created with its element name and species id number.  These
   tags are used to separate the atoms into separate files.

The -maxbo and -minbo options take one argument each that will
   make it so any printout will not contain any bonds with a bond
   order greater/less than the given argument.  This is often used
   in place of a bondControlFile, but it applies to all atoms
   independent of element so it is only useful if you want to
   apply conditions to all atoms.

The -maxbl and -minbl options take one argument each that will
   make it so any printout will not contain any bonds with a
   length greater/less than the given argument.  This is often
   used in place of a bondControlFile, but it applies to all
   atoms independent of element so it is only useful if you want
   to apply conditions to all atoms.

The default output style is to make a scatter plot using the
   species of the system as the grouping criterion for both the
   Q* and bond order.

The -scatter option will create multiple files that can be
   imported by an associated Origin script and made into a scatter
   plot.  The files will contain specific groups of atoms based on
   either the element and species id or the group control file.

The -profile option will produce a simple data file that can be
   plotted as a curve indicating the bond profile along a specific
   axis.  The bond profile at a given point along a profile_axis
   is simply the sum of the bond orders for all the bonds that
   cross that point divided by the number of bonds.  The
   acceptable values for profile_axis are "a", "b", or "c".  At
   present the -profile option will only make one file for the
   whole system and is not separated according to any criteria.

The -model option will produce multiple openDX format files that
   describe a model of the system containing specific sets of
   atoms and bonds that can be drawn using colour coding to
   illustrate their relative bond order and charge transfer.  The
   openDX files contain info to draw the atoms as spheres with
   the colour statically coded according to the values in the
   database.  The size of the sphere is given as the covalent
   radius by default, but this can be adjusted by the -rf
   suboption.  The radius_factor that is given for the -rf
   suboption will multiply against the default radius.  The -grey
   suboption will make the spheres appear in grey-scale instead
   of colour.  The -bo3c suboption will alter the structure to
   show three-centre bonds where they occur and to show
   two-centre bonds in the places they occur.  The most important
   suboptions for -model are the -atom and -bond suboptions.  If
   the -atom suboption is given, then the openDX models will
   contain specific sets of atoms and the bonds in the model must
   contain at least one atom from that set.  If the -bond
   suboption is given, then the openDX models will contain
   specific sets of bonds and the atoms in the model must be
   present as one of the ends of the bond.

The -vtk option will do the same as the preceding option, but
   will produce only one file and in the VTK legacy format.  The
   colour coding has yet to be implemented.

The -mesh option will produce a 3D mesh of data points for
   display in openDX.  The values at each mesh point are computed
   by averaging the values of select data (e.g. BO values at the
   bond midpoint) from nearby data points.  The suboption -atom
   is used to request that the properties of a specific atom are
   subtracted from the mesh point values to create a difference
   mesh.  The -step option is used to define the mesh grid.  The
   -fwhm option is used to specify how to weight the nearby
   values.  This is the full width at half the maximum height of
   the Gaussian function that is used for applying the weights.
   The -limit suboption specifies the maximum distance of
   separation between a bond point and a mesh point.

Other notes:

There are four important distinctions that must be made
   concerning the atoms in the calculation.  First, there is a
   total number of atoms in the system.  These atoms are divided
   into three groups.  The first is type, the second is species,
   and the third is element.

Element refers to the elemental type.  An atom is either a
   carbon or a silicon or a nitrogen etc.  The set of elements is
   the smallest set.  It may contain only two elements (Si and N
   for g-Si3N4).  If we make a substitution of a carbon atom for
   a silicon in the tetrahedral configuration, then the set of
   elements increases by one.  It now contains (Si, N, and C).
   If we make a second substitution in the same system of another
   carbon for a silicon in the octahedral configuration, the
   elemental set stays the same.  This scheme is pretty obviously
   inflexible and set in stone.  The element set only deals with
   the elements in the system, and not the system relations
   itself.

Species refers to the atoms that share the same global symmetry
   and/or a similar local configuration prior to any system
   modification.  So, in crystalline g-Si3N4 there are two
   species of Si.  One in a tetrahedral configuration and one in
   an octahedral configuration.  There is still only one species
   of N.  For complex systems with little or no symmetry the
   species are defined by their local configurational
   similarities.  Those atoms that have similar neighbour atoms
   at similar distances will be treated as the same species.

Type refers to the atoms that should be considered equivalent by
   the calculation.  Often the assignment of types is the exact
   same as the assignment of species.  In some cases it is desired
   to distinguish atoms of the same species and make them have
   different potential coefficients (i.e. a different potential
   type or different type).  This is true for XANES calculations
   or cases where we just assign all atoms in a region to have
   different types because they cannot be easily grouped in any
   other way.

USAGE:
  makeBOND.py [-data DATA_FILE] [-pos POS_FILE]
              [-bc BOND_CONTROL_FILE] [-gc GROUP_CONTROL_FILE]
              [-maxbo MAX_BO] [-minbo MIN_BO]
              [-maxbl MAX_BL] [-minbl MIN_BL]
              [-scatter]
              [-dipole]
              [-profile PROFILE_AXIS]
              [-model [-atom | -bond] [-rf RADIUS_FACTOR]
                      [-grey] [-bo3c]]
              [-vtk]
              [-mesh [-atom COMP_ATOM] [-steps A B C]
                     [-fwhm FWHM] [-limit LIMIT_DIST]]
              [-help]
"""

import argparse as ap
import math
import os
import sys
from datetime import datetime

from element_data import ElementData


# ============================================================
#  Physical constant
# ============================================================
BOHR_RAD = 0.5291772180  # Bohr radius in Angstroms


# ============================================================
#  ScriptSettings -- command-line / rc-file parameters
# ============================================================
class ScriptSettings():
    """The instance variables of this object are the user settings
    that control the program.  The variable values are pulled from
    a list that is created within a resource control file and that
    are then reconciled with command-line parameters."""

    def __init__(self):
        """Define default values for the program parameters by
        pulling them from the resource control file in the default
        location: $OLCAO_RC/makeBONDrc.py or from the current
        working directory if a local copy of makeBONDrc.py is
        present."""

        # Read default variables from the resource control file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit(
                "Error: $OLCAO_RC is not set. "
                "See instructions."
            )
        sys.path.insert(1, rc_dir)
        from makeBONDrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile CL arguments with rc-file defaults.
        self.reconcile(args)

        # Record the command line that was used.
        self.record_clp()

    def assign_rc_defaults(self, rc):
        """Assign all default values from the rc dictionary."""

        # Input file names.
        self.data_file = rc["data_file"]
        self.pos_file = rc["pos_file"]

        # Control files.
        self.bond_control_file = rc["bond_control_file"]
        self.group_control_file = rc["group_control_file"]
        self.bond_profile_file = rc["bond_profile_file"]

        # Output mode flags.
        self.scatter_plot = rc["scatter_plot"]
        self.dipole = rc["dipole"]
        self.system_model_dx = rc["system_model_dx"]
        self.bo3c = rc["bo3c"]
        self.system_model_vtk = rc["system_model_vtk"]
        self.mesh_dx = rc["mesh_dx"]
        self.atom_based_dx = rc["atom_based_dx"]
        self.bond_based_dx = rc["bond_based_dx"]
        self.bond_profile_plot = rc["bond_profile_plot"]
        self.profile_axis = rc["profile_axis"]

        # Bond order / bond length thresholds.
        self.min_bo = rc["min_bo"]
        self.max_bo = rc["max_bo"]
        self.min_bl = rc["min_bl"]
        self.max_bl = rc["max_bl"]

        # Model tuning parameters.
        self.radius_factor = rc["radius_factor"]
        self.grey_scale = rc["grey_scale"]

        # Mesh parameters.
        self.comp_atom = rc["comp_atom"]
        self.fwhm = rc["fwhm"]
        self.num_mesh_points = list(rc["num_mesh_points"])
        self.limit_dist = rc["limit_dist"]

    def parse_command_line(self):
        """Build the argparse parser and return parsed args."""

        prog_name = "makeBOND.py"

        description_text = """
Parse the bond order output file produced by the OLCAO bond
program and create files suitable for plotting by openDX,
Paraview/VisIT, Origin, or similar tools.

The input data file (default: gs_bond-mb.raw) contains raw
bond order and Q* data.  The position file (default:
structure.dat) provides lattice vectors and atom positions.

Output modes (mutually exclusive -- last one wins if
multiple are given; scatter is the default):
  -scatter   Scatter-plot files grouped by atom/bond tag.
  -model     openDX structural model with coloured bonds.
  -vtk       VTK Legacy structural model.
  -mesh      3-D openDX mesh of bond order values.
  -profile   Bond-profile curve along a lattice axis.
"""

        epilog_text = """
Defaults are given in ./makeBONDrc.py or $OLCAO_RC/makeBONDrc.py.
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

        # ---- Input files ----
        parser.add_argument(
            '-data', dest='data_file', type=str,
            default=self.data_file,
            help=(
                'Raw bond order data file. '
                f'Default: {self.data_file}'
            ),
        )
        parser.add_argument(
            '-pos', dest='pos_file', type=str,
            default=self.pos_file,
            help=(
                'Structure file with lattice and positions. '
                f'Default: {self.pos_file}'
            ),
        )

        # ---- Control files ----
        parser.add_argument(
            '-bc', dest='bond_control_file', type=str,
            default=self.bond_control_file,
            help=(
                'Bond control file for element-pair BO/BL '
                'filtering.'
            ),
        )
        parser.add_argument(
            '-gc', dest='group_control_file', type=str,
            default=self.group_control_file,
            help=(
                'Group control file defining output grouping '
                'tags.'
            ),
        )

        # ---- Output modes ----
        parser.add_argument(
            '-scatter', dest='scatter_plot',
            action='store_true',
            default=self.scatter_plot,
            help='Produce scatter-plot output files.',
        )
        parser.add_argument(
            '-dipole', dest='dipole',
            action='store_true',
            default=self.dipole,
            help='Compute and print dipole information.',
        )
        parser.add_argument(
            '-model', dest='system_model_dx',
            action='store_true',
            default=self.system_model_dx,
            help='Produce openDX structural model files.',
        )
        parser.add_argument(
            '-vtk', dest='system_model_vtk',
            action='store_true',
            default=self.system_model_vtk,
            help='Produce VTK Legacy structural model.',
        )
        parser.add_argument(
            '-profile', dest='profile_axis', type=str,
            default=self.profile_axis,
            help=(
                'Produce bond profile along the given '
                'axis (a, b, or c).'
            ),
        )
        parser.add_argument(
            '-mesh', dest='mesh_dx',
            action='store_true',
            default=self.mesh_dx,
            help='Produce 3-D openDX bond order mesh.',
        )

        # ---- Model suboptions ----
        parser.add_argument(
            '-atom', dest='atom_based_dx',
            action='store_true',
            default=self.atom_based_dx,
            help=(
                'Make openDX/mesh model atom-based '
                '(default).'
            ),
        )
        parser.add_argument(
            '-bond', dest='bond_based_dx',
            action='store_true',
            default=self.bond_based_dx,
            help='Make openDX model bond-based.',
        )
        parser.add_argument(
            '-rf', dest='radius_factor', type=float,
            default=self.radius_factor,
            help=(
                'Multiplicative factor for atom sphere '
                f'radii. Default: {self.radius_factor}'
            ),
        )
        parser.add_argument(
            '-grey', dest='grey_scale',
            action='store_true',
            default=self.grey_scale,
            help='Use grey-scale atom colours.',
        )
        parser.add_argument(
            '-bo3c', dest='bo3c',
            action='store_true',
            default=self.bo3c,
            help='Show three-centre bonds in the model.',
        )

        # ---- BO/BL thresholds ----
        parser.add_argument(
            '-maxbo', dest='max_bo', type=float,
            default=self.max_bo,
            help=(
                'Maximum bond order to include. '
                f'Default: {self.max_bo}'
            ),
        )
        parser.add_argument(
            '-minbo', dest='min_bo', type=float,
            default=self.min_bo,
            help=(
                'Minimum bond order to include. '
                f'Default: {self.min_bo}'
            ),
        )
        parser.add_argument(
            '-maxbl', dest='max_bl', type=float,
            default=self.max_bl,
            help=(
                'Maximum bond length to include. '
                f'Default: {self.max_bl}'
            ),
        )
        parser.add_argument(
            '-minbl', dest='min_bl', type=float,
            default=self.min_bl,
            help=(
                'Minimum bond length to include. '
                f'Default: {self.min_bl}'
            ),
        )

        # ---- Mesh suboptions ----
        parser.add_argument(
            '-comp_atom', dest='comp_atom', type=int,
            default=self.comp_atom,
            help=(
                'Atom index for difference mesh. '
                f'Default: {self.comp_atom}'
            ),
        )
        parser.add_argument(
            '-steps', dest='num_mesh_points',
            nargs=3, type=int,
            default=self.num_mesh_points,
            help=(
                'Number of mesh points along a, b, c. '
                f'Default: {self.num_mesh_points}'
            ),
        )
        parser.add_argument(
            '-fwhm', dest='fwhm', type=float,
            default=self.fwhm,
            help=(
                'Gaussian FWHM for mesh weighting. '
                f'Default: {self.fwhm}'
            ),
        )
        parser.add_argument(
            '-limit', dest='limit_dist', type=float,
            default=self.limit_dist,
            help=(
                'Max bond-to-mesh-point distance. '
                f'Default: {self.limit_dist}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command-line arguments with rc defaults.

        Command-line values override rc defaults.  Additionally,
        some mutual-exclusivity logic is applied:
        - If -bond is given, atom_based_dx is turned off.
        - If -profile was given a value, bond_profile_plot is set.
        - If no output mode was explicitly selected, scatter_plot
          becomes the default.
        """

        self.data_file = args.data_file
        self.pos_file = args.pos_file
        self.bond_control_file = args.bond_control_file
        self.group_control_file = args.group_control_file

        self.scatter_plot = args.scatter_plot
        self.dipole = args.dipole
        self.system_model_dx = args.system_model_dx
        self.system_model_vtk = args.system_model_vtk
        self.mesh_dx = args.mesh_dx

        self.atom_based_dx = args.atom_based_dx
        self.bond_based_dx = args.bond_based_dx
        self.radius_factor = args.radius_factor
        self.grey_scale = args.grey_scale
        self.bo3c = args.bo3c

        self.min_bo = args.min_bo
        self.max_bo = args.max_bo
        self.min_bl = args.min_bl
        self.max_bl = args.max_bl

        self.comp_atom = args.comp_atom
        self.fwhm = args.fwhm
        self.num_mesh_points = list(args.num_mesh_points)
        self.limit_dist = args.limit_dist

        # Handle -profile: if a value was given, enable the
        #   bond profile plot mode.
        if args.profile_axis:
            self.bond_profile_plot = True
            self.profile_axis = args.profile_axis
        else:
            self.bond_profile_plot = False
            self.profile_axis = ""

        # If -bond is given, turn off atom-based mode.
        if self.bond_based_dx:
            self.atom_based_dx = False

        # Assign a default output of scatter plot if no
        #   output mode was given.
        if (not self.system_model_dx
                and not self.bond_profile_plot
                and not self.mesh_dx
                and not self.system_model_vtk):
            self.scatter_plot = True

        # Compute derived Gaussian parameters from fwhm.
        sigma = 1.0 / (
            2.0 * math.sqrt(2.0 * math.log(2.0))
            / self.fwhm
        )
        self.alpha = 1.0 / (2.0 * sigma * sigma)

    def record_clp(self):
        """Append the command line used to the 'command' file."""
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
#  BondData -- all data arrays for the bond analysis
# ============================================================
class BondData():
    """Container for all data structures used during the bond
    order analysis.

    This class holds the atom positions, bond information, tags,
    element data, extended-cell data, and DX/VTK output lists.
    It is populated progressively by the functions below.
    """

    def __init__(self, settings):
        """Initialise the BondData object with the given
        ScriptSettings and prepare the element database."""

        self.settings = settings

        # Initialise the element data from the database.
        self.elem_data = ElementData()
        self.elem_data.init_element_data()

        # ---- Core atom/bond arrays (1-indexed) ----
        #   Index 0 is a placeholder and must not be used.

        # Number of atoms in the system.
        self.num_atoms = 0

        # Total number of bonds in the system (after
        #   filtering by the control file requirements).
        self.num_system_bonds = 0

        # Element name for each atom.
        #   element_name[i] = string, e.g. "SI", "N".
        self.element_name = [None]

        # Atomic number Z for each atom.
        self.element_z = [None]

        # Effective charge (Q*) for each atom.
        self.atom_q_star = [None]

        # Filter names and values read from the data file
        #   for each atom.
        #   num_filters[i] = int count of filters for atom i.
        #   filter_name[i][k] = string filter name.
        #   filter_value[i][k] = string filter value.
        self.num_filters = [None]
        self.filter_name = [None]
        self.filter_value = [None]

        # Bond information indexed by system bond number.
        #   bond_atom[b] = [atom1, atom2] (the two atoms).
        #   bond_length[b] = float bond length.
        #   bond_order[b] = float bond order value.
        self.bond_atom = [None]
        self.bond_length = [None]
        self.bond_order = [None]

        # Bond information indexed by atom.
        #   num_bonds_by_atom[i] = int count.
        #   bond_atom_by_atom[i][j] = other atom index.
        #   bond_length_by_atom[i][j] = float.
        #   bond_order_by_atom[i][j] = float.
        self.num_bonds_by_atom = [None]
        self.bond_atom_by_atom = [None]
        self.bond_length_by_atom = [None]
        self.bond_order_by_atom = [None]

        # Atom and bond tags (strings built from filters).
        self.atom_tag = [None]
        self.bond_tag = [None]
        self.unique_atom_tags = []
        self.unique_bond_tags = []

        # Atom intrinsic properties from the element database.
        #   atom_radius[i] = covalent radius * radius_factor.
        #   atom_color[i] = DX colour index.
        #   atom_charge_transfer[i] = Q* minus valence charge.
        self.atom_radius = [None]
        self.atom_color = [None]
        self.atom_charge_transfer = [None]

        # ---- Lattice and position data ----

        # Real-space lattice vectors (3x3, 1-indexed).
        #   real_lattice[axis][component]
        #   axis: 1=a, 2=b, 3=c;  component: 1=x, 2=y, 3=z.
        self.real_lattice = [
            None,
            [None, 0.0, 0.0, 0.0],
            [None, 0.0, 0.0, 0.0],
            [None, 0.0, 0.0, 0.0],
        ]

        # Lattice vector components (a, b, c in x, y, z).
        self.a = [None, 0.0, 0.0, 0.0]
        self.b = [None, 0.0, 0.0, 0.0]
        self.c = [None, 0.0, 0.0, 0.0]

        # Magnitudes of lattice vectors (Angstroms).
        self.mag_a = 0.0
        self.mag_b = 0.0
        self.mag_c = 0.0

        # Cell determinant and inverse-lattice factors for
        #   converting x,y,z to a,b,c coordinates.
        self.cell_determinant = 0.0
        self.x_factor = [None, 0.0, 0.0, 0.0]
        self.y_factor = [None, 0.0, 0.0, 0.0]
        self.z_factor = [None, 0.0, 0.0, 0.0]

        # Atom positions in Cartesian (Angstroms).
        self.x_pos = [None]
        self.y_pos = [None]
        self.z_pos = [None]

        # Atom positions in fractional*magnitude coords.
        self.a_pos = [None]
        self.b_pos = [None]
        self.c_pos = [None]

        # Bond midpoint positions in a, b, c coords.
        self.a_bond_pos = [None]
        self.b_bond_pos = [None]
        self.c_bond_pos = [None]

        # Bond midpoint positions in x, y, z coords.
        self.x_bond_pos = [None]
        self.y_bond_pos = [None]
        self.z_bond_pos = [None]

        # Per-atom statistics.
        self.atom_min_bo = [None]
        self.atom_max_bo = [None]
        self.atom_avg_bo = [None]
        self.atom_min_bl = [None]
        self.atom_max_bl = [None]
        self.atom_avg_bl = [None]

        # ---- Extended-cell data ----
        #   Used when bonds cross cell boundaries.

        self.num_bonds_ext = 0
        self.num_atoms_ext = 0

        # Extended positions (Cartesian and a,b,c).
        self.x_pos_ext = [None]
        self.y_pos_ext = [None]
        self.z_pos_ext = [None]
        self.a_pos_ext = [None]
        self.b_pos_ext = [None]
        self.c_pos_ext = [None]

        # Extended bond definitions.
        self.bond_atom_ext = [None]
        self.bond_order_ext = [None]
        self.bond_tag_ext = [None]
        self.atom_tag_ext = [None]

        # ---- DX output lists ----
        self.unique_tag_dx = ""
        self.num_atoms_dx = 0
        self.num_atoms_ext_dx = 0
        self.num_bonds_dx = 0
        self.x_pos_ext_dx = [None]
        self.y_pos_ext_dx = [None]
        self.z_pos_ext_dx = [None]
        self.atom_map_2_dx = [None]
        self.atom_color_dx = [None]
        self.atom_radius_dx = [None]
        self.atom_charge_transfer_dx = [None]
        self.bond_atom_dx = [None]
        self.bond_order_dx = [None]

        # Bond positions for DX (abc and xyz).
        self.abc_bond_pos = [None]
        self.xyz_bond_pos = [None]

        # ---- Mesh data ----
        self.num_scan_points = 0
        self.delta = [
            None,
            [None, 0.0, 0.0, 0.0],
            [None, 0.0, 0.0, 0.0],
            [None, 0.0, 0.0, 0.0],
        ]
        self.scan_points = [None, [None], [None], [None]]
        self.mesh_bo_dx = [None]

        # ---- Bond control limits ----
        self.num_bo_limits = 0
        self.num_bl_limits = 0
        self.bo_limits = [None]
        self.bl_limits = [None]

        # ---- Group control ----
        self.group_name = [None]
        self.num_groups = 0


# ============================================================
#  Core processing functions
# ============================================================

def get_system_structure(bd):
    """Read lattice vectors and atom positions from the
    structure file.

    This function reads the system structure from the position
    file (default: structure.dat).  The key variables defined
    when this is done are:

    1) a[], b[], c[]:  Range = 1-3.
       Meaning = Cell's a, b, c vectors in x, y, z coordinate
                 components (Angstroms).

    2) mag_a, mag_b, mag_c:  Magnitude of the cell's a, b, c
       vectors in Angstroms.

    3) num_atoms:  The number of atoms in the system.

    4) x_pos[], y_pos[], z_pos[]:  Range = 1..num_atoms.
       Meaning = Location of each atom in x, y, z orthogonal
                 coordinates (Angstroms).

    5) a_pos[], b_pos[], c_pos[]:  Range = 1..num_atoms.
       Meaning = Location of each atom in a, b, c coordinates.
    """
    s = bd.settings

    with open(s.pos_file, "r") as pos:
        for line in pos:

            # ---- Read cell vectors ----
            if "CELL_VECTORS" in line:
                # Read the a vector (in a.u.).
                vals = pos.readline().split()
                for i in range(3):
                    bd.a[i + 1] = (
                        float(vals[i]) * BOHR_RAD
                    )
                bd.mag_a = math.sqrt(
                    bd.a[1]**2 + bd.a[2]**2 + bd.a[3]**2
                )

                # Read the b vector (in a.u.).
                vals = pos.readline().split()
                for i in range(3):
                    bd.b[i + 1] = (
                        float(vals[i]) * BOHR_RAD
                    )
                bd.mag_b = math.sqrt(
                    bd.b[1]**2 + bd.b[2]**2 + bd.b[3]**2
                )

                # Read the c vector (in a.u.).
                vals = pos.readline().split()
                for i in range(3):
                    bd.c[i + 1] = (
                        float(vals[i]) * BOHR_RAD
                    )
                bd.mag_c = math.sqrt(
                    bd.c[1]**2 + bd.c[2]**2 + bd.c[3]**2
                )

                # Compute the cell determinant.
                bd.cell_determinant = (
                    bd.c[1] * bd.a[2] * bd.b[3]
                    - bd.c[1] * bd.a[3] * bd.b[2]
                    + bd.a[3] * bd.c[2] * bd.b[1]
                    + bd.c[3] * bd.a[1] * bd.b[2]
                    - bd.a[2] * bd.c[3] * bd.b[1]
                    + bd.c[2] * bd.a[1] * bd.b[3]
                )

                # Compute the factors for converting
                #   x,y,z -> a,b,c coordinates.
                bd.x_factor[1] = (
                    bd.c[3] * bd.b[2]
                    - bd.c[2] * bd.b[3]
                )
                bd.x_factor[2] = (
                    bd.a[3] * bd.c[2]
                    - bd.a[2] * bd.c[3]
                )
                bd.x_factor[3] = (
                    bd.a[2] * bd.b[3]
                    - bd.a[3] * bd.b[2]
                )
                bd.y_factor[1] = (
                    bd.c[1] * bd.b[3]
                    - bd.c[3] * bd.b[1]
                )
                bd.y_factor[2] = (
                    bd.c[3] * bd.a[1]
                    - bd.c[1] * bd.a[3]
                )
                bd.y_factor[3] = (
                    bd.b[1] * bd.a[3]
                    - bd.a[1] * bd.b[3]
                )
                bd.z_factor[1] = (
                    bd.b[1] * bd.c[2]
                    - bd.c[1] * bd.b[2]
                )
                bd.z_factor[2] = (
                    bd.c[1] * bd.a[2]
                    - bd.c[2] * bd.a[1]
                )
                bd.z_factor[3] = (
                    bd.a[1] * bd.b[2]
                    - bd.b[1] * bd.a[2]
                )

                # Store the real lattice matrix.
                bd.real_lattice[1][1] = bd.a[1]
                bd.real_lattice[1][2] = bd.a[2]
                bd.real_lattice[1][3] = bd.a[3]
                bd.real_lattice[2][1] = bd.b[1]
                bd.real_lattice[2][2] = bd.b[2]
                bd.real_lattice[2][3] = bd.b[3]
                bd.real_lattice[3][1] = bd.c[1]
                bd.real_lattice[3][2] = bd.c[2]
                bd.real_lattice[3][3] = bd.c[3]

            # ---- Read atom positions ----
            if "NUM_ATOM_SITES" in line:
                # Get the number of atoms to read.
                vals = pos.readline().split()
                bd.num_atoms = int(vals[0])

                # Position the read marker past the label.
                pos.readline()

                # Read the position input sorted by type.
                for i in range(1, bd.num_atoms + 1):
                    vals = pos.readline().split()

                    # Positions in Angstroms (converted
                    #   from a.u.).
                    x = float(vals[3]) * BOHR_RAD
                    y = float(vals[4]) * BOHR_RAD
                    z = float(vals[5]) * BOHR_RAD
                    bd.x_pos.append(x)
                    bd.y_pos.append(y)
                    bd.z_pos.append(z)

                    # Convert to a, b, c coordinates.
                    a_val = (
                        (bd.x_factor[1] * x
                         + bd.y_factor[1] * y
                         + bd.z_factor[1] * z)
                        / bd.cell_determinant
                    ) * bd.mag_a
                    b_val = (
                        (bd.x_factor[2] * x
                         + bd.y_factor[2] * y
                         + bd.z_factor[2] * z)
                        / bd.cell_determinant
                    ) * bd.mag_b
                    c_val = (
                        (bd.x_factor[3] * x
                         + bd.y_factor[3] * y
                         + bd.z_factor[3] * z)
                        / bd.cell_determinant
                    ) * bd.mag_c
                    bd.a_pos.append(a_val)
                    bd.b_pos.append(b_val)
                    bd.c_pos.append(c_val)


def read_control_files(bd):
    """Read the bond control file and/or group control file.

    This function will read the control file if one is requested.
    There are two types of control files depending on the type
    of output desired.

    The first (bond control) is used to limit the bond order and
    bond length of individual atomic elements.  The file format
    is lines of: TAG ELEM1 ELEM2 MIN MAX, where TAG is "BO" or
    "BL".

    The second (group control) is used to divide up the output
    based on certain group definitions such as those of a HIN
    file.  Each name on a single line corresponds to a filter
    name available in the raw data file.  Default grouping is
    by ELEMENT_NAME and SPECIES_ID.
    """
    s = bd.settings

    # If a bond control file was given then open and parse it.
    if s.bond_control_file:
        # Initialise counters for bond order and bond length
        #   thresholds.
        bd.num_bo_limits = 0
        bd.num_bl_limits = 0

        with open(s.bond_control_file, "r") as ctrl:
            for line in ctrl:
                vals = line.split()
                if not vals:
                    continue

                if vals[0] == "BO":
                    bd.num_bo_limits += 1
                    elem1 = vals[1].upper()
                    elem2 = vals[2].upper()
                    # Make sure element name 1 comes
                    #   before element name 2
                    #   alphabetically.
                    if elem1 > elem2:
                        elem1, elem2 = elem2, elem1
                    bd.bo_limits.append([
                        None,
                        elem1,
                        elem2,
                        float(vals[3]),
                        float(vals[4]),
                    ])

                elif vals[0] == "BL":
                    bd.num_bl_limits += 1
                    elem1 = vals[1].upper()
                    elem2 = vals[2].upper()
                    # Make sure element name 1 comes
                    #   before element name 2
                    #   alphabetically.
                    if elem1 > elem2:
                        elem1, elem2 = elem2, elem1
                    bd.bl_limits.append([
                        None,
                        elem1,
                        elem2,
                        float(vals[3]),
                        float(vals[4]),
                    ])

                else:
                    sys.exit(
                        f"Control File Error: "
                        f"{vals[0]} must be BO or BL"
                    )
    else:
        # Initialise counters for bond order and bond length
        #   thresholds.
        bd.num_bo_limits = 0
        bd.num_bl_limits = 0

    # If a grouping control file was given, then open and
    #   parse it.
    if s.group_control_file:
        with open(s.group_control_file, "r") as ctrl:
            # Obtain the line that defines how to create the
            #   bond and atom tags.
            line = ctrl.readline()
            vals = line.split()
            # Build the 1-indexed group_name list.
            for i, name in enumerate(vals):
                bd.group_name.append(name)
            bd.num_groups = len(vals)
    else:
        bd.group_name.append("ELEMENT_NAME")
        bd.group_name.append("SPECIES_ID")
        bd.num_groups = 2


def read_data(bd):
    """Read raw bond order data from the data file.

    This function reads the raw data produced by the OLCAO bond
    program.  It populates the bond arrays (bond_atom,
    bond_length, bond_order) and per-atom bond arrays
    (bond_atom_by_atom, etc.).

    The data file contains, for each atom:
    - A set of filter name/value pairs (ELEMENT_NAME, SPECIES_ID,
      TYPE_ID, etc.) terminated by ATOM_CHARGE.
    - Orbital charges (read and skipped).
    - A bond count followed by lines of: bonded_atom bond_length
      bond_order.
    - A bond angle count followed by angle data (skipped).

    Each bond is checked against the control file requirements
    via check_requirements().
    """
    s = bd.settings
    ed = bd.elem_data

    # At this point we have to perform an ugly kludge until
    #   this is fixed properly.  We will need to open the data
    #   file, read it all, and determine the set of
    #   ELEMENT_NAMEs.

    # FIX
    # Open the raw data file for the first time to collect
    #   element names.
    with open(s.data_file, "r") as data:
        for line in data:
            vals = line.split()
            if vals and vals[0] == "ELEMENT_NAME":
                bd.element_name.append(vals[1])

    # Open the raw data file for the second time to read
    #   all bond data.
    with open(s.data_file, "r") as data:

        # Initialise the counter for total bonds.
        bd.num_system_bonds = 0

        # Read the header to get num_atoms.
        header_vals = data.readline().split()
        bd.num_atoms = int(header_vals[-1])

        # Initialise the number of bonds for each atom to
        #   zero.  Note that the number of bonds by atom is
        #   accumulated by looking at both atoms listed in
        #   the bonds.  This is because the reverse bonds
        #   are not listed.  (e.g. if a Si to C bond is
        #   listed, then that C to Si bond is NOT listed.)
        for i in range(1, bd.num_atoms + 1):
            bd.num_bonds_by_atom.append(0)
            bd.num_filters.append(0)
            bd.filter_name.append([None])
            bd.filter_value.append([None])
            bd.bond_atom_by_atom.append([None])
            bd.bond_length_by_atom.append([None])
            bd.bond_order_by_atom.append([None])

        # Reset element_name to proper size (will be
        #   populated below, overwriting the kludge above).
        bd.element_name = [None] * (bd.num_atoms + 1)

        for i in range(1, bd.num_atoms + 1):

            # Read the filter names from the data file
            #   until the "ATOM_CHARGE" label is found.
            #   This indicates the end of the filter names.
            bd.num_filters[i] = 0
            while True:
                line = data.readline()
                if not line:
                    break
                vals = line.split()
                if not vals:
                    continue

                # Save the element name for this atom
                #   separately since it will come in handy
                #   to have the name readily accessible
                #   when working with the bonds later.
                if vals[0] == "ELEMENT_NAME":
                    bd.element_name[i] = vals[1]

                # Store the atom's effective charge when
                #   found.
                if vals[0] == "ATOM_CHARGE":
                    bd.atom_q_star.append(float(vals[1]))

                # Read past the orbital charges and quit.
                if vals[0] == "ATOM_ORBITAL_CHARGE":
                    num_orbitals = int(vals[1])
                    for _ in range(num_orbitals):
                        data.readline()
                    break

                # Increment the number of filters since we
                #   didn't abort.
                bd.num_filters[i] += 1

                # Save the filter name and value.
                bd.filter_name[i].append(vals[0])
                bd.filter_value[i].append(vals[1])

            # Get the intrinsic data for this element from
            #   the database.
            if bd.element_name[i].lower() == "3c":
                # Fake atom for a 3-centre bond.
                bd.atom_radius.append(
                    0.25 * s.radius_factor
                )
                bd.atom_charge_transfer.append(
                    bd.atom_q_star[i]
                )
                bd.atom_color.append(100)
            else:
                z = ed.get_element_z(
                    bd.element_name[i].lower()
                )
                bd.element_z.append(z)
                bd.atom_radius.append(
                    ed.coval_radii[z] * s.radius_factor
                )
                ct = bd.atom_q_star[i]
                for qn_l in range(ed.max_qn_l + 1):
                    ct -= ed.vale_charge[z][qn_l]
                bd.atom_charge_transfer.append(ct)
                if not s.grey_scale:
                    bd.atom_color.append(ed.color_dx[z])
                else:
                    bd.atom_color.append(ed.grey_dx[z])

            # Obtain the number of bonds that this atom has
            #   with other atoms that have a higher index
            #   number.  Consider bonds between atom i and
            #   atom j; we only count bonds where i < j.
            #   Otherwise we are double counting.
            bond_header = data.readline().split()
            num_bonds_temp = int(bond_header[-1])

            # Add to the number of bonds that this atom
            #   has.  It can be decremented later if
            #   certain bonds do not fit the requested
            #   criteria.
            bd.num_bonds_by_atom[i] += num_bonds_temp

            # Loop through those bonds to obtain the bonded
            #   atom's index number, the bond order, and
            #   the bond length.
            for j in range(1, num_bonds_temp + 1):
                vals = data.readline().split()

                bonded_atom = int(vals[0])
                b_length = float(vals[1])
                b_order = float(vals[2])

                # Store the bond information associated
                #   with this atom.  Note that the second
                #   index for the bond number of the
                #   current atom must be expressed
                #   carefully.  It is possible for earlier
                #   atoms to have had bonds to this (i)
                #   atom so that the number of bonds by
                #   atom to the current (i) atom is
                #   non-zero.  So, we must start indexing
                #   after those already recorded bonds.
                idx = (
                    j + bd.num_bonds_by_atom[i]
                    - num_bonds_temp
                )
                # Extend the by-atom lists as needed.
                while len(bd.bond_atom_by_atom[i]) <= idx:
                    bd.bond_atom_by_atom[i].append(None)
                    bd.bond_length_by_atom[i].append(None)
                    bd.bond_order_by_atom[i].append(None)
                bd.bond_atom_by_atom[i][idx] = bonded_atom
                bd.bond_length_by_atom[i][idx] = b_length
                bd.bond_order_by_atom[i][idx] = b_order

                # Store the bond information for the atom
                #   at the other end of the bond.
                bd.num_bonds_by_atom[bonded_atom] += 1
                nba = bd.num_bonds_by_atom[bonded_atom]
                while (len(bd.bond_atom_by_atom[bonded_atom])
                       <= nba):
                    bd.bond_atom_by_atom[
                        bonded_atom].append(None)
                    bd.bond_length_by_atom[
                        bonded_atom].append(None)
                    bd.bond_order_by_atom[
                        bonded_atom].append(None)
                bd.bond_atom_by_atom[bonded_atom][nba] = i
                bd.bond_length_by_atom[
                    bonded_atom][nba] = b_length
                bd.bond_order_by_atom[
                    bonded_atom][nba] = b_order

                # Store the bond information for the
                #   system.
                bd.num_system_bonds += 1
                bd.bond_atom.append(
                    [None, i, bonded_atom]
                )
                bd.bond_length.append(b_length)
                bd.bond_order.append(b_order)

                # Check that this bond fits the
                #   requirements laid out in the bond
                #   control file.
                if not check_requirements(bd):
                    # Decrement the number of bonds for
                    #   this atom in the 'ByAtom' array
                    #   and the 'system' array, and the
                    #   'ByAtom' array for the atom on
                    #   the other end of the bond.
                    bd.num_bonds_by_atom[i] -= 1
                    bd.num_system_bonds -= 1
                    bd.bond_atom.pop()
                    bd.bond_length.pop()
                    bd.bond_order.pop()
                    bd.num_bonds_by_atom[
                        bonded_atom] -= 1

            # Obtain the number of bond angles for this
            #   atom and skip past them.
            angle_vals = data.readline().split()
            num_angles = int(angle_vals[1])
            for _ in range(num_angles):
                data.readline()


def check_requirements(bd):
    """Check that the most recently added bond passes all filters.

    This function checks the current bond information to make
    sure that it conforms to the requests on the command line
    and in the control file.  It checks:

    1) The bond length and bond order against the global
       min/max thresholds from the command line.
    2) The bond order against any element-pair-specific BO
       limits from the bond control file.
    3) The bond length against any element-pair-specific BL
       limits from the bond control file.

    Returns True if the bond passes, False otherwise.
    """
    s = bd.settings
    n = bd.num_system_bonds  # Current bond index.

    # Compare to the min and max BO and BL command line
    #   requirements first.
    if (bd.bond_length[n] < s.min_bl
            or bd.bond_length[n] > s.max_bl
            or bd.bond_order[n] < s.min_bo
            or bd.bond_order[n] > s.max_bo):
        return False

    # Identify the current elements (upper-cased).
    elem1 = bd.element_name[bd.bond_atom[n][1]].upper()
    elem2 = bd.element_name[bd.bond_atom[n][2]].upper()

    # Ensure alphabetical order for comparison.
    if elem1 > elem2:
        elem1, elem2 = elem2, elem1

    # Search for the current atom element pair in the BO
    #   limits list.
    for i in range(1, bd.num_bo_limits + 1):
        if (bd.bo_limits[i][1] == elem1
                and bd.bo_limits[i][2] == elem2):
            if (bd.bond_order[n] < bd.bo_limits[i][3]
                    or bd.bond_order[n]
                    > bd.bo_limits[i][4]):
                return False
            break

    # Search for the current atom element pair in the BL
    #   limits list.
    for i in range(1, bd.num_bl_limits + 1):
        if (bd.bl_limits[i][1] == elem1
                and bd.bl_limits[i][2] == elem2):
            if (bd.bond_length[n] < bd.bl_limits[i][3]
                    or bd.bond_length[n]
                    > bd.bl_limits[i][4]):
                return False
            break

    return True


def implicit_data(bd):
    """Compute implicit data for each atom and bond.

    This function builds tags for atoms and bonds using the
    filter values and group names, identifies unique atom tags
    and unique bond tags, and ensures that the tag for atom 1
    of each bond is alphabetically earlier than the tag for
    atom 2.
    """
    # Build up the tag for each atom using the values for each
    #   filter that matches the requested group names.
    for i in range(1, bd.num_atoms + 1):
        # Initialise the atom tag with an empty string.
        tag = ""

        # Compare each requested group to each available
        #   filter and build the tag with the values.
        for j in range(1, bd.num_groups + 1):
            for k in range(1, bd.num_filters[i] + 1):
                if bd.filter_name[i][k] == \
                        bd.group_name[j]:
                    tag += (
                        bd.group_name[j][0]
                        + bd.filter_value[i][k]
                        + "_"
                    )

        # Remove the trailing underscore from the last
        #   tag addition.
        if tag.endswith("_"):
            tag = tag[:-1]

        bd.atom_tag.append(tag)

    # Make sure that the tag for atom 1 of each bond is
    #   alphanumerically earlier than the tag name for atom 2
    #   for every bond.
    for i in range(1, bd.num_system_bonds + 1):
        if bd.atom_tag[bd.bond_atom[i][1]] > \
                bd.atom_tag[bd.bond_atom[i][2]]:
            bd.bond_atom[i][1], bd.bond_atom[i][2] = \
                bd.bond_atom[i][2], bd.bond_atom[i][1]

    # Build up the tag for each bond using the tags for each
    #   atom in the bond.
    for i in range(1, bd.num_system_bonds + 1):
        bd.bond_tag.append(
            bd.atom_tag[bd.bond_atom[i][1]]
            + "__"
            + bd.atom_tag[bd.bond_atom[i][2]]
        )

    # Make a list of the unique tags for all atoms.
    bd.unique_atom_tags = [None, bd.atom_tag[1]]
    for i in range(2, bd.num_atoms + 1):
        if bd.atom_tag[i] not in bd.unique_atom_tags:
            bd.unique_atom_tags.append(bd.atom_tag[i])

    # Make a list of the unique tags for all bonds.
    bd.unique_bond_tags = [None, bd.bond_tag[1]]
    for i in range(2, bd.num_system_bonds + 1):
        if bd.bond_tag[i] not in bd.unique_bond_tags:
            bd.unique_bond_tags.append(bd.bond_tag[i])


# ============================================================
#  Helper / geometry functions
# ============================================================

def get_bond_positions(bd):
    """Calculate bond midpoint positions in a, b, c coordinates.

    For each bond, the midpoint between the two atoms is
    computed along each of the three cell axes.  If the
    midpoint distance in a given direction exceeds the bond
    length (indicating a periodic-boundary crossing), the
    midpoint is shifted to account for the periodic image.
    """
    for i in range(1, bd.num_system_bonds + 1):
        a1 = bd.bond_atom[i][1]
        a2 = bd.bond_atom[i][2]

        # For each of the three coordinate directions
        #   (a, b, c) we check to see if the mid point
        #   between the bonds in that direction is less than
        #   the total bond length.  If it is then the mid
        #   point for that direction is within the cell.  If
        #   it is not, then the mid point for that direction
        #   is outside the cell and so we have to shift it
        #   back within the cell, but we must check which
        #   direction to shift.  For example consider atoms
        #   on opposite sides of a large cell.  The
        #   calculated a-axis mid point is the centre of the
        #   cell, but this is much larger than the bond
        #   length so the true bond mid point is through the
        #   periodic cell.  The question is now, what side?
        #   First we add the magnitude of the cell in the a
        #   direction, and then if it is still outside, we
        #   subtract to make it inside.

        # Do the a-axis first.
        if (abs(bd.a_pos[a1] - bd.a_pos[a2])
                < bd.bond_length[i] + 0.1):
            a_mid = (bd.a_pos[a1] + bd.a_pos[a2]) / 2.0
        else:
            a_mid = (
                bd.a_pos[a1] + bd.a_pos[a2] + bd.mag_a
            ) / 2.0
            if a_mid > bd.mag_a:
                a_mid -= bd.mag_a
        bd.a_bond_pos.append(a_mid)

        # Do the b-axis second.
        if (abs(bd.b_pos[a1] - bd.b_pos[a2])
                < bd.bond_length[i] + 0.1):
            b_mid = (bd.b_pos[a1] + bd.b_pos[a2]) / 2.0
        else:
            b_mid = (
                bd.b_pos[a1] + bd.b_pos[a2] + bd.mag_b
            ) / 2.0
            if b_mid > bd.mag_b:
                b_mid -= bd.mag_b
        bd.b_bond_pos.append(b_mid)

        # Do the c-axis last.
        if (abs(bd.c_pos[a1] - bd.c_pos[a2])
                < bd.bond_length[i] + 0.1):
            c_mid = (bd.c_pos[a1] + bd.c_pos[a2]) / 2.0
        else:
            c_mid = (
                bd.c_pos[a1] + bd.c_pos[a2] + bd.mag_c
            ) / 2.0
            if c_mid > bd.mag_c:
                c_mid -= bd.mag_c
        bd.c_bond_pos.append(c_mid)


def get_extra_atom_info(bd):
    """Compute min, max, and average BO and BL for each atom.

    This function loops through each atom and its bonds to
    calculate per-atom statistics: minimum, maximum, and
    average bond length and bond order.
    """
    # Loop through each atom and accumulate data about its
    #   bonds.
    for i in range(1, bd.num_atoms + 1):
        # Initialise the stats for this atom.
        if bd.num_bonds_by_atom[i] == 0:
            bd.atom_min_bl.append(0.0)
            bd.atom_max_bl.append(0.0)
            bd.atom_avg_bl.append(0.0)
            bd.atom_min_bo.append(0.0)
            bd.atom_max_bo.append(0.0)
            bd.atom_avg_bo.append(0.0)
        else:
            bd.atom_min_bl.append(1000.0)
            bd.atom_max_bl.append(0.0)
            bd.atom_avg_bl.append(0.0)
            bd.atom_min_bo.append(1000.0)
            bd.atom_max_bo.append(-1000.0)
            bd.atom_avg_bo.append(0.0)

        for j in range(1, bd.num_bonds_by_atom[i] + 1):
            bl = bd.bond_length_by_atom[i][j]
            bo = bd.bond_order_by_atom[i][j]

            if bl < bd.atom_min_bl[i]:
                bd.atom_min_bl[i] = bl
            if bl > bd.atom_max_bl[i]:
                bd.atom_max_bl[i] = bl
            bd.atom_avg_bl[i] += bl

            if bo < bd.atom_min_bo[i]:
                bd.atom_min_bo[i] = bo
            if bo > bd.atom_max_bo[i]:
                bd.atom_max_bo[i] = bo
            bd.atom_avg_bo[i] += bo

        # Compute the AVG values for the BO and BL for
        #   this atom.
        if bd.num_bonds_by_atom[i] > 0:
            bd.atom_avg_bl[i] /= bd.num_bonds_by_atom[i]
            bd.atom_avg_bo[i] /= bd.num_bonds_by_atom[i]


def get_extended_positions(bd):
    """Find atom/bond positions extending outside the cell.

    Sometimes it is necessary to represent bonds extending
    outside the cell connecting an atom on one side of a cell
    with the image of another atom from the other side.
    Periodic boundary conditions make atoms from opposite
    sides of a cell connect.

    This function:
    1) Initialises the extended arrays as copies of the
       in-cell arrays.
    2) For each bond whose atom-atom distance exceeds the
       known bond length, searches over periodic images
       (-2..+2 in each direction) to find the correct
       outside-the-cell atom positions.
    3) Creates new 'atoms' outside the cell and new bond
       definitions connecting them.
    4) Records the a, b, c and x, y, z midpoint positions
       for each extended bond.
    """
    # Initialise the counters.  num_atoms_ext tracks the
    #   number of end points for the bonds including those
    #   outside the cell.  num_bonds_ext tracks the number
    #   of bonds including those outside the cell.
    bd.num_bonds_ext = bd.num_system_bonds
    bd.num_atoms_ext = bd.num_atoms

    # Initialise values for the extended atomic positions
    #   in x, y, z coordinates.
    bd.x_pos_ext = list(bd.x_pos)
    bd.y_pos_ext = list(bd.y_pos)
    bd.z_pos_ext = list(bd.z_pos)

    # Initialise values for the extended atomic positions
    #   in a, b, c coordinates.
    bd.a_pos_ext = list(bd.a_pos)
    bd.b_pos_ext = list(bd.b_pos)
    bd.c_pos_ext = list(bd.c_pos)

    # Initialise the values for the extended bond atom
    #   definitions.
    bd.bond_atom_ext = [None]
    for i in range(1, bd.num_system_bonds + 1):
        bd.bond_atom_ext.append(
            [None, bd.bond_atom[i][1], bd.bond_atom[i][2]]
        )

    # Initialise the values for the bond order of each bond.
    bd.bond_order_ext = list(bd.bond_order)

    # Initialise the tags for each bond.
    bd.bond_tag_ext = list(bd.bond_tag)

    # Initialise the tags for each atom.
    bd.atom_tag_ext = list(bd.atom_tag)

    # Ensure bond position arrays are long enough.
    while len(bd.a_bond_pos) <= bd.num_system_bonds:
        bd.a_bond_pos.append(0.0)
    while len(bd.b_bond_pos) <= bd.num_system_bonds:
        bd.b_bond_pos.append(0.0)
    while len(bd.c_bond_pos) <= bd.num_system_bonds:
        bd.c_bond_pos.append(0.0)
    while len(bd.x_bond_pos) <= bd.num_system_bonds:
        bd.x_bond_pos.append(0.0)
    while len(bd.y_bond_pos) <= bd.num_system_bonds:
        bd.y_bond_pos.append(0.0)
    while len(bd.z_bond_pos) <= bd.num_system_bonds:
        bd.z_bond_pos.append(0.0)

    # Consider each bond in turn to determine if the bond
    #   connecting these atoms extends outside the cell.
    #   If so, define a new 'atom' location outside the cell
    #   that is used as an anchor for the bond.
    for i in range(1, bd.num_system_bonds + 1):
        # Determine if the distance between the atoms of
        #   this bond is greater than the bond length.  If
        #   so, find the position of the outside image atom.
        a1 = bd.bond_atom[i][1]
        a2 = bd.bond_atom[i][2]

        dx = bd.x_pos[a1] - bd.x_pos[a2]
        dy = bd.y_pos[a1] - bd.y_pos[a2]
        dz = bd.z_pos[a1] - bd.z_pos[a2]
        min_dist = [
            None,
            math.sqrt(dx*dx + dy*dy + dz*dz),
            0.0,
        ]
        min_dist[2] = min_dist[1]

        if min_dist[1] > bd.bond_length[i] + 0.1:
            # Find the locations of the two new 'atoms'.
            min_x = [None, 0.0, 0.0]
            min_y = [None, 0.0, 0.0]
            min_z = [None, 0.0, 0.0]

            for j in range(-2, 3):
                for k in range(-2, 3):
                    for l in range(-2, 3):
                        for m in (1, 2):
                            ba = bd.bond_atom[i][m]
                            # Shift the test position.
                            tx = (
                                bd.x_pos[ba]
                                + j * bd.a[1]
                                + k * bd.b[1]
                                + l * bd.c[1]
                            )
                            ty = (
                                bd.y_pos[ba]
                                + j * bd.a[2]
                                + k * bd.b[2]
                                + l * bd.c[2]
                            )
                            tz = (
                                bd.x_pos[ba]
                                + j * bd.a[3]
                                + k * bd.b[3]
                                + l * bd.c[3]
                            )
                            # Correctly compute tz using
                            #   z_pos, not x_pos.
                            tz = (
                                bd.z_pos[ba]
                                + j * bd.a[3]
                                + k * bd.b[3]
                                + l * bd.c[3]
                            )

                            # Determine the opposite atom.
                            other = 2 if m == 1 else 1
                            oa = bd.bond_atom[i][other]

                            ddx = tx - bd.x_pos[oa]
                            ddy = ty - bd.y_pos[oa]
                            ddz = tz - bd.z_pos[oa]
                            dist = math.sqrt(
                                ddx*ddx
                                + ddy*ddy
                                + ddz*ddz
                            )

                            if dist < min_dist[m]:
                                min_dist[m] = dist
                                min_x[m] = tx
                                min_y[m] = ty
                                min_z[m] = tz

            for j_idx in (1, 2):
                # One extra atom for each new bond end
                #   extending outside the cell.
                bd.num_atoms_ext += 1

                # Add the new atom locations to the
                #   extended atom list.
                bd.x_pos_ext.append(min_x[j_idx])
                bd.y_pos_ext.append(min_y[j_idx])
                bd.z_pos_ext.append(min_z[j_idx])

                na = bd.num_atoms_ext
                a_ext = (
                    (bd.x_factor[1]
                     * bd.x_pos_ext[na]
                     + bd.y_factor[1]
                     * bd.y_pos_ext[na]
                     + bd.z_factor[1]
                     * bd.z_pos_ext[na])
                    / bd.cell_determinant
                ) * bd.mag_a
                b_ext = (
                    (bd.x_factor[2]
                     * bd.x_pos_ext[na]
                     + bd.y_factor[2]
                     * bd.y_pos_ext[na]
                     + bd.z_factor[2]
                     * bd.z_pos_ext[na])
                    / bd.cell_determinant
                ) * bd.mag_b
                c_ext = (
                    (bd.x_factor[3]
                     * bd.x_pos_ext[na]
                     + bd.y_factor[3]
                     * bd.y_pos_ext[na]
                     + bd.z_factor[3]
                     * bd.z_pos_ext[na])
                    / bd.cell_determinant
                ) * bd.mag_c
                bd.a_pos_ext.append(a_ext)
                bd.b_pos_ext.append(b_ext)
                bd.c_pos_ext.append(c_ext)

                # Include the tags for this atom.
                bd.atom_tag_ext.append(
                    bd.atom_tag[bd.bond_atom[i][j_idx]]
                )

            # Make a new bond definition that includes one
            #   old atom from inside the cell, and one new
            #   image 'atom' outside the cell.  Also carry
            #   over the bond order value.
            bd.num_bonds_ext += 1
            bd.bond_atom_ext.append([
                None,
                bd.bond_atom_ext[i][2],
                bd.num_atoms_ext - 1,
            ])
            bd.bond_order_ext.append(
                bd.bond_order_ext[i]
            )
            bd.bond_tag_ext.append(bd.bond_tag[i])

            ne = bd.num_bonds_ext

            # Make sure that the tag for atom 1 of this
            #   bond is alphanumerically earlier than the
            #   tag for atom 2.
            if (bd.atom_tag_ext[bd.bond_atom_ext[ne][1]]
                    > bd.atom_tag_ext[
                        bd.bond_atom_ext[ne][2]]):
                (bd.bond_atom_ext[ne][1],
                 bd.bond_atom_ext[ne][2]) = (
                    bd.bond_atom_ext[ne][2],
                    bd.bond_atom_ext[ne][1],
                )

            # Assign the bond tag for this extended bond.
            bd.bond_tag_ext[ne] = (
                bd.atom_tag_ext[
                    bd.bond_atom_ext[ne][1]]
                + bd.atom_tag_ext[
                    bd.bond_atom_ext[ne][2]]
            )

            # Record the position of this extended bond.
            a_bp = (
                bd.a_pos_ext[bd.bond_atom_ext[ne][1]]
                + bd.a_pos_ext[bd.bond_atom_ext[ne][2]]
            ) / 2.0
            b_bp = (
                bd.b_pos_ext[bd.bond_atom_ext[ne][1]]
                + bd.b_pos_ext[bd.bond_atom_ext[ne][2]]
            ) / 2.0
            c_bp = (
                bd.c_pos_ext[bd.bond_atom_ext[ne][1]]
                + bd.c_pos_ext[bd.bond_atom_ext[ne][2]]
            ) / 2.0
            while len(bd.a_bond_pos) <= ne:
                bd.a_bond_pos.append(0.0)
                bd.b_bond_pos.append(0.0)
                bd.c_bond_pos.append(0.0)
                bd.x_bond_pos.append(0.0)
                bd.y_bond_pos.append(0.0)
                bd.z_bond_pos.append(0.0)
            bd.a_bond_pos[ne] = a_bp
            bd.b_bond_pos[ne] = b_bp
            bd.c_bond_pos[ne] = c_bp
            rl = bd.real_lattice
            bd.x_bond_pos[ne] = (
                a_bp * rl[1][1] / bd.mag_a
                + b_bp * rl[2][1] / bd.mag_b
                + c_bp * rl[3][1] / bd.mag_c
            )
            bd.y_bond_pos[ne] = (
                a_bp * rl[1][2] / bd.mag_a
                + b_bp * rl[2][2] / bd.mag_b
                + c_bp * rl[3][2] / bd.mag_c
            )
            bd.z_bond_pos[ne] = (
                a_bp * rl[1][3] / bd.mag_a
                + b_bp * rl[2][3] / bd.mag_b
                + c_bp * rl[3][3] / bd.mag_c
            )

            # Use the old bond definition and reassign one
            #   atom to be the other image 'atom' outside
            #   the cell.
            bd.bond_atom_ext[i][2] = bd.num_atoms_ext

            # Make sure tag ordering is correct.
            if (bd.atom_tag_ext[bd.bond_atom_ext[i][1]]
                    > bd.atom_tag_ext[
                        bd.bond_atom_ext[i][2]]):
                (bd.bond_atom_ext[i][1],
                 bd.bond_atom_ext[i][2]) = (
                    bd.bond_atom_ext[i][2],
                    bd.bond_atom_ext[i][1],
                )

            # Assign the bond tag for this extended bond.
            bd.bond_tag_ext[i] = (
                bd.atom_tag_ext[bd.bond_atom_ext[i][1]]
                + bd.atom_tag_ext[bd.bond_atom_ext[i][2]]
            )

        # Record the position of this bond (done for both
        #   extended and non-extended cases).
        a_bp = (
            bd.a_pos_ext[bd.bond_atom_ext[i][1]]
            + bd.a_pos_ext[bd.bond_atom_ext[i][2]]
        ) / 2.0
        b_bp = (
            bd.b_pos_ext[bd.bond_atom_ext[i][1]]
            + bd.b_pos_ext[bd.bond_atom_ext[i][2]]
        ) / 2.0
        c_bp = (
            bd.c_pos_ext[bd.bond_atom_ext[i][1]]
            + bd.c_pos_ext[bd.bond_atom_ext[i][2]]
        ) / 2.0
        bd.a_bond_pos[i] = a_bp
        bd.b_bond_pos[i] = b_bp
        bd.c_bond_pos[i] = c_bp
        rl = bd.real_lattice
        bd.x_bond_pos[i] = (
            a_bp * rl[1][1] / bd.mag_a
            + b_bp * rl[2][1] / bd.mag_b
            + c_bp * rl[3][1] / bd.mag_c
        )
        bd.y_bond_pos[i] = (
            a_bp * rl[1][2] / bd.mag_a
            + b_bp * rl[2][2] / bd.mag_b
            + c_bp * rl[3][2] / bd.mag_c
        )
        bd.z_bond_pos[i] = (
            a_bp * rl[1][3] / bd.mag_a
            + b_bp * rl[2][3] / bd.mag_b
            + c_bp * rl[3][3] / bd.mag_c
        )


# ============================================================
#  Scatter-plot output
# ============================================================

def print_scatter_data(bd):
    """Print scatter-plot data files for atoms and bonds.

    Produces:
    - ATOM_<tag>.dat files with per-atom Q*, BO/BL statistics,
      and a, b, c positions for each unique atom tag.
    - <bond_tag>.dat files with BO, BL, atom indices, and
      a, b, c midpoint positions for each unique bond tag.
    """
    # Obtain the bond position information.
    get_bond_positions(bd)

    # Obtain the extra atom information (MaxBO, MinBO, etc.).
    get_extra_atom_info(bd)

    num_unique_atom = len(bd.unique_atom_tags) - 1
    num_unique_bond = len(bd.unique_bond_tags) - 1

    # Print the atomic data first.
    for i in range(1, num_unique_atom + 1):
        tag = bd.unique_atom_tags[i]
        fname = f"ATOM_{tag}.dat"
        with open(fname, "w") as f:
            # Print the header for this file.
            f.write(
                "atom  QStar MaxBL MinBL AVGBL "
                "MaxBO MinBO AVGBO"
                "     a     b     c\n"
            )

            for j in range(1, bd.num_atoms + 1):
                if bd.atom_tag[j] == tag:
                    f.write(
                        f"{j:5d} "
                        f"{bd.atom_q_star[j]:5.3f} "
                        f"{bd.atom_max_bl[j]:5.3f} "
                        f"{bd.atom_min_bl[j]:5.3f} "
                        f"{bd.atom_avg_bl[j]:5.3f} "
                        f"{bd.atom_max_bo[j]:5.3f} "
                        f"{bd.atom_min_bo[j]:5.3f} "
                        f"{bd.atom_avg_bo[j]:5.3f} "
                        f"{bd.a_pos[j]:5.3f} "
                        f"{bd.b_pos[j]:5.3f} "
                        f"{bd.c_pos[j]:5.3f}\n"
                    )

    # Print the bond data second.
    for i in range(1, num_unique_bond + 1):
        tag = bd.unique_bond_tags[i]
        fname = f"{tag}.dat"
        with open(fname, "w") as f:
            # Print the header for this file.
            f.write(
                "BO     BL(au)    #1    #2"
                "    a      b      c\n"
            )

            for j in range(1, bd.num_system_bonds + 1):
                if bd.bond_tag[j] == tag:
                    f.write(
                        f"{bd.bond_order[j]:5.4f} "
                        f"{bd.bond_length[j]:5.4f} "
                        f"{bd.bond_atom[j][1]:5d} "
                        f"{bd.bond_atom[j][2]:5d} "
                        f"{bd.a_bond_pos[j]:5.4f} "
                        f"{bd.b_bond_pos[j]:5.4f} "
                        f"{bd.c_bond_pos[j]:5.4f}\n"
                    )


# ============================================================
#  openDX model output
# ============================================================

def print_model_data(bd):
    """Print openDX model files for structural visualisation.

    Produces, for each unique tag:
    - atomPos_<tag>.dx : atom positions, colours, and radii.
    - <tag>.dx : Q* data (positions, charge transfer, radii)
      and bond order data (extended positions, connections,
      bond order values).

    The model can be atom-based (files contain specific sets
    of atoms; bonds must include at least one atom from the
    set) or bond-based (files contain specific sets of bonds;
    atoms must be endpoints of included bonds).
    """
    s = bd.settings

    # Now we need to collect the bond data since some of the
    #   atoms may be bonded to atoms outside the current cell.
    #   Consider an atom near the far corner of a cubic cell.
    #   It will not bond to an atom in the near corner all the
    #   way across the cell.  Instead it will bond to that same
    #   atom in the neighbour cell.  We must determine the
    #   position of that atom in the neighbour cell, and we
    #   must define a bond to it.
    get_extended_positions(bd)

    # If the openDX file is to be based on the atoms, then the
    #   files that are made will contain a specific set of
    #   atoms, and the bonds must contain at least one atom
    #   from that specific set.  If the openDX file is to be
    #   based on bonds, then only atoms that are present in the
    #   unique bond sets will be included.
    if s.atom_based_dx:
        num_unique_tags = len(bd.unique_atom_tags) - 1
    else:
        num_unique_tags = len(bd.unique_bond_tags) - 1

    # Print the lattice information first.
    print_lattice_dx(bd)

    # Now, loop to print all the Q* and BO information of
    #   each requested kind, and the requested atom positions
    #   (in a separate file).
    for i in range(1, num_unique_tags + 1):

        # Create the atom list and bond list for this file.
        create_dx_list(bd, i)

        tag = bd.unique_tag_dx

        # ---- Print the atom positions file ----
        # This is necessary if we want to plot the atoms but
        #   don't want to use charge data to colour it.
        with open(f"atomPos_{tag}.dx", "w") as f:

            # Atom positions.
            f.write(
                "object 1 class array type float "
                "rank 1 shape 3 items "
                f"{bd.num_atoms_dx} data follows\n"
            )
            for j in range(1, bd.num_atoms_dx + 1):
                f.write(
                    f"{bd.x_pos_ext_dx[j]} "
                    f"{bd.y_pos_ext_dx[j]} "
                    f"{bd.z_pos_ext_dx[j]}\n"
                )
            f.write("\n")

            # Colour values.
            f.write(
                "object 2 class array type float "
                "rank 0 items "
                f"{bd.num_atoms_dx} data follows\n"
            )
            for j in range(1, bd.num_atoms_dx + 1):
                f.write(f"{bd.atom_color_dx[j]}\n")
            f.write(
                'attribute "dep" string '
                '"positions"\n\n'
            )

            # Atom sphere size.
            f.write(
                "object 3 class array type float "
                "rank 0 items "
                f"{bd.num_atoms_dx} data follows\n"
            )
            for j in range(1, bd.num_atoms_dx + 1):
                f.write(f"{bd.atom_radius_dx[j]}\n")
            f.write(
                'attribute "dep" string '
                '"positions"\n\n'
            )

            # Combine into the spheres field.
            f.write('object "atoms" class field\n')
            f.write('component "positions" value 1\n')
            f.write('component "data" value 2\n')
            f.write('component "sizes" value 3\n')
            f.write("end\n\n")

        # ---- Print the main DX file ----
        with open(f"{tag}.dx", "w") as f:

            # Q* data: atom positions.
            f.write(
                "object 1 class array type float "
                "rank 1 shape 3 items "
                f"{bd.num_atoms_dx} data follows\n"
            )
            for j in range(1, bd.num_atoms_dx + 1):
                f.write(
                    f"{bd.x_pos_ext_dx[j]} "
                    f"{bd.y_pos_ext_dx[j]} "
                    f"{bd.z_pos_ext_dx[j]}\n"
                )
            f.write("\n")

            # Charge transfer data.
            f.write(
                "object 2 class array type float "
                "rank 0 items "
                f"{bd.num_atoms_dx} data follows\n"
            )
            for j in range(1, bd.num_atoms_dx + 1):
                f.write(
                    f"{bd.atom_charge_transfer_dx[j]}"
                    "\n"
                )
            f.write(
                'attribute "dep" string '
                '"positions"\n\n'
            )

            # Atom sphere size.
            f.write(
                "object 3 class array type float "
                "rank 0 items "
                f"{bd.num_atoms_dx} data follows\n"
            )
            for j in range(1, bd.num_atoms_dx + 1):
                f.write(f"{bd.atom_radius_dx[j]}\n")
            f.write(
                'attribute "dep" string '
                '"positions"\n\n'
            )

            # Combine into the Q* field.
            f.write('object "qstar" class field\n')
            f.write('component "positions" value 1\n')
            f.write('component "data" value 2\n')
            f.write('component "sizes" value 3\n')
            f.write("\n")

            # ---- Bond data ----

            # All atom positions (including extended).
            f.write(
                "object 4 class array type float "
                "rank 1 shape 3 items "
                f"{bd.num_atoms_ext_dx} data follows\n"
            )
            for j in range(1, bd.num_atoms_ext_dx + 1):
                f.write(
                    f"{bd.x_pos_ext_dx[j]} "
                    f"{bd.y_pos_ext_dx[j]} "
                    f"{bd.z_pos_ext_dx[j]}\n"
                )
            f.write("\n")

            # Connections between data points.
            f.write(
                "object 5 class array type int "
                "rank 1 shape 2 items "
                f"{bd.num_bonds_dx} data follows\n"
            )
            for j in range(1, bd.num_bonds_dx + 1):
                # Subtract 1 because openDX is C-based
                #   and array refs start from 0.
                fa = bd.bond_atom_dx[j][1] - 1
                sa = bd.bond_atom_dx[j][2] - 1
                f.write(f"{fa}  {sa}\n")
            f.write(
                'attribute "ref" string '
                '"positions"\n'
            )
            f.write(
                'attribute "element type" string '
                '"lines"\n\n'
            )

            # Bond order values.
            f.write(
                "object 6 class array type float "
                "rank 0 items "
                f"{bd.num_bonds_dx} data follows\n"
            )
            for j in range(1, bd.num_bonds_dx + 1):
                f.write(f"{bd.bond_order_dx[j]}\n")
            f.write(
                'attribute "dep" string '
                '"connections"\n\n'
            )

            # Combine into the BO field.
            f.write('object "bo" class field\n')
            f.write('component "positions" value 4\n')
            f.write(
                'component "connections" value 5\n'
            )
            f.write('component "data" value 6\n')
            f.write("\n")

            # End the DX file.  The DX file format
            #   requires a blank line at the end.
            f.write("end\n\n")

        # Release memory that is not necessary any more.
        bd.bond_order_dx = [None]
        bd.abc_bond_pos = [None]
        bd.xyz_bond_pos = [None]
        bd.bond_atom_dx = [None]
        bd.atom_radius_dx = [None]
        bd.atom_color_dx = [None]
        bd.atom_charge_transfer_dx = [None]
        bd.x_pos_ext_dx = [None]
        bd.y_pos_ext_dx = [None]
        bd.z_pos_ext_dx = [None]


def print_lattice_dx(bd):
    """Print lattice information in openDX format.

    Writes lattice.dx containing the unit cell as a
    2x2x2 grid field for openDX visualisation.
    """
    # Print the lattice information first.  It does not need
    #   to be included in every file.
    with open("lattice.dx", "w") as f:
        a = bd.a
        b = bd.b
        c = bd.c
        f.write("object 1 class gridpositions "
                "counts 2 2 2\n")
        f.write("origin 0 0 0\n")
        f.write(f"delta {a[1]} {a[2]} {a[3]}\n")
        f.write(f"delta {b[1]} {b[2]} {b[3]}\n")
        f.write(f"delta {c[1]} {c[2]} {c[3]}\n")
        f.write("\n")
        f.write("object 2 class gridconnections "
                "counts 2 2 2\n")
        f.write("object 3 class array type float "
                "rank 0 items 8 data follows\n")
        f.write("1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n")
        f.write('attribute "dep" string '
                '"positions"\n')
        f.write("\n")
        f.write('object "lattice" class field\n')
        f.write('component "positions" value 1\n')
        f.write('component "connections" value 2\n')
        f.write('component "data" value 3\n')
        f.write("end\n\n")


def create_dx_list(bd, current_iter):
    """Create atom and bond lists for a specific DX file.

    Depending on whether the model is atom-based or bond-based,
    this function collects the relevant atoms and bonds for the
    unique tag corresponding to current_iter.  It also produces
    a mapping between old atom numbers and new DX atom numbers.
    """
    s = bd.settings

    if s.atom_based_dx:
        # Obtain the tag to use in collecting atoms.
        bd.unique_tag_dx = bd.unique_atom_tags[current_iter]

        # Loop through all the atoms and make a list of the
        #   ones that belong to this group.  Also, make a
        #   mapping between the old atom index number and
        #   the new one.
        bd.num_atoms_dx = 0
        bd.num_atoms_ext_dx = 0
        bd.atom_map_2_dx = [None] * (bd.num_atoms_ext + 1)
        bd.x_pos_ext_dx = [None]
        bd.y_pos_ext_dx = [None]
        bd.z_pos_ext_dx = [None]
        bd.atom_color_dx = [None]
        bd.atom_radius_dx = [None]
        bd.atom_charge_transfer_dx = [None]

        for i in range(1, bd.num_atoms_ext + 1):
            if bd.atom_tag_ext[i] == bd.unique_tag_dx:
                add_dx_atom(bd, i)

        # Loop through all the bonds and include the atoms
        #   in the bonds that were not already collected.
        bd.num_bonds_dx = 0
        bd.bond_atom_dx = [None]
        bd.bond_order_dx = [None]
        bd.abc_bond_pos = [None]
        bd.xyz_bond_pos = [None]

        for i in range(1, bd.num_bonds_ext + 1):
            # Check to see if one of the atoms in the bond
            #   contains an atom that was already included.
            #   If so we may need to add the other atom.
            if bd.unique_tag_dx in bd.bond_tag_ext[i]:
                # Add the atom that was not already
                #   accounted for, if there is one.
                a1 = bd.bond_atom_ext[i][1]
                a2 = bd.bond_atom_ext[i][2]
                if bd.atom_tag_ext[a1] != \
                        bd.unique_tag_dx:
                    add_dx_atom(bd, a1)
                elif bd.atom_tag_ext[a2] != \
                        bd.unique_tag_dx:
                    add_dx_atom(bd, a2)

                # Add this bond to the openDX list.
                add_dx_bond(bd, i)

    else:
        # Obtain the tag to use in collecting bonds.
        bd.unique_tag_dx = bd.unique_bond_tags[current_iter]

        # Initialise the count of atoms.
        bd.num_atoms_dx = 0
        bd.num_atoms_ext_dx = 0
        bd.atom_map_2_dx = [None] * (bd.num_atoms_ext + 1)
        bd.x_pos_ext_dx = [None]
        bd.y_pos_ext_dx = [None]
        bd.z_pos_ext_dx = [None]
        bd.atom_color_dx = [None]
        bd.atom_radius_dx = [None]
        bd.atom_charge_transfer_dx = [None]

        # Loop through all the bonds (including the
        #   extended set) and include only those that match
        #   the requested bond type.  Also include the
        #   atoms associated with those included bonds.
        bd.num_bonds_dx = 0
        bd.bond_atom_dx = [None]
        bd.bond_order_dx = [None]
        bd.abc_bond_pos = [None]
        bd.xyz_bond_pos = [None]

        for i in range(1, bd.num_bonds_ext + 1):
            if bd.unique_tag_dx == bd.bond_tag_ext[i]:
                add_dx_atom(bd, bd.bond_atom_ext[i][1])
                add_dx_atom(bd, bd.bond_atom_ext[i][2])
                add_dx_bond(bd, i)


def add_dx_atom(bd, atom_number):
    """Add an atom to the current openDX atom list.

    Accumulates both the regular and extended atom count,
    copies the orthogonal coordinate positions, and stores
    the atom mapping.  For atoms inside the cell with matching
    tags, also copies colour, radius, and charge transfer data.
    """
    s = bd.settings

    # Accumulate both the regular and extended atom count.
    bd.num_atoms_ext_dx += 1

    # Copy the orthogonal coordinate positions.
    bd.x_pos_ext_dx.append(bd.x_pos_ext[atom_number])
    bd.y_pos_ext_dx.append(bd.y_pos_ext[atom_number])
    bd.z_pos_ext_dx.append(bd.z_pos_ext[atom_number])

    # Store the atom associations between the original and
    #   the DX.  The index number is the old atom number,
    #   and the value is the new DX atom number.
    # Extend the map if needed.
    while len(bd.atom_map_2_dx) <= atom_number:
        bd.atom_map_2_dx.append(None)
    bd.atom_map_2_dx[atom_number] = bd.num_atoms_ext_dx

    # Copy important information for the atoms that will be
    #   shown which includes only those atoms where the given
    #   atom number is less than the number of atoms in the
    #   cell (excludes the extended set).  In the case of the
    #   atom_based_dx we add the further restriction that the
    #   atom tag must exactly match the current unique DX tag.
    if (((s.atom_based_dx
          and bd.atom_tag_ext[atom_number]
          == bd.unique_tag_dx)
         or s.bond_based_dx)
            and atom_number <= bd.num_atoms):
        bd.num_atoms_dx += 1
        bd.atom_color_dx.append(
            bd.atom_color[atom_number]
        )
        bd.atom_radius_dx.append(
            bd.atom_radius[atom_number]
        )
        bd.atom_charge_transfer_dx.append(
            bd.atom_charge_transfer[atom_number]
        )


def add_dx_bond(bd, bond_number):
    """Add a bond to the current openDX bond list.

    Maps the bond atom indices from the original numbering to
    the DX numbering, and copies the bond order and midpoint
    position data.
    """
    bd.num_bonds_dx += 1
    bd.bond_atom_dx.append([
        None,
        bd.atom_map_2_dx[
            bd.bond_atom_ext[bond_number][1]],
        bd.atom_map_2_dx[
            bd.bond_atom_ext[bond_number][2]],
    ])
    bd.bond_order_dx.append(
        bd.bond_order_ext[bond_number]
    )

    bd.abc_bond_pos.append([
        None,
        bd.a_bond_pos[bond_number],
        bd.b_bond_pos[bond_number],
        bd.c_bond_pos[bond_number],
    ])
    bd.xyz_bond_pos.append([
        None,
        bd.x_bond_pos[bond_number],
        bd.y_bond_pos[bond_number],
        bd.z_bond_pos[bond_number],
    ])


# ============================================================
#  3-D mesh output
# ============================================================

def print_3d_mesh_data(bd):
    """Print a 3-D openDX mesh of bond order values.

    For each unique bond tag, creates a regular 3-D grid,
    evaluates the weighted bond order at each grid point,
    and writes the result as an openDX field.
    """
    s = bd.settings

    # Now we need to collect the bond data and atom data for
    #   the extended lattice because we need to know about all
    #   contributions from bonds that extend into neighbouring
    #   cells.
    s.bond_based_dx = True
    get_extended_positions(bd)

    # Print the lattice information since it is only needed
    #   once.
    print_lattice_dx(bd)

    # Create the openDX mesh.
    make_dx_mesh(bd)

    num_unique_bond = len(bd.unique_bond_tags) - 1

    # Now loop to print the mesh for each requested tag.
    for unique_tag in range(1, num_unique_bond + 1):

        # Create the atom list and bond list for this file.
        create_dx_list(bd, unique_tag)

        # Define the 3D mesh of points.
        eval_dx_mesh(bd)

        # Print the 3D mesh of points.
        print_dx_mesh(bd, unique_tag)

        # Release unused memory.
        bd.bond_order_dx = [None]
        bd.abc_bond_pos = [None]
        bd.xyz_bond_pos = [None]
        bd.bond_atom_dx = [None]
        bd.atom_radius_dx = [None]
        bd.atom_color_dx = [None]
        bd.atom_charge_transfer_dx = [None]
        bd.x_pos_ext_dx = [None]
        bd.y_pos_ext_dx = [None]
        bd.z_pos_ext_dx = [None]


def make_dx_mesh(bd):
    """Create the 3-D mesh grid of scan points.

    Computes the total number of mesh points and their x, y, z
    positions based on the lattice vectors and the requested
    number of grid points along each axis.
    """
    s = bd.settings
    nmp = s.num_mesh_points

    # Compute the total number of points in this scan.
    bd.num_scan_points = nmp[0] * nmp[1] * nmp[2]

    # Compute the x, y, z delta necessary along each a, b, c
    #   axis.  delta[a,b,c][x,y,z]
    for axis in range(1, 4):  # Loop over x, y, z
        # a axis contribution to x, y, z delta.
        bd.delta[1][axis] = (
            bd.real_lattice[1][axis] / (nmp[0] - 1)
        )
        # b axis contribution to x, y, z delta.
        bd.delta[2][axis] = (
            bd.real_lattice[2][axis] / (nmp[1] - 1)
        )
        # c axis contribution to x, y, z delta.
        bd.delta[3][axis] = (
            bd.real_lattice[3][axis] / (nmp[2] - 1)
        )

    # Compute the position of each point along the path.
    # Reinitialise scan_points as 1-indexed lists for x, y, z.
    bd.scan_points = [None, [None], [None], [None]]
    point = 0

    for a_pt in range(nmp[0]):
        for b_pt in range(nmp[1]):
            for c_pt in range(nmp[2]):
                point += 1
                for axis in range(1, 4):
                    bd.scan_points[axis].append(
                        bd.delta[1][axis] * a_pt
                        + bd.delta[2][axis] * b_pt
                        + bd.delta[3][axis] * c_pt
                    )


def eval_dx_mesh(bd):
    """Evaluate bond order values at each mesh point.

    For each mesh point, considers each bond in the system
    that is within the limit distance.  The contribution of
    each nearby bond is weighted by a Gaussian function of
    the bond-to-mesh-point distance.  The weights are
    normalised so they sum to one, and the weighted average
    bond order is stored.
    """
    s = bd.settings

    # Reinitialise the mesh BO array.
    bd.mesh_bo_dx = [None] * (bd.num_scan_points + 1)

    # For each mesh point consider each bond in the system
    #   that is sufficiently close to consider contributing
    #   to the bond order mesh plot.
    for point in range(1, bd.num_scan_points + 1):

        # Initialise a count of the number of bonds that
        #   will contribute to the current mesh point.
        num_contrib_bonds = 0
        gaussian_dist_factor = [None]
        contrib_bond_number = [None]

        # Loop over all bonds to check which ones should
        #   contribute and how much.
        for bond in range(1, bd.num_bonds_dx + 1):
            # Compute the separation between this bond
            #   location and this point.
            diff_sq = 0.0
            for axis in range(1, 4):
                diff_sq += (
                    bd.scan_points[axis][point]
                    - bd.xyz_bond_pos[bond][axis]
                ) ** 2
            distance = math.sqrt(diff_sq)

            # Abort if the distance is too great.
            if distance > s.limit_dist:
                continue

            # Increment the number of bonds contributing.
            num_contrib_bonds += 1

            # Determine the distance weighting factor.
            gaussian_dist_factor.append(
                math.exp(-1.0 * s.alpha * distance**2)
            )

            # Save the associated bond number.
            contrib_bond_number.append(bond)

        # Compute the normalised weighting factor.  (The
        #   sum of all weighting factors should equal one.)
        total_weight = 0.0
        for cb in range(1, num_contrib_bonds + 1):
            total_weight += gaussian_dist_factor[cb]
        if total_weight > 0.0:
            for cb in range(1, num_contrib_bonds + 1):
                gaussian_dist_factor[cb] /= total_weight

        # Accumulate the total bond order for this mesh
        #   point.
        bd.mesh_bo_dx[point] = 0.0
        if num_contrib_bonds == 0:
            print(f"NO BONDS {point}")
        print(num_contrib_bonds)
        for cb in range(1, num_contrib_bonds + 1):
            bd.mesh_bo_dx[point] += (
                gaussian_dist_factor[cb]
                * bd.bond_order_dx[
                    contrib_bond_number[cb]]
            )


def print_dx_mesh(bd, unique_tag):
    """Print the evaluated mesh to an openDX file.

    Writes the mesh as a gridpositions/gridconnections/data
    field in the openDX format.
    """
    s = bd.settings
    nmp = s.num_mesh_points

    out_file = f"{unique_tag}.dx"

    with open(out_file, "w") as f:
        f.write(
            f"object 1 class gridpositions counts "
            f"{nmp[0]} {nmp[1]} {nmp[2]}\n"
        )
        f.write("origin 0 0 0\n")
        f.write(
            f"delta {bd.delta[1][1]} "
            f"{bd.delta[1][2]} {bd.delta[1][3]}\n"
        )
        f.write(
            f"delta {bd.delta[2][1]} "
            f"{bd.delta[2][2]} {bd.delta[2][3]}\n"
        )
        f.write(
            f"delta {bd.delta[3][1]} "
            f"{bd.delta[3][2]} {bd.delta[3][3]}\n\n"
        )

        f.write(
            f"object 2 class gridconnections counts "
            f"{nmp[0]} {nmp[1]} {nmp[2]}\n"
        )
        f.write(
            f"object 3 class array type float rank 0 "
            f"items {bd.num_scan_points} data follows\n"
        )

        point_count = 0
        for _ in range(nmp[0]):
            for _ in range(nmp[1]):
                for _ in range(nmp[2]):
                    point_count += 1
                    val = bd.mesh_bo_dx[point_count]
                    f.write(f"{val:12.8f} ")
                    if point_count % 5 == 0:
                        f.write("\n")

        # Add a final newline if necessary.
        if point_count != 0:
            f.write("\n")

        # Tack on the ending information.
        f.write(
            'attribute "dep" string '
            '"positions"\n\n'
        )
        f.write('object "dataField" class field\n')
        f.write('component "positions" value 1\n')
        f.write('component "connections" value 2\n')
        f.write('component "data" value 3\n')
        f.write("end\n\n")


# ============================================================
#  Bond profile output
# ============================================================

def print_bond_profile(bd):
    """Print bond profile data along a lattice axis.

    This function produces a simple data file that can be
    plotted as a curve indicating the bond profile along a
    specific axis.  The bond profile at a given point along
    the profile_axis is simply the sum of the bond orders
    for all the bonds that cross that point divided by the
    number of bonds.

    Steps:
    1) Create arrays where each index represents a point
       on the requested axis at 0.01 Angstrom spacing.
    2) For each bond, increment the cumulative bond order
       and bond count arrays for all points between the
       two bond endpoints.
    3) Divide to get the average bond order at each point.
    """
    s = bd.settings

    # To compute the bond profile we:
    #   1) Create a pair of arrays where each index represents
    #      a point on the requested axis.  The number of
    #      indices must be sufficiently high to present a clean
    #      looking graph.
    #   2) The values in one array will represent the
    #      cumulative number of bonds at that point.  The
    #      values in the other array will represent the
    #      cumulative bond order at that point.
    #   3) Consider each bond in turn.  Increment the array
    #      values for both arrays for the points between the
    #      bond endpoints.

    # Obtain the positions of atoms outside the cell that
    #   immediately bond to the atoms inside the cell.
    get_extended_positions(bd)

    # Copy the axis of choice for easy reference.
    if s.profile_axis == "a":
        axis_pos_ext = bd.a_pos_ext
        axis_mag = bd.mag_a
    elif s.profile_axis == "b":
        axis_pos_ext = bd.b_pos_ext
        axis_mag = bd.mag_b
    elif s.profile_axis == "c":
        axis_pos_ext = bd.c_pos_ext
        axis_mag = bd.mag_c
    else:
        print(
            f"Invalid profile axis choice "
            f"{s.profile_axis}.  Aborting."
        )
        return

    # Create the axis_points, cumul_bond_order, and
    #   cumul_num_bonds arrays.
    num_points = int(axis_mag / 0.01) + 1

    # axis_points[0] = -0.01 (sentinel).  The strange start,
    #   end, and reference scheme is for numerical accuracy.
    axis_points = [-0.01]
    cumul_bond_order = [0.0]
    cumul_num_bonds = [0]
    for i in range(num_points):
        axis_points.append(i * 0.01)
        cumul_bond_order.append(0.0)
        cumul_num_bonds.append(0)

    # Determine the maximum number of binary searches that
    #   will be necessary for each bond to find its array
    #   index.
    num_searches = math.log(num_points) / math.log(2)
    if num_searches - int(num_searches) > 0:
        num_searches = int(num_searches) + 1
    else:
        num_searches = int(num_searches)

    # Consider each atom in the extended system to determine
    #   its index.
    atom_index = [None]
    for i in range(1, bd.num_atoms_ext + 1):
        if axis_pos_ext[i] < 0:
            atom_index.append(0)
        elif axis_pos_ext[i] > axis_mag:
            atom_index.append(num_points)
        else:
            # Perform a binary search for the most
            #   appropriate index number.
            search_point = num_points // 2
            last_point = 0
            for _ in range(num_searches):
                if (axis_points[search_point]
                        < axis_pos_ext[i]):
                    current = search_point
                    search_point += (
                        abs(last_point - search_point)
                        // 2
                    )
                    last_point = current
                else:
                    current = search_point
                    search_point -= (
                        abs(last_point - search_point)
                        // 2
                    )
                    last_point = current
            atom_index.append(search_point)

    for i in range(1, bd.num_bonds_ext + 1):
        # Obtain the indices for the atoms of this bond.
        left = atom_index[bd.bond_atom_ext[i][1]]
        right = atom_index[bd.bond_atom_ext[i][2]]

        if left > right:
            left, right = right, left

        # Increment the cumulative arrays for all points
        #   between the left end point and the right end.
        while left < right:
            cumul_bond_order[left] += (
                bd.bond_order_ext[i]
            )
            cumul_num_bonds[left] += 1
            left += 1

    # Print out the results.
    with open(s.bond_profile_file, "w") as f:
        axis_label = s.profile_axis + "-axis"
        f.write(
            f"{axis_label}  bondProfile  bondOrder"
            f"  numBonds\n"
        )
        for i in range(1, num_points):
            if cumul_bond_order[i] != 0:
                bp = (
                    cumul_bond_order[i]
                    / cumul_num_bonds[i]
                )
            else:
                bp = 0.0
            f.write(
                f"{axis_points[i]:5.4f} "
                f"{bp:5.4f} "
                f"{cumul_bond_order[i]:5.4f} "
                f"{cumul_num_bonds[i]:5d}\n"
            )


# ============================================================
#  VTK Legacy output
# ============================================================

def print_model_vtk(bd):
    """Print a structural model in VTK Legacy format.

    This function takes the data initially displayed in openDX
    format in print_model_data() and instead displays it in
    the legacy VTK format for use in programs such as Paraview
    or VisIT.

    Only the model shape is implemented for now; colour tables
    representing bond strength and charge transfer will be
    included later.

    This function reuses many of the same helper functions as
    the openDX output (getExtendedPositions, createDXList,
    etc.) so as not to change the program too much.

    Notes on VTK format:
    - For Paraview, set the first line to the newest VTK
      DataFile version.
    - For VisIT, change the version to 2.0 as VisIT may not
      support newer versions.
    - VERTICES are used so that Paraview can apply spherical
      glyphs to atom positions.
    - If spherical glyphs are not applying to vertices, try
      changing the vertex point count as a quick fix.
    """
    s = bd.settings

    # Now we need to collect the bond data since some of the
    #   atoms may be bonded to atoms outside the current cell.
    #   Consider an atom near the far corner of a cubic cell.
    #   It will not bond to an atom in the near corner all the
    #   way across the cell.  Instead it will bond to that same
    #   atom in the neighbour cell.  We must determine the
    #   position of that atom in the neighbour cell, and we
    #   must define a bond to it.
    get_extended_positions(bd)

    # Since this function will only be bond based for the time
    #   being, then the tags used for it need only specify
    #   atoms in the unique bonds.
    num_unique_tags = len(bd.unique_bond_tags) - 1

    # Now, loop to print all the Q* and BO information of each
    #   requested kind, and the requested atom positions.  For
    #   early implementation only the atoms and bonds will be
    #   printed.  The Q* data will be included either by the
    #   original author or the next person to touch this code.
    #   Please don't forget to include it.
    for i in range(1, num_unique_tags + 1):

        # Create the atom list and bond list for this file.
        #   Since openDX and VTK both list the positions of
        #   atoms/points in x y z format, simply using this
        #   function as is saves time.
        create_dx_list(bd, i)

        tag = bd.unique_tag_dx

        # First we open a VTK file and print the header.
        #   For Paraview visualisation, or to take advantage
        #   of newer updates to VTK, change the first line to
        #   the newest release of VTK DataFile version, but
        #   for VisIT visualisation you must change the
        #   version number to 2.0 as VisIT has not been
        #   updated in some time.
        with open(f"{tag}.vtk", "w") as f:

            f.write("# vtk DataFile Version 4.2\n")
            f.write(
                "This is a model produced by the "
                "OLCAO program suite\n"
            )
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n\n")

            # First print all atom positions (including
            #   extended).
            f.write(
                f"POINTS {bd.num_atoms_ext_dx} float\n"
            )
            for j in range(1, bd.num_atoms_ext_dx + 1):
                f.write(
                    f"{bd.x_pos_ext_dx[j]} "
                    f"{bd.y_pos_ext_dx[j]} "
                    f"{bd.z_pos_ext_dx[j]}\n"
                )
            f.write("\n")

            # Now list the connections between data points
            #   as VERTICES.  The header must include the
            #   number of vertices followed by the total
            #   number of data points used to express them.
            total_vert = 2 * bd.num_atoms_ext_dx
            f.write(
                f"VERTICES {bd.num_atoms_ext_dx} "
                f"{total_vert}\n"
            )
            for j in range(1, bd.num_atoms_ext_dx + 1):
                # We need to attach a vertex to each
                #   point/atom.  Since C arrays start at
                #   zero, we must start from there.
                #   Attempts to use in Paraview can
                #   sometimes be "hit or miss" when
                #   applying glyphs to vertices.  If
                #   spherical glyphs are not applying to
                #   one or more vertices, change this next
                #   print line to a number higher than "1"
                #   as a quick fix.
                k = j - 1
                f.write(f"1 {k}\n")
            f.write("\n")

            # Now list the LINES (bonds).  The header
            #   must include the number of lines, followed
            #   by the total number of data points used.
            #   Since we only use start and endpoints, that
            #   is 3x the number of bonds.
            total_bond_pts = 3 * bd.num_bonds_dx
            f.write(
                f"LINES {bd.num_bonds_dx} "
                f"{total_bond_pts}\n"
            )
            for j in range(1, bd.num_bonds_dx + 1):
                # Subtract 1 because, like openDX, VTK is
                #   C-based and array refs start from zero.
                #   Since we only use start and end points,
                #   the first number is always 2.  If more
                #   points are needed (e.g. a midpoint for
                #   a cylinder), the number would be 3.
                fa = bd.bond_atom_dx[j][1] - 1
                sa = bd.bond_atom_dx[j][2] - 1
                f.write(f"2 {fa} {sa}\n")
            f.write("\n")

            # ---- POINT_DATA section ----
            # This section adds the specified colour table
            #   to the bonds and atoms.  The colour tables
            #   for VTK may be left as default, or you may
            #   construct a colour table and specify how
            #   it's used.  Each entry in the lookup table
            #   is a rgba array (red green blue alpha)
            #   where alpha is transparency.  To save time,
            #   the default colour table is used.  In the
            #   future a colour table should be generated
            #   in each file, specifying what colour is
            #   associated with what type of atom or what
            #   strength of bond.
            f.write(
                f"POINT_DATA {bd.num_atoms_ext_dx}\n"
            )
            f.write("SCALARS sample_scalars float\n")
            f.write("LOOKUP_TABLE default\n")

            # Scaling factor for spherical glyphs based
            #   on atomic radii.
            f.write("SCALARS atomic_radius float\n")
            f.write("LOOKUP_TABLE default\n")
            for j in range(1, bd.num_atoms_dx + 1):
                f.write(f"{bd.atom_radius_dx[j]}\n")

            # Charge transfer data as point data.
            f.write(
                "SCALARS charge_transfer float\n"
            )
            f.write("LOOKUP_TABLE default\n")
            for j in range(1, bd.num_atoms_dx + 1):
                f.write(
                    f"{bd.atom_charge_transfer_dx[j]}"
                    "\n"
                )

        # Release memory.
        bd.bond_order_dx = [None]
        bd.abc_bond_pos = [None]
        bd.xyz_bond_pos = [None]
        bd.bond_atom_dx = [None]
        bd.atom_radius_dx = [None]
        bd.atom_color_dx = [None]
        bd.atom_charge_transfer_dx = [None]
        bd.x_pos_ext_dx = [None]
        bd.y_pos_ext_dx = [None]
        bd.z_pos_ext_dx = [None]


# ============================================================
#  Main program
# ============================================================

def main():
    """Main execution flow for makeBOND.py.

    The overall flow is:
    1. Initialise environment (element data, defaults).
    2. Parse command-line parameters.
    3. Read system structure (lattice, atom positions).
    4. Read control files (bond filtering, grouping).
    5. Read raw bond order data.
    6. Compute implicit data (atom/bond tags).
    7. Produce the requested output (scatter, model, mesh,
       profile, or VTK).
    """

    print("\n\nScript Executing.\n")

    # Get script settings from the resource control file and
    #   command-line parameters.  This also initialises the
    #   environment and records the command line used.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(f"\nInitialising the environment at..........{ts}.")
    settings = ScriptSettings()

    # Create the data container.
    bd = BondData(settings)

    # Get the system structure: Lattice constants, atom types,
    #   and atom positions.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(
        f"\nGetting system data at.................."
        f".{ts}."
    )
    get_system_structure(bd)

    # Read the control files if any were provided.  If not,
    #   create the default control definitions.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(
        f"\nReading control files at................"
        f".{ts}."
    )
    read_control_files(bd)

    # Read in the raw bond order data and store the results
    #   in convenient arrays.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(
        f"\nReading in raw data at.................."
        f".{ts}."
    )
    read_data(bd)

    # Compute the implicit data for each atom and bond.
    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(
        f"\nComputing the implicit data at.........."
        f".{ts}."
    )
    implicit_data(bd)

    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")

    # Print the results in the command-line-requested form.
    if settings.scatter_plot:
        # Print the results into multiple files in a form that
        #   can be imported by the associated Origin script (or
        #   some other program) for plotting of the bond order
        #   and charge transfer results as a scatter plot.
        print(
            f"\nPrinting scatter plot data at..........."
            f".{ts}."
        )
        print_scatter_data(bd)

    elif settings.system_model_dx:
        # Print the results in a form usable by the openDX
        #   program for drawing a structure model with bonds
        #   coloured according to bond order value, and atoms
        #   drawn with colours corresponding to charge transfer
        #   value, etc.
        print(
            f"\nPrinting openDX model data at..........."
            f".{ts}."
        )
        print_model_data(bd)

    elif settings.mesh_dx:
        # Print the results in the form of an openDX 3D mesh.
        #   Each mesh point contains the weighted sum of
        #   contributions from nearby bonds.  The results are
        #   the bond order results.
        print(
            f"\nPrinting openDX bond order mesh at......"
            f".{ts}."
        )
        print_3d_mesh_data(bd)

    elif settings.bond_profile_plot:
        # Print the results into one file to be plotted as a
        #   curve indicating the bond profile along the command
        #   line specified axis.
        print(
            f"\nPrinting bond profile data at..........."
            f".{ts}."
        )
        print_bond_profile(bd)

    elif settings.system_model_vtk:
        # Print the results in a form usable by Paraview/VisIT
        #   for drawing a structural model with (in the future)
        #   bonds coloured according to bond order value, and
        #   atoms drawn with colours corresponding to charge
        #   transfer value, etc.
        print(
            f"\nPrinting VTK Legacy model data at......."
            f".{ts}."
        )
        print_model_vtk(bd)

    ts = datetime.now().strftime("%b %d, %Y: %H:%M:%S")
    print(
        f"\nScript complete at......................"
        f".{ts}.\n"
    )


if __name__ == '__main__':
    # Everything before this point was a class/function
    #   definition or an import.  Only now do we actually
    #   start running the program.  This structure allows
    #   another Python program to import this script and
    #   call its functions internally.
    main()
