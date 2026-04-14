#!/usr/bin/env python3

"""make_reactions.py -- Create reaction templates for LAMMPS bond/react.

PROGRAM: make_reactions.py
PURPOSE: To create a set of pre- and post-reaction molecule files
   along with a mapping file for a bond/react based LAMMPS
   calculation.

ALGORITHM:
   This program requires two directories as input that each
   identify a precursor molecule type. It will then create a new
   directory based on the combination of their names. Inside that
   new directory, it will copy OLCAO skeleton input files from the
   originally identified directories and will then process them to
   create pre- and post-reaction files along with a mapping file.

   Each input skeleton file is parsed to find and record all of a
   "surface" atom (a.k.a., an S atom, not to be confused with
   sulfur) defined by a command line parameter. The script then
   identifies the set of atoms attached to each S atom that should
   be included in the list of the atoms to be deleted upon
   molecule bonding. That list is called the "delete list". Unless
   specified otherwise it is assumed that only the individual S
   atoms or perhaps just one of the S atoms will be part of the
   delete list.

   If it is desired to have more atoms in the delete list then the
   -d1 and -d2 options can be used to delete more atoms from
   molecules 1 and 2 respectively. In the default assumption the
   so-called "trigger atoms" that are bonded to the S atoms will
   be the points between the two molecules at which new bonds will
   be identified and formed. The physical vision is that the S
   atoms will break off to form S_2 (or whatever) and the
   remaining bare (trigger) atoms will bond. (In the simulation,
   the S_2 is never actually created. Instead the two S atoms are
   just deleted when the molecules bond.) For larger tree depths,
   the whole tree is removed. That story is modified as expected
   if the option to keep one S atom and have it serve as the
   bridge between the two trigger atoms is used.

   Once all the atoms in both molecules have been cataloged, then
   the second molecule will be systematically translated and
   oriented so that each of its S atoms would be in a proper
   position for reacting with each of the S atoms in the first
   molecule. Once so translated and rotated, for a given S, then
   the skeleton files for molecule1 and the modified molecule2
   are merged into a single file.

   At this point two product files are created that are suitable
   "molecule" files for use in LAMMPS bond/react calculations. The
   first is the pre-reaction molecule file that contains relevant
   segments of both molecules (where "relevant" is determined by
   the chain-length variable which is 4 by default so as to
   include the atom bonded to an S atom, the next nearest neighbor
   atoms, and the next-next nearest neighbor atoms). The second is
   the post-reaction molecule. We use a default chain length of 4
   because we will change angles to the atom bonded to the S atom.
   Those S-bonded atoms might be part of an angle that goes "two
   bonds deeper" into the molecule. (E.g., Consider the molecule:
   S-B-C-X. Because we want to modify the B atom and we need to
   account for angles, we must include the X atom in the template
   because of a possible B-C-X bond angle.)

   The atomic structure (number of atoms and their coordinates)
   are the same in both files. The difference is in the
   information about their bonding patterns (topology). When it is
   time for both the pre- and post-reaction LAMMPS files to be
   created, the OLCAO skeleton file containing the merged pair of
   molecules is read, bonds and bond angles are analyzed, and a
   subset of the original atoms from the skeleton file are taken
   according to a chain-length variable.

REQUIREMENTS:
   - Python 3.6+
   - structure_control.py (StructureControl class)
   - element_data.py (ElementData class)
   - bond_analysis.py (BondAnalysis class)
   - mod_struct.py (for molecular transformations)
   - $OLCAO_DATA/angles.dat (Hooke angle database)
   - $OLCAO_RC/make_reactionsrc.py (default parameters)

USAGE:
   make_reactions.py -m1 DIR1 -m2 DIR2
                    [-c1 CHAIN_LEN1] [-c2 CHAIN_LEN2]
                    [-d1 DEPTH TARGET] [-d2 DEPTH TARGET]
                    [-s ELEMENT DIST DEL_BOTH]
                    [-uptypes STEP]
                    [-consttypes]
                    [-help]

This is the Python port of the original Perl makeReactions script.
"""

import argparse as ap
import math
import os
import shutil
import subprocess
import sys
from collections import deque
from datetime import datetime


# ================================================================
# AngleData: reads the Hooke angle database from angles.dat.
#
# This is a Python port of the Perl AngleData.pm module. It reads
# the file $OLCAO_DATA/angles.dat and stores the Hooke angle
# coefficients used for matching computed bond angles to known
# angle types.
# ================================================================
class AngleData:
    """Read and provide access to the Hooke angle coefficient
    database.

    The angles.dat file contains entries of the form:
        Z1  Z_vertex  Z2  k_spring  rest_angle  tolerance
    where Z values are atomic numbers, k_spring is the spring
    constant, rest_angle is the equilibrium angle in degrees, and
    tolerance is the allowed deviation in degrees for matching.

    After init_angle_data() is called, the data is available
    through the hooke_angle_coeffs list (1-indexed) and
    num_hooke_angles.
    """

    def __init__(self):
        """Initialize the AngleData object."""
        self.num_hooke_angles = 0
        # 1-indexed: hooke_angle_coeffs[i] is a list of 7 entries:
        #   [None, Z1, Z_vertex, Z2, k_spring, rest_angle, tol]
        self.hooke_angle_coeffs = [None]
        self._init_done = False

    def init_angle_data(self):
        """Read the angle data file and populate the
        hooke_angle_coeffs array.

        The data file is located at $OLCAO_DATA/angles.dat and
        has the following format:

            NUM_HOOKE_ANGLES
            <count>
            HOOKE_ANGLE_COEFFS
            Z1  Z_vertex  Z2  k  rest_angle  tolerance
            ...

        Each entry describes one type of bond angle:
        - Z1, Z_vertex, Z2: atomic numbers of the three atoms
          forming the angle (Z1 <= Z2 by convention).
        - k: Hooke spring constant for this angle type.
        - rest_angle: equilibrium angle in degrees.
        - tolerance: allowed deviation (degrees) for matching
          a computed angle to this type.
        """
        if self._init_done:
            return
        self._init_done = True

        data_dir = os.environ.get('OLCAO_DATA', '')
        angle_file = os.path.join(data_dir, 'angles.dat')

        with open(angle_file) as f:
            # Read the NUM_HOOKE_ANGLES tag and count.
            line = self._read_next(f)
            if line[0] != 'NUM_HOOKE_ANGLES':
                sys.exit("Expecting NUM_HOOKE_ANGLES tag "
                         f"in {angle_file}")
            line = self._read_next(f)
            self.num_hooke_angles = int(line[0])

            # Read the HOOKE_ANGLE_COEFFS tag.
            line = self._read_next(f)
            if line[0] != 'HOOKE_ANGLE_COEFFS':
                sys.exit("Expecting HOOKE_ANGLE_COEFFS tag "
                         f"in {angle_file}")

            # Read each angle coefficient entry.
            for _ in range(self.num_hooke_angles):
                line = self._read_next(f)
                # Store as 1-indexed list:
                #   [None, Z1, Z_vertex, Z2, k, angle, tol]
                entry = [
                    None,
                    int(line[0]),    # Z1
                    int(line[1]),    # Z_vertex
                    int(line[2]),    # Z2
                    float(line[3]),  # k spring constant
                    float(line[4]),  # rest angle (degrees)
                    float(line[5]),  # tolerance (degrees)
                ]
                self.hooke_angle_coeffs.append(entry)

    @staticmethod
    def _read_next(fh):
        """Read and split the next non-empty line from a file
        handle, stripping whitespace."""
        while True:
            raw = fh.readline()
            if not raw:
                sys.exit("Unexpected end of angle data file.")
            parts = raw.strip().split()
            if parts:
                return parts


# ================================================================
# ScriptSettings: command line parsing and defaults.
# ================================================================
class ScriptSettings:
    """The instance variables of this object are the user settings
    that control the program. The variable values are pulled from
    a list that is created within a resource control file and that
    are then reconciled with command line parameters.

    Following the XYZ.py pattern, defaults come from the rc file
    ($OLCAO_RC/make_reactionsrc.py or ./make_reactionsrc.py) and
    are then overridden by any command line arguments.
    """

    def __init__(self):
        """Define default values by pulling them from the resource
        control file in the default location:
        $OLCAO_RC/make_reactionsrc.py or from the current working
        directory if a local copy is present."""

        # Read default variables from the resource control file.
        rc_dir = os.environ.get('OLCAO_RC')
        if not rc_dir:
            sys.exit("Error: $OLCAO_RC is not set. "
                     "See installation instructions.")
        sys.path.insert(1, rc_dir)
        from make_reactionsrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values from the rc defaults file.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc file.
        self.reconcile(args)

        # Record the command line parameters that were used.
        self.record_clp()

    def assign_rc_defaults(self, rc):
        """Assign default values from the resource control
        dictionary to instance variables.

        Parameters
        ----------
        rc : dict
            Dictionary from make_reactionsrc.parameters_and_defaults
        """

        # Chain length parameters: number of bond-hops from the
        #   S atom to include in the reaction template for each
        #   molecule.
        self.chain_len1 = rc["chain_len1"]
        self.chain_len2 = rc["chain_len2"]

        # Delete tree parameters: depth and target element for
        #   atoms to remove upon binding.
        self.del_tree_depth1 = rc["del_tree_depth1"]
        self.del_tree_depth2 = rc["del_tree_depth2"]
        self.del_tree_target1 = rc["del_tree_target1"]
        self.del_tree_target2 = rc["del_tree_target2"]

        # Surface atom parameters.
        self.s_element = rc["s_element"]
        self.ss_dist = rc["ss_dist"]
        self.del_both_s = rc["del_both_s"]

        # Type numbering parameters.
        self.type_step = rc["type_step"]
        self.const_types = rc["const_types"]

        # Input directories (required; no defaults).
        self.in_dir1 = None
        self.in_dir2 = None

    def parse_command_line(self):
        """Create the argparse parser and parse the command line.

        Returns
        -------
        argparse.Namespace
            Parsed command line arguments.
        """
        prog_name = "make_reactions.py"

        description_text = """\
make_reactions.py -- Create reaction templates for LAMMPS
bond/react.

Creates pre- and post-reaction molecule template files and a
mapping file from two input molecule directories that each contain
an OLCAO skeleton (.skl) file. The templates are suitable for use
with the bond/react fix in LAMMPS.

IMPORTANT: Each molecule skeleton file must use space group 1_a
(all atomic coordinates specified explicitly) and must use
Cartesian coordinates.
"""

        epilog_text = """\
Defaults are given in ./make_reactionsrc.py or
$OLCAO_RC/make_reactionsrc.py.
"""

        parser = ap.ArgumentParser(
            prog=prog_name,
            formatter_class=ap.RawDescriptionHelpFormatter,
            description=description_text,
            epilog=epilog_text,
        )

        self._add_parser_arguments(parser)
        return parser.parse_args()

    def _add_parser_arguments(self, parser):
        """Add all command line arguments to the parser.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            The parser to add arguments to.
        """
        # -m1 and -m2: molecule directories (required).
        parser.add_argument(
            '-m1', dest='in_dir1', type=str,
            required=True,
            help=(
                "Directory containing the first molecule's "
                "OLCAO skeleton file. Required."
            ),
        )

        parser.add_argument(
            '-m2', dest='in_dir2', type=str,
            required=True,
            help=(
                "Directory containing the second molecule's "
                "OLCAO skeleton file. Required."
            ),
        )

        # -c1 and -c2: chain lengths.
        # The chain length variable defines the number of atoms
        #   (starting at the point of "contact" between the two
        #   molecules) to include in the reaction template from
        #   each molecule. Consider the molecules CH4 and
        #   B10C2H12 with hydrogen serving as the S atoms. The
        #   contact points are an S atom bound to the C atom from
        #   CH4 and an S atom bound to a B or C from B10C2H12.
        #   Then -c1 of 3 for CH4 would include the S atom bound
        #   to the central C, the central C itself, and all three
        #   other H (S) atoms. A -c2 of 2 would include the S
        #   atom bound to B or C and then the B or C atom itself.
        parser.add_argument(
            '-c1', dest='chain_len1', type=int,
            default=self.chain_len1,
            help=(
                "Chain length for molecule 1: number of "
                "bond-hops from the S atom to include in the "
                "template. "
                f"Default: {self.chain_len1}"
            ),
        )

        parser.add_argument(
            '-c2', dest='chain_len2', type=int,
            default=self.chain_len2,
            help=(
                "Chain length for molecule 2. "
                f"Default: {self.chain_len2}"
            ),
        )

        # -d1 and -d2: delete tree depth and target.
        # The delete tree depth variable is used to remove atoms
        #   from the molecules upon binding. Consider a molecule:
        #
        #       H       H
        #       |       |
        #  H  H-C-H   H-C-H  H
        #  |    |       |    |
        #H-C--Si---N---Si---C-H
        #  |    |       |    |
        #  H  H-C-H   H-C-H  H
        #       |       |
        #       H       H
        #
        #   If "-d1 1 si_1" is used then the far left CH3 group
        #   will be added to the delete tree and the left Si will
        #   be the binding point. If "-d1 2 n_1" is used then all
        #   atoms to the left of the N will be in the delete tree
        #   and the N will be the binding point.
        parser.add_argument(
            '-d1', dest='d1', nargs=2,
            metavar=('DEPTH', 'TARGET'),
            default=None,
            help=(
                "Delete tree for molecule 1: DEPTH (int) and "
                "TARGET (element_species, e.g. 'si_1'). "
                "Default: depth=0 (only S atom deleted)"
            ),
        )

        parser.add_argument(
            '-d2', dest='d2', nargs=2,
            metavar=('DEPTH', 'TARGET'),
            default=None,
            help=(
                "Delete tree for molecule 2: DEPTH (int) and "
                "TARGET (element_species, e.g. 'n_1'). "
                "Default: depth=0 (only S atom deleted)"
            ),
        )

        # -s: surface atom element, S-S distance, and del_both_s.
        # The surface atom is the element from the periodic table
        #   that will serve as the "surface" atom. The separation
        #   distance is the distance between the two S atoms in
        #   the templates. del_both_s controls whether one or both
        #   S atoms will be deleted upon reacting: 1 = delete
        #   both; 0 = delete only one (the one from molecule 2).
        parser.add_argument(
            '-s', dest='s_params', nargs=3,
            metavar=('ELEMENT', 'DIST', 'DEL_BOTH'),
            default=None,
            help=(
                "Surface atom parameters: element name, S-S "
                "distance, and whether to delete both S atoms "
                "(1) or only one (0). "
                f"Default: {self.s_element} {self.ss_dist} "
                f"{self.del_both_s}"
            ),
        )

        # -uptypes: increment type numbers for bonded atoms.
        # By default type numbers assigned to atoms in molecule
        #   skeleton files remain constant throughout the
        #   simulation (type_step=0). When this option is given,
        #   atoms that bind are assigned a type number that is
        #   higher by type_step than the highest existing type for
        #   that element.
        parser.add_argument(
            '-uptypes', dest='type_step', type=int,
            default=self.type_step,
            help=(
                "Increment in type number for bonded atoms "
                "and their neighbours. "
                f"Default: {self.type_step}"
            ),
        )

        # -consttypes: keep type numbers constant when merging.
        # When given, the type numbers of each molecule remain
        #   unchanged when combined. E.g., C with type 1 on
        #   molecule 1 and C with type 1 on molecule 2 will both
        #   remain type 1 in the merged system.
        parser.add_argument(
            '-consttypes', dest='const_types',
            action='store_true',
            default=self.const_types,
            help=(
                "Keep atom types constant from precursor "
                "files when merging. "
                f"Default: {self.const_types}"
            ),
        )

    def reconcile(self, args):
        """Reconcile the command line arguments with the rc file
        defaults. Command line values take precedence.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed command line arguments.
        """
        self.in_dir1 = args.in_dir1
        self.in_dir2 = args.in_dir2
        self.chain_len1 = args.chain_len1
        self.chain_len2 = args.chain_len2
        self.type_step = args.type_step
        self.const_types = args.const_types

        # Handle -d1 and -d2 (each takes depth + target).
        if args.d1 is not None:
            self.del_tree_depth1 = int(args.d1[0])
            self.del_tree_target1 = args.d1[1]
        if args.d2 is not None:
            self.del_tree_depth2 = int(args.d2[0])
            self.del_tree_target2 = args.d2[1]

        # Handle -s (element, distance, del_both_s).
        if args.s_params is not None:
            self.s_element = args.s_params[0].lower()
            self.ss_dist = float(args.s_params[1])
            self.del_both_s = int(args.s_params[2])

    def record_clp(self):
        """Record the command line parameters to a 'command'
        file for reproducibility."""
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


# ================================================================
# MakeReactions: the main workhorse class.
# ================================================================
class MakeReactions:
    """Orchestrates the creation of LAMMPS bond/react template
    files from two precursor molecule directories.

    This class encapsulates all the data structures and
    subroutines needed to:
    1. Prepare molecules (center them and unify cell sizes).
    2. Catalog surface (S) atoms and their bonding neighbours.
    3. Generate merged molecule files by translating and
       rotating molecule 2 to align with molecule 1.
    4. Catalog element names, species IDs, and molecule
       identities.
    5. Prune the merged molecules to the chain-length depth.
    6. Build pre- and post-reaction LAMMPS template files
       and the associated mapping file.
    """

    def __init__(self, settings):
        """Initialize the MakeReactions object with user
        settings and load required databases.

        Parameters
        ----------
        settings : ScriptSettings
            The reconciled user settings from the command line
            and the rc file.
        """
        self.s = settings  # Shorthand for settings.

        # Physical constants.
        self.PI = math.pi
        self.EPSILON = 1e-5

        # Import the structure_control module.
        from structure_control import StructureControl
        self.StructureControl = StructureControl

        # Import the element_data module.
        from element_data import ElementData
        self.ed = ElementData()
        self.ed.init_element_data()

        # Initialize the angle database.
        self.angle_data = AngleData()
        self.angle_data.init_angle_data()

        # ---- Molecule-level data ----
        self.in_file1 = ""
        self.in_file2 = ""
        self.mol_name1 = ""
        self.mol_name2 = ""
        self.in_file_cent1 = ""
        self.in_file_cent2 = ""
        self.rxn_dir = ""
        self.num_atoms1 = 0
        self.num_atoms2 = 0

        # ---- S-atom data (per-molecule) ----
        # These are 1-indexed by S-atom number within each mol.
        self.num_s_atoms1 = 0
        self.num_s_atoms2 = 0
        # Coordinates: [s_atom][axis] (axis 1..3).
        self.s_atom_coords1 = [None]
        self.s_atom_coords2 = [None]
        self.trigger_coords1 = [None]
        self.trigger_coords2 = [None]
        # Atom numbers within the individual molecules.
        self.s_atom_number1 = [None]
        self.s_atom_number2 = [None]
        self.trigger_atom_number1 = [None]
        self.trigger_atom_number2 = [None]
        # Binding atom info.
        self.bonded_atom_number1 = [None]
        self.bonded_atom_number2 = [None]
        self.bonded_elem_name1 = [None]
        self.bonded_elem_name2 = [None]
        self.bonded_species1 = [None]
        self.bonded_species2 = [None]
        # Delete trees.
        self.num_del_atoms1 = [None]
        self.num_del_atoms2 = [None]
        self.del_tree1 = [None]
        self.del_tree2 = [None]

        # ---- Merged-molecule data ----
        self.merged_s_atom_number1 = [None]
        self.merged_s_atom_number2 = [None]
        # [s1][s2][atom] element name.
        self.merged_atom_elem_name = {}
        # [s1][s2][atom][phase] species ID (phase 1=pre,2=post).
        self.merged_atom_species_id = {}
        # [s1][s2][atom] molecule name.
        self.merged_atom_mol_name = {}
        # [s1][s2][del_idx] atom numbers in delete tree.
        self.merged_del_tree = {}
        # [s1][s2] count of atoms to delete.
        self.merged_num_del_atoms = {}

        # ---- Pruned-molecule data ----
        self.pruned_s_atom_number1 = {}
        self.pruned_s_atom_number2 = {}
        self.pruned_trigger_atom_number1 = {}
        self.pruned_trigger_atom_number2 = {}
        self.pruned_bonded_atom_number1 = {}
        self.pruned_bonded_atom_number2 = {}
        self.pruned_bonded_elem_name1 = {}
        self.pruned_bonded_elem_name2 = {}
        self.pruned_bonded_species1 = {}
        self.pruned_bonded_species2 = {}
        self.pruned_atom_elem_name = {}
        self.pruned_atom_species_id = {}
        self.pruned_atom_mol_name = {}
        self.pruned_del_tree = {}
        self.pruned_num_atoms = {}
        self.pruned_num_atoms1 = {}
        self.pruned_num_atoms2 = {}

        # ---- Edge data ----
        # edge_id[id][s1][s2] = atom number of edge atom.
        self.edge_id = {}
        self.pruned_edge_id = {}
        self.num_edge_ids = {}

        # ---- Keep-atom mask for pruning ----
        self.keep_atom = []

    def run(self):
        """Execute the full makeReactions pipeline.

        This is the main entry point that orchestrates the
        entire workflow:

        1. Finalize the environment (arrange molecules,
           create directories).
        2. Prepare molecules (center and unify cell sizes).
        3. Catalog S atoms and delete trees for both
           molecules.
        4. Generate merged molecule files (translate/rotate
           molecule 2 to align with molecule 1).
        5. Catalog element/species/molecule IDs for merged
           molecules.
        6. Prune merged molecules to the chain-length depth.
        7. Build pre- and post-reaction LAMMPS templates and
           the mapping file.
        """
        # Finalize the environment on the basis of required
        #   input parameters.
        self.finalize_env()

        # Prepare the molecules by creating skeleton files
        #   with the molecules centered in each of their
        #   respective cells and by ensuring that each cell
        #   is sufficiently large to accommodate both
        #   molecules at their fullest possible extension in
        #   all directions and orientations.
        self.prep_molecules()

        # Catalog the positions of the S atoms, their bound
        #   neighbours for alignment and trigger purposes upon
        #   merging and binding, and the atoms in the delete
        #   trees in both molecules which includes the specific
        #   binding atom and its key information.
        self.catalog_s(
            1, self.in_file_cent1,
            self.s_atom_coords1, self.trigger_coords1,
            self.s_atom_number1, self.trigger_atom_number1,
            self.bonded_atom_number1, self.bonded_elem_name1,
            self.bonded_species1,
            self.s.del_tree_depth1,
            self.s.del_tree_target1,
            self.num_del_atoms1, self.del_tree1,
        )
        self.catalog_s(
            2, self.in_file_cent2,
            self.s_atom_coords2, self.trigger_coords2,
            self.s_atom_number2, self.trigger_atom_number2,
            self.bonded_atom_number2, self.bonded_elem_name2,
            self.bonded_species2,
            self.s.del_tree_depth2,
            self.s.del_tree_target2,
            self.num_del_atoms2, self.del_tree2,
        )

        # Generate a set of merged molecule files with
        #   different orientations.
        self.generate_merged_mols()

        # Catalog the element names, species numbers, and
        #   molecule names for each member of the set of
        #   merged molecules. Additionally, identify the
        #   species number of the binding atoms and make
        #   them higher than the highest species number for
        #   that element so that they are distinguishable
        #   when combined into the full LAMMPS simulation
        #   with all other molecules.
        self.catalog_elem_spec_mol_id()

        # Remove atoms that are far from the bond reaction
        #   point.
        self.prune_merged_molecule()

        # Build LAMMPS bond/react templates from the pruned
        #   olcao.skl files. There is one phase to a complete
        #   reaction: the atoms identified as binding atoms
        #   create a bond between themselves and the atoms to
        #   be deleted are removed. Only one pre-/post-reaction
        #   template pair is needed.
        self.make_reaction_templates(1)  # Pre-reaction.
        self.make_reaction_templates(2)  # Post-reaction.

    # ============================================================
    # finalize_env
    # ============================================================
    def finalize_env(self):
        """Finalize the environment based on the required input
        parameters.

        Activities:
        - Sort the two molecule directories alphanumerically
          so that the -m1 molecule is always before -m2.
        - Create the reaction directory.
        - Establish input skeleton file names.
        - Copy skeleton files into the reaction directory.
        - Change into the reaction directory.
        """
        # Arrange the order of the molecules alphabetically.
        if self.s.in_dir1 >= self.s.in_dir2:
            self.s.in_dir1, self.s.in_dir2 = (
                self.s.in_dir2, self.s.in_dir1
            )

        # Define the reaction directory name.
        self.rxn_dir = (
            self.s.in_dir1 + "__" + self.s.in_dir2
        )

        # Make the reaction directory (where we will stay for
        #   the remainder of the script execution).
        os.makedirs(self.rxn_dir, exist_ok=True)

        # Establish variables that refer to the input skeleton
        #   files.
        self.in_file1 = self.s.in_dir1 + ".skl"
        self.in_file2 = self.s.in_dir2 + ".skl"

        # The molecule names are the same as the directory
        #   names.
        self.mol_name1 = self.s.in_dir1
        self.mol_name2 = self.s.in_dir2

        # Copy the molecule skeleton files into the reaction
        #   directory. (Note that if in_file1 and in_file2 are
        #   actually the same, then the second one will be
        #   copied over the first. This is okay.)
        src1 = os.path.join(
            self.s.in_dir1, self.in_file1
        )
        src2 = os.path.join(
            self.s.in_dir2, self.in_file2
        )
        shutil.copy2(src1, self.rxn_dir)
        shutil.copy2(src2, self.rxn_dir)

        # Enter into the reaction template directory.
        os.chdir(self.rxn_dir)

    # ============================================================
    # prep_molecules
    # ============================================================
    def prep_molecules(self):
        """Prepare molecules by centering them and unifying their
        simulation cell sizes.

        The sequence of activity here is a bit tricky and
        convoluted because there are limits to what kinds of
        information we have access to at each stage. Also, we
        don't want to write a lot of custom code if we can get
        away with just using existing subroutines from
        StructureControl.

        At the end of this subroutine, we will want each molecule
        to be in a new skeleton file. The molecules will be
        centered in the simulation cells. Also, the sizes of the
        simulation cells will be such that the full breadth of
        both molecules in any direction can be contained.

        We proceed by reading each molecule, getting the maximum
        cell magnitudes from each, and then making the final cell
        cubic with sides that are a sum of those two maxima.
        """
        # Create a directory to store the centered results in.
        os.makedirs("centered", exist_ok=True)

        # Molecule 1 first. Read the file, get the magnitudes
        #   of the lattice vectors, and find the largest.
        sc1 = self.StructureControl()
        sc1.read_input_file(self.in_file1, use_file_species=True)
        max_side1 = max(sc1.mag[1], sc1.mag[2], sc1.mag[3])

        # Molecule 2 second. Read the file, get the magnitudes
        #   of the lattice vectors, and find the largest.
        sc2 = self.StructureControl()
        sc2.read_input_file(self.in_file2, use_file_species=True)
        max_side2 = max(sc2.mag[1], sc2.mag[2], sc2.mag[3])

        total_buffer = max_side1 + max_side2

        # Now we can do some actual work on the molecules.

        # Read the first olcao skeleton file and then compute
        #   the "crystal" parameters which (as a side effect)
        #   will center the molecule. (Note that the original
        #   purpose of computeCrystalParameters was to deal
        #   with input files that do not contain lattice
        #   parameters.)
        sc1_work = self.StructureControl()
        sc1_work.read_input_file(
            self.in_file1, use_file_species=True
        )
        sc1_work.set_buffer(total_buffer)
        sc1_work.compute_crystal_parameters()

        # Print the first molecule to a center-shifted skeleton
        #   file. The new name is the old one with "_cent" before
        #   the .skl extension.
        base1 = self.in_file1.rsplit('.', 1)[0]
        self.in_file_cent1 = base1 + "_cent.skl"
        title1 = ''.join(sc1_work.system_title[1:])
        sc1_work.print_olcao(
            filename=os.path.join(
                "centered", self.in_file_cent1
            ),
            title=title1,
            style="cartesian",
        )

        # Repeat the same procedure for the second molecule.
        sc2_work = self.StructureControl()
        sc2_work.read_input_file(
            self.in_file2, use_file_species=True
        )
        sc2_work.set_buffer(total_buffer)
        sc2_work.compute_crystal_parameters()

        # Print the center-shifted second molecule.
        base2 = self.in_file2.rsplit('.', 1)[0]
        self.in_file_cent2 = base2 + "_cent.skl"
        title2 = ''.join(sc2_work.system_title[1:])
        sc2_work.print_olcao(
            filename=os.path.join(
                "centered", self.in_file_cent2
            ),
            title=title2,
            style="cartesian",
        )

    # ============================================================
    # catalog_s
    # ============================================================
    def catalog_s(
        self, mol_number, in_file,
        s_atom_coords, trigger_coords,
        s_atom_number, trigger_atom_number,
        bonded_atom_number, bonded_elem_name,
        bonded_species,
        del_tree_depth, del_tree_target,
        num_del_atoms, del_tree,
    ):
        """Catalog the positions of S atoms and their bonding
        neighbours in a molecule.

        This subroutine reads the centered skeleton file and uses
        bond_analysis.py to compute bond lengths. It then
        identifies all surface (S) atoms, records their positions
        and the positions of the atoms they are bonded to (trigger
        atoms), and builds the delete tree for each S atom.

        Parameters
        ----------
        mol_number : int
            1 or 2, indicating which molecule is being processed.
        in_file : str
            Filename of the centered skeleton file.
        s_atom_coords : list
            (Output) Coordinates of S atoms; [s_idx][axis].
        trigger_coords : list
            (Output) Coordinates of trigger atoms; [s_idx][axis].
        s_atom_number : list
            (Output) Atom numbers of S atoms.
        trigger_atom_number : list
            (Output) Atom numbers of trigger atoms.
        bonded_atom_number : list
            (Output) Atom numbers of binding atoms.
        bonded_elem_name : list
            (Output) Element names of binding atoms.
        bonded_species : list
            (Output) Species numbers of binding atoms.
        del_tree_depth : int
            Depth of the delete tree.
        del_tree_target : str
            Element_species target (e.g. "si_1").
        num_del_atoms : list
            (Output) Number of atoms to delete for each S atom.
        del_tree : list
            (Output) Atom numbers in the delete tree; [s][idx].
        """
        # Read the centered olcao skeleton file and extract
        #   key data from it.
        sc = self.StructureControl()
        centered_file = os.path.join("centered", in_file)
        sc.read_input_file(
            centered_file, use_file_species=True
        )
        num_atoms = sc.num_atoms

        # Request that the bond lengths between all atoms be
        #   computed using bond_analysis.py (Python replacement
        #   for the Perl bondAnalysis program).
        bond_file = f"bondAnalysis{mol_number}.bl"
        subprocess.run(
            [
                "bond_analysis.py",
                "-bf", "1.1",
                "-bl",
                "-i", centered_file,
                "-o", os.path.join("centered", bond_file),
            ],
            check=True,
        )

        # Read the bondAnalysis results to get bonding info.
        sc.read_bond_analysis_bl(
            os.path.join("centered", bond_file)
        )

        # Initialize the number of S atoms to zero. Then parse
        #   through the molecule data to get the number of S
        #   atoms and build the bonded atom data and delete list.
        num_s_atoms = 0

        for atom in range(1, num_atoms + 1):
            # Check if we found an S atom.
            if sc.atom_element_name[atom] != self.s.s_element:
                continue

            # If this S atom has more than one bond, then it is
            #   not permitted to be a surface atom.
            if sc.num_bonds[atom] > 1:
                continue

            # Increment count of the number of S atoms.
            num_s_atoms += 1

            # Store the atom number of this S atom.
            while len(s_atom_number) <= num_s_atoms:
                s_atom_number.append(None)
            s_atom_number[num_s_atoms] = atom

            # Initialize the count of the number of atoms to
            #   delete from this molecule for this S atom.
            while len(num_del_atoms) <= num_s_atoms:
                num_del_atoms.append(0)
            num_del_atoms[num_s_atoms] = 0

            # Ensure del_tree has space for this S atom index.
            while len(del_tree) <= num_s_atoms:
                del_tree.append([None])

            # Get the atom number of the trigger atom (also used
            #   for alignment). The trigger atom is the atom
            #   bonded to the S atom.
            while len(trigger_atom_number) <= num_s_atoms:
                trigger_atom_number.append(None)
            trigger_atom_number[num_s_atoms] = (
                sc.bonded[atom][1]
            )
            temp_bonded_atom_num = (
                trigger_atom_number[num_s_atoms]
            )

            # Store additional relevant information for this
            #   S-OtherAtom pair.
            while len(s_atom_coords) <= num_s_atoms:
                s_atom_coords.append([None, 0, 0, 0])
            while len(trigger_coords) <= num_s_atoms:
                trigger_coords.append([None, 0, 0, 0])
            for axis in range(1, 4):
                s_atom_coords[num_s_atoms][axis] = (
                    sc.direct_xyz[atom][axis]
                )
                trigger_coords[num_s_atoms][axis] = (
                    sc.direct_xyz[
                        trigger_atom_number[num_s_atoms]
                    ][axis]
                )

            # Ensure bonded_atom_number has space.
            while len(bonded_atom_number) <= num_s_atoms:
                bonded_atom_number.append(None)
            while len(bonded_elem_name) <= num_s_atoms:
                bonded_elem_name.append(None)
            while len(bonded_species) <= num_s_atoms:
                bonded_species.append(None)

            # Find all the atoms that need to be deleted from
            #   this molecule and record their atom numbers and
            #   associated information.
            #
            #   Note that presently, this is a weak algorithm
            #   that only deals with the depth == 0,1 cases.
            #   In the future, it would be good to have the
            #   depth and the target work together so that if
            #   a molecule had multiple atoms of the same target
            #   type we could plan to skip one to "get to" the
            #   other. For now though, that can't happen.
            if del_tree_depth == 0:
                # Add the S atom to the delete tree only if we
                #   are deleting both S atoms, OR, if we are
                #   only deleting one of the S atoms, then we
                #   only delete them from molecule 2.
                if self.s.del_both_s == 1:
                    num_del_atoms[num_s_atoms] += 1
                    idx = num_del_atoms[num_s_atoms]
                    while len(del_tree[num_s_atoms]) <= idx:
                        del_tree[num_s_atoms].append(None)
                    del_tree[num_s_atoms][idx] = atom

                    # In this case, the alignment atom and the
                    #   binding atom are the same atom.
                    bonded_atom_number[num_s_atoms] = (
                        temp_bonded_atom_num
                    )
                else:
                    # Only delete the S atom from molecule 2.
                    #   The binding atom for mol 2 is the same
                    #   as the trigger atom. The binding atom
                    #   for mol 1 is the S atom itself while
                    #   the alignment atom from mol 1 is the
                    #   trigger atom.
                    if mol_number == 2:
                        num_del_atoms[num_s_atoms] += 1
                        idx = num_del_atoms[num_s_atoms]
                        while (len(del_tree[num_s_atoms])
                               <= idx):
                            del_tree[num_s_atoms].append(None)
                        del_tree[num_s_atoms][idx] = atom
                        bonded_atom_number[num_s_atoms] = (
                            temp_bonded_atom_num
                        )
                    else:
                        bonded_atom_number[num_s_atoms] = atom
            else:
                # In this case we need to execute a traversal
                #   of the molecule in such a way that we add
                #   any atom we meet into the delete list until
                #   we find the target atom. At that point we
                #   do not add any atoms bonded to the target
                #   to the delete tree (except for the atom we
                #   "came from" of course). Then we continue
                #   the traversal until there are no more
                #   accessible points on the molecule.

                # Get the element name and species number of
                #   the target.
                parts = del_tree_target.lower().split('_')
                target_element = parts[0]
                target_species = int(parts[1])

                # Initialize: assume no atoms will be deleted.
                marked_atoms = [0] * (num_atoms + 1)

                # Add atoms to delete via recursive traversal.
                self._add_atoms_to_delete(
                    target_element, target_species,
                    num_del_atoms, del_tree,
                    sc.num_bonds, sc.bonded,
                    atom,
                    sc.atom_element_name,
                    sc.atom_species_id,
                    bonded_atom_number,
                    num_s_atoms, marked_atoms,
                )

                # Leave with temp_bonded_atom_num referring to
                #   the binding atom.
                temp_bonded_atom_num = (
                    bonded_atom_number[num_s_atoms]
                )

            # Get the element name of the binding atom.
            bonded_elem_name[num_s_atoms] = (
                sc.atom_element_name[temp_bonded_atom_num]
            )

            # Get the species number of the binding atom.
            bonded_species[num_s_atoms] = (
                sc.atom_species_id[temp_bonded_atom_num]
            )

        # Save the total number of atoms and number of S atoms
        #   for this molecule.
        if mol_number == 1:
            self.num_atoms1 = num_atoms
            self.num_s_atoms1 = num_s_atoms
        else:
            self.num_atoms2 = num_atoms
            self.num_s_atoms2 = num_s_atoms

    # ============================================================
    # _add_atoms_to_delete (recursive helper)
    # ============================================================
    def _add_atoms_to_delete(
        self, target_element, target_species,
        num_del_atoms, del_tree,
        num_bonds, bonded,
        current_atom,
        atom_element_name, atom_species_id,
        bonded_atom_number,
        num_s_atoms, marked_atoms,
    ):
        """Recursively traverse the molecular graph, adding atoms
        to the delete tree until the target element/species is
        found.

        When the target atom is reached, it is NOT added to the
        delete tree. Instead, it is recorded as the binding atom
        and no further traversal beyond it occurs. All other
        visited atoms are added to the delete tree.

        Parameters
        ----------
        target_element : str
            Lower-case element name of the target.
        target_species : int
            Species number of the target.
        num_del_atoms : list
            (In/Out) Count of delete-tree atoms for each S atom.
        del_tree : list
            (In/Out) Delete-tree atom numbers [s][idx].
        num_bonds : list
            Number of bonds for each atom.
        bonded : list
            Bonding list: bonded[atom][bond_idx] = partner.
        current_atom : int
            The atom currently being visited.
        atom_element_name : list
            Element names of all atoms.
        atom_species_id : list
            Species IDs of all atoms.
        bonded_atom_number : list
            (In/Out) Binding atom number for each S atom.
        num_s_atoms : int
            Current S atom index.
        marked_atoms : list
            (In/Out) Visited flags for each atom.
        """
        # Check if the current atom is the target. If so, then
        #   progress no further along the molecule, record this
        #   atom as the binding atom, and then permit the
        #   traversal to continue.
        if (atom_element_name[current_atom] == target_element
                and atom_species_id[current_atom]
                == target_species):
            bonded_atom_number[num_s_atoms] = current_atom
        else:
            # Increment the number of atoms to delete to include
            #   the current atom.
            num_del_atoms[num_s_atoms] += 1

            # Record this atom into the delete list.
            idx = num_del_atoms[num_s_atoms]
            while len(del_tree[num_s_atoms]) <= idx:
                del_tree[num_s_atoms].append(None)
            del_tree[num_s_atoms][idx] = current_atom

            # Record that this atom is going to be deleted.
            marked_atoms[current_atom] = 1

            # Visit each atom bonded to the current one and add
            #   it to the delete list as needed.
            for bond in range(1, num_bonds[current_atom] + 1):
                partner = bonded[current_atom][bond]
                if marked_atoms[partner] == 0:
                    self._add_atoms_to_delete(
                        target_element, target_species,
                        num_del_atoms, del_tree,
                        num_bonds, bonded,
                        partner,
                        atom_element_name, atom_species_id,
                        bonded_atom_number,
                        num_s_atoms, marked_atoms,
                    )

    # ============================================================
    # generate_merged_mols
    # ============================================================
    def generate_merged_mols(self):
        """Generate merged molecule files by translating and
        rotating molecule 2 to align with molecule 1.

        Repeatedly transforms the position of molecule 2 such
        that every combination of the four atoms comprising:
          - a bonded atom of molecule 1,
          - the associated S atom of molecule 1,
          - the S atom of molecule 2,
          - the bonded atom of molecule 2
        are all in a line in that listed order.

        The transformation consists of:
        1. A translation making the two S atoms coincident.
        2. A rotation bringing the bonded-to-S vectors of both
           molecules into anti-alignment.
        3. A final translation separating the S atoms by ss_dist.

        The two molecules are then printed into a single "merged"
        file that contains all atoms.

        This process is only applied to unique combinations of
        binding atom types. Duplicate pairs are skipped.

        A merged delete tree is also produced.
        """
        # Create a directory to store the merged results in.
        os.makedirs("merged", exist_ok=True)

        num_unique_bonded_names = 0
        unique_bonded_names = [None]

        for s1 in range(1, self.num_s_atoms1 + 1):
            # Some care must be taken if the S atoms and
            #   alignment atoms are co-linear with each other.
            #   That may happen easily if the two molecules are
            #   the same because we then naturally grab identical
            #   atoms.
            # If the atoms from the two molecules are colinear,
            #   then the axis of rotation cannot be defined by
            #   taking three of the atomic coordinates and using
            #   them to define a plane because they can only
            #   define a line. Thus, if such a case is detected,
            #   we will need to select a third point such that the
            #   angle of rotation is properly computed too.

            # The translation vector that will bring the two S
            #   atoms ss_dist apart from each other. Note that
            #   this is a constant for each s1 and its alignment
            #   pair. Note that this is the *second* translation
            #   that will be applied although we compute this
            #   translation vector first. That is the reason that
            #   we ultimately produce trans_vector2 instead of
            #   trans_vector1.
            # Store diff_vector1 in the 1-indexed [None, dx, dy, dz]
            # layout so it can be passed directly to the 1-indexed
            # vector helpers (get_vector_angle, etc.) further below.
            diff_vector1 = [None, 0.0, 0.0, 0.0]
            for axis in range(1, 4):
                diff_vector1[axis] = (
                    self.s_atom_coords1[s1][axis]
                    - self.trigger_coords1[s1][axis]
                )
            diff1_mag = math.sqrt(
                diff_vector1[1] ** 2 +
                diff_vector1[2] ** 2 +
                diff_vector1[3] ** 2
            )
            for axis in range(1, 4):
                diff_vector1[axis] /= diff1_mag
            # trans_vector2 stays 0-indexed: it is consumed immediately
            # below as a Cartesian shift by atom-building code that
            # indexes with ``ax + 1`` into other 1-indexed arrays.
            trans_vector2 = [
                diff_vector1[ax] * self.s.ss_dist
                for ax in range(1, 4)
            ]

            # Iterate over all S atoms of the second molecule.
            for s2 in range(1, self.num_s_atoms2 + 1):
                # Produce a merged delete tree.
                self.merged_del_tree[(s1, s2)] = [None]
                for i in range(
                    1, self.num_del_atoms1[s1] + 1
                ):
                    self.merged_del_tree[(s1, s2)].append(
                        self.del_tree1[s1][i]
                    )
                for i in range(
                    1, self.num_del_atoms2[s2] + 1
                ):
                    self.merged_del_tree[(s1, s2)].append(
                        self.del_tree2[s2][i]
                        + self.num_atoms1
                    )
                self.merged_num_del_atoms[(s1, s2)] = (
                    self.num_del_atoms1[s1]
                    + self.num_del_atoms2[s2]
                )

                # Produce a name identifying the binding pair.
                curr_bonded_name = (
                    f"{self.bonded_elem_name1[s1]} "
                    f"{self.bonded_species1[s1]} "
                    f"{self.bonded_elem_name2[s2]} "
                    f"{self.bonded_species2[s2]}"
                )

                # Determine if this bonded name is unique so far.
                if curr_bonded_name in unique_bonded_names:
                    continue
                num_unique_bonded_names += 1
                unique_bonded_names.append(curr_bonded_name)

                # The translation vector that will bring the two
                #   S atoms coincident with each other.  trans_vector1
                #   is consumed by the command-line builder further
                #   down (which reads ``trans_vector1[0..2]``), so it
                #   stays 0-indexed.
                trans_vector1 = [
                    self.s_atom_coords1[s1][ax + 1]
                    - self.s_atom_coords2[s2][ax + 1]
                    for ax in range(3)
                ]

                # The coordinates of the points that define a
                #   plane for the rotation are given by:
                #   - the alignment atom of molecule 1,
                #   - the S atom of molecule 1, and
                #   - the alignment atom of molecule 2 *after*
                #     it is translated.
                # There is a tricky part: the position of the
                #   second alignment atom is only going to be
                #   right *after* the previous translation. We
                #   don't have that value yet, so we compute the
                #   translation of that atom manually here.
                # Another tricky part: if the three atoms are
                #   colinear, we need to pick an arbitrary but
                #   non-colinear position.

                # Compute the final coordinate assuming non-
                #   colinear atoms.  Stored 1-indexed so it can
                #   later be fed to the 1-indexed cross-product
                #   helper without any reshaping.
                final_coord = [None] + [
                    trans_vector1[ax]
                    + self.trigger_coords2[s2][ax + 1]
                    for ax in range(3)
                ]

                # Check if the atoms are colinear by comparing
                #   normalised difference vectors.  Both diff vectors
                #   use the 1-indexed [None, dx, dy, dz] layout.
                diff_vector2 = [None] + [
                    self.s_atom_coords1[s1][axis]
                    - final_coord[axis]
                    for axis in range(1, 4)
                ]
                diff2_mag = math.sqrt(
                    diff_vector2[1] ** 2 +
                    diff_vector2[2] ** 2 +
                    diff_vector2[3] ** 2
                )
                diff_vector2 = [None] + [
                    diff_vector2[axis] / diff2_mag
                    for axis in range(1, 4)
                ]

                # Assume colinear. If any direction differs by
                #   more than epsilon, they are not colinear.
                not_colinear = False
                for axis in range(1, 4):
                    if (abs(abs(diff_vector1[axis])
                            - abs(diff_vector2[axis]))
                            > self.EPSILON):
                        not_colinear = True

                if not not_colinear:
                    # Make an arbitrary non-colinear coordinate.
                    # Create a vector not parallel to diff2 by
                    #   scaling components differently.
                    temp_vector = [
                        None,
                        diff_vector2[1] * 2.0,
                        diff_vector2[2] * 3.0,
                        diff_vector2[3] * 4.0,
                    ]
                    # Use a temporary SC for the cross product; it
                    # returns a 1-indexed vector directly.
                    sc_tmp = self.StructureControl()
                    final_coord = (
                        sc_tmp.normalized_cross_product(
                            diff_vector2, temp_vector
                        )
                    )

                # Build the plane coords for the -rotP option
                #   of mod_struct.py.
                plane_coords = [
                    self.s_atom_coords1[s1][1],
                    self.s_atom_coords1[s1][2],
                    self.s_atom_coords1[s1][3],
                    self.trigger_coords1[s1][1],
                    self.trigger_coords1[s1][2],
                    self.trigger_coords1[s1][3],
                    final_coord[1],
                    final_coord[2],
                    final_coord[3],
                ]

                # Compute the rotation angle.
                sc_tmp = self.StructureControl()
                rot_angle = sc_tmp.get_vector_angle(
                    diff_vector2, diff_vector1
                )
                rot_angle = (
                    (self.PI - rot_angle) * 180.0 / self.PI
                )

                # Origin of rotation (the S atom of molecule 1).
                origin = [
                    self.s_atom_coords1[s1][1],
                    self.s_atom_coords1[s1][2],
                    self.s_atom_coords1[s1][3],
                ]

                # Apply the translation, rotation, translation
                #   sequence using mod_struct.py (Python
                #   replacement for the Perl modStruct program).
                cmd = [
                    "mod_struct.py",
                    "-i",
                    os.path.join(
                        "centered", self.in_file_cent2
                    ),
                    "-o", os.path.join(
                        "centered", "temp.skl"
                    ),
                    "-trans",
                    str(trans_vector1[0]),
                    str(trans_vector1[1]),
                    str(trans_vector1[2]),
                    "-rotP",
                    str(plane_coords[0]),
                    str(plane_coords[1]),
                    str(plane_coords[2]),
                    str(plane_coords[3]),
                    str(plane_coords[4]),
                    str(plane_coords[5]),
                    str(plane_coords[6]),
                    str(plane_coords[7]),
                    str(plane_coords[8]),
                    "-angle", str(rot_angle),
                    "-orig",
                    str(origin[0]),
                    str(origin[1]),
                    str(origin[2]),
                    "-trans",
                    str(trans_vector2[0]),
                    str(trans_vector2[1]),
                    str(trans_vector2[2]),
                ]
                subprocess.run(cmd, check=True)

                # The just-created skeleton file can now be
                #   merged with the skeleton file of the first
                #   molecule.
                out_file = f"olcao_{s1}_{s2}.skl"

                # Open both files and merge them into one. As
                #   we merge them, we also create a catalog of
                #   the element-species assignments for each
                #   atom in each molecule.
                mol1_path = os.path.join(
                    "centered", self.in_file_cent1
                )
                mol2_path = os.path.join(
                    "centered", "temp.skl"
                )
                merge_path = os.path.join(
                    "merged", out_file
                )

                with (
                    open(mol1_path) as f1,
                    open(mol2_path) as f2,
                    open(merge_path, 'w') as fm,
                ):
                    self._merge_skl_files(
                        f1, f2, fm, s1, s2
                    )

        # Record the mapping of S atom numbers from the original
        #   molecules to the merged single molecule. For mol 1,
        #   the numbers are the same.
        for s1 in range(1, self.num_s_atoms1 + 1):
            while len(self.merged_s_atom_number1) <= s1:
                self.merged_s_atom_number1.append(None)
            self.merged_s_atom_number1[s1] = (
                self.s_atom_number1[s1]
            )

        # For mol 2, shift by the number of atoms in mol 1.
        for s2 in range(1, self.num_s_atoms2 + 1):
            while len(self.merged_s_atom_number2) <= s2:
                self.merged_s_atom_number2.append(None)
            self.merged_s_atom_number2[s2] = (
                self.s_atom_number2[s2] + self.num_atoms1
            )
            self.bonded_atom_number2[s2] += self.num_atoms1
            self.trigger_atom_number2[s2] += self.num_atoms1

    def _merge_skl_files(self, f1, f2, fm, s1, s2):
        """Merge two skeleton files into one merged file.

        This helper reads both molecule skeleton files
        line-by-line and writes a combined output with:
        - Merged title blocks
        - Unified cell parameters (from molecule 2)
        - Combined atom list
        - Standard footer (space 1_a, supercell 1 1 1, full)

        Parameters
        ----------
        f1, f2 : file objects
            Open handles for molecule 1 and 2 skeleton files.
        fm : file object
            Open handle for the merged output file.
        s1, s2 : int
            S-atom indices for the current pair.
        """
        # Put the titles for both molecules into the merged file.
        fm.write("title\n")
        fm.write("MOLECULE 1\n")
        next(f1)  # Read past the "title" keyword line.
        for line in f1:
            parts = line.strip().split()
            if parts and parts[0] == 'end':
                break
            fm.write(line)
        fm.write("MOLECULE 2\n")
        next(f2)  # Read past the "title" keyword line.
        for line in f2:
            parts = line.strip().split()
            if parts and parts[0] == 'end':
                break
            fm.write(line)
        fm.write("end\n")

        # Copy the lattice information. We arbitrarily use the
        #   cell parameters of the second molecule to define the
        #   cell of the merged molecules (both should be large
        #   enough).
        next(f1)  # Read past "cell" of molecule 1.
        next(f1)  # Read past the parameters.
        next(f2)  # Read past "cell" of molecule 2.
        cell_line = next(f2)  # Read the parameters of mol 2.
        fm.write("cell\n")
        fm.write(cell_line)
        next(f1)  # Read past "cartesian" in molecule 1.
        next(f2)  # Read past "cartesian" in molecule 2.
        total = self.num_atoms1 + self.num_atoms2
        fm.write(f"cartesian {total}\n")

        # Print the atoms from each molecule.
        for _ in range(self.num_atoms1):
            line = next(f1)
            fm.write(line)
        for _ in range(self.num_atoms2):
            line = next(f2)
            fm.write(line)
        fm.write("space 1_a\n")
        fm.write("supercell 1 1 1\n")
        fm.write("full\n")

    # ============================================================
    # catalog_elem_spec_mol_id
    # ============================================================
    def catalog_elem_spec_mol_id(self):
        """Catalog element names, species IDs, and molecule names
        for each atom in each merged molecule.

        Populate the mergedAtomSpeciesID array that lists the
        species number for each atom in the merged molecule,
        taking into account that the molecules might be different.
        Thus, atoms that originally have the same element + species
        names (e.g. both "c1") might now need to be specified as
        different types (e.g., one is now "c1" and the other is
        now "c2" if they came from different molecules but they
        might need to remain as both "c1" if the molecules are
        the same). Additionally, the type numbers for the atoms
        that participate in bonding between the molecules may need
        to change in the post-reaction template. That is also
        computed here.
        """
        num_unique_bonded_names = 0
        unique_bonded_names = [None]

        for s1 in range(1, self.num_s_atoms1 + 1):
            for s2 in range(1, self.num_s_atoms2 + 1):
                # Produce a name identifying the binding pair.
                curr_bonded_name = (
                    f"{self.bonded_elem_name1[s1]} "
                    f"{self.bonded_species1[s1]} "
                    f"{self.bonded_elem_name2[s2]} "
                    f"{self.bonded_species2[s2]}"
                )

                # Skip if this pair has been seen before.
                if curr_bonded_name in unique_bonded_names:
                    continue
                num_unique_bonded_names += 1
                unique_bonded_names.append(curr_bonded_name)

                # Read the merged skeleton file.
                merged_file = f"olcao_{s1}_{s2}.skl"
                sc = self.StructureControl()
                sc.read_input_file(
                    os.path.join("merged", merged_file),
                    use_file_species=True,
                )
                num_atoms = sc.num_atoms

                # Perform a bond analysis to get bonded atoms.
                subprocess.run(
                    [
                        "bond_analysis.py",
                        "-bl", "-bf", "1.1",
                        "-i",
                        os.path.join("merged", merged_file),
                    ],
                    check=True,
                )

                # Read the bondAnalysis results.
                sc.read_bond_analysis_bl("bondAnalysis.bl")

                # Assign the element name, species ID, and
                #   molecule name for each atom in the merged
                #   pre- and post-reaction templates.
                for atom in range(1, num_atoms + 1):
                    # Assign element name.
                    self.merged_atom_elem_name[
                        (s1, s2, atom)
                    ] = sc.atom_element_name[atom]

                    # Assign initial pre- and post-reaction
                    #   species numbers. (To be corrected below
                    #   according to the type_step.)
                    self.merged_atom_species_id[
                        (s1, s2, atom, 1)
                    ] = sc.atom_species_id[atom]
                    self.merged_atom_species_id[
                        (s1, s2, atom, 2)
                    ] = sc.atom_species_id[atom]

                    # Assign the molecule name.
                    if atom <= self.num_atoms1:
                        self.merged_atom_mol_name[
                            (s1, s2, atom)
                        ] = self.mol_name1
                    else:
                        self.merged_atom_mol_name[
                            (s1, s2, atom)
                        ] = self.mol_name2

                # If type_step is non-zero, correct the species
                #   ID numbers of the bonded atoms and the
                #   neighbours of the bonded atoms in the
                #   post-reaction merged molecule.
                if self.s.type_step > 0:
                    tb1 = self.bonded_atom_number1[s1]
                    tb2 = self.bonded_atom_number2[s2]

                    # Correct the bonded atoms first.
                    self.merged_atom_species_id[
                        (s1, s2, tb1, 2)
                    ] = (
                        sc.num_species[sc.atom_element_id[tb1]]
                        + self.s.type_step
                    )
                    self.merged_atom_species_id[
                        (s1, s2, tb2, 2)
                    ] = (
                        sc.num_species[sc.atom_element_id[tb2]]
                        + self.s.type_step
                    )

                    # Correct neighbours of the bonded atoms.
                    for bond in range(
                        1, sc.num_bonds[tb1] + 1
                    ):
                        partner = sc.bonded[tb1][bond]
                        self.merged_atom_species_id[
                            (s1, s2, partner, 2)
                        ] = (
                            sc.num_species[
                                sc.atom_element_id[tb1]
                            ]
                            + self.s.type_step + 1
                        )
                    for bond in range(
                        1, sc.num_bonds[tb2] + 1
                    ):
                        partner = sc.bonded[tb2][bond]
                        self.merged_atom_species_id[
                            (s1, s2, partner, 2)
                        ] = (
                            sc.num_species[
                                sc.atom_element_id[tb2]
                            ]
                            + self.s.type_step + 1
                        )

    # ============================================================
    # prune_merged_molecule
    # ============================================================
    def prune_merged_molecule(self):
        """Remove atoms that are far from the bond reaction point.

        This subroutine reads each merged molecule file in the
        "merged" directory and produces a new file for each in the
        "pruned" directory. The new file is a pruned version of
        the merged molecule. Pruning refers to the act of removing
        atoms that are a sufficient number of hops (defined by the
        chain length variable) away from the reacting S atoms.

        One of the tasks executed here includes reading a
        bondAnalysis.bl file and constructing a representation of
        the bonding data. In the future, this subroutine should be
        moved into the StructureControl module and modified to
        correctly deal with periodic boundary conditions. At
        present, this subroutine ignores PBC.
        """
        os.makedirs("pruned", exist_ok=True)

        num_atoms = self.num_atoms1 + self.num_atoms2

        num_unique_bonded_names = 0
        unique_bonded_names = [None]

        for s1 in range(1, self.num_s_atoms1 + 1):
            for s2 in range(1, self.num_s_atoms2 + 1):
                # Produce a name identifying the binding pair.
                curr_bonded_name = (
                    f"{self.bonded_elem_name1[s1]} "
                    f"{self.bonded_species1[s1]} "
                    f"{self.bonded_elem_name2[s2]} "
                    f"{self.bonded_species2[s2]}"
                )

                # Skip if this pair has been seen before.
                if curr_bonded_name in unique_bonded_names:
                    continue
                num_unique_bonded_names += 1
                unique_bonded_names.append(curr_bonded_name)

                # Compute the bond analysis for this merged mol.
                merged_mol_file = f"olcao_{s1}_{s2}.skl"
                subprocess.run(
                    [
                        "bond_analysis.py",
                        "-i",
                        os.path.join(
                            "merged", merged_mol_file
                        ),
                        "-bl", "-bf", "1.1",
                    ],
                    check=True,
                )

                # Read the bondAnalysis.bl file using SC.
                sc = self.StructureControl()
                sc.set_num_atoms(num_atoms)
                sc.read_bond_analysis_bl("bondAnalysis.bl")

                # Now we need to prune this molecule. Start by
                #   assuming that no atoms will be kept and then
                #   traverse from the current S atom.

                # Assume no atoms are retained except for the
                #   current S atoms. Initialize the edge count.
                self.num_edge_ids[(s1, s2)] = 0
                self.keep_atom = [0] * (num_atoms + 1)
                self.keep_atom[
                    self.merged_s_atom_number1[s1]
                ] = 1
                self.keep_atom[
                    self.merged_s_atom_number2[s2]
                ] = 1

                # Perform a BFS of the bonding network. Mark all
                #   atoms within the chain length as keepers.
                # We call it twice: once for each S atom with
                #   different chain lengths so that the two
                #   original molecules can be treated differently.

                # Initialize pruned atom count to 2 (the two S
                #   atoms). This is the number of atoms *left*
                #   after pruning.
                self.pruned_num_atoms[(s1, s2)] = 2

                # Probe the first molecule.
                self._mark_atoms_to_keep(
                    self.s.chain_len1,
                    self.merged_s_atom_number1[s1],
                    s1, s2,
                    sc.num_bonds, sc.bonded,
                )

                # Number of atoms in pruned first molecule.
                self.pruned_num_atoms1[(s1, s2)] = (
                    self.pruned_num_atoms[(s1, s2)] - 1
                )  # Subtract mol2 S atom.

                # Probe the second molecule.
                self._mark_atoms_to_keep(
                    self.s.chain_len2,
                    self.merged_s_atom_number2[s2],
                    s1, s2,
                    sc.num_bonds, sc.bonded,
                )

                # Number of atoms in pruned second molecule.
                self.pruned_num_atoms2[(s1, s2)] = (
                    self.pruned_num_atoms[(s1, s2)]
                    - self.pruned_num_atoms1[(s1, s2)]
                )

                # Create the pruned skeleton file by copying
                #   the merged file and omitting pruned atoms.
                merged_path = os.path.join(
                    "merged", merged_mol_file
                )
                pruned_path = os.path.join(
                    "pruned", merged_mol_file
                )

                with (
                    open(merged_path) as fin,
                    open(pruned_path, 'w') as fout,
                ):
                    self._write_pruned_skl(
                        fin, fout, s1, s2, num_atoms
                    )

    def _mark_atoms_to_keep(
        self, max_chain_len, current_atom,
        s1, s2, num_bonds, bonded,
    ):
        """BFS traversal of the bonding network to identify atoms
        that should not be pruned away.

        This is a breadth-first search starting from the given
        S atom (trigger atom). Atoms within max_chain_len hops
        are marked as keepers. Atoms at the exact chain-length
        boundary are also marked as edge atoms.

        Parameters
        ----------
        max_chain_len : int
            Maximum number of hops to include.
        current_atom : int
            Starting atom (an S atom).
        s1, s2 : int
            S-atom pair indices.
        num_bonds : list
            Bond counts per atom (from bond analysis).
        bonded : list
            Bonding partners per atom.
        """
        # Only the trigger atoms will have a non-zero keep
        #   value when this subroutine is called. (They have a
        #   value of 1.) The neighbours bonded to the trigger
        #   atoms with keep==0 need to be added to the queue
        #   along with a record of their chain level (2).
        queue = deque()
        chain_level_queue = deque()
        current_chain_level = 2

        for bond in range(1, num_bonds[current_atom] + 1):
            partner = bonded[current_atom][bond]
            if self.keep_atom[partner] == 0:
                chain_level_queue.append(
                    current_chain_level
                )
                queue.append(partner)

        # While the queue is not empty:
        #   Shift the current atom and chain level from the
        #   queue. If keep==0 and within chain length, mark
        #   atom to keep. Otherwise, skip to the next item.
        #
        #   If within the chain length, add each bonded atom
        #   (that isn't already a keeper and isn't already in
        #   the queue) to the queue.
        #   Otherwise, mark this atom as an edge and record it.
        while queue:
            current_atom = queue.popleft()
            current_chain_level = chain_level_queue.popleft()

            if (self.keep_atom[current_atom] == 0
                    and current_chain_level <= max_chain_len):
                # Mark this atom as a keeper and increment the
                #   count of retained atoms.
                self.keep_atom[current_atom] = 1
                self.pruned_num_atoms[(s1, s2)] += 1
            else:
                continue

            if current_chain_level < max_chain_len:
                for bond in range(
                    1, num_bonds[current_atom] + 1
                ):
                    partner = bonded[current_atom][bond]
                    # Check if already a keeper or in the queue.
                    if (self.keep_atom[partner] == 0
                            and partner not in queue):
                        chain_level_queue.append(
                            current_chain_level + 1
                        )
                        queue.append(partner)
            else:
                # Identify the edge atom.
                self.num_edge_ids[(s1, s2)] += 1
                eid = self.num_edge_ids[(s1, s2)]
                self.edge_id[(eid, s1, s2)] = current_atom

    def _write_pruned_skl(
        self, fin, fout, s1, s2, num_atoms
    ):
        """Write a pruned skeleton file from a merged file.

        Copies all content from the merged skeleton file to the
        pruned file, but replaces the cartesian line with the
        pruned atom count and omits atoms that were pruned.
        Along the way, maps specialty index numbers from the
        merged to the pruned atom numbering.

        Parameters
        ----------
        fin, fout : file objects
            Input (merged) and output (pruned) file handles.
        s1, s2 : int
            S-atom pair indices.
        num_atoms : int
            Total atoms in the merged file.
        """
        for line in fin:
            parts = line.strip().split()
            if parts and parts[0].startswith('cartesian'):
                # Print new cartesian line with pruned count.
                pruned_count = self.pruned_num_atoms[(s1, s2)]
                fout.write(f"cartesian {pruned_count}\n")

                # Go through all atoms, omitting pruned ones.
                pruned_atom_count = 0
                for atom in range(1, num_atoms + 1):
                    atom_line = next(fin)
                    if self.keep_atom[atom] != 1:
                        continue

                    fout.write(atom_line)
                    pruned_atom_count += 1

                    # Map specialty index numbers.
                    if (atom
                            == self.merged_s_atom_number1[s1]):
                        self.pruned_s_atom_number1[
                            (s1, s2)
                        ] = pruned_atom_count
                    elif (atom
                          == self.merged_s_atom_number2[s2]):
                        self.pruned_s_atom_number2[
                            (s1, s2)
                        ] = pruned_atom_count
                    if (atom
                            == self.trigger_atom_number1[s1]):
                        self.pruned_trigger_atom_number1[
                            (s1, s2)
                        ] = pruned_atom_count
                    elif (atom
                          == self.trigger_atom_number2[s2]):
                        self.pruned_trigger_atom_number2[
                            (s1, s2)
                        ] = pruned_atom_count
                    if (atom
                            == self.bonded_atom_number1[s1]):
                        self.pruned_bonded_atom_number1[
                            (s1, s2)
                        ] = pruned_atom_count
                    elif (atom
                          == self.bonded_atom_number2[s2]):
                        self.pruned_bonded_atom_number2[
                            (s1, s2)
                        ] = pruned_atom_count

                    # Map edge IDs.
                    for eid in range(
                        1,
                        self.num_edge_ids[(s1, s2)] + 1,
                    ):
                        if (atom
                                == self.edge_id[
                                    (eid, s1, s2)
                                ]):
                            self.pruned_edge_id[
                                (eid, s1, s2)
                            ] = pruned_atom_count
                            break

                    # Map delete-tree IDs.
                    ndel = self.merged_num_del_atoms[(s1, s2)]
                    for did in range(1, ndel + 1):
                        dt = self.merged_del_tree[(s1, s2)]
                        if atom == dt[did]:
                            self.pruned_del_tree[
                                (s1, s2, did)
                            ] = pruned_atom_count

                    # Map element names and species IDs.
                    self.pruned_atom_elem_name[
                        (s1, s2, pruned_atom_count)
                    ] = self.merged_atom_elem_name[
                        (s1, s2, atom)
                    ]
                    self.pruned_atom_species_id[
                        (s1, s2, pruned_atom_count, 1)
                    ] = self.merged_atom_species_id[
                        (s1, s2, atom, 1)
                    ]
                    self.pruned_atom_species_id[
                        (s1, s2, pruned_atom_count, 2)
                    ] = self.merged_atom_species_id[
                        (s1, s2, atom, 2)
                    ]
                    self.pruned_atom_mol_name[
                        (s1, s2, pruned_atom_count)
                    ] = self.merged_atom_mol_name[
                        (s1, s2, atom)
                    ]
            else:
                fout.write(line)

    # ============================================================
    # make_reaction_templates
    # ============================================================
    def make_reaction_templates(self, rxn_phase):
        """Build LAMMPS bond/react templates from the pruned
        olcao.skl files.

        The reaction requires a pre- and a post-reaction template
        as well as a map file that relates the two.

        Pre-reaction template (rxn_phase == 1):
          Contains all atoms in the pruned file and all bonds
          among atoms of the originally distinct molecules, but
          no bonds connecting the two molecules. All elements
          (masses) and types are assigned per the pruned skeleton
          file. For the reaction to occur, an S atom from one
          molecule must be near an S atom from the other.

        Post-reaction template (rxn_phase == 2):
          Contains the same atoms but all bonds and angles that
          used to include either S atom are removed. One new bond
          between the now-bare atoms that were bonded to the S
          atoms is added. The types of the atoms and bonds may
          change to re-define the topology.

        Template map file (written with the pre-reaction):
          All atoms in the pre-reaction template are the same as
          in the post. Equivalences are 1-to-1 (1 1, 2 2, etc.).
          Bonding IDs are the atoms that bind. Delete IDs are the
          S atoms. Edge IDs are from the pruning step.

        Parameters
        ----------
        rxn_phase : int
            1 for pre-reaction, 2 for post-reaction.
        """
        # Create a directory to store the templates.
        os.makedirs("rxnTemplates", exist_ok=True)

        num_unique_bonded_names = 0
        unique_bonded_names = [None]

        for s1 in range(1, self.num_s_atoms1 + 1):
            for s2 in range(1, self.num_s_atoms2 + 1):
                # Produce a name identifying the binding pair.
                curr_bonded_name = (
                    f"{self.bonded_elem_name1[s1]} "
                    f"{self.bonded_species1[s1]} "
                    f"{self.bonded_elem_name2[s2]} "
                    f"{self.bonded_species2[s2]}"
                )

                # Skip if this pair has been seen before.
                if curr_bonded_name in unique_bonded_names:
                    continue
                num_unique_bonded_names += 1
                unique_bonded_names.append(curr_bonded_name)

                # Read the pruned skeleton file.
                pruned_file = f"olcao_{s1}_{s2}.skl"
                sc = self.StructureControl()
                sc.read_input_file(
                    os.path.join("pruned", pruned_file),
                    use_file_species=True,
                )
                num_atoms = sc.num_atoms

                # ---- Bond length analysis ----
                subprocess.run(
                    [
                        "bond_analysis.py",
                        "-i",
                        os.path.join(
                            "pruned", pruned_file
                        ),
                        "-bl", "-bf", "1.1",
                    ],
                    check=True,
                )

                # Read bond data manually (same as Perl).
                (
                    num_bonds, bonded, bond_tag_id,
                    unique_bond_tags, num_unique_bond_tags,
                ) = self._read_bond_data(
                    "bondAnalysis.bl", sc, s1, s2
                )

                # ---- Bond angle analysis ----
                subprocess.run(
                    [
                        "bond_analysis.py",
                        "-i",
                        os.path.join(
                            "pruned", pruned_file
                        ),
                        "-ba", "-bf", "1.1",
                    ],
                    check=True,
                )

                (
                    num_bond_angles, angle_bonded,
                    angle_tag_id, unique_angle_tags,
                    num_unique_angle_tags, num_angles_total,
                    bond_angles_ext,
                ) = self._read_angle_data(
                    "bondAnalysis.ba", sc, s1, s2
                )

                # ---- Compute ordered species IDs ----
                (
                    ordered_species_id,
                    ordered_species_names,
                ) = self._compute_ordered_species(
                    sc, num_atoms, s1, s2
                )

                # ---- Construct angle list for this phase ----
                phase_angle_list = self._build_phase_angles(
                    rxn_phase, num_atoms,
                    num_bond_angles, angle_bonded,
                    s1, s2,
                )

                # ---- Construct bond list for this phase ----
                (
                    phase_bond_list, phase_num_bonds,
                    num_bonds, bonded,
                    bond_tag_id,
                    unique_bond_tags,
                    num_unique_bond_tags,
                ) = self._build_phase_bonds(
                    rxn_phase, num_atoms,
                    num_bonds, bonded,
                    bond_tag_id,
                    unique_bond_tags,
                    num_unique_bond_tags,
                    s1, s2, sc,
                )

                # ---- Define bonding and edge IDs ----
                bonding_id1 = 0
                bonding_id2 = 0
                if rxn_phase == 1:
                    bonding_id1 = (
                        self.pruned_trigger_atom_number1[
                            (s1, s2)
                        ]
                    )
                    bonding_id2 = (
                        self.pruned_trigger_atom_number2[
                            (s1, s2)
                        ]
                    )

                # Construct the edge IDs for this phase.
                phase_num_edge_ids = 0
                phase_pruned_edge_id = [None]
                if rxn_phase == 1:
                    for eid in range(
                        1,
                        self.num_edge_ids[(s1, s2)] + 1,
                    ):
                        pe = self.pruned_edge_id[
                            (eid, s1, s2)
                        ]
                        ps1 = self.pruned_s_atom_number1[
                            (s1, s2)
                        ]
                        ps2 = self.pruned_s_atom_number2[
                            (s1, s2)
                        ]
                        if pe != ps1 and pe != ps2:
                            phase_num_edge_ids += 1
                            phase_pruned_edge_id.append(pe)

                # ---- Define file names ----
                pb1 = self.pruned_bonded_atom_number1[
                    (s1, s2)
                ]
                pb2 = self.pruned_bonded_atom_number2[
                    (s1, s2)
                ]
                rxn_binding_pair = (
                    f"{sc.atom_element_name[pb1]}-"
                    f"{sc.atom_species_id[pb1]}_"
                    f"{sc.atom_element_name[pb2]}-"
                    f"{sc.atom_species_id[pb2]}"
                )

                if rxn_phase % 2 != 0:
                    pre_or_post = "preRxn"
                else:
                    pre_or_post = "postRxn"

                map_file = f"{rxn_binding_pair}.map"
                template = (
                    f"{pre_or_post}.{rxn_binding_pair}.data"
                )

                # ---- Write the template file ----
                self._write_template(
                    template, num_atoms, rxn_phase,
                    phase_bond_list, len(phase_bond_list) - 1,
                    phase_angle_list, len(phase_angle_list) - 1,
                    sc, s1, s2,
                    ordered_species_id,
                    num_bonds, bonded,
                    bond_tag_id, unique_bond_tags,
                    num_bond_angles, angle_bonded,
                    angle_tag_id, unique_angle_tags,
                    bond_angles_ext,
                )

                # ---- Write the map file (pre-reaction only)
                if rxn_phase % 2 != 0:
                    self._write_map_file(
                        map_file, num_atoms,
                        phase_num_edge_ids,
                        phase_pruned_edge_id,
                        bonding_id1, bonding_id2,
                        s1, s2,
                    )

    def _read_bond_data(self, filename, sc, s1, s2):
        """Read bond data from a bondAnalysis.bl file and build
        local bonding arrays with cross-molecule bond rejection.

        Returns
        -------
        tuple
            (num_bonds, bonded, bond_tag_id,
             unique_bond_tags, num_unique_bond_tags)
        """
        num_unique_bond_tags = 0
        unique_bond_tags = [None]
        num_atoms = sc.num_atoms
        num_bonds = [0] * (num_atoms + 1)
        bonded = [[None] for _ in range(num_atoms + 1)]
        bond_tag_id = [[None] for _ in range(num_atoms + 1)]

        with open(filename) as f:
            while True:
                line = f.readline()
                if not line:
                    break
                values = line.split()
                if not values:
                    continue

                parts = values[0].split('_')
                curr_atom = int(parts[1])
                atom_tag = (
                    f"{sc.atom_element_name[curr_atom]} "
                    f"{sc.atom_species_id[curr_atom]}"
                )
                pn1 = self.pruned_num_atoms1[(s1, s2)]
                if curr_atom <= pn1:
                    atom_tag += f" {self.mol_name1}"
                else:
                    atom_tag += f" {self.mol_name2}"

                n_bonds = int(values[2])
                num_bonds[curr_atom] = n_bonds
                num_bond_rows = math.ceil(n_bonds / 4.0)
                bond_index = 0

                for _ in range(num_bond_rows):
                    row_vals = f.readline().split()
                    n_in_row = len(row_vals) // 2
                    for b in range(n_in_row):
                        bp = row_vals[2 * b].split('_')
                        partner = int(bp[1])
                        bond_index += 1

                        # Temporarily add to check bonding.
                        while (len(bonded[curr_atom])
                               <= bond_index):
                            bonded[curr_atom].append(None)
                        bonded[curr_atom][bond_index] = partner

                        # Reject cross-molecule bonds.
                        if not self._check_bonding(
                            curr_atom, bond_index,
                            s1, s2, bonded,
                        ):
                            bond_index -= 1
                            num_bonds[curr_atom] -= 1
                            bonded[curr_atom].pop()
                            continue

                        # Build the bond tag.
                        bonded_tag = (
                            f"{sc.atom_element_name[partner]}"
                            f" {sc.atom_species_id[partner]}"
                        )
                        if partner <= pn1:
                            bonded_tag += (
                                f" {self.mol_name1}"
                            )
                        else:
                            bonded_tag += (
                                f" {self.mol_name2}"
                            )

                        # Make the tag ordered.
                        if atom_tag > bonded_tag:
                            bond_tag = (
                                f"{bonded_tag} {atom_tag}"
                            )
                        else:
                            bond_tag = (
                                f"{atom_tag} {bonded_tag}"
                            )

                        found = 0
                        for t in range(
                            1, num_unique_bond_tags + 1
                        ):
                            if (bond_tag
                                    == unique_bond_tags[t]):
                                found = t
                                break
                        if found == 0:
                            num_unique_bond_tags += 1
                            unique_bond_tags.append(bond_tag)
                            found = num_unique_bond_tags

                        while (len(bond_tag_id[curr_atom])
                               <= bond_index):
                            bond_tag_id[curr_atom].append(None)
                        bond_tag_id[curr_atom][
                            bond_index
                        ] = found

        return (
            num_bonds, bonded, bond_tag_id,
            unique_bond_tags, num_unique_bond_tags,
        )

    def _read_angle_data(self, filename, sc, s1, s2):
        """Read bond angle data from a bondAnalysis.ba file.

        Returns
        -------
        tuple
            (num_bond_angles, angle_bonded, angle_tag_id,
             unique_angle_tags, num_unique_angle_tags,
             num_angles_total, bond_angles_ext)
        """
        num_atoms = sc.num_atoms
        num_bond_angles = [0] * (num_atoms + 1)
        angle_bonded = [
            [None] for _ in range(num_atoms + 1)
        ]
        angle_tag_id = [
            [None] for _ in range(num_atoms + 1)
        ]
        bond_angles_ext = [
            [None] for _ in range(num_atoms + 1)
        ]

        num_unique_angle_tags = 0
        unique_angle_tags = [None]
        num_angles_total = 0

        pn1 = self.pruned_num_atoms1[(s1, s2)]

        with open(filename) as f:
            while True:
                line = f.readline()
                if not line:
                    break
                values = line.split()
                if not values:
                    continue

                curr_atom = int(values[0])
                n_angles = int(values[4])
                num_bond_angles[curr_atom] = n_angles
                num_angles_total += n_angles

                if curr_atom <= pn1:
                    curr_mol = self.mol_name1
                else:
                    curr_mol = self.mol_name2

                for angle in range(1, n_angles + 1):
                    vals = f.readline().split()

                    # Round the bond angle.
                    raw_angle = float(vals[-1])
                    decimal = raw_angle - int(raw_angle)
                    if decimal > 0.75:
                        curr_bond_angle = int(raw_angle) + 1.0
                    elif decimal > 0.25:
                        curr_bond_angle = int(raw_angle) + 0.5
                    else:
                        curr_bond_angle = float(int(raw_angle))

                    while (len(bond_angles_ext[curr_atom])
                           <= angle):
                        bond_angles_ext[curr_atom].append(None)
                    bond_angles_ext[curr_atom][angle] = (
                        curr_bond_angle
                    )

                    # Create a tag for this angle.
                    # vals contains:
                    #   0=tag1, 1=vertex_tag, 2=tag2,
                    #   3=atom1, 4=vertex, 5=atom2, last=angle
                    a1 = int(vals[3])
                    v = int(vals[4])
                    a2 = int(vals[5])

                    while (len(angle_bonded[curr_atom])
                           <= angle):
                        angle_bonded[curr_atom].append(
                            [None, None, None]
                        )

                    if vals[0] <= vals[2]:
                        angle_tag = (
                            f"{sc.atom_element_name[a1]} "
                            f"{sc.atom_species_id[a1]} "
                            f"{curr_mol} "
                            f"{sc.atom_element_name[v]} "
                            f"{sc.atom_species_id[v]} "
                            f"{curr_mol} "
                            f"{sc.atom_element_name[a2]} "
                            f"{sc.atom_species_id[a2]} "
                            f"{curr_mol}"
                        )
                        angle_bonded[curr_atom][angle] = (
                            [None, a1, a2]
                        )
                    else:
                        angle_tag = (
                            f"{sc.atom_element_name[a2]} "
                            f"{sc.atom_species_id[a2]} "
                            f"{curr_mol} "
                            f"{sc.atom_element_name[v]} "
                            f"{sc.atom_species_id[v]} "
                            f"{curr_mol} "
                            f"{sc.atom_element_name[a1]} "
                            f"{sc.atom_species_id[a1]} "
                            f"{curr_mol}"
                        )
                        angle_bonded[curr_atom][angle] = (
                            [None, a2, a1]
                        )

                    # Match to the Hooke angle database.
                    z1 = self.ed.get_element_z(
                        angle_tag.split()[0]
                    )
                    z2 = self.ed.get_element_z(
                        angle_tag.split()[3]
                    )
                    z3 = self.ed.get_element_z(
                        angle_tag.split()[6]
                    )
                    if z1 > z3:
                        z1, z3 = z3, z1

                    found_hooke = 0
                    ad = self.angle_data
                    for h in range(
                        1, ad.num_hooke_angles + 1
                    ):
                        hc = ad.hooke_angle_coeffs[h]
                        if (z1 == hc[1] and z2 == hc[2]
                                and z3 == hc[3]):
                            if (abs(curr_bond_angle - hc[5])
                                    <= hc[6]):
                                found_hooke = h
                                break

                    if found_hooke == 0:
                        print(
                            f"Angle number = {angle} "
                            f"Angle = {curr_bond_angle}"
                        )
                        print(f"atomZSet = {z1} {z2} {z3}")
                        sys.exit(
                            "Cannot find angle in the "
                            "database"
                        )

                    angle_tag += (
                        f" {ad.hooke_angle_coeffs[found_hooke][5]}"
                        f" {found_hooke}"
                    )

                    # Check if this angle tag is unique.
                    found = 0
                    for t in range(
                        1, num_unique_angle_tags + 1
                    ):
                        if angle_tag == unique_angle_tags[t]:
                            found = t
                            break
                    if found == 0:
                        num_unique_angle_tags += 1
                        found = num_unique_angle_tags
                        unique_angle_tags.append(angle_tag)

                    while (len(angle_tag_id[curr_atom])
                           <= angle):
                        angle_tag_id[curr_atom].append(None)
                    angle_tag_id[curr_atom][angle] = found

        return (
            num_bond_angles, angle_bonded,
            angle_tag_id, unique_angle_tags,
            num_unique_angle_tags, num_angles_total,
            bond_angles_ext,
        )

    def _compute_ordered_species(
        self, sc, num_atoms, s1, s2,
    ):
        """Compute an ordered list of species ID numbers for
        atoms in the current reaction template.

        It is crucially important to realize that the ordered
        species ID numbers here are really only valid for the
        current reaction template and not the entire simulation.
        To make things better, it will be necessary to merge all
        types from all templates and from the LAMMPS data file
        into one unified listing. That task is reserved for the
        "condense" script.

        Returns
        -------
        tuple
            (ordered_species_id, ordered_species_names) lists.
        """
        num_unique = 0
        unique_names = [None]
        ordered_species_id = [None] * (num_atoms + 1)
        ordered_species_names = [None] * (num_atoms + 1)

        for atom in range(1, num_atoms + 1):
            # Create an identifier for the current atom.
            if not self.s.const_types:
                curr_name = (
                    f"{sc.atom_element_name[atom]}"
                    f"{sc.atom_species_id[atom]} "
                    f"{self.pruned_atom_mol_name[(s1, s2, atom)]}"
                )
            else:
                curr_name = (
                    f"{sc.atom_element_name[atom]}"
                    f"{sc.atom_species_id[atom]}"
                )

            # Determine if this identifier has been seen.
            found = 0
            for i in range(1, num_unique + 1):
                if curr_name == unique_names[i]:
                    found = i
                    break

            if found == 0:
                num_unique += 1
                found = num_unique
                unique_names.append(curr_name)

            ordered_species_id[atom] = found
            ordered_species_names[atom] = curr_name

        return ordered_species_id, ordered_species_names

    def _build_phase_angles(
        self, rxn_phase, num_atoms,
        num_bond_angles, angle_bonded,
        s1, s2,
    ):
        """Build the list of bond angles for the current reaction
        phase.

        Include all bond angles as defined for the molecules
        before bonding, and include the same set post-bonding
        minus any bond angles that touched the S atoms in some
        way.

        In the future it might be necessary to *add* bond angles
        through the bonding atoms (after the S are removed), but
        presently we do not do that. (Hence the bonded molecules
        may be too floppy.)

        Returns
        -------
        list
            1-indexed list of [None, atom, angle] entries.
        """
        phase_angle_list = [None]

        for atom in range(1, num_atoms + 1):
            for angle in range(
                1, num_bond_angles[atom] + 1
            ):
                angle_set = [
                    angle_bonded[atom][angle][1],
                    atom,
                    angle_bonded[atom][angle][2],
                ]

                # Exclude angles that cross between molecules.
                if not self._check_angle_bonding(
                    atom, angle, s1, s2, angle_bonded
                ):
                    continue

                if rxn_phase == 1:
                    # Keep all bond angles.
                    phase_angle_list.append(
                        [None, atom, angle]
                    )
                elif rxn_phase == 2:
                    # Exclude angles including delete-tree atoms.
                    ndel = self.merged_num_del_atoms[(s1, s2)]
                    found = False
                    for did in range(1, ndel + 1):
                        dt = self.pruned_del_tree.get(
                            (s1, s2, did), None
                        )
                        if dt is None:
                            continue
                        if (atom == dt
                                or angle_set[0] == dt
                                or angle_set[2] == dt):
                            found = True
                            break
                    if found:
                        continue
                    phase_angle_list.append(
                        [None, atom, angle]
                    )

        return phase_angle_list

    def _build_phase_bonds(
        self, rxn_phase, num_atoms,
        num_bonds, bonded,
        bond_tag_id,
        unique_bond_tags, num_unique_bond_tags,
        s1, s2, sc,
    ):
        """Build the list of bonds for the current reaction phase.

        This is done in two parts. In the first part we go through
        all known bonds and include only those that are actually
        needed. The second part explicitly adds bonds that would
        not otherwise be present (the new bond between the
        previously S-bonded atoms in the post-reaction case).

        Returns
        -------
        tuple
            Updated (phase_bond_list, phase_num_bonds,
                     num_bonds, bonded,
                     bond_tag_id,
                     unique_bond_tags, num_unique_bond_tags)
        """
        phase_bond_list = [None]
        phase_num_bonds = 0

        for atom in range(1, num_atoms + 1):
            for bond in range(1, num_bonds[atom] + 1):
                # Avoid double-counting: only include when the
                #   bonded partner has a higher index.
                if bonded[atom][bond] < atom:
                    continue

                # Only permit bonds within the same molecule.
                if not self._check_bonding(
                    atom, bond, s1, s2, bonded
                ):
                    continue

                if rxn_phase == 1:
                    # Include all bonds.
                    phase_num_bonds += 1
                    phase_bond_list.append(
                        [None, atom, bond]
                    )
                elif rxn_phase == 2:
                    # Exclude bonds involving delete-tree atoms.
                    ndel = self.merged_num_del_atoms[(s1, s2)]
                    found = False
                    for did in range(1, ndel + 1):
                        dt = self.pruned_del_tree.get(
                            (s1, s2, did), None
                        )
                        if dt is None:
                            continue
                        if (atom == dt
                                or bonded[atom][bond] == dt):
                            found = True
                            break
                    if found:
                        continue
                    phase_num_bonds += 1
                    phase_bond_list.append(
                        [None, atom, bond]
                    )

        # Add special bonds not in bondAnalysis results.
        new_bond_atom1 = 0
        new_bond_atom2 = 0

        # The bond between atoms that used to be bonded to an
        #   S atom is present in the post-reaction template.
        if rxn_phase == 2:
            new_bond_atom1 = self.pruned_bonded_atom_number1[
                (s1, s2)
            ]
            new_bond_atom2 = self.pruned_bonded_atom_number2[
                (s1, s2)
            ]

        if new_bond_atom1 > 0:
            # Increment bonds for the first atom.
            num_bonds[new_bond_atom1] += 1

            # Record the second atom as the bonded partner.
            while (len(bonded[new_bond_atom1])
                   <= num_bonds[new_bond_atom1]):
                bonded[new_bond_atom1].append(None)
            bonded[new_bond_atom1][
                num_bonds[new_bond_atom1]
            ] = new_bond_atom2

            # Record into the phase bond list.
            phase_num_bonds += 1
            phase_bond_list.append(
                [None, new_bond_atom1,
                 num_bonds[new_bond_atom1]]
            )

            # Make a tag for this bond pair.
            sp1 = self.pruned_atom_species_id[
                (s1, s2, new_bond_atom1, rxn_phase)
            ]
            tag1 = (
                f"{self.pruned_atom_elem_name[(s1, s2, new_bond_atom1)]}"
                f" {sp1}"
                f" {self.pruned_atom_mol_name[(s1, s2, new_bond_atom1)]}"
            )
            sp2 = self.pruned_atom_species_id[
                (s1, s2, new_bond_atom2, rxn_phase)
            ]
            tag2 = (
                f"{self.pruned_atom_elem_name[(s1, s2, new_bond_atom2)]}"
                f" {sp2}"
                f" {self.pruned_atom_mol_name[(s1, s2, new_bond_atom2)]}"
            )
            if tag1 > tag2:
                bond_tag = f"{tag2} {tag1}"
            else:
                bond_tag = f"{tag1} {tag2}"

            found = 0
            for t in range(1, num_unique_bond_tags + 1):
                if bond_tag == unique_bond_tags[t]:
                    found = t
                    break
            if found == 0:
                num_unique_bond_tags += 1
                unique_bond_tags.append(bond_tag)
                found = num_unique_bond_tags

            while (len(bond_tag_id[new_bond_atom1])
                   <= num_bonds[new_bond_atom1]):
                bond_tag_id[new_bond_atom1].append(None)
            bond_tag_id[new_bond_atom1][
                num_bonds[new_bond_atom1]
            ] = found

        return (
            phase_bond_list, phase_num_bonds,
            num_bonds, bonded,
            bond_tag_id,
            unique_bond_tags, num_unique_bond_tags,
        )

    def _check_bonding(self, atom, bond, s1, s2, bonded):
        """Determine whether a bond is valid (not crossing
        between the two molecules).

        Make sure that we never include any bond between
        molecule 1 and molecule 2. The two molecules should be
        placed "far enough" apart that they never have bonds
        between them. However, just to be certain, we do this
        check.

        Parameters
        ----------
        atom : int
            Atom index.
        bond : int
            Bond index for this atom.
        s1, s2 : int
            S-atom pair indices.
        bonded : list
            Bonding list.

        Returns
        -------
        bool
            True if the bond is within the same molecule.
        """
        pn1 = self.pruned_num_atoms1[(s1, s2)]
        partner = bonded[atom][bond]
        if ((atom <= pn1 and partner > pn1)
                or (atom > pn1 and partner <= pn1)):
            return False
        return True

    def _check_angle_bonding(
        self, atom, angle, s1, s2, angle_bonded
    ):
        """Determine whether an angle is valid (not crossing
        between the two molecules).

        Make sure that we never include any angle between
        molecule 1 and molecule 2.

        Parameters
        ----------
        atom : int
            Vertex atom index.
        angle : int
            Angle index for this vertex atom.
        s1, s2 : int
            S-atom pair indices.
        angle_bonded : list
            Angle bonding list.

        Returns
        -------
        bool
            True if the angle is within the same molecule.
        """
        pn1 = self.pruned_num_atoms1[(s1, s2)]
        a1 = angle_bonded[atom][angle][1]
        a2 = angle_bonded[atom][angle][2]

        if atom <= pn1:
            if a1 > pn1 or a2 > pn1:
                return False
        else:
            if a1 <= pn1 or a2 <= pn1:
                return False
        return True

    def _write_template(
        self, template_name, num_atoms, rxn_phase,
        phase_bond_list, phase_num_bonds,
        phase_angle_list, phase_num_angles,
        sc, s1, s2,
        ordered_species_id,
        num_bonds, bonded,
        bond_tag_id, unique_bond_tags,
        num_bond_angles, angle_bonded,
        angle_tag_id, unique_angle_tags,
        bond_angles_ext,
    ):
        """Write a pre- or post-reaction LAMMPS template file.

        Parameters
        ----------
        template_name : str
            Output file name (e.g. "preRxn.c-1_b-1.data").
        num_atoms : int
            Number of atoms in the pruned template.
        rxn_phase : int
            1 for pre-reaction, 2 for post-reaction.
        phase_bond_list : list
            1-indexed list of [None, atom, bond] entries.
        phase_num_bonds : int
            Number of bonds in this phase.
        phase_angle_list : list
            1-indexed list of [None, atom, angle] entries.
        phase_num_angles : int
            Number of angles in this phase.
        sc : StructureControl
            Structure control instance with pruned data loaded.
        s1, s2 : int
            S-atom pair indices.
        ordered_species_id : list
            Species ID for each atom in the template.
        num_bonds, bonded : lists
            Bond data arrays.
        bond_tag_id, unique_bond_tags : lists
            Bond type data.
        num_bond_angles, angle_bonded : lists
            Angle data arrays.
        angle_tag_id, unique_angle_tags : lists
            Angle type data.
        bond_angles_ext : list
            Computed bond angles.
        """
        path = os.path.join("rxnTemplates", template_name)
        with open(path, 'w') as out:
            out.write(f"{template_name}\n\n")

            out.write(f"{num_atoms} atoms\n")
            out.write(f"{phase_num_bonds} bonds\n")
            out.write(f"{phase_num_angles} angles\n")
            out.write("0 dihedrals\n")
            out.write("0 impropers\n")

            # Masses section.
            out.write("\nMasses\n\n")
            for atom in range(1, num_atoms + 1):
                z = sc.atom_element_id[atom]
                mass = self.ed.atomic_masses[z]
                en = self.pruned_atom_elem_name[
                    (s1, s2, atom)
                ]
                sp = self.pruned_atom_species_id[
                    (s1, s2, atom, rxn_phase)
                ]
                mn = self.pruned_atom_mol_name[
                    (s1, s2, atom)
                ]
                out.write(
                    f"{atom} {mass} "
                    f"# {en} {sp} {mn}\n"
                )

            # Types section.
            out.write("\nTypes\n\n")
            for atom in range(1, num_atoms + 1):
                en = self.pruned_atom_elem_name[
                    (s1, s2, atom)
                ]
                sp = self.pruned_atom_species_id[
                    (s1, s2, atom, rxn_phase)
                ]
                mn = self.pruned_atom_mol_name[
                    (s1, s2, atom)
                ]
                out.write(
                    f"{atom} {ordered_species_id[atom]} "
                    f"# {en} {sp} {mn}\n"
                )

            # Bonds section.
            out.write("\nBonds\n\n")
            for b in range(1, phase_num_bonds + 1):
                ca = phase_bond_list[b][1]
                cb = phase_bond_list[b][2]
                bt = bond_tag_id[ca][cb]
                out.write(
                    f"{b} {bt} {ca} "
                    f"{bonded[ca][cb]} "
                    f"# {unique_bond_tags[bt]}\n"
                )

            # Angles section.
            out.write("\nAngles\n\n")
            for a in range(1, phase_num_angles + 1):
                ca = phase_angle_list[a][1]
                cang = phase_angle_list[a][2]
                at = angle_tag_id[ca][cang]
                a1 = angle_bonded[ca][cang][1]
                a2 = angle_bonded[ca][cang][2]
                out.write(
                    f"{a} {at} {a1} {ca} {a2} "
                    f" # {unique_angle_tags[at]}\n"
                )

    def _write_map_file(
        self, map_file_name, num_atoms,
        phase_num_edge_ids, phase_pruned_edge_id,
        bonding_id1, bonding_id2,
        s1, s2,
    ):
        """Write the reaction map file.

        All atoms in the pre-reaction template are the same as
        those in the post reaction template. Therefore, there
        will be one equivalency for each atom (1 1, 2 2, 3 3,
        etc.). The bonding IDs will be the atom numbers of the
        two atoms that were bonded to the S atoms. The S atoms
        will be classified as DeleteIDs. The edge IDs will be
        all edges that were found when pruning the molecule.

        Parameters
        ----------
        map_file_name : str
            Output file name.
        num_atoms : int
            Total atoms in the template.
        phase_num_edge_ids : int
            Number of edge atoms.
        phase_pruned_edge_id : list
            1-indexed edge atom IDs.
        bonding_id1, bonding_id2 : int
            Atom IDs of the two bonding atoms.
        s1, s2 : int
            S-atom pair indices.
        """
        ndel = self.merged_num_del_atoms[(s1, s2)]
        path = os.path.join("rxnTemplates", map_file_name)
        with open(path, 'w') as out:
            out.write(f"# Map file: {map_file_name}\n\n")
            out.write(f"{phase_num_edge_ids} edgeIDs\n")
            out.write(f"{ndel} deleteIDs\n")
            out.write(f"{num_atoms} equivalences\n\n")

            out.write("InitiatorIDs\n\n")
            out.write(f"{bonding_id1}\n")
            out.write(f"{bonding_id2}\n\n")

            out.write("DeleteIDs\n\n")
            for did in range(1, ndel + 1):
                dt = self.pruned_del_tree[(s1, s2, did)]
                out.write(f"{dt}\n")
            out.write("\n")

            out.write("EdgeIDs\n\n")
            for eid in range(1, phase_num_edge_ids + 1):
                out.write(
                    f"{phase_pruned_edge_id[eid]}\n"
                )

            out.write("\nEquivalences\n\n")
            for i in range(1, num_atoms + 1):
                out.write(f"{i} {i}\n")


# ================================================================
# Main entry point
# ================================================================
def main():
    """Main entry point for make_reactions.py.

    Get script settings from a combination of the resource
    control file and parameters given by the user on the
    command line. Then execute the full reaction template
    pipeline.
    """
    # Get script settings from a combination of the resource
    #   control file and parameters given by the user on the
    #   command line.
    settings = ScriptSettings()

    # Create the MakeReactions object and run the pipeline.
    reactor = MakeReactions(settings)
    reactor.run()


if __name__ == '__main__':
    # Everything before this point was a class/function
    #   definition or a request to import information from
    #   external modules. Only now do we actually start running
    #   the program. The purpose of this is to allow another
    #   python program to import *this* script and call its
    #   functions internally.
    main()
