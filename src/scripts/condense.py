#!/usr/bin/env python3

"""condense.py -- Create LAMMPS condensation input files.

PROGRAM: condense.py
PURPOSE: To create all the necessary input files for running a LAMMPS
   condensation on a set of molecules with defined reaction types.

INPUT:
Expected input: Keywords that define sections.  Sections follow the
Keywords.  The Keywords may have an accompanying Value that defines the number
of lines in the rest of the section.  The Keywords and
their sections are defined as follows:

   (1) Composition:
   Keyword: "composition";
   Required Value: Following the word "composition" is an integer
      that signifies the number of different types of molecules that will be
      used in the simulation.  Each type of molecule must be present in the
      precursor database.
   Section Definition: One line for each type of molecule.  Each
      line contains a molecule name followed by a family name followed by an
      integer.
      IMPORTANT NOTE: The lines MUST be alphabetized by
      <MOLECULENAME>.
   Molecule Name: A molecule name is specified in two parts and
      referred to in total with <MOLECULENAME>.  The parts of <MOLECULENAME> are
      separated by a "_" and are given by (1) The chemical formula; and (2) an
      integer that signifies a particular type assignment for atoms in that
      molecule and/or different possible conformations or atomic arrangements
      for the same chemical formula.  In general, the precursor database is
      expected to contain a directory called <MOLECULENAME>.  Inside that
      directory *at a minimum* is a file with the name <MOLECULENAME>.skl.  So,
      for a specific example, if <MOLECULENAME> is "h2o_1" then the precursor
      database should have a "h2o_1/" directory with a file in it with the name
      "h2o_1.skl".
   Family Name: A family name is any simple string.  All
      molecules with the same family name will share the same element-species
      types.  For example, if molecules ch4 and c2h6 both have the same family
      name, then a c1 from ch4 and a c1 from c2h6 will *not* be recognized as
      different types in LAMMPS.  On the other hand, if they have different
      family names, then the two c1 atoms *will be* recognized as different
      types.  If you want a molecule to not be a part of a broader family then
      just make the molecule name and the family name equal to each other.
   Integer: The <MOLECULENAME> part of the line must be followed
      by an integer that indicates how many of that type of molecule should be
      placed in the simulation box.
   NOTE on <MOLECULENAME>: the integer part of the
      <MOLECULENAME> is used to distinguish different possible type assignments
      for atoms in the molecule so that, for example, in a simplified simulation
      all of the atoms could have type = 1.  Or, if the molecule was somewhat
      complicated, different atoms of the same element could be given different
      types to distinguish their different electronic environments.  An example
      would be b10c2h12.  This icosahedral molecule will have different types of
      B atoms because some will have C neighbors and some will not.  As another
      complicated example using b10c2h12, the positions of the 2 C atoms in the
      icosahedra can vary with three different possible options: (1) para; (2)
      meta; and (3) ortho.  Each of those molecules can be specified using a
      different integer number as part of the <MOLECULENAME>.

   (2) Cell:
   Keyword: "cell_size"
   Required Value: Floating point size of the cell in Angstroms
   Section Definition: Absent. No section information needed.

   (3) Max molecule speed:
   Keyword: "max_speed"
   Required Value: Floating point maximum speed of each molecule
   Section Definition: Absent. No section information needed.

   (4) Stages:
   Keyword: "stage"
   Required Value: Seven space-separated values:
      (a) squish_factor -- final fractional cell size after
          compression (0 = no deformation);
      (b) squish_step_size -- apply deformation every this
          many timesteps;
      (c) ensemble_type -- "nve" (Newtonian, no thermostat)
          or "nvt" (Nose-Hoover thermostat);
      (d) temp_start -- starting temperature (K);
      (e) temp_end -- ending temperature (K);
      (f) temp_damp -- thermostat damping constant (fs);
          for NVT this is the Tdamp parameter in the LAMMPS "fix nvt" command;
          ignored for NVE;
      (g) run_steps -- total number of timesteps for the
          stage.
   The "stage" keyword may be repeated multiple times in the
      input file.

   (5) Reactions:
   Keyword: "reactions"
   Required Value: Integer number of different kinds of reactions.
   Section Definition: One line with 5 components for each type
      of reaction.  Each line contains a pair of <MOLECULENAME> + element-type
      designations followed by a real number between 0.0 and 1.0 inclusive.  The
      <MOLECULENAME> values must be consistent with those specified in the
      "Components" section and the element-type values must be elements and
      types that are present in the given molecules.  The real number is the
      probability that the reaction will occur once the molecular fragments are
      sufficiently close.  So, for example a reaction leading with some
      probability (say 0.86) to the binding of two ch4 molecules would be
      expressed on one line as: "ch4_1 c-1 ch4_1 c-1 0.86".
   NOTE: The working expectation in all cases is that every
      reaction will require that an H atom is discharged from each molecule
      before the binding between the requested underlying atoms can take place.

Example input:

   composition 3
   b10c2h12 Family2 30
   c2h6 Family1 3
   ch4 Family1 8

   cell_size 100.0
   max_speed 5.0

   stage 0.33 500 nvt 800 800 500 100000

   reactions 6
   b10c2h12 b-1 b10c2h12 b-1 0.85
   b10c2h12 b-1 ch4 c-1 0.85
   b10c2h12 b-1 c2h6 c-1 0.75
   ch4 c-1 ch4 c-1 0.5
   ch4 c-1 c2h6 c-1 0.7
   c2h6 c1 c2h6 c-1 0.8


Prepare molecules and reactions.

The above input file is read and, if the precursor database does not already
contain the desired information, the "makeReactions" script is executed for each
of the molecule reaction pairs.  The makeReactions script produces a number of
output files for each pair.  The most important are the preRxn and postRxn
molecule template files that will be used to drive the molecule binding process.


Run Packmol

The skeleton files associated with the molecules listed in the condense input
file will be grabbed, converted to PDB format (if the PDB is not already
present) and sent as input to the packmol program.


Prepare LAMMPS input

   The output PDB from packmol that includes all the molecules will be converted
   into a skeleton file and then into a LAMMPS data file.  At the same time, the
   reaction templates will be gathered and compared to the LAMMPS data file to
   ensure consistency for all bond types and bond-angle types across all files.
   The results of that analysis will be used to construct a LAMMPS "in" file.  A
   slurm submission file is also created.

USAGE: condense.py [-i $inputFile] [-f] | [-help]
"""

import argparse as ap
import glob
import os
import random
import re
import shutil
import subprocess
import sys
from datetime import datetime


# ================================================================
# BondData: Read Hooke bond force constants from bonds.dat
# ================================================================

class BondData:
    """Stores Hooke bond force constants from $OLCAO_DATA/bonds.dat.

    This class reads the bond data file which contains spring constant and rest
    bond length information for pairs of elements identified by their atomic Z
    numbers.

    Attributes
    ----------
    num_hooke_bonds : int
        The number of Hooke bond entries in the database.
    hooke_bond_coeffs : list of list
        Doubly 1-indexed.  ``hooke_bond_coeffs[bond]`` is itself a
        1-indexed list ``[None, Z1, Z2, k, r0]`` where slot 1 is the
        atomic number of the first element, slot 2 the second, slot 3
        the spring constant, and slot 4 the rest bond length in
        Angstroms.  Slot 0 on every level is the unused sentinel
        matching Perl's ``$hookeBondCoeffs[$bond][1..4]`` layout.
    """

    def __init__(self):
        self.num_hooke_bonds = 0
        # 1-indexed: hooke_bond_coeffs[bond] = [None, Z1, Z2, k, r0]
        self.hooke_bond_coeffs = [None]
        self._init_data = False

    def init_bond_data(self):
        """Read the bond data file.

        The file is located at $OLCAO_DATA/bonds.dat and is
        structured with tagged sections:
          NUM_HOOKE_BONDS
          <count>
          HOOKE_BOND_COEFFS
          <Z1> <Z2> <k> <r0>  (one line per bond)

        This method is safe to call multiple times; subsequent calls return
        immediately.
        """

        # Skip if already initialized.
        if self._init_data:
            return
        self._init_data = True

        olcao_data = os.environ.get('OLCAO_DATA')
        if not olcao_data:
            sys.exit("Error: $OLCAO_DATA is not set.")
        fn = os.path.join(olcao_data, "bonds.dat")

        with open(fn, 'r') as f:
            # Read the number of bonds.
            values = self._prep_line(f)
            if values[0] != "NUM_HOOKE_BONDS":
                sys.exit(f"Expecting NUM_HOOKE_BONDS tag in {fn}")
            values = self._prep_line(f)
            self.num_hooke_bonds = int(values[0])

            # Read the Hooke bond coefficients.
            values = self._prep_line(f)
            if values[0] != "HOOKE_BOND_COEFFS":
                sys.exit(f"Expecting HOOKE_BOND_COEFFS tag in {fn}")
            for bond in range(1, self.num_hooke_bonds + 1):
                values = self._prep_line(f)
                # 1-indexed inner: slot 0 is the None sentinel, slots
                # 1..4 carry [Z1, Z2, k, r0] to match Perl's
                # ``$hookeBondCoeffs[$bond][1..4]`` convention.
                self.hooke_bond_coeffs.append([
                    None,
                    int(values[0]),    # Z of first element (slot 1)
                    int(values[1]),    # Z of second element (slot 2)
                    float(values[2]),  # k spring constant (slot 3)
                    float(values[3]),  # Rest bond length (slot 4)
                ])

    @staticmethod
    def _prep_line(f):
        """Read and split a non-empty line from file handle f.

        Skips blank lines and strips leading/trailing whitespace.
        """
        while True:
            line = f.readline()
            if not line:
                return[]
            line = line.strip()
            if line:
                return line.split()


# ================================================================
# AngleData: Read Hooke angle force constants from angles.dat
# ================================================================

class AngleData:
    """Stores Hooke angle force constants from
    $OLCAO_DATA/angles.dat.

    This class reads the angle data file which contains spring constant, rest
    angle, and tolerance information for triplets of elements identified by
    their atomic Z numbers.

    Attributes
    ----------
    num_hooke_angles : int
        The number of Hooke angle entries in the database.
    hooke_angle_coeffs : list of list
        Doubly 1-indexed.  ``hooke_angle_coeffs[angle]`` is itself a
        1-indexed list ``[None, Z1, Z_vertex, Z2, k, angle_deg,
        tolerance]``.  Slot 1 is Z of the first end atom, slot 2 is
        Z of the vertex atom, slot 3 is Z of the second end atom,
        slot 4 is the angular spring constant, slot 5 is the rest
        angle in degrees, and slot 6 is the allowed relaxation in
        degrees.  Slot 0 on every level is the unused sentinel
        matching Perl's ``$hookeAngleCoeffs[$angle][1..6]`` layout.
    """

    def __init__(self):
        self.num_hooke_angles = 0
        # 1-indexed: hooke_angle_coeffs[angle] =
        #   [None, Z1, Zv, Z2, k, angle_deg, tolerance]
        self.hooke_angle_coeffs = [None]
        self._init_data = False

    def init_angle_data(self):
        """Read the angle data file.

        The file is located at $OLCAO_DATA/angles.dat and is
        structured with tagged sections:
          NUM_HOOKE_ANGLES
          <count>
          HOOKE_ANGLE_COEFFS
          <Z1> <Zv> <Z2> <k> <angle> <tol>  (one per angle)

        This method is safe to call multiple times; subsequent calls return
        immediately.
        """

        # Skip if already initialized.
        if self._init_data:
            return
        self._init_data = True

        olcao_data = os.environ.get('OLCAO_DATA')
        if not olcao_data:
            sys.exit("Error: $OLCAO_DATA is not set.")
        fn = os.path.join(olcao_data, "angles.dat")

        with open(fn, 'r') as f:
            # Read the number of angles.
            values = self._prep_line(f)
            if values[0] != "NUM_HOOKE_ANGLES":
                sys.exit(f"Expecting NUM_HOOKE_ANGLES tag in {fn}")
            values = self._prep_line(f)
            self.num_hooke_angles = int(values[0])

            # Read the Hooke angle coefficients.
            values = self._prep_line(f)
            if values[0] != "HOOKE_ANGLE_COEFFS":
                sys.exit(f"Expecting HOOKE_ANGLE_COEFFS tag in {fn}")
            for angle in range(1, self.num_hooke_angles + 1):
                values = self._prep_line(f)
                # 1-indexed inner: slot 0 is the None sentinel, slots
                # 1..6 carry [Z1, Zv, Z2, k, angle_deg, tolerance] to
                # match Perl's ``$hookeAngleCoeffs[$angle][1..6]``
                # convention.
                self.hooke_angle_coeffs.append([
                    None,
                    int(values[0]),    # Z first element (slot 1)
                    int(values[1]),    # Z vertex element (slot 2)
                    int(values[2]),    # Z second element (slot 3)
                    float(values[3]),  # k angle spring const (slot 4)
                    float(values[4]),  # Rest angle, degrees (slot 5)
                    float(values[5]),  # Tolerance, degrees (slot 6)
                ])

    @staticmethod
    def _prep_line(f):
        """Read and split a non-empty line from file handle f."""
        while True:
            line = f.readline()
            if not line:
                return[]
            line = line.strip()
            if line:
                return line.split()


# ================================================================
# ScriptSettings: Parse rc defaults and command-line arguments
# ================================================================

class ScriptSettings:
    """Holds all user-controllable settings for condense.py.

    The instance variables of this object are the user settings that control the
    program.  The variable values are pulled from a resource control file and
    then reconciled with command line parameters.
    """

    def __init__(self):
        """Initialize settings from the rc file and command line.

        Default values come from the resource control file at
        $OLCAO_RC/condenserc.py (or a local copy in the current working
        directory).  Command line arguments override rc file defaults.
        """

        # Read default variables from the resource control file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit("Error: $OLCAO_RC is not set. " "See instructions.")
        sys.path.insert(1, rc_dir)
        from condenserc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc file.
        self.reconcile(args)

        # Record the command line parameters that were used.
        self.record_clp()

    def assign_rc_defaults(self, default_rc):
        """Assign default values from the resource control file.

        Parameters
        ----------
        default_rc : dict
            The dictionary returned by condenserc.py's parameters_and_defaults()
            function.
        """

        self.input_file = default_rc["input_file"]
        self.force_collision = default_rc["force_collision"]
        self.cell_size = default_rc["cell_size"]
        self.max_speed = default_rc["max_speed"]
        self.rxn_template_dir = default_rc["rxn_template_dir"]

    def parse_command_line(self):
        """Build and run the argument parser.

        Returns
        -------
        argparse.Namespace
            Parsed command line arguments.
        """

        prog_name = "condense.py"

        description_text = """
condense.py -- Create LAMMPS condensation input files.

Version 1.0 (Python port of the Perl condense script).

Requires: structure_control.py, element_data.py,
          bond_analysis.py, pdb2skl.py, packmol.

Creates all necessary input files for running a LAMMPS
condensation on a set of molecules with defined reaction
types.  This includes the LAMMPS data file, the LAMMPS
input file, reaction templates with normalized types,
and a slurm submission file.
"""

        epilog_text = """
Defaults are given in ./condenserc.py or \
$OLCAO_RC/condenserc.py.
"""

        parser = ap.ArgumentParser(prog=prog_name,
            formatter_class=(ap.RawDescriptionHelpFormatter),
            description=description_text, epilog=epilog_text)

        self.add_parser_arguments(parser)

        return parser.parse_args()

    def add_parser_arguments(self, parser):
        """Add all command line arguments to the parser.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            The parser to add arguments to.
        """

        # The -i option specifies the name of the input file. If not given, the
        # default value from the rc file is used (typically "condense.in").
        parser.add_argument('-i', '--input', dest='input_file', type=str,
            default=self.input_file, help=("Name of the condense input file. "
                f"Default: {self.input_file}"))

        # The -f option forces the lammps.dat input file to arrange the
        # molecules such that they are put on a forced collision course near the
        # center of the cell.  If the number of molecules is greater than three,
        # this option will not work.
        parser.add_argument('-f', '--force', dest='force_collision',
            action='store_true', default=self.force_collision,
            help=("Force collision mode: arrange molecules "
                    "on a collision course near cell center. "
                    "Only works for <= 3 molecules."))

    def reconcile(self, args):
        """Reconcile command line arguments with rc defaults.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed command line arguments.
        """
        self.input_file = args.input_file
        self.force_collision = args.force_collision

    def record_clp(self):
        """Record the command line parameters to the command file.

        Appends the date and full command line to a file named "command" in the
        current working directory.
        """
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write("Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


# ================================================================
# Condense: Main workhorse class
# ================================================================

class Condense:
    """Orchestrates the LAMMPS condensation file creation.

    This class holds all data structures needed for the condense
    workflow and implements each major step as a method:

    1. init_env()               -- Initialize databases.
    2. parse_input_file()       -- Read the condense input.
    3. compute_implicit_input() -- Derive atom/molecule counts.
    4. copy_reaction_templates()-- Copy templates from DB.
    5. run_packmol()            -- Pack molecules into the cell.
    6. create_lammps_files()    -- Write LAMMPS data/input/slurm.
    7. normalize_types()        -- Unify types across all files.
    """

    def __init__(self, settings):
        """Initialize the Condense object.

        Parameters
        ----------
        settings : ScriptSettings
            User-controllable settings parsed from the rc file and command line.
        """

        self.settings = settings

        # ----- Simulation parameters -----
        self.cell_size = settings.cell_size
        self.max_speed = settings.max_speed
        self.rxn_template_dir = settings.rxn_template_dir

        # ----- Stage data (populated by parse_input_file) -----
        self.num_stages = 0
        # All stage arrays are 1-indexed.
        self.squish_factor = [None]
        self.squish_step_size = [None]
        self.ensemble_type = [None]
        self.ensemble_temp_start = [None]
        self.ensemble_temp_end = [None]
        self.ensemble_temp_damp = [None]
        self.run_steps = [None]

        # ----- Composition data -----
        self.num_molecule_types = 0
        # 1-indexed arrays for molecule data.
        self.molecule_name = [None]
        self.family_name = [None]
        self.num_molecules = [None]
        self.molecule_to_family_map = {}
        self.num_mol_atoms = [None]

        # ----- Atom-molecule mapping -----
        # These are 1-indexed and populated by compute_implicit_input.
        self.atom_molecule_id = [None]
        self.atom_molecule_name = [None]

        # ----- Reaction data -----
        self.num_reaction_types = 0
        # rxn_mol_name[side][rxn]: side is 1 or 2.
        self.rxn_mol_name = [None,[None],[None]]
        # rxn_binding[side][rxn]: side is 1 or 2.
        self.rxn_binding = [None,[None],[None]]
        self.rxn_probability = [None]

        # ----- Database objects -----
        self.element_data = None
        self.bond_data = None
        self.angle_data = None
        self.struct = None

        # ----- Precursor database path -----
        olcao_data = os.environ.get('OLCAO_DATA', '')
        self.precursor_db = os.path.join(olcao_data, "precursorDB")

    # ============================================================
    # Step 1: Initialize the environment
    # ============================================================

    def init_env(self):
        """Initialize the element, bond, and angle databases.

        Reads the periodic table data from $OLCAO_DATA/elements.dat, the Hooke
        bond force constants from $OLCAO_DATA/bonds.dat, and the Hooke angle
        force constants from $OLCAO_DATA/angles.dat.
        """

        # Add the installation bin directory to the Python path so that we can
        # import structure_control and element_data.
        olcao_bin = os.environ.get('OLCAO_BIN', '')
        if olcao_bin and olcao_bin not in sys.path:
            sys.path.insert(0, olcao_bin)

        from element_data import ElementData
        from structure_control import StructureControl

        # Initialize the element database.
        self.element_data = ElementData()
        self.element_data.init_element_data()

        # Initialize the bond database.
        self.bond_data = BondData()
        self.bond_data.init_bond_data()

        # Initialize the angle database.
        self.angle_data = AngleData()
        self.angle_data.init_angle_data()

        # Store a reference to the StructureControl class for later
        # instantiation.
        self._sc_class = StructureControl

    # ============================================================
    # Step 2: Parse the input file
    # ============================================================

    def parse_input_file(self):
        """Read and parse the condense input file.

        The input file contains keyword-driven sections that define the
        composition, cell size, maximum speed, simulation stages, and reactions.
        See the module docstring for a complete description of the file format.
        """

        input_file = self.settings.input_file

        with open(input_file, 'r') as f:
            for line in f:
                values = line.strip().split()
                if not values:
                    continue

                keyword = values[0].lower()

                # ----- Composition -----
                # Defines the types and counts of molecules.
                if keyword == "composition":
                    self.num_molecule_types = int(values[1])

                    for mol in range(1, self.num_molecule_types + 1):
                        vals = f.readline().strip().split()
                        mol_name = vals[0].lower()
                        fam_name = vals[1].lower()
                        count = int(vals[2])

                        self.molecule_name.append(mol_name)
                        self.family_name.append(fam_name)
                        self.molecule_to_family_map[mol_name] = fam_name
                        self.num_molecules.append(count)

                # ----- Cell size -----
                elif keyword == "cell_size":
                    self.cell_size = float(values[1])

                # ----- Max speed -----
                elif keyword == "max_speed":
                    self.max_speed = float(values[1])

                # ----- Stage -----
                # Each "stage" line defines one compression/ thermostat stage.
                # Multiple stages are supported.
                elif keyword == "stage":
                    self.num_stages += 1
                    self.squish_factor.append(float(values[1]))
                    self.squish_step_size.append(int(values[2]))
                    self.ensemble_type.append(values[3])
                    self.ensemble_temp_start.append(float(values[4]))
                    self.ensemble_temp_end.append(float(values[5]))
                    self.ensemble_temp_damp.append(float(values[6]))
                    self.run_steps.append(int(values[7]))

                # ----- Reactions -----
                # Defines the types of reactions and their probabilities.
                elif keyword == "reactions":
                    self.num_reaction_types = int(values[1])

                    # Read the descriptor for each reaction that includes the
                    # name of each of the participating molecules and the type
                    # of binding that each molecule makes available for that
                    # reaction to occur.
                    for rxn in range(1, self.num_reaction_types + 1):
                        vals = f.readline().strip().split()
                        self.rxn_mol_name[1].append(vals[0].lower())
                        self.rxn_binding[1].append(vals[1].lower())
                        self.rxn_mol_name[2].append(vals[2].lower())
                        self.rxn_binding[2].append(vals[3].lower())
                        self.rxn_probability.append(float(vals[4]))

    # ============================================================
    # Step 3: Compute implicit input
    # ============================================================

    def compute_implicit_input(self):
        """Derive data not explicitly given in the input file.

        This method:
        1. Counts the number of atoms in each molecule by
           reading the corresponding skeleton (.skl) file from the precursor
           database.
        2. Assigns a molecule ID and molecule (family) name
           to every atom in the full model.
        """

        # Compute the number of atoms in each molecule type by reading each
        # molecule's skeleton file from the precursor database.
        for mol in range(1, self.num_molecule_types + 1):
            mol_name_lc = self.molecule_name[mol].lower()
            mol_file = os.path.join(self.precursor_db, mol_name_lc,
                f"{mol_name_lc}.skl")

            curr_num_atoms = 0
            with open(mol_file, 'r') as skl:
                for skl_line in skl:
                    if (skl_line.strip().startswith("cart")
                            or skl_line.strip().startswith(
                                "frac")):
                        vals = skl_line.strip().split()
                        curr_num_atoms = int(vals[1])
                        break

            self.num_mol_atoms.append(curr_num_atoms)

        # Now, assign a molecule number and family name to each atom in the
        # whole model.  The assumption is that atoms are ordered by molecule
        # type and then by repeat molecule within each type.  Also build a
        # 1-indexed map from molecule instance number to molecule type index so
        # that later code (e.g., velocity re-assignment after minimization) can
        # look up the atom count per molecule.
        atom_count = 0
        mol_count = 0
        self.molecule_type_of_mol = [None]
        for mol in range(1, self.num_molecule_types + 1):
            for repeat_mol in range(1, self.num_molecules[mol] + 1):
                mol_count += 1
                self.molecule_type_of_mol.append(mol)
                for atom in range(1, self.num_mol_atoms[mol] + 1):
                    atom_count += 1
                    self.atom_molecule_id.append(mol_count)
                    # Use the family name rather than the molecule name so that
                    # species types are shared within a family.
                    self.atom_molecule_name.append(self.family_name[mol])
        self.total_num_molecules = mol_count

    # ============================================================
    # Step 4: Copy reaction templates
    # ============================================================

    def copy_reaction_templates(self):
        """Copy reaction template files from the precursor DB.

        A database (called the precursorDB) should already contain certain very
        useful input files for the condensation process.  Therefore, the desired
        files can be copied from the database into the reaction templates
        directory, the packmol directory, and the lammps directory.  None of the
        files should need to be generated "on the fly" for this specific run of
        the condense script.

        The specific files that are needed can be divided into
        two groups: molecule files and reaction files.
        - Molecule Files: A PDB format file of each molecule is
          needed to build the input for the packmol program.  An OLCAO-style
          skeleton file is not needed explicitly, but it is expected to be
          present because it is used to generate all other files.
        - Reaction Files: For each pair of named molecules in a
          reaction a set of pre- and post-reaction template files are needed.
          Similarly, a map file that connects them is also needed.

        IMPORTANT NOTE! The reaction templates are *NOT* exactly the same as the
        files that will be used by LAMMPS.  The difference is that the type
        assignments in the template files and the type assignments in the LAMMPS
        input/data files may be (likely are) different.  The types from each of
        the reaction template pairs need to be assigned in a uniform consistent
        way across all files.  Hence the need for calling the normalize_types
        method, etc.
        """

        rxn_dir = self.rxn_template_dir

        # Create a directory for the combined set of reaction templates if it
        # does not already exist.  (It might already exist if this condense
        # script is being called as a part of an iterative refinement process.)
        os.makedirs(rxn_dir, exist_ok=True)

        # For each reaction, copy the reaction templates and map files from the
        # precursor database.
        print(f"numReactionTypes {self.num_reaction_types}")
        for rxn in range(1, self.num_reaction_types + 1):
            # Construct the name for the current reaction molecule pair ensuring
            # that the name will be uniquely ordered (alphabetically).
            name1 = self.rxn_mol_name[1][rxn]
            name2 = self.rxn_mol_name[2][rxn]
            if name1 < name2:
                curr_rxn_mol_pair = f"{name1}__{name2}"
            else:
                curr_rxn_mol_pair = f"{name2}__{name1}"

            # Check that the currMolRxnPair actually exists in the database.
            pair_path = os.path.join(self.precursor_db, curr_rxn_mol_pair)
            if not os.path.isdir(pair_path):
                sys.exit(f"Unable to find {pair_path}.")

            # For the current reaction, grab the appropriate reaction templates
            # and store them in the collective rxn_template_dir directory with
            # appropriate names to distinguish them.
            bind1 = self.rxn_binding[1][rxn]
            bind2 = self.rxn_binding[2][rxn]
            rxn_binding_pair = f"{bind1}_{bind2}"
            rxn_mol_binding_pair = f"{name1}_{bind1}_{name2}_{bind2}"

            # Copy all the pre- and post- template files.
            rxn_tmpl_src = os.path.join(pair_path, "rxnTemplates")

            shutil.copy2(
                os.path.join(rxn_tmpl_src, f"preRxn.{rxn_binding_pair}.data"),
                os.path.join(rxn_dir, f"preRxn.{rxn_mol_binding_pair}.data"))
            shutil.copy2(
                os.path.join(rxn_tmpl_src, f"postRxn.{rxn_binding_pair}.data"),
                os.path.join(rxn_dir, f"postRxn.{rxn_mol_binding_pair}.data"))

            # Copy the map files.
            shutil.copy2(os.path.join(rxn_tmpl_src, f"{rxn_binding_pair}.map"),
                os.path.join(rxn_dir, f"{rxn_mol_binding_pair}.map"))

    # ============================================================
    # Step 5: Run packmol
    # ============================================================

    def run_packmol(self):
        """Assemble the packmol input and run packmol.

        This method:
        1. Creates a "packmol" subdirectory.
        2. Writes a packmol input file that places the requested
           number of each molecule type inside a cubic cell.
        3. Copies the necessary PDB files from the precursor
           database.
        4. Executes the packmol program.
        5. Fixes the packmol PDB output to comply with the PDB
           standard as specified at
           http://www.wwpdb.org/documentation/file-format.

        The packmol program itself is an external dependency that must be
        installed and available on the system PATH.  There are lots of options
        for packmol.  It may be necessary to enforce a distance between
        molecules so that no inter-molecular "bonds" or "bond angles" are
        defined between the molecules before the simulation starts.
        """

        # Create a directory to run packmol in.
        os.makedirs("packmol", exist_ok=True)

        packmol_file = os.path.join("packmol", "packmol.in")

        # Open the packmol input file to be written.
        with open(packmol_file, 'w') as pack:
            # Add header comments.
            pack.write("#Packmol input file from condense.\n\n")

            # Insert the random seed and tolerance.
            rand_seed = random.randint(0, 99999)
            pack.write(f"seed {rand_seed}\n")
            pack.write("tolerance 5.0\n")

            # The input and output file types for the molecules are PDB.
            pack.write("filetype pdb\n")
            pack.write("output packmol.raw.pdb\n\n")

            for mol in range(1, self.num_molecule_types + 1):
                # Do not put a molecule into the packmol input file if there are
                # zero of them requested.  This might occur if a molecule can be
                # *assembled* from other components, but should not be present
                # in the initial system.
                if self.num_molecules[mol] == 0:
                    continue

                name = self.molecule_name[mol]
                pack.write(f"structure {name}.pdb\n")
                pack.write(f"  number {self.num_molecules[mol]}\n")
                # Pad the placement region so that no molecule straddles the
                # periodic boundary. Without padding, a molecule whose center is
                # near the box edge can have atoms that wrap around, causing
                # LAMMPS "Inconsistent image flags" errors and breaking fix
                # bond/react.
                pad = 5.0
                lo = pad
                hi = self.cell_size - pad
                pack.write(f"  inside cube {lo} {lo} {lo} {hi}\n")
                pack.write("end structure\n\n")

        # Copy in the necessary PDB files.
        for mol in range(1, self.num_molecule_types + 1):
            # As above, we may have certain molecules that are not present in
            # the initial structure but which can be made during the lammps
            # condensation process itself.
            if self.num_molecules[mol] == 0:
                continue

            name = self.molecule_name[mol].lower()
            mol_pdb = os.path.join(self.precursor_db, name, f"{name}.pdb")
            shutil.copy2(mol_pdb, os.path.join("packmol", f"{name}.pdb"))

        # Execute packmol from inside the packmol directory so that relative PDB
        # filenames resolve correctly.
        subprocess.run("packmol < packmol.in", shell=True, check=True,
            cwd="packmol")

        # Unfortunately, the packmol program seems to not follow the PDB
        # standard as specified here:
        #   http://www.wwpdb.org/documentation/file-format.
        # Therefore, we will now correct the PDB file that was output by
        # packmol.
        raw_pdb = os.path.join("packmol", "packmol.raw.pdb")
        fixed_pdb = os.path.join("packmol", "packmol.fixed.pdb")

        with open(raw_pdb, 'r') as raw, \
                open(fixed_pdb, 'w') as fixed:
            for line in raw:
                values = line.strip().split()
                if not values:
                    fixed.write(line)
                    continue

                if values[0] == "ATOM":
                    # Strip non-alphabetic characters from the element name
                    # (last field).
                    element = re.sub(r'[^a-zA-Z]', '', values[-1])
                    # Write first 54 characters of the original line, then
                    # occupancy and temperature factor, then element symbol
                    # right-justified in 12 chars.
                    fixed.write(line[:54])
                    fixed.write("  1.00  0.00")
                    fixed.write(f"{element:>12s}")
                    fixed.write("\n")
                else:
                    fixed.write(line)

    # ============================================================
    # Step 6: Create LAMMPS files
    # ============================================================

    def create_lammps_files(self):
        """Create the LAMMPS data, input, and slurm files.

        This is the largest and most complex step.  It:
        1. Converts the packmol PDB into a skeleton file using
           the pdb2skl.py script.
        2. Reads the skeleton file via StructureControl.
        3. Sets the lattice to the requested cubic cell size.
        4. Runs bond_analysis.py to compute bond lengths and
           bond angles.
        5. Computes ordered species assignments where each
           unique combination of element + species number + molecule family is a
           distinct LAMMPS atom type.
        6. Writes the LAMMPS data file (lammps.dat) with atoms,
           bonds, angles, masses, pair coefficients, bond coefficients, and
           angle coefficients.
        7. Writes the LAMMPS input file (lammps.in) with
           initialization, molecule definitions, bond reactions, simulation
           stages, and energy minimization.
        8. Writes a slurm submission file.
        """

        # Create a directory for the LAMMPS data files and enter it.  The Perl
        # script did chdir("lammps") for this entire block so that all byproduct
        # files (sginput, sgoutput, skl2PDB.map, bondAnalysis.*) land inside the
        # lammps directory.
        os.makedirs("lammps", exist_ok=True)
        saved_cwd = os.getcwd()
        # Absolutize the reaction template directory before changing into
        # lammps/ so the relative path still resolves correctly.
        self.rxn_template_dir = os.path.abspath(self.rxn_template_dir)
        os.chdir("lammps")

        # Convert the packmol PDB into a skeleton file.
        shutil.copy2(os.path.join(saved_cwd, "packmol", "packmol.fixed.pdb"),
            "packmol.fixed.pdb")
        subprocess.run(
            [
                "pdb2skl.py",
                "-i", "packmol.fixed.pdb",
                "-o", "packmol.skl",
                "--pdbtypes",
            ],
            check=True,
        )

        # Read the skeleton file.
        self.struct = self._sc_class()
        self.struct.read_input_file("packmol.skl", use_file_species=True)

        # Adjust the StructureControl real and reciprocal cell vector magnitudes
        # to the requested cubic cell.
        self.struct.set_lattice_from_mag_angle(
            [None, self.cell_size, self.cell_size, self.cell_size],
            [None, 90.0, 90.0, 90.0])

        # Extract necessary information.
        num_atoms = self.struct.num_atoms

        # Perform a bond analysis to get bond length and bond angle information.
        # Use the Python bond_analysis.py script instead of the Perl
        # bondAnalysis.
        subprocess.run(
            ["bond_analysis.py", "-bl", "-bf", "1.1", "-i", "packmol.skl",],
            check=True)
        subprocess.run(
            ["bond_analysis.py", "-ba", "-bf", "1.1", "-i", "packmol.skl",],
            check=True)

        # Read the bondAnalysis results and extract necessary information.
        self.struct.read_bond_analysis_bl("bondAnalysis.bl")
        self.struct.read_bond_analysis_ba("bondAnalysis.ba")

        # Convenient aliases for the structure data.
        bonded = self.struct.bonded
        bond_length_ext = self.struct.bond_length_ext
        bond_tag_id_struct = self.struct.bond_tag_id
        num_bonds = self.struct.num_bonds
        num_bonds_total = self.struct.num_bonds_total
        num_bond_angles = self.struct.num_bond_angles
        angle_bonded = self.struct.angle_bonded
        bond_angles_ext = self.struct.bond_angles_ext
        num_angles_total = self.struct.num_angles_total
        atom_element_name = self.struct.atom_element_name
        atom_species_id = self.struct.atom_species_id
        atomic_z = self.struct.atomic_z
        direct_xyz = self.struct.direct_xyz

        ed = self.element_data
        bd = self.bond_data
        ad = self.angle_data

        # ------------------------------------------------
        # Assign an *ordered* species number to each atom
        # ------------------------------------------------
        # This is computed according to the element, species number, and
        # molecule to which the atom belongs.  Note that this is not a species
        # number that is nested within the element ID number.  For example, if a
        # system has one type of molecule that contains 3 Si species and 2 O
        # species and 90 atoms then the ordered_species_id for each of the 90
        # atoms will be a number between 1 and 5 inclusive.  As another example
        # consider a system that contains H2O and CO2.  The O from H2O and the O
        # from CO2 will need to be identified as different species.

        num_ordered_species = 0
        num_unique_atom_tags = 0
        unique_atom_tags = [None]   # 1-indexed
        ordered_species_id = [None]  # 1-indexed

        for atom in range(1, num_atoms + 1):
            # Assemble an element-species-molecule tag for this atom.
            atom_tag = (
                f"{atom_element_name[atom]} "
                f"{atom_species_id[atom]} "
                f"{self.atom_molecule_name[atom]}"
            )

            # Determine if this atom tag has been seen before.
            found = 0
            for tag_idx in range(1, num_unique_atom_tags + 1):
                if atom_tag == unique_atom_tags[tag_idx]:
                    found = tag_idx
                    break

            if found == 0:
                num_unique_atom_tags += 1
                unique_atom_tags.append(atom_tag)
                ordered_species_id.append(num_unique_atom_tags)
            else:
                ordered_species_id.append(found)

        num_ordered_species = num_unique_atom_tags

        # ------------------------------------------------
        # Establish unique ordered species and their info
        # ------------------------------------------------
        ordered_species_tag = [None] * (num_ordered_species + 1)
        ordered_species_element = [None] * (num_ordered_species + 1)
        ordered_species_masses = [None] * (num_ordered_species + 1)
        # ordered_species_pair_coeffs[sp] = [eps, sigma]
        ordered_species_pair_coeffs = [None] * (num_ordered_species + 1)

        bond_count = 0
        angle_count = 0
        num_unique_bond_tags = 0
        num_unique_angle_tags = 0
        unique_bond_tags_local = [None]
        unique_angle_tags_local = [None]
        # bond_tag_id_local[atom][bond_idx] = tag id
        bond_tag_id_local = [None] * (num_atoms + 1)
        # angle_tag_id_local[atom][angle_idx] = tag id
        angle_tag_id_local = [None] * (num_atoms + 1)
        ordered_bond_type = [None]
        ordered_bonded_atoms = [None]
        ordered_angle_type = [None]
        angle_bonded_atoms = [None]
        unique_bond_coeffs = [None]
        unique_angle_coeffs = [None]

        for atom in range(1, num_atoms + 1):
            # Initialize per-atom bond/angle tag arrays.
            if bond_tag_id_local[atom] is None:
                bond_tag_id_local[atom] = [None]
            if angle_tag_id_local[atom] is None:
                angle_tag_id_local[atom] = [None]

            # --- Populate species info on first encounter
            sp = ordered_species_id[atom]
            if ordered_species_tag[sp] is None:
                # Element-species-molecule tag.
                ordered_species_tag[sp] = (
                    f"{atom_element_name[atom]} "
                    f"{atom_species_id[atom]} "
                    f"{self.atom_molecule_name[atom]}"
                )

                # Element name for OVITO visualization.
                ordered_species_element[sp] = (
                    atom_element_name[atom].capitalize()
                )

                # Atomic mass.
                z = atomic_z[atom]
                ordered_species_masses[sp] = ed.atomic_masses[z]

                # LJ pair coefficients.  The true LJ interaction is computed
                # from the combination of coefficients from different elements.
                ordered_species_pair_coeffs[sp] = [
                    ed.lj_pair_coeffs[z][0],
                    ed.lj_pair_coeffs[z][1],
                ]

            # --- Process bonds ---
            atom_num_bonds = (
                num_bonds[atom]
                if num_bonds[atom] is not None
                else 0
            )
            for bond_idx in range(1, atom_num_bonds + 1):
                bonded_atom = bonded[atom][bond_idx]

                # Only count each bond once (atom < bonded).
                if atom > bonded_atom:
                    bond_tag_id_local[atom].append(None)
                    continue

                bond_count += 1

                # Create a unique element-species-molecule tag for this bond.
                tag1 = (
                    f"{atom_element_name[atom]} "
                    f"{atom_species_id[atom]} "
                    f"{self.atom_molecule_name[atom]}"
                )
                tag2 = (
                    f"{atom_element_name[bonded_atom]} "
                    f"{atom_species_id[bonded_atom]} "
                    f"{self.atom_molecule_name[bonded_atom]}"
                )
                if tag1 < tag2:
                    current_tag = f"{tag1} {tag2}"
                else:
                    current_tag = f"{tag2} {tag1}"

                # Determine if this tag is unique.
                found = 0
                for t in range(1, num_unique_bond_tags + 1):
                    if (current_tag
                            == unique_bond_tags_local[t]):
                        found = t
                        break

                if found == 0:
                    num_unique_bond_tags += 1
                    unique_bond_tags_local.append(current_tag)
                    found = num_unique_bond_tags

                bond_tag_id_local[atom].append(found)

                # Extract the Z numbers by search because the tag was
                # alphabetized so there is no guarantee that the first part of
                # the tag corresponds to atom.
                parts = current_tag.split()
                atom1_z = ed.get_element_z(parts[0])
                atom2_z = ed.get_element_z(parts[3])

                # Make sure atom1_z <= atom2_z.
                if atom1_z > atom2_z:
                    atom1_z, atom2_z = atom2_z, atom1_z

                # Record the bonded atoms and bond types.  The inner
                # pair is 1-indexed ``[None, atom, bonded_atom]`` to
                # match Perl's
                # ``$orderedBondedAtoms[$bondCount][1..2]``.
                ordered_bonded_atoms.append([None, atom, bonded_atom])
                ordered_bond_type.append(found)

                # Collect the coefficients for each unique bond type.  This
                # should work if bonds.dat always has the lower Z atom listed
                # first.  hbc is 1-indexed: slots 1..4 are Z1, Z2, k, r0.
                # The stored entry is 1-indexed ``[None, k, r0]`` to
                # match Perl's ``$uniqueBondCoeffs[$bond][1..2]``.
                for hb in range(1, bd.num_hooke_bonds + 1):
                    hbc = bd.hooke_bond_coeffs[hb]
                    if (atom1_z == hbc[1]
                            and atom2_z == hbc[2]):
                        # Ensure the unique_bond_coeffs list is large enough.
                        while len(unique_bond_coeffs) <= (found):
                            unique_bond_coeffs.append(None)
                        unique_bond_coeffs[found] = [None, hbc[3], hbc[4]]

            # --- Process angles ---
            atom_num_angles = (
                num_bond_angles[atom]
                if (atom < len(num_bond_angles)
                    and num_bond_angles[atom] is not None)
                else 0
            )
            for angle_idx in range(1, atom_num_angles + 1):
                # Note that bond angles were not double listed in the
                # bondAnalysis.ba file (unlike the double listed bonds in
                # bondAnalysis.bl).  Thus, we don't need to check for double
                # counted bond angles here.

                angle_count += 1

                # Get the current tag for this bond angle and then extract the
                # atomic Z numbers for the vertex atom and the non-vertex atoms.
                # The vertex atom is $atom.
                ab = angle_bonded[atom]
                tag_part2 = (
                    f"{atom_element_name[atom]} "
                    f"{atom_species_id[atom]} "
                    f"{self.atom_molecule_name[atom]}"
                )
                a1 = ab[angle_idx][1]
                a2 = ab[angle_idx][2]
                tag_part1 = (
                    f"{atom_element_name[a1]} "
                    f"{atom_species_id[a1]} "
                    f"{self.atom_molecule_name[a1]}"
                )
                tag_part3 = (
                    f"{atom_element_name[a2]} "
                    f"{atom_species_id[a2]} "
                    f"{self.atom_molecule_name[a2]}"
                )
                if tag_part1 < tag_part3:
                    current_tag = f"{tag_part1} {tag_part2} {tag_part3}"
                else:
                    current_tag = f"{tag_part3} {tag_part2} {tag_part1}"

                # Extract the Z numbers.
                parts = current_tag.split()
                atom1_z = ed.get_element_z(parts[0])
                atom2_z = ed.get_element_z(parts[3])
                atom3_z = ed.get_element_z(parts[6])

                # Make sure atom1_z <= atom3_z.
                if atom1_z > atom3_z:
                    atom1_z, atom3_z = atom3_z, atom1_z

                # Find the angle in the database.  hac is 1-indexed:
                # slots 1..6 are Z1, Zv, Z2, k, angle_deg, tolerance.
                found_angle = 0
                for ha in range(1, ad.num_hooke_angles + 1):
                    hac = ad.hooke_angle_coeffs[ha]
                    if (atom1_z == hac[1]
                            and atom2_z == hac[2]
                            and atom3_z == hac[3]):
                        angle_val = bond_angles_ext[atom][angle_idx]
                        if (abs(angle_val - hac[5])
                                <= hac[6]):
                            found_angle = ha
                            break

                if found_angle == 0:
                    angle_val = bond_angles_ext[atom][angle_idx]
                    print(f"Angle number = {angle_idx} Angle = {angle_val}")
                    print(f"atomZSet = {atom1_z} {atom2_z} {atom3_z}")
                    sys.exit("Cannot find angle in the " "database")
                else:
                    hac = ad.hooke_angle_coeffs[found_angle]
                    # Append the rest angle (slot 5) and the Hooke ID.
                    current_tag = f"{current_tag} {hac[5]} {found_angle}"

                # Determine if this tag is unique.
                found = 0
                for t in range(1, num_unique_angle_tags + 1):
                    if current_tag == unique_angle_tags_local[ t]:
                        found = t
                        break

                if found == 0:
                    num_unique_angle_tags += 1
                    unique_angle_tags_local.append(current_tag)
                    found = num_unique_angle_tags

                # Initialize per-atom angle tag if needed.
                if angle_tag_id_local[atom] is None:
                    angle_tag_id_local[atom] = [None]
                angle_tag_id_local[atom].append(found)

                # Collect the coefficients for each unique bond angle type.
                # hac is 1-indexed: slots 1..6 are Z1, Zv, Z2, k,
                # angle_deg, tolerance.  The stored entry is 1-indexed
                # ``[None, k, angle_deg]`` to match Perl's
                # ``$uniqueAngleCoeffs[$angle][1..2]``.
                for ha in range(1, ad.num_hooke_angles + 1):
                    hac = ad.hooke_angle_coeffs[ha]
                    if (((atom1_z == hac[1]
                          and atom2_z == hac[2]
                          and atom3_z == hac[3])
                         or (atom3_z == hac[1]
                             and atom2_z == hac[2]
                             and atom1_z == hac[3]))):
                        while len(unique_angle_coeffs) <= (found):
                            unique_angle_coeffs.append(None)
                        unique_angle_coeffs[found] = [None, hac[4], hac[5]]

                # The inner triple is 1-indexed
                # ``[None, a1, atom, a2]`` to match Perl's
                # ``$angleBondedAtoms[$angleCount][1..3]`` layout.
                # Slot 1 is the first end atom, slot 2 is the
                # *vertex* atom, slot 3 is the second end atom —
                # same role ordering as Perl.
                angle_bonded_atoms.append([None, a1, atom, a2])
                ordered_angle_type.append(found)

        # ------------------------------------------------
        # Write the LAMMPS data file (lammps.dat)
        # ------------------------------------------------
        lammps_dat = "lammps.dat"
        with open(lammps_dat, 'w') as lmp:
            # Print the header.
            lmp.write("lammps.dat\n\n")

            # Print number of items and number of types.
            lmp.write(f"{num_atoms} atoms\n")
            lmp.write(f"{bond_count} bonds\n")
            lmp.write(f"{angle_count} angles\n")
            lmp.write("0 dihedrals\n")
            lmp.write("0 impropers\n\n")
            lmp.write(f"{num_ordered_species} atom types\n")
            lmp.write(f"{num_unique_bond_tags} bond types\n")
            lmp.write(f"{num_unique_angle_tags} angle types\n\n")
            cs = self.cell_size
            lmp.write(f"0.000 {cs} xlo xhi\n")
            lmp.write(f"0.000 {cs} ylo yhi\n")
            lmp.write(f"0.000 {cs} zlo zhi\n\n")

            # Print the Mass information for each atom type.
            lmp.write("Masses\n\n")
            for sp in range(1, num_ordered_species + 1):
                lmp.write(
                    f"{sp} "
                    f"{ordered_species_masses[sp]} "
                    f"# {ordered_species_tag[sp]}\n"
                )

            # Print the pair coefficients for the LJ interaction of each atom
            # type.
            lmp.write("\nPair Coeffs\n\n")
            for sp in range(1, num_ordered_species + 1):
                pc = ordered_species_pair_coeffs[sp]
                lmp.write(
                    f"{sp} {pc[0]} {pc[1]} "
                    f"# {ordered_species_tag[sp]}\n"
                )

            # Print the bond coefficients for the spring-like interaction of
            # each bond type.  bc is 1-indexed: slot 1 is k, slot 2 is r0.
            lmp.write("\nBond Coeffs\n\n")
            for b in range(1, num_unique_bond_tags + 1):
                bc = unique_bond_coeffs[b]
                lmp.write(
                    f"{b} {bc[1]} {bc[2]} "
                    f"# {unique_bond_tags_local[b]}\n"
                )

            # Print the bond angle coefficients for the spring-like
            # interaction.  ac is 1-indexed: slot 1 is k, slot 2 is
            # the rest angle in degrees.
            lmp.write("\nAngle Coeffs\n\n")
            for a in range(1, num_unique_angle_tags + 1):
                ac = unique_angle_coeffs[a]
                lmp.write(
                    f"{a} {ac[1]} {ac[2]} "
                    f"# {unique_angle_tags_local[a]}\n"
                )

            # If we are going to force a collision, we first translate the
            # molecules into position before printing.
            if self.settings.force_collision:
                current_mol = 0
                for atom in range(1, num_atoms + 1):
                    if (self.atom_molecule_id[atom]
                            != current_mol):
                        current_mol = self.atom_molecule_id[atom]

                        # Find the atom that starts the next molecule after this
                        # one.
                        next_mol_atom = atom
                        for na in range(atom + 1, num_atoms + 1):
                            next_mol_atom = na
                            if (self.atom_molecule_id[na]
                                    != current_mol):
                                next_mol_atom = na - 1
                                break

                        # Shift the molecule into the center.
                        self.struct.get_min_max_xyz(atom, next_mol_atom)
                        self.struct.shift_xyz_center(atom, next_mol_atom)

                        # Compute the x-axis distance to shift.
                        if current_mol == 1:
                            trans_dist = [None, 0.0, 0.0, 0.0]
                        elif current_mol == 2:
                            trans_dist = [None, cs * -0.25, 0.0, 0.0,]
                        elif current_mol == 3:
                            trans_dist = [None, cs * 0.25, 0.0, 0.0,]
                        else:
                            print(
                                "Force doesn't work for "
                                "more than 3 molecules."
                            )
                            sys.exit(1)

                        self.struct.translate_atoms(atom, next_mol_atom,
                            trans_dist)

            # Print the atom information.
            lmp.write("\nAtoms\n\n")
            for atom in range(1, num_atoms + 1):
                xyz = direct_xyz[atom]
                lmp.write(
                    f"{atom} "
                    f"{self.atom_molecule_id[atom]} "
                    f"{ordered_species_id[atom]} 0.00 "
                    f" {xyz[1]} {xyz[2]} {xyz[3]}"
                    f" # {atom_element_name[atom]} "
                    f"{atom_species_id[atom]} "
                    f"{self.atom_molecule_name[atom]}\n"
                )

            # Print initial velocities for all atoms as zero.  The actual
            # per-molecule translational velocities are assigned later in
            # lammps.in after a force relaxation (minimize) step.
            lmp.write("\nVelocities\n\n")
            for atom in range(1, num_atoms + 1):
                lmp.write(
                    f"{atom} 0.0 0.0 0.0"
                    f"  # {atom_element_name[atom]} "
                    f"{atom_species_id[atom]} "
                    f"{self.atom_molecule_name[atom]}\n"
                )

            # Print the bond information.  ba is 1-indexed:
            # slot 1 is the atom, slot 2 is the bonded atom.
            lmp.write("\nBonds\n\n")
            for bond in range(1, bond_count + 1):
                bt = ordered_bond_type[bond]
                ba = ordered_bonded_atoms[bond]
                lmp.write(
                    f"{bond} {bt} "
                    f"{ba[1]} {ba[2]} # "
                    f"{unique_bond_tags_local[bt]}\n"
                )

            # Print the bond angle information.  aa is 1-indexed:
            # slot 1 is the first end atom, slot 2 is the vertex,
            # slot 3 is the second end atom (matching Perl).
            lmp.write("\nAngles\n\n")
            for angle in range(1, angle_count + 1):
                at = ordered_angle_type[angle]
                aa = angle_bonded_atoms[angle]
                lmp.write(
                    f"{angle} {at} "
                    f"{aa[1]} {aa[2]} {aa[3]} # "
                    f"{unique_angle_tags_local[at]}\n"
                )

        # ------------------------------------------------
        # Write the LAMMPS input file (lammps.in)
        # ------------------------------------------------
        lammps_in = "lammps.in"
        with open(lammps_in, 'w') as lmpin:
            # Header with initialization settings.
            lmpin.write("""
# A lammps command file for condensing a periodic cell.

# Initialization

units real
dimension 3
boundary p p p
atom_style full
pair_style lj/cut 10.0
bond_style harmonic
angle_style harmonic
neigh_modify delay 100 every 50 check yes
#comm_modify mode single cutoff 20.0
pair_modify shift yes mix sixthpower
newton on

# Atom definition
read_data lammps.dat extra/bond/per/atom 25 extra/special/per/atom 25

special_bonds lj/coul 0 1 1

########################
# Initialization
########################

timestep 0.5

# Set up output

#restart 50000 restart.fast1.0 restart.fast1.1

dump coord_dump all custom 500 dump.coarse id type element x y z
""")

            # dump_modify element list.
            lmpin.write("dump_modify coord_dump element ")
            for sp in range(1, num_ordered_species + 1):
                lmpin.write(f" {ordered_species_element[sp]}")
            lmpin.write("\n")

            lmpin.write("""
dump_modify coord_dump sort id


compute     bond all property/local batom1 batom2 btype
compute     dist all bond/local dist
dump        bond_dump all local 500 dump.bond.coarse c_bond[1] c_bond[2] c_bond[3] c_dist
dump_modify bond_dump &
              label BONDS &
              colname 1 batom1 &
              colname 2 batom2 &
              colname 3 btype &
              colname 4 dist

#compute ang all property/local aatom1 aatom2 aatom3 atype
#dump angles all local 500 dump.angle.coarse index c_ang[1] c_ang[2] c_ang[3] c_ang[4]


########################
# Coarse run
########################

region simcell block EDGE EDGE EDGE EDGE EDGE EDGE units box

""")

            # Insert the molecule definitions.  There will be eight for every
            # reaction because there are four phases to a reaction and a pre and
            # post aspect to each phase.
            lmpin.write("# Molecules\n")
            for rxn in range(1, self.num_reaction_types + 1):
                n1 = self.rxn_mol_name[1][rxn]
                b1 = self.rxn_binding[1][rxn]
                n2 = self.rxn_mol_name[2][rxn]
                b2 = self.rxn_binding[2][rxn]
                pre_file = f"preRxn.{n1}_{b1}_{n2}_{b2}.data"
                post_file = f"postRxn.{n1}_{b1}_{n2}_{b2}.data"
                rxn_dir = self.rxn_template_dir
                shutil.copy2(os.path.join(rxn_dir, pre_file), pre_file,)
                lmpin.write(f"molecule MOLpre{rxn} {pre_file}\n")
                shutil.copy2(os.path.join(rxn_dir, post_file), post_file,)
                lmpin.write(f"molecule MOLpost{rxn} {post_file}\n")

            # Insert the bond reactions.
            lmpin.write("\n# Bond Reactions\n")
            lmpin.write("fix reaction all bond/react " "stabilization no &\n")
            for rxn in range(1, self.num_reaction_types + 1):
                # Prepare the map file and print the react command for each
                # reaction.
                rand_seed = random.randint(0, 99999)
                n1 = self.rxn_mol_name[1][rxn]
                b1 = self.rxn_binding[1][rxn]
                n2 = self.rxn_mol_name[2][rxn]
                b2 = self.rxn_binding[2][rxn]
                map_file = f"{n1}_{b1}_{n2}_{b2}.map"
                rxn_dir = self.rxn_template_dir
                shutil.copy2(os.path.join(rxn_dir, map_file), map_file,)
                prob = self.rxn_probability[rxn]
                line = (
                    f"  react RXN{rxn} all 100 0.0 "
                    f"3.0 MOLpre{rxn} MOLpost{rxn} "
                    f"{map_file} prob {prob} "
                    f"{rand_seed}"
                )
                if rxn < self.num_reaction_types:
                    lmpin.write(f"{line} &\n")
                else:
                    lmpin.write(f"{line}\n")

            # Minimize molecular geometry first to relax the enormous E_mol from
            # packmol coordinates. Without this, the spring forces are so large
            # that atoms fly apart on the first step. Minimize zeros all
            # velocities, so afterward we re-assign per-molecule rigid-body
            # translations using LAMMPS velocity/set commands grouped by
            # molecule ID.
            lmpin.write("\n# Relax molecular geometry\n"
                "minimize 1e-6 1e-8 10000 100000\n\n")

            # Re-assign per-molecule rigid-body velocities.  We define a
            # temporary group per molecule, assign its velocity, then delete the
            # group.  In normal mode each molecule gets a uniform random
            # velocity with components drawn from [-max_speed, +max_speed].  In
            # force_collision mode (-f flag), hardcoded velocities aim the
            # molecules toward each other (max 3 molecules).
            lmpin.write("# Per-molecule translational " "velocities\n")
            ms = self.max_speed
            current_atom = 1
            for mol in range(1, self.total_num_molecules + 1):
                mol_type_idx = self.molecule_type_of_mol[mol]
                n_atoms = self.num_mol_atoms[mol_type_idx]

                if not self.settings.force_collision:
                    # Normal mode: random velocity for each molecule.
                    vx = random.uniform(-ms, ms)
                    vy = random.uniform(-ms, ms)
                    vz = random.uniform(-ms, ms)
                else:
                    # Force collision mode: aim molecules toward each other
                    # along x-axis.
                    if mol == 1:
                        vx, vy, vz = 0.0, 0.0, 0.0
                    elif mol == 2:
                        vx, vy, vz = 1.0, 0.0, 0.0
                    elif mol == 3:
                        vx, vy, vz = -0.7, 0.0, 0.0
                    else:
                        print("Force collision mode does not work for more "
                            "than three molecules.")
                        sys.exit(1)

                lo = current_atom
                hi = current_atom + n_atoms - 1
                lmpin.write(
                    f"group mol_tmp id {lo}:{hi}\n"
                    f"velocity mol_tmp set {vx} {vy} {vz} units box\n"
                    f"group mol_tmp delete\n"
                )
                current_atom += n_atoms
            lmpin.write("\n")

            # Fill the LAMMPS input file with the requested stages.  Each stage
            # specifies an ensemble type (nve or nvt) and optional box
            # deformation.  NVE uses pure Newtonian dynamics (no thermostat) so
            # translational velocities are conserved.  NVT applies a Nose-Hoover
            # thermostat to control temperature between temp_start and temp_end
            # with the given damping constant.
            for stage in range(1, self.num_stages + 1):
                lmpin.write(f"# Stage {stage}\n")

                # Turn on deformation if requested.
                sf = self.squish_factor[stage]
                if sf > 0:
                    ss = self.squish_step_size[stage]
                    lmpin.write(
                        f"fix squish all deform {ss}"
                        f" x scale {sf}"
                        f" y scale {sf}"
                        f" z scale {sf} remap x\n"
                    )

                # Write the ensemble fix.
                ens = self.ensemble_type[stage].lower()
                if ens == "nve":
                    lmpin.write("fix ensemble all nve\n")
                elif ens == "nvt":
                    t_start = self.ensemble_temp_start[stage]
                    t_end = self.ensemble_temp_end[stage]
                    t_damp = self.ensemble_temp_damp[stage]
                    lmpin.write(
                        f"fix ensemble all nvt temp"
                        f" {t_start} {t_end}"
                        f" {t_damp}\n"
                    )
                else:
                    print(
                        f"Unknown ensemble type "
                        f"'{ens}' in stage {stage}."
                    )
                    sys.exit(1)

                # Run the simulation stage.
                rs = self.run_steps[stage]
                lmpin.write(f"run {rs}\n")

                # Turn off the fixes.
                if sf > 0:
                    lmpin.write("unfix squish\n")
                lmpin.write("unfix ensemble\n\n")

            # Final tail of the LAMMPS input file.
            lmpin.write("""
undump coord_dump
undump bond_dump

# Final density calculation
variable final_density equal density
print "${final_density}" append density_evolve

# Total energy calculation
###--------------------Define Settings------------------------------
compute eng all pe/atom
compute eatoms all reduce sum c_eng
#
###---------------------Run Minimization-----------------------------
reset_timestep 0
timestep 0.1
thermo 1
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms density
#dump 2 all custom 1 md.lammpstrj id type x y z
min_style cg
minimize 1e-18 1e-18 50000 100000
variable natoms equal "count(all)"
variable teng equal "c_eatoms"
print "${teng}" append totE_evolve
variable length equal "lx"
#variable ecoh equal "v_teng/v_natoms"
#
print "system_total_energy ${teng};"
print "Number of atoms is ${natoms};"
print "Lattice constant (Angstrom) is ${length};"
""")

        # ------------------------------------------------
        # Write the slurm submission file
        # ------------------------------------------------
        slurm_file = "slurm"
        with open(slurm_file, 'w') as slurm:
            slurm.write("""
#!/bin/bash
#SBATCH -p rulisp-lab,general
#SBATCH -A rulisp-lab
#SBATCH -J lmp
#SBATCH -o lmp.o%J
#SBATCH -e lmp.e%J
#SBATCH -N 1
#SBATCH -n 125
#SBATCH -t 00:15:00
#SBATCH --mem=2G
#
export OMP_NUM_THREADS=1
cd $SLURM_SUBMIT_DIR
mpirun lmp -in lammps.in
""")

        # Store data needed by normalize_types.
        self._bond_count = bond_count
        self._angle_count = angle_count
        self._num_unique_bond_tags = num_unique_bond_tags
        self._num_unique_angle_tags = num_unique_angle_tags
        self._unique_bond_tags_local = unique_bond_tags_local
        self._unique_angle_tags_local = unique_angle_tags_local
        self._num_ordered_species = num_ordered_species

        # Return to the original working directory.
        os.chdir(saved_cwd)

    # ============================================================
    # Step 7: Normalize types
    # ============================================================

    def normalize_types(self):
        """Unify atom, bond, and angle types across all files.

        We need to open each molecule file and also the LAMMPS data file and
        collect information about all the bond types that are present.  Those
        bond types will be analyzed and assembled into a list of the unique bond
        types.  Then each file will be opened again and all the bond types will
        be re-assigned according to the unique list of *all* bond types.

        The same process is applied to atom types and angle types to ensure
        complete consistency across the LAMMPS data file and all pre/post
        reaction template files.
        """

        ed = self.element_data
        bd = self.bond_data
        ad = self.angle_data

        # Initialize unique type trackers.  ``unique_bond_types`` and
        # ``unique_angle_types`` are doubly 1-indexed: the outer
        # dimension iterates over bond/angle type indices (slot 0 is
        # the None sentinel), and every stored entry is itself a
        # 1-indexed ``[None, bt1, bt2]`` / ``[None, at1, at2, at3,
        # at4]`` list mirroring Perl's
        # ``@uniqueBondTypes[$ub][1..2]`` /
        # ``@uniqueAngleTypes[$ua][1..4]`` layout.
        # ``unique_atom_types[ua]`` is a flat string token, so its
        # inner axis has no indexing to choose.
        num_unique_bonds = 0
        num_unique_angles = 0
        num_unique_atom_types = 0
        unique_bond_types = [None]   # 1-indexed
        unique_angle_types = [None]  # 1-indexed
        unique_atom_types = [None]   # 1-indexed

        # ------------------------------------------------
        # Read and parse the lammps.dat file into memory
        # ------------------------------------------------
        lammps_dat = os.path.join("lammps", "lammps.dat")
        with open(lammps_dat, 'r') as f:
            lammps_lines = f.read().splitlines()

        num_atom_types = 0
        num_bond_types_lmp = 0
        num_angle_types_lmp = 0

        line_num = 0
        while line_num < len(lammps_lines):
            line = lammps_lines[line_num]

            if "atom types" in line:
                vals = line.split()
                num_atom_types = int(vals[0])

            elif "bond types" in line:
                vals = line.split()
                num_bond_types_lmp = int(vals[0])

            elif "angle types" in line:
                vals = line.split()
                num_angle_types_lmp = int(vals[0])

            elif line.strip() == "Masses":
                line_num += 1  # skip blank line
                for _a in range(num_atom_types):
                    line_num += 1
                    vals = lammps_lines[line_num].split()
                    # Construct the atom type tag from the comment fields:
                    # element species molecule. Format: "idx mass # elem sp mol"
                    atom_type = f"{vals[3]} {vals[4]} {vals[5]}"

                    found = 0
                    for ua in range(1, num_unique_atom_types + 1):
                        if (unique_atom_types[ua]
                                == atom_type):
                            found = ua
                            break
                    if found == 0:
                        num_unique_atom_types += 1
                        unique_atom_types.append(atom_type)

            elif line.strip() == "Bond Coeffs":
                line_num += 1  # skip blank line
                for _b in range(num_bond_types_lmp):
                    line_num += 1
                    vals = lammps_lines[line_num].split()
                    # Format: "idx k r0 # e1 s1 m1 e2 s2 m2"
                    bt1 = f"{vals[4]} {vals[5]} {vals[6]}"
                    bt2 = f"{vals[7]} {vals[8]} {vals[9]}"

                    # Each unique_bond_types entry is [None, bt1, bt2]
                    # so slots 1..2 carry the two end-atom tags.
                    found = 0
                    for ub in range(1, num_unique_bonds + 1):
                        if (unique_bond_types[ub][1]
                                == bt1
                                and unique_bond_types[
                                    ub][2] == bt2):
                            found = ub
                            break
                    if found == 0:
                        num_unique_bonds += 1
                        unique_bond_types.append([None, bt1, bt2])

            elif line.strip() == "Angle Coeffs":
                line_num += 1  # skip blank line
                for _a in range(num_angle_types_lmp):
                    line_num += 1
                    vals = lammps_lines[line_num].split()
                    # Format: "idx k angle # e1 s1 m1 e2 s2 m2 e3 s3 m3 angle
                    #   hookeid"
                    at1 = f"{vals[4]} {vals[5]} {vals[6]}"
                    at2 = f"{vals[7]} {vals[8]} {vals[9]}"
                    at3 = f"{vals[10]} {vals[11]} {vals[12]}"
                    at4 = f"{vals[13]} {vals[14]}"

                    # Each unique_angle_types entry is 1-indexed
                    # [None, at1, at2, at3, at4] so slots 1..4 carry
                    # the two end-atom, vertex, and rest-angle tags.
                    found = 0
                    for ua in range(1, num_unique_angles + 1):
                        if (unique_angle_types[ua][1]
                                == at1
                                and unique_angle_types[
                                    ua][2] == at2
                                and unique_angle_types[
                                    ua][3] == at3
                                and unique_angle_types[
                                    ua][4] == at4):
                            found = ua
                            break
                    if found == 0:
                        num_unique_angles += 1
                        unique_angle_types.append([None, at1, at2, at3, at4])

            line_num += 1

        # ------------------------------------------------
        # Get the molecule file names
        # ------------------------------------------------
        lammps_dir = os.path.join("lammps")
        mol_files_pre = sorted(glob.glob(os.path.join(lammps_dir, "preRxn*")))
        mol_files_post = sorted(glob.glob(os.path.join(lammps_dir, "postRxn*")))
        mol_files = mol_files_pre + mol_files_post

        # ------------------------------------------------
        # Open each molecule file and collect all unique atom and bond and angle
        # types.
        # ------------------------------------------------
        for mol_path in mol_files:
            with open(mol_path, 'r') as mf:
                num_mol_atoms = 0
                num_mol_bonds = 0
                num_mol_angles = 0

                for mline in mf:
                    mline_s = mline.strip()

                    if "atoms" in mline_s:
                        vals = mline_s.split()
                        num_mol_atoms = int(vals[0])

                    elif "bonds" in mline_s:
                        vals = mline_s.split()
                        num_mol_bonds = int(vals[0])

                    elif "angles" in mline_s:
                        vals = mline_s.split()
                        num_mol_angles = int(vals[0])

                    elif "Types" in mline_s:
                        mf.readline()  # skip blank
                        for _a in range(num_mol_atoms):
                            vals = mf.readline().split()
                            # Adjust the molecule name to the family name.
                            mol_nm = vals[5]
                            if mol_nm in (self.molecule_to_family_map):
                                vals[5] = self.molecule_to_family_map[mol_nm]

                            atom_type = f"{vals[3]} {vals[4]} {vals[5]}"
                            found = 0
                            for ua in range(1, num_unique_atom_types + 1):
                                if (unique_atom_types[ua]
                                        == atom_type):
                                    found = ua
                                    break
                            if found == 0:
                                num_unique_atom_types += 1
                                unique_atom_types.append(atom_type)

                    elif "Bonds" in mline_s:
                        mf.readline()  # skip blank
                        for _b in range(num_mol_bonds):
                            vals = mf.readline().split()
                            # Adjust molecule names to family names.
                            for idx in[7, 10]:
                                mn = vals[idx]
                                if mn in (self.molecule_to_family_map):
                                    vals[idx] = self.molecule_to_family_map[mn]
                            bt1 = f"{vals[5]} {vals[6]} {vals[7]}"
                            bt2 = f"{vals[8]} {vals[9]} {vals[10]}"
                            # Slots 1..2 of each entry carry the two
                            # end-atom tags.
                            found = 0
                            for ub in range(1, num_unique_bonds + 1):
                                if (unique_bond_types
                                        [ub][1] == bt1
                                        and
                                        unique_bond_types
                                        [ub][2] == bt2):
                                    found = ub
                                    break
                            if found == 0:
                                num_unique_bonds += 1
                                unique_bond_types.append([None, bt1, bt2])

                    elif "Angles" in mline_s:
                        mf.readline()  # skip blank
                        for _a in range(num_mol_angles):
                            vals = mf.readline().split()
                            # Adjust molecule names.
                            for idx in[8, 11, 14]:
                                mn = vals[idx]
                                if mn in (self.molecule_to_family_map):
                                    vals[idx] = self.molecule_to_family_map[mn]
                            at1 = f"{vals[6]} {vals[7]} {vals[8]}"
                            at2 = f"{vals[9]} {vals[10]} {vals[11]}"
                            at3 = f"{vals[12]} {vals[13]} {vals[14]}"
                            at4 = f"{vals[15]} {vals[16]}"
                            # Slots 1..4 carry the two end-atom, vertex,
                            # and rest-angle tags (1-indexed layout).
                            found = 0
                            for ua in range(1, num_unique_angles + 1):
                                if (unique_angle_types
                                        [ua][1] == at1
                                        and
                                        unique_angle_types
                                        [ua][2] == at2
                                        and
                                        unique_angle_types
                                        [ua][3] == at3
                                        and
                                        unique_angle_types
                                        [ua][4] == at4):
                                    found = ua
                                    break
                            if found == 0:
                                num_unique_angles += 1
                                unique_angle_types.append(
                                    [None, at1, at2, at3, at4])

        # ------------------------------------------------
        # Compute coefficients for each unique bond type
        # ------------------------------------------------
        unique_bond_coeffs = [None] * (num_unique_bonds + 1)
        for ub in range(1, num_unique_bonds + 1):
            # unique_bond_types[ub] is 1-indexed [None, bt1, bt2].
            vals1 = unique_bond_types[ub][1].split()
            atom1_z = ed.get_element_z(vals1[0])
            vals2 = unique_bond_types[ub][2].split()
            atom2_z = ed.get_element_z(vals2[0])

            # Make sure atom1_z <= atom2_z.
            if atom1_z > atom2_z:
                atom1_z, atom2_z = atom2_z, atom1_z

            # hbc is 1-indexed: slots 1..4 are Z1, Z2, k, r0.
            # Store a 1-indexed [None, k, r0] pair.
            for hb in range(1, bd.num_hooke_bonds + 1):
                hbc = bd.hooke_bond_coeffs[hb]
                if (((atom1_z == hbc[1]
                      and atom2_z == hbc[2])
                     or (atom2_z == hbc[1]
                         and atom1_z == hbc[2]))):
                    unique_bond_coeffs[ub] = [None, hbc[3], hbc[4]]

        # ------------------------------------------------
        # Compute coefficients for each unique angle type
        # ------------------------------------------------
        unique_angle_coeffs = [None] * (num_unique_angles + 1)
        for ua in range(1, num_unique_angles + 1):
            # unique_angle_types[ua] is 1-indexed
            # [None, at1, at2, at3, at4].
            vals1 = unique_angle_types[ua][1].split()
            atom1_z = ed.get_element_z(vals1[0])

            vals2 = unique_angle_types[ua][2].split()
            atom2_z = ed.get_element_z(vals2[0])

            vals3 = unique_angle_types[ua][3].split()
            atom3_z = ed.get_element_z(vals3[0])

            # Make sure atom1_z <= atom3_z.
            if atom1_z > atom3_z:
                atom1_z, atom3_z = atom3_z, atom1_z

            vals4 = unique_angle_types[ua][4].split()
            curr_bond_angle = float(vals4[0])

            # hac is 1-indexed: slots 1..6 are Z1, Zv, Z2, k,
            # angle_deg, tolerance.  Store a 1-indexed
            # [None, k, angle_deg] pair — same shape as the
            # create_lammps_files producer.  The historical Perl
            # port carried an extra hooke-angle-id tail here, but
            # it was never read downstream; dropping it keeps both
            # producers structurally identical.
            for ha in range(1, ad.num_hooke_angles + 1):
                hac = ad.hooke_angle_coeffs[ha]
                if (atom1_z == hac[1]
                        and atom2_z == hac[2]
                        and atom3_z == hac[3]
                        and abs(curr_bond_angle - hac[5])
                        <= hac[6]):
                    unique_angle_coeffs[ua] = [None, hac[4], hac[5]]

        # ------------------------------------------------
        # Compute pair coefficients and masses for each unique atom type
        # ------------------------------------------------
        unique_lj_pair_coeffs = [None] * (num_unique_atom_types + 1)
        unique_masses = [None] * (num_unique_atom_types + 1)
        for at in range(1, num_unique_atom_types + 1):
            vals = unique_atom_types[at].split()
            atom1_z = ed.get_element_z(vals[0])
            lj = ed.lj_pair_coeffs[atom1_z]
            unique_lj_pair_coeffs[at] = f"{lj[0]} {lj[1]}"
            unique_masses[at] = ed.atomic_masses[atom1_z]

        # ------------------------------------------------
        # Rewrite the lammps.dat file with unified types
        # ------------------------------------------------
        with open(lammps_dat, 'w') as lmp:
            num_atoms_lmp = 0
            num_bonds_lmp = 0
            num_angles_lmp = 0

            line_num = 0
            while line_num < len(lammps_lines):
                line = lammps_lines[line_num]

                if "atom types" in line:
                    lmp.write(f"{num_unique_atom_types} atom types\n")

                elif re.match(r'^\d+\s+atoms\s*$', line.strip()):
                    vals = line.split()
                    num_atoms_lmp = int(vals[0])
                    lmp.write(f"{num_atoms_lmp} atoms\n")

                elif line.strip() == "Masses":
                    lmp.write("Masses\n\n")
                    for ua in range(1, num_unique_atom_types + 1):
                        lmp.write(
                            f"{ua} "
                            f"{unique_masses[ua]} "
                            f"# {unique_atom_types[ua]}"
                            f"\n"
                        )
                    # Skip original mass lines.
                    line_num += num_atom_types + 1

                elif line.strip() == "Pair Coeffs":
                    lmp.write("Pair Coeffs\n\n")
                    for ua in range(1, num_unique_atom_types + 1):
                        lmp.write(
                            f"{ua} "
                            f"{unique_lj_pair_coeffs[ua]}"
                            f" # "
                            f"{unique_atom_types[ua]}\n"
                        )
                    # Skip original pair coeff lines.
                    line_num += num_atom_types + 1

                elif "bond types" in line:
                    lmp.write(f"{num_unique_bonds} bond types\n")

                elif re.match(r'^\d+\s+bonds\s*$', line.strip()):
                    vals = line.split()
                    num_bonds_lmp = int(vals[0])
                    lmp.write(f"{num_bonds_lmp} bonds\n")

                elif line.strip() == "Bond Coeffs":
                    # bc is 1-indexed: slot 1 is k, slot 2 is r0.
                    # unique_bond_types[ub] is 1-indexed [None, bt1, bt2].
                    lmp.write("Bond Coeffs\n\n")
                    for ub in range(1, num_unique_bonds + 1):
                        bc = unique_bond_coeffs[ub]
                        lmp.write(
                            f"{ub} {bc[1]} {bc[2]} "
                            f"# "
                            f"{unique_bond_types[ub][1]}"
                            f" "
                            f"{unique_bond_types[ub][2]}"
                            f"\n"
                        )
                    # Skip original bond coeff lines.
                    line_num += num_bond_types_lmp + 1

                elif line.strip() == "Bonds":
                    lmp.write("Bonds\n\n")
                    line_num += 1  # skip blank line
                    for _b in range(num_bonds_lmp):
                        line_num += 1
                        vals = lammps_lines[line_num].split()
                        bt1 = f"{vals[5]} {vals[6]} {vals[7]}"
                        bt2 = f"{vals[8]} {vals[9]} {vals[10]}"

                        # Alphabetize the bond tag.
                        if bt1 > bt2:
                            bond_tag = f"{bt2} {bt1}"
                        else:
                            bond_tag = f"{bt1} {bt2}"

                        # Find the matching unique bond.  Each entry
                        # is 1-indexed [None, bt1, bt2] so the two
                        # end-atom tags live at slots 1 and 2.
                        for ub in range(1, num_unique_bonds + 1):
                            if (((bt1 ==
                                  unique_bond_types
                                  [ub][1]
                                  and bt2 ==
                                  unique_bond_types
                                  [ub][2])
                                 or (bt1 ==
                                     unique_bond_types
                                     [ub][2]
                                     and bt2 ==
                                     unique_bond_types
                                     [ub][1]))):
                                bond_num = int(vals[0])
                                lmp.write(
                                    f"{bond_num} {ub} "
                                    f"{vals[2]} "
                                    f"{vals[3]} "
                                    f"{vals[4]} "
                                    f"{bond_tag}\n"
                                )
                                break
                    line_num += 1  # skip trailing blank

                elif "angle types" in line:
                    lmp.write(f"{num_unique_angles} angle types\n")

                elif re.match(r'^\d+\s+angles\s*$', line.strip()):
                    vals = line.split()
                    num_angles_lmp = int(vals[0])
                    lmp.write(f"{num_angles_lmp} angles\n")

                elif line.strip() == "Angle Coeffs":
                    # ac is 1-indexed: slot 1 is k, slot 2 is the
                    # rest angle in degrees.  unique_angle_types[ua]
                    # is 1-indexed [None, at1, at2, at3, at4].
                    lmp.write("Angle Coeffs\n\n")
                    for ua in range(1, num_unique_angles + 1):
                        ac = unique_angle_coeffs[ua]
                        lmp.write(
                            f"{ua} {ac[1]} {ac[2]} "
                            f"# "
                            f"{unique_angle_types[ua][1]}"
                            f" "
                            f"{unique_angle_types[ua][2]}"
                            f" "
                            f"{unique_angle_types[ua][3]}"
                            f" "
                            f"{unique_angle_types[ua][4]}"
                            f"\n"
                        )
                    # Skip original angle coeff lines.
                    line_num += num_angle_types_lmp + 1

                elif line.strip() == "Angles":
                    lmp.write("Angles\n\n")
                    line_num += 1  # skip blank line
                    for _a in range(num_angles_lmp):
                        line_num += 1
                        vals = lammps_lines[line_num].split()
                        at1 = f"{vals[6]} {vals[7]} {vals[8]}"
                        at2 = f"{vals[9]} {vals[10]} {vals[11]}"
                        at3 = f"{vals[12]} {vals[13]} {vals[14]}"
                        at4 = f"{vals[15]} {vals[16]}"

                        # Alphabetize the angle tag.
                        if at1 > at3:
                            angle_tag = f"{at3} {at2} {at1} {at4}"
                        else:
                            angle_tag = f"{at1} {at2} {at3} {at4}"

                        # Find the matching unique angle.  uat is
                        # 1-indexed [None, at1, at2, at3, at4].
                        for ua in range(1, num_unique_angles + 1):
                            uat = unique_angle_types[ua]
                            if (((at1 == uat[1]
                                  and at2 == uat[2]
                                  and at3 == uat[3]
                                  and at4 == uat[4])
                                 or (at3 == uat[1]
                                     and at2 == uat[2]
                                     and at1 == uat[3]
                                     and at4 == uat[4]))):
                                angle_num = int(vals[0])
                                lmp.write(
                                    f"{angle_num} "
                                    f"{ua} "
                                    f"{vals[2]} "
                                    f"{vals[3]} "
                                    f"{vals[4]} # "
                                    f"{angle_tag}\n"
                                )
                                break
                    line_num += 1  # skip trailing blank

                else:
                    lmp.write(f"{line}\n")

                line_num += 1

        # ------------------------------------------------
        # Adjust the contents of each molecule file
        # ------------------------------------------------
        for mol_path in mol_files:
            with open(mol_path, 'r') as mf:
                mol_lines = mf.read().splitlines()

            with open(mol_path, 'w') as mf:
                num_mol_atoms = 0
                num_mol_bonds = 0
                num_mol_angles = 0
                line_num = 0

                while line_num < len(mol_lines):
                    line = mol_lines[line_num]

                    if "atoms" in line:
                        mf.write(f"{line}\n")
                        vals = line.split()
                        num_mol_atoms = int(vals[0])

                    elif "bonds" in line:
                        mf.write(f"{line}\n")
                        vals = line.split()
                        num_mol_bonds = int(vals[0])

                    elif "angles" in line:
                        mf.write(f"{line}\n")
                        vals = line.split()
                        num_mol_angles = int(vals[0])

                    elif "Types" in line:
                        mf.write(f"{line}\n\n")
                        line_num += 1  # skip blank
                        for _a in range(num_mol_atoms):
                            line_num += 1
                            vals = mol_lines[line_num].split()

                            # Adjust molecule name to family name.
                            mn = vals[5]
                            if mn in (self.molecule_to_family_map):
                                vals[5] = self.molecule_to_family_map[mn]

                            atom_type = f"{vals[3]} {vals[4]} {vals[5]}"

                            # Find the matching unique atom type.
                            for ua in range(1, num_unique_atom_types + 1):
                                if (unique_atom_types[ua]
                                        == atom_type):
                                    atom_num = int(vals[0])
                                    mf.write(
                                        f"{atom_num} "
                                        f"{ua} # "
                                        f"{atom_type}\n"
                                    )
                                    break

                    elif "Bonds" in line:
                        mf.write(f"{line}\n\n")
                        line_num += 1  # skip blank
                        for _b in range(num_mol_bonds):
                            line_num += 1
                            vals = mol_lines[line_num].split()

                            # Adjust molecule names.
                            for idx in[7, 10]:
                                mn = vals[idx]
                                if mn in (self.molecule_to_family_map):
                                    vals[idx] = self.molecule_to_family_map[mn]

                            bt1 = f"{vals[5]} {vals[6]} {vals[7]}"
                            bt2 = f"{vals[8]} {vals[9]} {vals[10]}"

                            # Alphabetize the bond tag.
                            if bt1 > bt2:
                                bond_tag = f"{bt2} {bt1}"
                            else:
                                bond_tag = f"{bt1} {bt2}"

                            # Find the matching bond.  Each entry is
                            # 1-indexed [None, bt1, bt2].
                            for ub in range(1, num_unique_bonds + 1):
                                if (((unique_bond_types
                                      [ub][1] == bt1
                                      and
                                      unique_bond_types
                                      [ub][2] == bt2)
                                     or
                                     (unique_bond_types
                                      [ub][2] == bt1
                                      and
                                      unique_bond_types
                                      [ub][1] == bt2))):
                                    bond_num = int(vals[0])
                                    mf.write(
                                        f"{bond_num} "
                                        f"{ub} "
                                        f"{vals[2]} "
                                        f"{vals[3]} "
                                        f"{vals[4]} "
                                        f"{bond_tag}\n"
                                    )
                                    break

                    elif "Angles" in line:
                        mf.write(f"{line}\n\n")
                        line_num += 1  # skip blank
                        for _a in range(num_mol_angles):
                            line_num += 1
                            vals = mol_lines[line_num].split()

                            # Adjust molecule names.
                            for idx in[8, 11, 14]:
                                mn = vals[idx]
                                if mn in (self.molecule_to_family_map):
                                    vals[idx] = self.molecule_to_family_map[mn]

                            at1 = f"{vals[6]} {vals[7]} {vals[8]}"
                            at2 = f"{vals[9]} {vals[10]} {vals[11]}"
                            at3 = f"{vals[12]} {vals[13]} {vals[14]}"
                            at4 = f"{vals[15]} {vals[16]}"

                            # Alphabetize.
                            if at1 > at3:
                                angle_tag = f"{at3} {at2} {at1} {at4}"
                            else:
                                angle_tag = f"{at1} {at2} {at3} {at4}"

                            # Find matching angle.  uat is 1-indexed
                            # [None, at1, at2, at3, at4].
                            for ua in range(1, num_unique_angles + 1):
                                uat = unique_angle_types[ua]
                                if (((at1 == uat[1]
                                      and at2 == uat[2]
                                      and at3 == uat[3]
                                      and at4 == uat[4])
                                     or
                                     (at3 == uat[1]
                                      and at2 == uat[2]
                                      and at1 == uat[3]
                                      and at4 == uat[4]))):
                                    angle_num = int(vals[0])
                                    mf.write(
                                        f"{angle_num} "
                                        f"{ua} "
                                        f"{vals[2]} "
                                        f"{vals[3]} "
                                        f"{vals[4]} # "
                                        f"{angle_tag}"
                                        f"\n"
                                    )
                                    break

                    else:
                        mf.write(f"{line}\n")

                    line_num += 1


# ================================================================
# Main entry point
# ================================================================

def main():
    """Main function: orchestrate the condense workflow.

    The overall flow is:
    1. Parse settings from the rc file and command line.
    2. Initialize the element, bond, and angle databases.
    3. Read the condense input file.
    4. Compute implicit information (atom counts, molecule
       assignments).
    5. Copy reaction templates from the precursor database.
    6. Run packmol to place molecules in the simulation cell.
    7. Create LAMMPS data, input, and slurm files.
    8. Normalize atom/bond/angle types across all files.
    """

    # Get script settings from a combination of the resource control file and
    # parameters given by the user on the command line.
    settings = ScriptSettings()

    # Create the Condense object and execute each step.
    condense = Condense(settings)

    # Initialize the environment (element/bond/angle DBs).
    condense.init_env()

    # Read the main input file.
    condense.parse_input_file()

    # Compute implicit information not explicitly given in the input file.
    condense.compute_implicit_input()

    # Obtain the reaction templates from the database.
    condense.copy_reaction_templates()

    # Assemble the packmol input file and run packmol.
    condense.run_packmol()

    # Create LAMMPS data_file, LAMMPS in_file, and the slurm submission file.
    condense.create_lammps_files()

    # Make all types and bond types consistent across all files.
    condense.normalize_types()


if __name__ == '__main__':
    # Everything before this point was a class/function definition or a request
    # to import information from external modules.  Only now do we actually
    # start running the program.  The purpose of this is to allow another python
    # program to import *this* script and call its functions internally.
    main()
