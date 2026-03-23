#!/usr/bin/env python3

"""runIsoAtoms.py -- Run isolated atom calculations and build
a potential function database.

PROGRAM:  runIsoAtoms
PURPOSE:  This program will run an isolated atom calculation
          for every atom in the element database.  The
          resultant converged potentials will be used to
          create an atomic potential function database.

          For each element (or a single user-specified
          element), the script:
            1. Creates a directory named after the element.
            2. Writes a minimal olcao.skl file with a single
               atom in a 20x20x20 Angstrom cell.
            3. Runs makeinput to generate OLCAO input files.
            4. Runs uolcao to perform the SCF calculation.
          After all elements are computed, the script:
            5. Creates an "atomPDB" directory.
            6. For each element, copies the converged SCF
               potential coefficients (gs_scfV-fb.dat) into
               "coeff1" and "coeff.isolated" files.
            7. Copies the existing potential definition file
               (pot1) from $OLCAO_DATA/atomicPDB/.

          The -nocore option can be used to request an
          all-electron calculation.  It is somewhat counter-
          intuitive that an option called -nocore should
          yield an all-electron result.  But, the idea is
          that -nocore *will not* include any core orbitals
          in the list of orbitals that are orthogonalized
          away.  Therefore, all of the core orbitals *will
          be* included in the list of orbitals to compute.

USAGE:    runIsoAtoms.py [-e ELEMENT_Z] [--nocore] [--help]

OPTIONS:

  -e/--element  Atomic Z number of a single element to
                process.  When not given (or 0), all
                elements in the database are processed.

  --nocore      Pass the -nocore option to makeinput so
                that the all-electron electronic structure
                is computed.

  --help        Print usage information and exit.
"""

import argparse as ap
import os
import shutil
import subprocess
import sys
from datetime import datetime

from element_data import ElementData


# ----------------------------------------------------------------
#  ScriptSettings
# ----------------------------------------------------------------

class ScriptSettings():
    """User settings that control the program.

    Variable values are pulled from a resource control file
    (runIsoAtomsrc.py) and reconciled with command line
    parameters.  The resource control file is loaded from
    $OLCAO_RC (or from the current working directory if a
    local copy is present).
    """

    def __init__(self):
        """Initialize settings from the rc file and command
        line."""

        # Read default variables from the resource control
        #   file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit(
                "Error: $OLCAO_RC is not set. "
                "See instructions."
            )
        sys.path.insert(1, rc_dir)
        from runIsoAtomsrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc
        #   defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc
        #   file.
        self.reconcile(args)

        # Record the command line parameters that were used.
        self.recordCLP()

    def assign_rc_defaults(self, default_rc):
        """Pull each default value from the rc dictionary
        and assign it to an instance variable."""

        # Atomic Z number of a single element to process.
        #   0 means do all elements.
        self.target_z = default_rc["target_z"]

        # Whether to use the -nocore option when running
        #   makeinput.
        self.nocore = default_rc["nocore"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv."""

        prog_name = "runIsoAtoms.py"

        description_text = """\
PROGRAM:  runIsoAtoms
PURPOSE:  Run an isolated atom calculation for every element
          (or a single specified element) and build an atomic
          potential function database from the converged SCF
          potentials.
"""

        epilog_text = """\
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults: ./runIsoAtomsrc.py or $OLCAO_RC/runIsoAtomsrc.py.
"""

        parser = ap.ArgumentParser(
            prog=prog_name,
            formatter_class=ap.RawDescriptionHelpFormatter,
            description=description_text,
            epilog=epilog_text,
        )

        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()

    def add_parser_arguments(self, parser):
        """Define all command line arguments."""

        # The -e option takes one atomic Z number as an
        #   argument and will run an isolated atom
        #   calculation for that element only.  When not
        #   given (or 0), all elements are processed.
        parser.add_argument(
            '-e', '--element',
            dest='target_z',
            type=int,
            default=self.target_z,
            help=(
                'Atomic Z number of a single element to '
                'process. 0 = all elements. '
                f'Default: {self.target_z}'
            ),
        )

        # The --nocore option will pass on the -nocore
        #   option to makeinput so that the all-electron
        #   electronic structure is computed.
        parser.add_argument(
            '--nocore',
            dest='nocore',
            action='store_true',
            default=self.nocore,
            help=(
                'Pass -nocore to makeinput for '
                'all-electron calculation. '
                f'Default: {self.nocore}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command line values with rc defaults."""
        self.target_z = args.target_z
        self.nocore = args.nocore

    def recordCLP(self):
        """Append the command line invocation to the
        'command' file for provenance tracking."""
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


# ----------------------------------------------------------------
#  Helper: determine the element range
# ----------------------------------------------------------------

def get_element_range(target_z, num_elements):
    """Determine the range of elements to process.

    If target_z is 0, all elements (1..num_elements) are
    processed.  Otherwise, only the single element with
    Z = target_z is processed.

    Args:
        target_z (int): Atomic Z number (0 = all).
        num_elements (int): Total number of elements in
            the database.

    Returns:
        tuple of int: (element_init, element_fin) defining
            the inclusive range of elements to process.
    """
    if target_z == 0:
        return 1, num_elements
    else:
        return target_z, target_z


# ----------------------------------------------------------------
#  Isolated atom SCF calculations
# ----------------------------------------------------------------

def make_atom_pot(settings, element_names, num_elements):
    """Create and run isolated atom SCF calculations.

    For each element in the specified range, this function:
      1. Creates a directory named after the element.
      2. Writes a minimal olcao.skl file containing a single
         atom of that element placed at the origin of a
         20x20x20 Angstrom cubic cell with 90-degree angles.
      3. Runs makeinput to generate the OLCAO input files.
         If --nocore was requested, the -nocore flag is
         passed to makeinput.
      4. Runs uolcao (the unified OLCAO SCF solver) to
         perform the self-consistent field calculation.

    The resulting converged potential is stored in the
    element's directory as gs_scfV-fb.dat and will be
    harvested later by make_pdb().

    Args:
        settings (ScriptSettings): Program settings.
        element_names (list of str): 1-indexed list of
            lowercase element names from the database.
        num_elements (int): Total number of elements.
    """
    olcao_bin = os.getenv('OLCAO_BIN', '')

    element_init, element_fin = get_element_range(
        settings.target_z, num_elements
    )

    # Save the starting directory so we can return to it
    #   after processing each element.
    start_dir = os.getcwd()

    for element in range(element_init, element_fin + 1):
        elem_name = element_names[element]

        # Create and enter the element directory.
        os.makedirs(elem_name, exist_ok=True)
        os.chdir(elem_name)

        # Write the minimal skeleton file for an isolated
        #   atom.  The 20x20x20 cell ensures the atom is
        #   well-separated from its periodic images.
        with open("olcao.skl", "w") as skl:
            skl.write("title\n")
            skl.write(f"Isolated {elem_name}\n")
            skl.write("end\n")
            skl.write("cell\n")
            skl.write("   20 20 20 90 90 90\n")
            skl.write("fract 1\n")
            skl.write(
                f"{elem_name} 0.000 0.000 0.000\n"
            )
            skl.write("space 1_a\n")
            skl.write("supercell 1 1 1\n")
            skl.write("prim\n")

        # Run makeinput to generate OLCAO input files.
        #   If --nocore was requested, pass -nocore to
        #   makeinput.
        makeinput_cmd = "makeinput.py"
        if settings.nocore:
            makeinput_cmd += " -nocore"
        subprocess.run(
            makeinput_cmd, shell=True, check=True
        )

        # Run uolcao to perform the SCF calculation.
        uolcao_path = os.path.join(olcao_bin, "uolcao")
        subprocess.run(
            uolcao_path, shell=True, check=True
        )

        # Return to the starting directory.
        os.chdir(start_dir)


# ----------------------------------------------------------------
#  Potential database builder
# ----------------------------------------------------------------

def make_pdb(settings, element_names, num_elements):
    """Create a new atomic potential function database.

    For each element in the specified range, this function:
      1. Creates a sub-directory under "atomPDB/" named
         after the element.
      2. Opens the just-computed SCF converged potential
         file (gs_scfV-fb.dat) from the element's
         calculation directory.
      3. Skips the header line of the SCF file.
      4. Copies the remaining potential coefficients into
         two files:
           - coeff1: the standard coefficient file used by
             OLCAO for subsequent calculations.
           - coeff.isolated: a duplicate that explicitly
             marks these coefficients as coming from an
             isolated atom calculation.  (Even though
             presently both files contain the same
             information, the distinction is maintained for
             clarity and potential future use.)
      5. Copies the existing potential definition file
         (pot1) from $OLCAO_DATA/atomicPDB/<element>/pot1
         into the new database directory.

    Args:
        settings (ScriptSettings): Program settings.
        element_names (list of str): 1-indexed list of
            lowercase element names from the database.
        num_elements (int): Total number of elements.
    """
    olcao_data = os.getenv('OLCAO_DATA', '')

    element_init, element_fin = get_element_range(
        settings.target_z, num_elements
    )

    # Create the directory that will contain the potential
    #   function database.
    start_dir = os.getcwd()
    os.makedirs("atomPDB", exist_ok=True)
    os.chdir("atomPDB")

    for element in range(element_init, element_fin + 1):
        elem_name = element_names[element]

        # Create and enter the directory for this element.
        os.makedirs(elem_name, exist_ok=True)
        os.chdir(elem_name)

        # Path to the just-computed SCF converged potential
        #   coefficients file.
        scf_path = os.path.join(
            start_dir, elem_name, "gs_scfV-fb.dat"
        )

        # Open the SCF coefficient file and read past its
        #   header line.
        with open(scf_path, 'r') as scf_in:
            scf_in.readline()  # Skip the header.

            # Copy the rest of the SCF data into the new
            #   potential coefficient files.
            remaining = scf_in.read()

        # Write the standard coefficient file.
        with open("coeff1", 'w') as coeff_out:
            coeff_out.write(remaining)

        # Write the isolated-atom coefficient file.
        #   This is a duplicate that explicitly marks
        #   these coefficients as from an isolated atom
        #   calculation.
        with open("coeff.isolated", 'w') as iso_out:
            iso_out.write(remaining)

        # Copy the existing potential definition file
        #   from the OLCAO data directory.
        pot_src = os.path.join(
            olcao_data, "atomicPDB",
            elem_name, "pot1"
        )
        shutil.copy2(pot_src, "pot1")

        # Return to the atomPDB directory.
        os.chdir("..")

    # Return to the original starting directory.
    os.chdir(start_dir)


# ----------------------------------------------------------------
#  Main
# ----------------------------------------------------------------

def main():
    """Main entry point for runIsoAtoms.

    Runs isolated atom SCF calculations for every element
    (or a single specified element) and then builds an
    atomic potential function database from the converged
    potentials.

    The workflow mirrors the original Perl script:
      1. Parse command line and set defaults.
      2. Initialize the element database.
      3. Run isolated atom SCF calculations.
      4. Build the potential function database.
    """

    # Get script settings from a combination of the
    #   resource control file and parameters given by the
    #   user on the command line.
    settings = ScriptSettings()

    # Read the element data into the Python module.
    #   This provides the element names, number of
    #   elements, and other periodic table data.
    ed = ElementData()
    ed.init_element_data()

    # Get the number of possible elements to compute and
    #   the lowercase names of the elements from the
    #   periodic table.
    num_elements = ed.get_num_elements()
    element_names = ed.element_names

    # Create input for every atom and run the SCF
    #   calculation.
    make_atom_pot(settings, element_names, num_elements)

    # Create a new atomic potential function database
    #   from the converged SCF potentials.
    make_pdb(settings, element_names, num_elements)


if __name__ == '__main__':
    # Everything before this point was a subroutine
    #   definition or a request to import information from
    #   external modules.  Only now do we actually start
    #   running the program.  The purpose of this is to
    #   allow another python program to import *this* script
    #   and call its functions internally.
    main()
