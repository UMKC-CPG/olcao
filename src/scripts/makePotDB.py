#!/usr/bin/env python3

# PROGRAM:  makePotDB
# PURPOSE:  This program will create a default potential database
#           using information about each element stored in the
#           element_data.py module. This program will also create
#           an initial set of potential coefficients based on the
#           numerical potential created when the atomic wave
#           function basis set was created (makeBDB). The program
#           that actually created the numerical potential is the
#           atomSCF program.
#
# USAGE:    makePotDB.py [-e ELEMENT] [-f FORK] [-h]
#
# OPTIONS:
#
# The -e/--element option takes one atomic Z number as an argument
#   and will create the atomic potential files for that element
#   only.  When set to 0 (the default), all elements are
#   processed.
# The -f/--fork option allows the user to request that the
#   Gaussian fitting of the numerical potential be performed in
#   parallel with the given number of processes. (The fitting can
#   take a little while, but not forever.)
# The -h/--help option prints usage information.

import argparse as ap
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
from datetime import datetime

from element_data import ElementData


# ----------------------------------------------------------------
#  Element data singleton — initialized once, used throughout.
# ----------------------------------------------------------------
_ed = ElementData()
_ed.init_element_data()


# ----------------------------------------------------------------
#  ScriptSettings
# ----------------------------------------------------------------

class ScriptSettings():
    """The instance variables of this object are the user settings
       that control the program. The variable values are pulled
       from a list that is created within a resource control file
       and that are then reconciled with command line
       parameters."""


    def __init__(self):
        """Define default values for the parameters by pulling
        them from the resource control file in the default
        location: $OLCAO_RC/makePotDBrc.py or from the current
        working directory if a local copy of makePotDBrc.py is
        present."""

        # Read default variables from the resource control file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit(
                "Error: $OLCAO_RC is not set. See instructions."
            )
        sys.path.insert(1, rc_dir)
        from makePotDBrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc file.
        self.reconcile(args)

        # Record the command line parameters that were used.
        self.recordCLP()


    def assign_rc_defaults(self, default_rc):
        """Pull each default value from the rc dictionary and
        assign it to an instance variable."""

        # Target element Z number (0 = all elements).
        self.target_z = default_rc["element"]

        # Number of parallel processes for Gaussian fitting.
        self.num_forks = default_rc["fork"]


    def parse_command_line(self):
        """Build the argument parser and parse sys.argv."""

        prog_name = "makePotDB.py"

        description_text = """\
PROGRAM:  makePotDB
PURPOSE:  Create a default potential database using information
          about each element stored in element_data.py. Also
          create an initial set of potential coefficients based
          on the numerical potential created when the atomic
          wave function basis set was created (makeBDB). The
          program that actually created the numerical potential
          is the atomSCF program.
"""

        epilog_text = """\
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults are given in ./makePotDBrc.py or \
$OLCAO_RC/makePotDBrc.py.
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

        # The -element option takes one atomic Z number as an
        #   argument and will create the atomic potential files
        #   for that element only.
        parser.add_argument(
            '-e', '--element',
            dest='target_z',
            type=int,
            default=self.target_z,
            help=(
                'Atomic Z number of the single element to '
                'process. 0 means all elements. '
                f'Default: {self.target_z}'
            ),
        )

        # The -fork option allows the user to request that the
        #   Gaussian fitting of the numerical potential be
        #   performed in parallel with the given number of
        #   processes.
        parser.add_argument(
            '-f', '--fork',
            dest='num_forks',
            type=int,
            default=self.num_forks,
            help=(
                'Number of parallel processes for Gaussian '
                'fitting. '
                f'Default: {self.num_forks}'
            ),
        )


    def reconcile(self, args):
        """Reconcile command line values with rc defaults."""
        self.target_z = args.target_z
        self.num_forks = args.num_forks


    def recordCLP(self):
        """Append the command line invocation to the 'command'
        file for provenance tracking."""
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
#  Helper functions
# ----------------------------------------------------------------

def get_element_range(target_z, num_elements):
    """Return the (inclusive) range of element indices to process.

    If target_z is 0, process all elements (1 .. num_elements).
    Otherwise, process only the single requested element.

    Args:
        target_z (int): Target element Z number (0 = all).
        num_elements (int): Total number of elements in the
            database.

    Returns:
        range: A Python range object over the element indices.
    """
    if target_z == 0:
        return range(1, num_elements + 1)
    else:
        return range(target_z, target_z + 1)


def get_element_data(element):
    """Extract potential definition data for a given element.

    Retrieves the number of Gaussian terms and the minimum and
    maximum exponential alpha values that define the potential
    expansion for the specified element.

    Args:
        element (int): Atomic Z number of the element.

    Returns:
        tuple: (num_pot_alphas, min_alpha, max_alpha)
            num_pot_alphas (int): Number of Gaussian terms for
                the potential function.
            min_alpha (float): Minimum exponential Gaussian
                coefficient for the potential.
            max_alpha (float): Maximum exponential Gaussian
                coefficient for the potential.
    """
    num_pot_alphas = _ed.num_terms_pot[element]
    min_alpha = _ed.min_term_pot[element]
    max_alpha = _ed.max_term_pot[element]
    return (num_pot_alphas, min_alpha, max_alpha)


def get_database_paths():
    """Determine the atomic database paths from $OLCAO_DATA.

    The atomicPDB directory holds the potential function database.
    The atomicBDB directory holds the basis set database (which
    contains the numerical potentials created by atomSCF).

    Returns:
        tuple: (atomic_pdb, atomic_bdb)
            atomic_pdb (str): Path to the atomic potential
                database.
            atomic_bdb (str): Path to the atomic basis set
                database.
    """
    data_dir = os.getenv('OLCAO_DATA')
    if not data_dir:
        sys.exit(
            "Error: $OLCAO_DATA is not set. See instructions."
        )
    atomic_pdb = os.path.join(data_dir, 'atomicPDB')
    atomic_bdb = os.path.join(data_dir, 'atomicBDB')
    return (atomic_pdb, atomic_bdb)


# ----------------------------------------------------------------
#  Core subroutines
# ----------------------------------------------------------------

def make_dirs(target_z, num_elements, atomic_pdb):
    """Create the element directories in the potential database.

    For each element in the target range, create a subdirectory
    (named by element symbol) inside the atomicPDB directory.
    The database directory itself is also created if it does not
    already exist.

    Args:
        target_z (int): Target element Z number (0 = all).
        num_elements (int): Total number of elements.
        atomic_pdb (str): Path to the atomicPDB directory.
    """
    elem_range = get_element_range(target_z, num_elements)

    # Create the database directory if it does not already exist.
    os.makedirs(atomic_pdb, exist_ok=True)

    # Make the directory for each element only if it does not
    #   already exist.
    for element in elem_range:
        elem_dir = os.path.join(
            atomic_pdb, _ed.element_names[element]
        )
        os.makedirs(elem_dir, exist_ok=True)


def make_pot_definition(target_z, num_elements, atomic_pdb):
    """Write the potential definition file (pot1) for each element.

    The pot1 file specifies the nuclear charge, a default alpha
    value, a default covalent radius, and the number and range of
    Gaussian exponents used to expand the potential.

    Each pot1 file has the format:
        NUCLEAR_CHARGE__ALPHA
        <Z> 20.000000
        COVALENT_RADIUS
        1.000000
        NUM_ALPHAS
        <num_pot_alphas>
        ALPHAS
        <min_alpha> <max_alpha>

    Args:
        target_z (int): Target element Z number (0 = all).
        num_elements (int): Total number of elements.
        atomic_pdb (str): Path to the atomicPDB directory.
    """
    elem_range = get_element_range(target_z, num_elements)

    # Print the potential definition file for each requested
    #   element.
    for element in elem_range:
        elem_name = _ed.element_names[element]
        pot_file = os.path.join(
            atomic_pdb, elem_name, "pot1"
        )

        # Do not touch this element if it already is done.
        if os.path.exists(pot_file):
            continue

        # Extract information needed for this element.
        num_pot_alphas, min_alpha, max_alpha = \
            get_element_data(element)

        # Write the potential definition.
        with open(pot_file, 'w') as pot:
            pot.write("NUCLEAR_CHARGE__ALPHA\n")
            pot.write(f"{element:f} {20.0:f}\n")
            pot.write("COVALENT_RADIUS\n")
            pot.write(f"{1.0:f}\n")
            pot.write("NUM_ALPHAS\n")
            pot.write(f"{num_pot_alphas:d}\n")
            pot.write("ALPHAS\n")
            pot.write(f"{min_alpha:e} {max_alpha:e}\n")


def _fit_one_element(element, atomic_pdb, atomic_bdb):
    """Fit the potential for a single element (worker function).

    This function is designed to be called in a child process
    (via multiprocessing). It:
      1. Reads the numerical potential from the atomicBDB
         directory (created by atomSCF during makeBDB).
      2. Modifies the header parameters to reflect the potential
         definition for this atom.
      3. Runs the gaussFit program to perform a linear Gaussian
         fit.
      4. Appends three zero columns (total charge, valence
         charge, spin polarization) to each coefficient row.
      5. Moves the result to coeff1 and cleans up.

    Args:
        element (int): Atomic Z number.
        atomic_pdb (str): Path to the atomicPDB directory.
        atomic_bdb (str): Path to the atomicBDB directory.
    """
    elem_name = _ed.element_names[element]
    pot_data = "numericalPot.dat"
    elem_dir = os.path.join(atomic_pdb, elem_name)

    # Do not touch this element if it already is done.
    coeff_file = os.path.join(elem_dir, "coeff1")
    if os.path.exists(coeff_file):
        return

    # Extract information needed for this element.
    num_pot_alphas, min_alpha, max_alpha = \
        get_element_data(element)

    # Modify the parameters in the numerical potential file to
    #   reflect the definition of the potential for this atom.

    # Open the numerical potential file that was fitted
    #   non-linearly for the contracting process. (READING)
    #   Note that this potential file is from the atomic
    #   contraction.
    src_pot = os.path.join(
        atomic_bdb, elem_name, pot_data
    )
    dst_pot = os.path.join(elem_dir, pot_data)

    with open(src_pot, 'r') as numpot, \
         open(dst_pot, 'w') as fitpot:

        # Read and modify the header information.
        # Line 1: copy as-is.
        fitpot.write(numpot.readline())
        # Line 2: replace with 0.
        numpot.readline()
        fitpot.write("0\n")
        # Line 3: copy as-is.
        fitpot.write(numpot.readline())
        # Line 4: replace with num_pot_alphas.
        numpot.readline()
        fitpot.write(f"{num_pot_alphas}\n")
        # Line 5: copy as-is.
        fitpot.write(numpot.readline())
        # Line 6: replace with min_alpha.
        numpot.readline()
        fitpot.write(f"{min_alpha}\n")
        # Line 7: copy as-is.
        fitpot.write(numpot.readline())
        # Line 8: replace with max_alpha.
        numpot.readline()
        fitpot.write(f"{max_alpha}\n")
        # Remaining lines: copy as-is.
        for line in numpot:
            fitpot.write(line)

    # Execute the fitting program.
    with open(
        os.path.join(elem_dir, "progress.fit"), 'w'
    ) as prog_out:
        subprocess.run(
            ["gaussFit"],
            stdin=open(dst_pot, 'r'),
            stdout=prog_out,
            cwd=elem_dir,
            check=True,
        )

    # Append some zeros to each row representing the known
    #   total and valence charge, and the known spin
    #   polarization. These values are read by OLCAO, but not
    #   used.

    # Open the file, read the data to memory and strip
    #   trailing whitespace. Then close the file.
    gauss_fit = os.path.join(elem_dir, "gauss.fit")
    with open(gauss_fit, 'r') as f:
        coeff_data = [line.rstrip() for line in f]

    # Re-open the file for writing now. Print each line with
    #   a set of 3 zero values appended (except the first).
    #   Then close the file.
    with open(gauss_fit, 'w') as f:
        f.write(f"{coeff_data[0]}\n")
        for line in coeff_data[1:]:
            f.write(f"{line} 0.0 0.0 0.0\n")

    # Move the result to the correct file name and clean up.
    shutil.move(gauss_fit, coeff_file)
    for cleanup in [
        "gaussFit.out", "progress.fit", pot_data
    ]:
        path = os.path.join(elem_dir, cleanup)
        if os.path.exists(path):
            os.remove(path)


def make_pot_coeffs(
    target_z, num_elements, num_forks, atomic_pdb, atomic_bdb
):
    """Create the initial potential coefficients for each element.

    For each element in the target range, this function:
      1. Copies the numerical potential from atomicBDB with
         modified header parameters.
      2. Runs gaussFit to perform a linear Gaussian fit.
      3. Appends zero columns for total charge, valence charge,
         and spin polarization.
      4. Saves the result as coeff1 in the element directory.

    If num_forks > 1, the fitting is performed in parallel using
    a multiprocessing pool.

    Args:
        target_z (int): Target element Z number (0 = all).
        num_elements (int): Total number of elements.
        num_forks (int): Number of parallel worker processes.
        atomic_pdb (str): Path to the atomicPDB directory.
        atomic_bdb (str): Path to the atomicBDB directory.
    """
    elem_range = get_element_range(target_z, num_elements)

    if num_forks <= 1:
        # Serial execution.
        for element in elem_range:
            _fit_one_element(element, atomic_pdb, atomic_bdb)
    else:
        # Parallel execution using a process pool.
        with mp.Pool(processes=num_forks) as pool:
            pool.starmap(
                _fit_one_element,
                [
                    (elem, atomic_pdb, atomic_bdb)
                    for elem in elem_range
                ],
            )


# ----------------------------------------------------------------
#  Main
# ----------------------------------------------------------------

def main():

    # Get script settings from a combination of the resource
    #   control file and parameters given by the user on the
    #   command line.
    settings = ScriptSettings()

    # Get the number of elements from the element data module.
    num_elements = len(_ed.element_names) - 1  # 1-indexed

    # Get the database paths from the environment.
    atomic_pdb, atomic_bdb = get_database_paths()

    # For each element, create the directory.
    make_dirs(settings.target_z, num_elements, atomic_pdb)

    # Enter each directory and create the pot files for the
    #   potential definition.
    make_pot_definition(
        settings.target_z, num_elements, atomic_pdb
    )

    # Enter each directory and create an initial set of
    #   coefficients.
    make_pot_coeffs(
        settings.target_z,
        num_elements,
        settings.num_forks,
        atomic_pdb,
        atomic_bdb,
    )


if __name__ == '__main__':
    # Everything before this point was a subroutine definition
    #   or a request to import information from external modules.
    #   Only now do we actually start running the program. The
    #   purpose of this is to allow another python program to
    #   import *this* script and call its functions internally.
    main()
