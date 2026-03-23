#!/usr/bin/env python3

"""makeBDB.py -- Build the atomic basis set database for OLCAO.

PROGRAM:  makeBDB
PURPOSE:  This program will perform all the necessary operations
          to construct a default atomic basis set for OLCAO
          calculations using the atomic data stored in the
          element_data.py module.  There are three steps to the
          basis set construction.

          (1) An atomic SCF input file is created and a
          calculation is done using it to obtain a numerical
          free atom potential function.

          (2) The numerical potential is fit with Gaussians as
          accurately as possible.  (This may take some
          considerable time.)

          (3) A pair of input files for contracting a Gaussian
          basis set are constructed.  The first contains
          specific orbitals sequestered to the core where they
          can be orthogonalized out of the secular (eigenvalue)
          equation, the second includes those core orbitals in
          the valence section so that they are included in the
          secular equation.

          NOTE:  The contraction is not actually performed at
          this stage.  The wave functions are contracted at the
          time the 'makeinput' script is invoked.  This allows
          the user to modify the basis set used for any purpose
          if they desire.  It is recommended however, that
          because the basis sets are generally applicable to an
          exceptionally wide variety of systems that they not
          be modified without cause to 'improve' a result for
          a specific system.

          Advice and criticism on the basis set construction
          are always welcome.

USAGE:    makeBDB.py [--step {all,scf,fit,cont}]
                     [--fork NUMFORKS] [--element Z]
                     [--force] [--verbose] [--help]

OPTIONS:

  --step    Which step(s) to perform:
              "all"  = do all steps (default).
              "scf"  = do only the atomic SCF step.
              "fit"  = do only the Gaussian fitting step.
              "cont" = do only the contraction file step.
  --fork    Number of parallel processes to use for the
            Gaussian fitting of the numerical potential.
            (The non-linear fitting can take a very long time.)
  --element Atomic Z number of a single element to process.
            When set to 0 (the default), all elements are
            treated.
  --force   Overwrite any previously obtained results for the
            requested operations.
  --verbose Print progress reports of what the program is doing.
  --help    Print usage information and exit.
"""

import argparse as ap
import math
import multiprocessing as mp
import os
import subprocess
import sys
from datetime import datetime

from element_data import ElementData


# ----------------------------------------------------------------
#  Element data singleton -- initialized once, used throughout.
# ----------------------------------------------------------------
_ed = ElementData()
_ed.init_element_data()


# Tags for each orbital angular momentum quantum number.
#   s=0 -> 'S', p=1 -> 'P', d=2 -> 'D', f=3 -> 'F'
ORBITAL_TAGS = ['S', 'P', 'D', 'F']


# ----------------------------------------------------------------
#  ScriptSettings
# ----------------------------------------------------------------

class ScriptSettings():
    """User settings that control the program.

    Variable values are pulled from a resource control file
    (makeBDBrc.py) and reconciled with command line parameters.
    The resource control file is loaded from $OLCAO_RC (or from
    the current working directory if a local copy is present).
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
        from makeBDBrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc
        #   file.
        self.reconcile(args)

        # At this point, the command line parameters are set
        #   and accepted.  When this initialization subroutine
        #   returns the script will start running.  So, we use
        #   this as a good spot to record the command line
        #   parameters that were used.
        self.recordCLP()

    def assign_rc_defaults(self, default_rc):
        """Pull each default value from the rc dictionary and
        assign it to an instance variable."""

        # Which step(s) to perform: all, scf, fit, or cont.
        self.step = default_rc["step"]

        # Target element Z number (0 = all elements).
        self.target_z = default_rc["element"]

        # Number of parallel processes for Gaussian fitting.
        self.num_forks = default_rc["fork"]

        # Flag to overwrite previously obtained results.
        self.force = default_rc["force"]

        # Flag to print progress reports.
        self.verbose = default_rc["verbose"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv."""

        prog_name = "makeBDB.py"

        description_text = """\
PROGRAM:  makeBDB
PURPOSE:  Construct a default atomic basis set for OLCAO
          calculations using the atomic data stored in the
          element_data.py module.  Three steps:
          (1) Create an atomic SCF input file and run it
              to obtain a numerical free-atom potential.
          (2) Fit the numerical potential with Gaussians.
          (3) Create contraction input files (with-core
              and no-core variants).
          The contraction itself is deferred to makeinput.
"""

        epilog_text = """\
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults: ./makeBDBrc.py or $OLCAO_RC/makeBDBrc.py.
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

        # The --step option selects which step(s) to perform.
        #   "all"  = do all three steps (default).
        #   "scf"  = do only the atomic SCF step.
        #   "fit"  = do only the Gaussian fitting step.
        #   "cont" = do only the contraction file step.
        parser.add_argument(
            '-s', '--step',
            dest='step',
            type=str,
            choices=['all', 'scf', 'fit', 'cont'],
            default=self.step,
            help=(
                'Which step(s) to perform. '
                f'Default: {self.step}'
            ),
        )

        # The --element option takes one atomic Z number as
        #   its argument and will perform the requested
        #   operation(s) for that element only.  When set to 0
        #   (the default), all elements are treated.
        parser.add_argument(
            '-e', '--element',
            dest='target_z',
            type=int,
            default=self.target_z,
            help=(
                'Atomic Z number of a single element to '
                'process. 0 means all elements. '
                f'Default: {self.target_z}'
            ),
        )

        # The --fork option allows the user to request that
        #   the Gaussian fitting of the numerical potential
        #   be performed in parallel with the given number
        #   of processes.  (The non-linear fitting can take
        #   a very long time.)
        parser.add_argument(
            '-f', '--fork',
            dest='num_forks',
            type=int,
            default=self.num_forks,
            help=(
                'Number of parallel processes for '
                'Gaussian fitting. '
                f'Default: {self.num_forks}'
            ),
        )

        # The --force option will cause the program to
        #   overwrite any previously obtained results for
        #   the requested operations.
        parser.add_argument(
            '--force',
            dest='force',
            action='store_true',
            default=self.force,
            help=(
                'Overwrite previously obtained results. '
                f'Default: {self.force}'
            ),
        )

        # The --verbose option will print progress reports
        #   of what the program is doing.
        parser.add_argument(
            '-v', '--verbose',
            dest='verbose',
            action='store_true',
            default=self.verbose,
            help=(
                'Print progress reports. '
                f'Default: {self.verbose}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command line values with rc defaults."""
        self.step = args.step
        self.target_z = args.target_z
        self.num_forks = args.num_forks
        self.force = args.force
        self.verbose = args.verbose

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
    """Return the range of element indices to process.

    If target_z is 0, process all elements (1..num_elements).
    Otherwise, process only the single requested element.

    Args:
        target_z (int): Target element Z (0 = all).
        num_elements (int): Total elements in the database.

    Returns:
        range: Python range over the element indices.
    """
    if target_z == 0:
        return range(1, num_elements + 1)
    else:
        return range(target_z, target_z + 1)


def get_database_path():
    """Determine the atomicBDB path from $OLCAO_DATA.

    The atomicBDB directory holds the basis set database:
    atomic SCF results, Gaussian fits, and contraction
    files for every element.

    Returns:
        str: Absolute path to the atomicBDB directory.
    """
    data_dir = os.getenv('OLCAO_DATA')
    if not data_dir:
        sys.exit(
            "Error: $OLCAO_DATA is not set. "
            "See instructions."
        )
    return os.path.join(data_dir, 'atomicBDB')


def get_element_data(element):
    """Extract all element-specific data needed for SCF and
    contraction file generation.

    This function mirrors the Perl getElementData subroutine.
    It gathers Gaussian wave function parameters, orbital
    occupancies, and charge distributions from the element
    database, and computes derived quantities needed by both
    the atomic SCF input generator and the contraction file
    generator.

    The contraction files come in two variants:
      - "core" (index 0): core orbitals are sequestered and
        orthogonalized out of the secular equation.
      - "nocore" (index 1): all orbitals (including what
        would normally be core) are placed in the valence
        basis so that they are included in the secular
        equation.

    Args:
        element (int): Atomic Z number (1-indexed).

    Returns:
        dict: A dictionary with these keys:

        Gaussian expansion parameters:
            min_term_wf (float): Min Gaussian alpha for
                all orbital expansions.
            max_term_wf (float): Max Gaussian alpha for
                all orbital expansions.
            num_terms_wf (int): Number of Gaussians to
                use (taken from the s-type entry).

        Contraction file data:
            contract_core (list of list): Two [s,p,d,f]
                lists of core orbital counts.
                [0] = with-core, [1] = nocore (zeros).
            contract_vale (list of dict): Two dicts
                mapping basis level (1=MB, 2=FB, 3=EB)
                to [s,p,d,f] valence orbital counts.
                [0] = with-core, [1] = nocore.
            total_orbitals (list): [s,p,d,f] total
                orbital counts (core + MB + FB + EB).
            orbital_terms (list): Gaussian selection mask
                string per qn_l, or None if unused.

        Atomic SCF data:
            max_qn_l (int): Highest occupied angular
                momentum quantum number.
            num_total_core (int): Sum of core orbital
                counts across all qn_l.
            num_total_vale (int): Sum of occupied valence
                orbital counts across all qn_l.
            num_vale_orbitals (list): [s,p,d,f] counts
                of occupied valence orbitals.
            orbital_qn_n (list of list): Per-qn_l lists
                of principal quantum numbers for each
                occupied valence orbital.
            up_charge (list of list): Per-qn_l lists of
                up-spin electron counts per orbital.
            dn_charge (list of list): Per-qn_l lists of
                down-spin electron counts per orbital.
    """

    # Get the Gaussian wave function data.  The min and max
    #   define the range of Gaussian exponent (alpha) values.
    #   The number of terms is per angular momentum type
    #   [s, p, d, f]; we use the s-type count as the
    #   representative value.
    min_term_wf = _ed.min_term_wf[element]
    max_term_wf = _ed.max_term_wf[element]
    temp_num_terms_wf = _ed.get_num_terms_wf(element)
    num_terms_wf = temp_num_terms_wf[0]

    # Get the valence charge data for this element.
    #   vale_charge[qn_l] = number of valence electrons
    #   in orbital type qn_l.
    vale_charge = _ed.get_vale_charge(element)

    # Get the current element's core orbitals for both
    #   contract file variants.
    #   core variant (index 0): actual core orbital counts.
    #   nocore variant (index 1): all zeros (no core
    #   orbitals -- everything goes into the valence).
    core_orbitals_0 = list(
        _ed.get_core_orbitals(element)
    )
    core_orbitals_1 = [0, 0, 0, 0]

    # Get the current element's valence orbitals for the
    #   core case (variant 0).  Three basis levels:
    #   MB (minimal beyond core), FB (full beyond MB),
    #   EB (extended beyond FB).
    vale_orbitals_0 = {}
    for basis in range(1, 4):
        vale_orbitals_0[basis] = list(
            _ed.get_vale_orbitals(basis, element)
        )

    # Create the minimal basis for the nocore case
    #   (variant 1) as the sum of the core orbitals and
    #   the minimal basis from the core case.  This is
    #   because in the nocore variant, what would be core
    #   orbitals are folded into the valence MB level.
    vale_orbitals_1 = {}
    vale_orbitals_1[1] = [
        vale_orbitals_0[1][qn_l] + core_orbitals_0[qn_l]
        for qn_l in range(4)
    ]

    # Copy the full and extended basis definitions from
    #   the core variant to the nocore variant (these are
    #   the same in both cases).
    vale_orbitals_1[2] = list(vale_orbitals_0[2])
    vale_orbitals_1[3] = list(vale_orbitals_0[3])

    # Compute the total number of orbitals for each
    #   angular momentum quantum number.  This is the sum
    #   of core orbitals plus all three valence basis
    #   levels.  It is the same for both the core and
    #   nocore variants.
    total_orbitals = list(core_orbitals_0)
    for basis in range(1, 4):
        for qn_l in range(4):
            total_orbitals[qn_l] += (
                vale_orbitals_0[basis][qn_l]
            )

    # Obtain which Gaussian terms are used for which
    #   orbital types.  The selection mask is a compact
    #   boolean string (e.g. '11010') indicating which
    #   terms in the expansion are active.
    orbital_terms = [None] * 4
    for qn_l in range(4):
        if total_orbitals[qn_l] == 0:
            continue
        orbital_terms[qn_l] = (
            _ed.get_orbital_terms(qn_l, element)
        )

    # Compute the total number of core orbitals (sum over
    #   all angular momentum types).
    num_total_core = sum(core_orbitals_0)

    # Compute the number of occupied valence orbitals for
    #   each angular momentum type (s, p, d, f).  An
    #   orbital of type l holds at most (2l+1)*2 electrons.
    #   The number of occupied orbitals is determined by
    #   how many full + partially filled orbitals are
    #   needed to hold the valence charge.
    num_vale_orbitals = [0] * 4
    for qn_l in range(4):
        vc = int(vale_charge[qn_l])
        max_occ = (2 * qn_l + 1) * 2
        num_vale_orbitals[qn_l] = vc // max_occ
        if vc % max_occ > 0:
            num_vale_orbitals[qn_l] += 1

    # Compute the total number of occupied valence
    #   orbitals (sum over all angular momentum types).
    num_total_vale = sum(num_vale_orbitals)

    # Compute the maximum angular momentum quantum number
    #   (qn_l) present in this element's electronic
    #   structure.  This considers both core orbitals and
    #   valence charge.
    max_qn_l = 0
    for qn_l in range(4):
        if core_orbitals_0[qn_l] > 0:
            max_qn_l = qn_l
        if (vale_charge[qn_l] > 0
                and qn_l > max_qn_l):
            max_qn_l = qn_l

    # Compute the descriptive components for each of the
    #   occupied valence orbitals.  This includes the
    #   principal quantum number (QN_n), and the up-spin
    #   and down-spin electron counts for each occupied
    #   orbital.
    orbital_qn_n = [[] for _ in range(4)]
    up_charge = [[] for _ in range(4)]
    dn_charge = [[] for _ in range(4)]

    for qn_l in range(4):

        # Save the total valence charge in all orbitals
        #   of this angular momentum type.
        available = int(vale_charge[qn_l])

        # Maximum number of electrons that can be held by
        #   a single orbital of angular momentum qn_l.
        max_occ = (2 * qn_l + 1) * 2

        for orbital in range(
            1, num_vale_orbitals[qn_l] + 1
        ):
            # Determine the principal quantum number
            #   (QN_n) for this orbital.  It equals the
            #   number of core orbitals of this angular
            #   momentum type plus the orbital index plus
            #   the angular momentum quantum number.
            qn_n = (
                core_orbitals_0[qn_l]
                + orbital + qn_l
            )
            orbital_qn_n[qn_l].append(qn_n)

            # Determine the total charge for this orbital
            #   from the amount available for all orbitals
            #   of the current angular momentum type.
            if available - max_occ >= 0:
                this_charge = max_occ
                available -= this_charge
            else:
                this_charge = available

            # Separate the total charge for this orbital
            #   into up-spin and down-spin parts.
            up_charge[qn_l].append(
                math.floor(this_charge / 2.0)
            )
            dn_charge[qn_l].append(
                math.ceil(this_charge / 2.0)
            )

    return {
        # Gaussian expansion parameters.
        'min_term_wf': min_term_wf,
        'max_term_wf': max_term_wf,
        'num_terms_wf': num_terms_wf,

        # Contraction file data (two variants).
        'contract_core': [
            core_orbitals_0, core_orbitals_1
        ],
        'contract_vale': [
            vale_orbitals_0, vale_orbitals_1
        ],
        'total_orbitals': total_orbitals,
        'orbital_terms': orbital_terms,

        # Atomic SCF data.
        'max_qn_l': max_qn_l,
        'num_total_core': num_total_core,
        'num_total_vale': num_total_vale,
        'num_vale_orbitals': num_vale_orbitals,
        'orbital_qn_n': orbital_qn_n,
        'up_charge': up_charge,
        'dn_charge': dn_charge,
    }


# ----------------------------------------------------------------
#  Core subroutines
# ----------------------------------------------------------------

def make_dirs(target_z, num_elements, atomic_bdb):
    """Create the element directories in the basis database.

    For each element in the target range, create a
    subdirectory (named by element symbol) inside the
    atomicBDB directory.  The database directory itself is
    also created if it does not already exist.

    Args:
        target_z (int): Target element Z (0 = all).
        num_elements (int): Total elements in the database.
        atomic_bdb (str): Path to the atomicBDB directory.
    """
    elem_range = get_element_range(
        target_z, num_elements
    )

    # Create the database directory if it does not
    #   already exist.
    os.makedirs(atomic_bdb, exist_ok=True)

    # Make the directory for each element only if it does
    #   not already exist.
    for element in elem_range:
        elem_dir = os.path.join(
            atomic_bdb, _ed.element_names[element]
        )
        os.makedirs(elem_dir, exist_ok=True)


def make_atom_scf(
    target_z, num_elements, atomic_bdb, force, verbose
):
    """Create atomSCF input files and run the atomSCF program.

    For each element in the target range, this function:
      1. Creates the atomSCF.dat input file specifying the
         element, convergence parameters, and orbital
         configuration.
      2. Runs the atomSCF program to obtain a numerical
         free-atom potential function.

    The atomSCF program produces three output files:
      - atomSCF.out: main output log.
      - numericalPot.dat: numerical potential (used later
        by gaussFit in the fitting step).
      - numericalRho.dat: numerical charge density.

    Args:
        target_z (int): Target element Z (0 = all).
        num_elements (int): Total elements in the database.
        atomic_bdb (str): Path to the atomicBDB directory.
        force (bool): If True, overwrite previous results.
        verbose (bool): If True, print progress messages.
    """
    elem_range = get_element_range(
        target_z, num_elements
    )

    # Create the atomSCF input file and run the program
    #   for each requested element.
    for element in elem_range:
        elem_name = _ed.element_names[element]
        elem_dir = os.path.join(atomic_bdb, elem_name)
        scf_out = os.path.join(elem_dir, 'atomSCF.out')

        # If force is set, remove previous results so the
        #   element will be reprocessed.
        if force and os.path.exists(scf_out):
            for fname in [
                'atomSCF.out',
                'numericalPot.dat',
                'numericalRho.dat',
                'progress.scf',
            ]:
                path = os.path.join(elem_dir, fname)
                if os.path.exists(path):
                    os.remove(path)

        # Do not touch this element if it already is done.
        if os.path.exists(scf_out):
            continue

        # Extract information needed for this element.
        edata = get_element_data(element)

        # Open the atom SCF input file and write all of
        #   the input parameters.
        scf_path = os.path.join(
            elem_dir, 'atomSCF.dat'
        )
        with open(scf_path, 'w') as scf:

            # Title line.
            scf.write(
                "Default Atomic SCF Calculation\n"
            )

            # Element name.
            scf.write("ELEMENT_NAME\n")
            scf.write(f"{elem_name}\n")

            # Nuclear alpha parameter for the Gaussian
            #   representation of the nuclear charge.
            scf.write("NUCLEAR_ALPHA\n")
            scf.write("20.0\n")

            # Maximum number of SCF iterations.
            scf.write("MAX_ITERATION\n")
            scf.write("100\n")

            # SCF convergence tolerance.
            scf.write("CONVG_TOLERANCE\n")
            scf.write("1.0D-8\n")

            # Mixing factor for the SCF iteration
            #   (simple linear mixing of old and new
            #   charge densities).
            scf.write("MIXING_FACTOR\n")
            scf.write("0.4\n")

            # Exchange-correlation functional code.
            scf.write("EXCHANGE_CORRELATION_CODE\n")
            scf.write("1\n")

            # Spin polarization flag (0 = no spin
            #   polarization).
            scf.write("SPIN_POLARIZATION_FLAG\n")
            scf.write("0\n")

            # Relativistic flag (0 = non-relativistic).
            scf.write("RELATIVISTIC_FLAG\n")
            scf.write("0\n")

            # Atomic number.
            scf.write("ATOMIC_NUMBER\n")
            scf.write(f"{element}\n")

            # Shell charge for an external confining
            #   potential (0 = no shell).
            scf.write("SHELL_CHARGE\n")
            scf.write("0\n")

            # Shell radius for an external confining
            #   potential (0 = no shell).
            scf.write("SHELL_RADIUS\n")
            scf.write("0\n")

            # Radial grid parameters:
            #   max_radius  grid_density  core_cutoff
            scf.write("RADIAL_GRID_PARAMETERS\n")
            scf.write("80.0 5.0 50.0\n")

            # Maximum orbital angular momentum quantum
            #   number present in this element.
            scf.write("MAX_ORB_ANGMOM_QN\n")
            scf.write(f"{edata['max_qn_l']}\n")

            # Number of core orbitals.
            scf.write("NUM_CORE_ORBITALS\n")
            scf.write(
                f"{edata['num_total_core']}\n"
            )

            # Number of occupied valence orbitals,
            #   followed by a descriptor line for each:
            #   QN_n  QN_l  up_charge  dn_charge
            scf.write("NUM_VALE_ORBITALS\n")
            scf.write(
                f"{edata['num_total_vale']}\n"
            )
            for qn_l in range(4):
                n_orb = (
                    edata['num_vale_orbitals'][qn_l]
                )
                for idx in range(n_orb):
                    qn_n = (
                        edata['orbital_qn_n']
                        [qn_l][idx]
                    )
                    up = (
                        edata['up_charge']
                        [qn_l][idx]
                    )
                    dn = (
                        edata['dn_charge']
                        [qn_l][idx]
                    )
                    scf.write(
                        f"{qn_n} {qn_l}"
                        f" {up} {dn}\n"
                    )

        if verbose:
            print(
                f"Working on atomic SCF for element "
                f"{element} ({elem_name})"
            )

        # Execute the atomSCF program.  It reads
        #   atomSCF.dat from its working directory and
        #   produces atomSCF.out, numericalPot.dat, and
        #   numericalRho.dat.
        progress = os.path.join(
            elem_dir, 'progress.scf'
        )
        with open(progress, 'w') as prog:
            subprocess.run(
                ['atomSCF'],
                stdout=prog,
                cwd=elem_dir,
            )


def _fit_one_element(element, atomic_bdb, force):
    """Fit the potential for a single element (worker).

    This function is designed to be called either directly
    (serial mode) or via a multiprocessing pool (parallel
    mode).  It runs the gaussFit program which reads the
    numerical potential (from atomSCF) and performs a
    non-linear least-squares fit with Gaussian functions.

    The gaussFit program reads the numerical potential from
    stdin and writes:
      - gaussFit.out: main output log.
      - gauss.fit: fitted Gaussian coefficients.

    Args:
        element (int): Atomic Z number.
        atomic_bdb (str): Path to the atomicBDB directory.
        force (bool): If True, overwrite previous results.
    """
    elem_name = _ed.element_names[element]
    elem_dir = os.path.join(atomic_bdb, elem_name)
    done_file = os.path.join(elem_dir, 'gaussFit.out')

    # If force is set, remove previous results so the
    #   element will be reprocessed.
    if force and os.path.exists(done_file):
        for fname in [
            'gaussFit.out', 'gauss.fit', 'progress.fit'
        ]:
            path = os.path.join(elem_dir, fname)
            if os.path.exists(path):
                os.remove(path)

    # Do not touch this element if it already is done.
    if os.path.exists(done_file):
        return

    # Run the Gaussian fitting program.
    #   gaussFit reads the numerical potential from stdin
    #   and writes results to the current directory.
    pot_file = os.path.join(
        elem_dir, 'numericalPot.dat'
    )
    progress = os.path.join(
        elem_dir, 'progress.fit'
    )
    with open(pot_file, 'r') as infile, \
         open(progress, 'w') as outfile:
        subprocess.run(
            ['gaussFit'],
            stdin=infile,
            stdout=outfile,
            cwd=elem_dir,
        )


def make_gauss_fit(
    target_z, num_elements, num_forks, atomic_bdb, force
):
    """Fit the numerical potentials with Gaussian functions.

    For each element in the target range, this function runs
    the gaussFit program which performs a non-linear
    least-squares fit of the numerical potential (produced
    by atomSCF) with a set of Gaussian functions.  This step
    can be time-consuming for all 103 elements.

    If num_forks > 1, the fitting is performed in parallel
    using a multiprocessing pool.

    Args:
        target_z (int): Target element Z (0 = all).
        num_elements (int): Total elements in the database.
        num_forks (int): Number of parallel processes.
        atomic_bdb (str): Path to the atomicBDB directory.
        force (bool): If True, overwrite previous results.
    """
    elem_range = get_element_range(
        target_z, num_elements
    )

    if num_forks <= 1:
        # Serial execution.
        for element in elem_range:
            _fit_one_element(
                element, atomic_bdb, force
            )
    else:
        # Parallel execution using a process pool.
        with mp.Pool(processes=num_forks) as pool:
            pool.starmap(
                _fit_one_element,
                [
                    (elem, atomic_bdb, force)
                    for elem in elem_range
                ],
            )


def make_contracts(
    target_z, num_elements, atomic_bdb, force
):
    """Create the contraction input files for each element.

    For each element, two contraction files are created:

      1. contract1_<element>:
         Core orbitals are specified explicitly and will
         be orthogonalized out of the secular (eigenvalue)
         equation.  This is the standard OLCAO approach
         and reduces the size of the eigenvalue problem.

      2. nocore_contract1_<element>:
         All orbitals (including those normally treated as
         core) are placed in the valence section.  This
         means they will be included in the secular
         equation.  This variant is useful for comparing
         all-valence calculations.

    Each file specifies:
      - Element name
      - Core orbital counts [s, p, d, f]
      - Valence orbital counts at each basis level
        (MB = minimal, FB = full, EB = extended)
      - Gaussian expansion parameters
        (number of terms, min/max alpha)
      - Atomic number and nuclear alpha
      - Gaussian term selection masks per orbital type
      - A flag to skip charge density calculation

    NOTE: The contraction is not performed here.  These
    files are consumed later by the 'contract' program
    when 'makeinput' is invoked.

    Args:
        target_z (int): Target element Z (0 = all).
        num_elements (int): Total elements in the database.
        atomic_bdb (str): Path to the atomicBDB directory.
        force (bool): If True, overwrite previous results.
    """
    elem_range = get_element_range(
        target_z, num_elements
    )

    # Write the contract files for each requested element.
    for element in elem_range:
        elem_name = _ed.element_names[element]
        elem_dir = os.path.join(atomic_bdb, elem_name)

        # Define the file names for the two variants:
        #   [0] = with core orbitals sequestered.
        #   [1] = no core (all orbitals in valence).
        file_names = [
            os.path.join(
                elem_dir,
                f"contract1_{elem_name}"
            ),
            os.path.join(
                elem_dir,
                f"nocore_contract1_{elem_name}"
            ),
        ]

        # Clean up if we are forcing a recomputation.
        if force and os.path.exists(file_names[0]):
            for fname in file_names:
                if os.path.exists(fname):
                    os.remove(fname)

        # Extract information needed for this element.
        edata = get_element_data(element)

        # Write both contract file variants.
        #   file_idx 0 = with-core, 1 = nocore.
        for file_idx in range(2):
            core = edata['contract_core'][file_idx]
            vale = edata['contract_vale'][file_idx]

            with open(
                file_names[file_idx], 'w'
            ) as cf:

                # Write the element name.
                cf.write("ELEMENT_NAME\n")
                cf.write(f"{elem_name}\n")

                # Write the core orbital counts
                #   [s, p, d, f].
                cf.write("NUM_CORE_ORBITALS\n")
                cf.write(
                    ' '.join(
                        str(x) for x in core
                    ) + '\n'
                )

                # Write the valence orbital counts at
                #   each basis level beyond the previous
                #   level.
                #   MB = minimal beyond core.
                cf.write(
                    "NUM_VALE_ORBITALS_MB\n"
                )
                cf.write(
                    ' '.join(
                        str(x) for x in vale[1]
                    ) + '\n'
                )
                #   FB = full beyond MB.
                cf.write(
                    "NUM_VALE_ORBITALS_FB\n"
                )
                cf.write(
                    ' '.join(
                        str(x) for x in vale[2]
                    ) + '\n'
                )
                #   EB = extended beyond FB.
                cf.write(
                    "NUM_VALE_ORBITALS_EB\n"
                )
                cf.write(
                    ' '.join(
                        str(x) for x in vale[3]
                    ) + '\n'
                )

                # Write the maximum number of Gaussians
                #   for the orbital expansion.
                cf.write("MAX_NUM_GAUSSIANS\n")
                cf.write(
                    f"{edata['num_terms_wf']}\n"
                )

                # Write the min and max exponential
                #   alphas for the Gaussian orbital
                #   expansion.
                cf.write("MIN_MAX_ALPHAS\n")
                cf.write(
                    f"{edata['min_term_wf']} "
                    f"{edata['max_term_wf']}\n"
                )

                # Write the atomic Z number (as float,
                #   matching original Perl output format).
                cf.write("ATOMIC_NUMBER\n")
                cf.write(f"{element:f}\n")

                # Write the nuclear alpha parameter.
                cf.write("NUCLEAR_ALPHA\n")
                cf.write("20.0\n")

                # Write the list of which Gaussian terms
                #   to use for each orbital type.  Only
                #   orbital types with non-zero total
                #   orbitals are included.  The selection
                #   mask is a compact boolean string (e.g.
                #   '1111001111') where each character
                #   indicates whether that Gaussian term
                #   is active.
                for qn_l in range(4):
                    tot = edata['total_orbitals'][qn_l]
                    if tot == 0:
                        continue
                    tag = ORBITAL_TAGS[qn_l]
                    cf.write(
                        f"{tag}_GAUSSIAN_LIST\n"
                    )
                    mask = edata['orbital_terms'][qn_l]
                    cf.write(f"{mask}\n")

                # Make sure that we do not also compute
                #   the charge density at this stage.
                cf.write("CALCULATE_CHARGE\n")
                cf.write("0\n")


# ----------------------------------------------------------------
#  Main
# ----------------------------------------------------------------

def main():
    """Main entry point for makeBDB.

    Orchestrates the three steps of basis set construction:
      1. Atomic SCF (produces numerical potential).
      2. Gaussian fitting (fits potential with Gaussians).
      3. Contraction file generation (creates input for the
         contract program).
    """

    # Get script settings from a combination of the resource
    #   control file and parameters given by the user on the
    #   command line.
    settings = ScriptSettings()

    # Get the number of elements from the element data
    #   module.
    num_elements = _ed.num_elements

    # Get the database path from the environment.
    atomic_bdb = get_database_path()

    # For each element, create the directory.
    make_dirs(
        settings.target_z, num_elements, atomic_bdb
    )

    # Enter each directory and create the atomSCF.dat
    #   input files and run the atomic SCF calculations.
    if settings.step in ('all', 'scf'):
        make_atom_scf(
            settings.target_z,
            num_elements,
            atomic_bdb,
            settings.force,
            settings.verbose,
        )

    # Enter each directory and fit the potentials with
    #   Gaussian functions.
    if settings.step in ('all', 'fit'):
        make_gauss_fit(
            settings.target_z,
            num_elements,
            settings.num_forks,
            atomic_bdb,
            settings.force,
        )

    # Enter each directory and create the two contraction
    #   files (with-core and nocore variants).
    if settings.step in ('all', 'cont'):
        make_contracts(
            settings.target_z,
            num_elements,
            atomic_bdb,
            settings.force,
        )


if __name__ == '__main__':
    # Everything before this point was a subroutine
    #   definition or a request to import information from
    #   external modules.  Only now do we actually start
    #   running the program.  The purpose of this is to
    #   allow another python program to import *this* script
    #   and call its functions internally.
    main()
