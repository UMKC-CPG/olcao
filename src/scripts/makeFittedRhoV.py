#!/usr/bin/env python3

"""makeFittedRhoV — Fitted charge density and potential evaluator.

PROGRAM:  makeFittedRhoV.py
PURPOSE:  To evaluate computed OLCAO charge density and potential
   functions on a 3D real-space mesh for the purpose of
   visualization in a variety of forms. Data will be created
   to view the potential function, the total charge density,
   the valence charge density, the potential difference, and
   the valence charge difference. The differences are computed
   in comparison to the potential and charge from a single
   isolated atom calculation (as stored in the potential
   database). It is important to understand that the plot of
   the charge density is derived from the fitted charge.

AUTHOR:  Paul Rulis (Perl original); Python port by Claude.
LAST MODIFIED:  Mar. 20, 2026

USAGE:
   makeFittedRhoV.py -m numPointsA numPointsB numPointsC
                     [-i INPUT_FILE] [-spindn] [-noDX]
                     [-nc NUM_CELLS] [-h | --help]

The -m switch is used to define the number of mesh points
   along each a, b, c axis in the (possibly non-orthogonal)
   cell. Note that periodic boundary conditions are observed.
   This is a required switch; if it is not given, the script
   will stop.

The -i option is used to define which file should be used as
   the source of the SCF charge and potential data. If this
   option is not given then the value of "gs_scfV-fb.dat"
   will be used.

The -spindn switch is used specifically for spin-polarized
   calculations where the first and second sections of the
   potential function file are different. Recall that for
   non-spin-polarized calculations the first and second parts
   are identical and thus the second part can be ignored as
   duplicate data. When this switch is turned on, then the
   second section is specifically used as the data source.
   (Thus, if you want to see the spin up and spin down then
   you need to run this program twice, once for each spin
   direction.)

The -noDX switch is used to demand that the openDX plottable
   data *not* be created. By default, this data will be made
   because it does not cost much extra time and it will
   prevent the user from having to run the computation twice
   (once for the profiles and again for the openDX data).

The -nc option sets the number of levels of neighboring
   replicated cells to consider. 0 = original cell only,
   1 = first neighbors (27 total cells), 2 = two levels (125
   total cells). Default: 1.

PROGRAM NOTES:

This program will construct an input file for the OLCAOrhoV
   Fortran90 program and will then call that program with the
   appropriate command line parameters.

The input files for this program are:
   (1) A SCF potential file (defaults to "gs_scfV-fb.dat").
   (2) A structure file (always assumed to be
       "structure.dat").
   (3) A database of SCF potential files for single isolated
       atom calculations.

The output file created by this program is then used by the
   OLCAOrhoV Fortran90 program and it contains the following
   columns of data:
   Column #1: Exponential alphas for each Gaussian.
   Column #2: Potential coefficients for the SCF calculation.
   Column #3: Difference between potential coefficients of
              the SCF and isolated atom calculations.
   Column #4: Valence charge coefficients for the SCF
              calculation.
   Column #5: Difference between valence charge coefficients
              of the SCF and isolated atom calculations.
   Column #6: Total charge coefficients for the SCF
              calculation.

There is no column for the difference in total charge at
   this time although this could be added if such a need ever
   arose in the future.
"""

import argparse as ap
import os
import subprocess
import sys
from datetime import datetime


# ============================================================
# Utility: prep_line
# ============================================================

def prep_line(file_handle):
    """Read one line from an open file and return a list of
    tokens split on whitespace.

    This is the Python equivalent of the Perl
    StructureControl::prepLine subroutine. It reads a single
    line from *file_handle*, strips leading/trailing
    whitespace, and splits on whitespace. Python's
    str.split() with no argument handles leading whitespace
    automatically, so no explicit shift of an empty leading
    token is needed.

    Parameters
    ----------
    file_handle : file object
        An open file from which to read one line.

    Returns
    -------
    list of str
        The whitespace-separated tokens from the line.
        Returns an empty list if EOF is reached.
    """

    line = file_handle.readline()
    if not line:
        return []
    return line.strip().split()


# ============================================================
# ScriptSettings: command-line and rc-file handling
# ============================================================

class ScriptSettings():
    """Holds all user-configurable settings for makeFittedRhoV.

    The instance variables of this object are the user settings
    that control the program. The variable values are pulled
    from a list that is created within a resource control file
    and that are then reconciled with command line parameters.
    """

    def __init__(self):
        """Initialize settings from the rc file and CLI.

        Default values are loaded from the resource control
        file (makeFittedRhoVrc.py) found in $OLCAO_RC or in
        the current working directory. These are then
        overridden by any explicit command-line arguments.
        """

        # Read default variables from the resource control
        # file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit(
                "Error: $OLCAO_RC is not set. "
                "See instructions."
            )
        sys.path.insert(1, rc_dir)
        from makeFittedRhoVrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc
        # file.
        self.reconcile(args)

        # Record the command line parameters that were used.
        self.record_clp()

    def assign_rc_defaults(self, default_rc):
        """Copy rc-file defaults into instance attributes.

        Parameters
        ----------
        default_rc : dict
            Dictionary returned by
            makeFittedRhoVrc.parameters_and_defaults().
        """

        # File name defaults.
        self.struct_file = default_rc["struct_file"]
        self.in_file = default_rc["in_file"]
        self.out_file = default_rc["out_file"]
        self.atom_pot_db = default_rc["atom_pot_db"]

        # Behavioural flags and parameters.
        self.spin_dn = default_rc["spin_dn"]
        self.do_open_dx = default_rc["do_open_dx"]
        self.num_mesh_points = list(
            default_rc["num_mesh_points"]
        )
        self.num_cells = default_rc["num_cells"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv.

        Returns
        -------
        argparse.Namespace
            Parsed arguments.
        """

        prog_name = "makeFittedRhoV.py"

        description_text = """\
Fitted charge density and potential evaluator.

Evaluate computed OLCAO charge density and potential functions
on a 3D real-space mesh for visualization. Produces data for
the potential function, total/valence charge densities, and
their differences relative to isolated atom calculations.

The charge density plotted is derived from the fitted charge.
"""

        epilog_text = """\
Please contact Paul Rulis regarding questions.
Defaults are given in ./makeFittedRhoVrc.py or
$OLCAO_RC/makeFittedRhoVrc.py.
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
        """Register all CLI arguments with the parser.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            The parser to which arguments are added.
        """

        # The -m switch defines the number of mesh points
        # along each a, b, c axis in the (possibly
        # non-orthogonal) cell. Note that periodic boundary
        # conditions are observed. This is a required switch.
        parser.add_argument(
            '-m', '--mesh', nargs=3, dest='num_mesh_points',
            type=int, required=True,
            metavar=('A', 'B', 'C'),
            help=(
                "Number of mesh points along each a, b, c "
                "axis. Required. All values must be > 0."
            ),
        )

        # The -i option defines which file should be used as
        # the source of the SCF charge and potential data. If
        # not given, "gs_scfV-fb.dat" is used.
        parser.add_argument(
            '-i', '--input', dest='in_file', type=str,
            default=self.in_file,
            help=(
                "SCF potential/charge data file. "
                f"Default: {self.in_file}"
            ),
        )

        # The -spindn switch is used for spin-polarized
        # calculations to select the spin-down section of
        # the potential file instead of the total/spin-up
        # section. For non-spin-polarized calculations,
        # the two sections are identical so the spin-down
        # section is ignored by default.
        parser.add_argument(
            '-spindn', '--spindn', dest='spin_dn',
            action='store_true', default=self.spin_dn,
            help=(
                "Use the spin-down section of the potential "
                "file (for spin-polarized calculations)."
            ),
        )

        # The -noDX switch prevents the openDX plottable
        # data from being created. By default, openDX data
        # is made because it does not cost much extra time
        # and it prevents having to run the computation
        # twice (once for profiles, again for openDX data).
        parser.add_argument(
            '-noDX', '--noDX', dest='no_dx',
            action='store_true', default=False,
            help=(
                "Do NOT generate openDX data. By default "
                "openDX data is created alongside profiles."
            ),
        )

        # The -nc option sets the number of levels of
        # neighboring replicated cells to consider when
        # evaluating the charge density and potential
        # functions. 0 = original cell only (1 cell),
        # 1 = first neighbors (27 total cells),
        # 2 = two levels (125 total cells).
        parser.add_argument(
            '-nc', '--num-cells', dest='num_cells',
            type=int, default=self.num_cells,
            help=(
                "Number of levels of neighboring replicated "
                "cells to consider. 0=original only (1 "
                "cell), 1=first neighbors (27 cells), "
                "2=two levels (125 cells). "
                f"Default: {self.num_cells}"
            ),
        )

    def reconcile(self, args):
        """Reconcile parsed CLI args with rc-file defaults.

        After this method, instance attributes reflect the
        final settings: CLI values override rc defaults.

        Parameters
        ----------
        args : argparse.Namespace
            The parsed command-line arguments.
        """

        self.num_mesh_points = list(args.num_mesh_points)
        self.in_file = args.in_file
        self.spin_dn = args.spin_dn
        self.num_cells = args.num_cells

        # The -noDX flag inverts the do_open_dx default.
        if args.no_dx:
            self.do_open_dx = False

        # Validate mesh points: all must be positive
        # integers. (argparse already enforces int type.)
        for idx, val in enumerate(self.num_mesh_points):
            axis_name = ['a', 'b', 'c'][idx]
            if val <= 0:
                sys.exit(
                    f"Error: Mesh point along {axis_name} "
                    f"axis = {val}. "
                    "Mesh points must be > 0."
                )

    def record_clp(self):
        """Append the command line used to the 'command' file.

        This records the invocation for reproducibility,
        following the convention used across OLCAO scripts.
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


# ============================================================
# read_struct: Parse structure.dat for element/type info
# ============================================================

def read_struct(settings):
    """Read the structure file to get element types.

    This subroutine reads the structure.dat file to determine
    the number of atomic types and which chemical element is
    associated with each type. The information is needed to
    locate the correct isolated-atom coefficient files in the
    potential database.

    The structure file format is:
      - Lines 1-4: Cell parameters (skipped here).
      - Line 5: NUM_ATOM_SITES header (skipped).
      - Line 6: Number of atomic and potential sites.
      - Line 7: Header for atomic coordinate data (skipped).
      - Lines 8+: One line per atom site with fields:
        [atom_index, type_id, ..., element_name]

    This function extracts the type-to-element mapping by
    reading each atom line and noting when a new type number
    appears. The element name is always the last token on
    each atom line.

    Parameters
    ----------
    settings : ScriptSettings
        The program settings (provides struct_file path).

    Returns
    -------
    num_types : int
        The number of distinct atomic types.
    num_sites : int
        The total number of atomic/potential sites.
    type_elements : list of str
        Element name for each type, 1-indexed (index 0 is
        None). Element names are lowercase.
    """

    struct_file = settings.struct_file

    with open(struct_file, 'r') as f:

        # Read past the cell parameters (3 lines of lattice
        # vectors + 1 angle line) and the NUM_ATOM_SITES
        # header line. That is 5 lines total to skip.
        for _ in range(5):
            f.readline()

        # Read the number of atomic and potential sites in
        # the model.
        values = prep_line(f)
        num_sites = int(values[0])

        # Read past the header for the list of atomic
        # coordinate data.
        f.readline()

        # Initialize the type counter and the type-element
        # mapping list. Index 0 is unused (1-indexed).
        num_types = 0
        type_elements = [None]

        # Read each atom and collect the appropriate
        # information.
        for _ in range(num_sites):

            # Read and prepare the line for processing.
            values = prep_line(f)

            # Check if this is a new type and if so record
            # the element for this type. The type ID is in
            # the second column (index 1). The element name
            # is always the last token on the line.
            type_id = int(values[1])
            if type_id > num_types:
                num_types += 1
                type_elements.append(values[-1].lower())

    return num_types, num_sites, type_elements


# ============================================================
# compute_diff: SCF vs. isolated-atom difference calculation
# ============================================================

def compute_diff(settings, num_types, type_elements):
    """Compute the difference between SCF and isolated-atom
    coefficients.

    This is the core computation of the script. For each
    atomic type, it reads the corresponding isolated-atom
    coefficient file from the potential database and compares
    it term-by-term with the SCF potential/charge data from
    the converged calculation. The differences reveal how
    the potential and charge density change when atoms
    interact in the solid/molecule compared to when they are
    isolated.

    The SCF potential file has the following structure:
      Line 1: "NUM_TYPES <N>"
      Line 2: "TOTAL__OR__SPIN_UP" label
      For each type:
        Line: <num_terms> (number of Gaussian terms)
        num_terms lines of data with columns:
          col 0: potential coefficient
          col 1: exponential alpha (Gaussian exponent)
          col 2: total charge coefficient
          col 3: valence charge coefficient
          col 4: (additional column, often zero)
      Then: "SPIN_DN" label
      Same structure repeated for spin-down data.

    The isolated-atom coefficient file (coeff.isolated) in
    the database has the same term-by-term structure, with
    the same exponential alphas (Gaussian exponents) as the
    SCF file. This allows a direct term-by-term subtraction.

    The output data for each term contains 6 values:
      [0] Exponential alpha (Gaussian exponent)
      [1] SCF potential coefficient
      [2] Potential difference (SCF - isolated)
      [3] SCF valence charge coefficient
      [4] Valence charge difference (SCF - isolated)
      [5] SCF total charge coefficient

    Parameters
    ----------
    settings : ScriptSettings
        Program settings (provides in_file, atom_pot_db,
        and spin_dn).
    num_types : int
        Number of distinct atomic types.
    type_elements : list of str
        Element name for each type (1-indexed).

    Returns
    -------
    out_values : list of list of float
        The computed output data. Each inner list has 6
        floats. The outer list is 1-indexed (index 0 is
        None).
    num_terms : int
        Total number of terms across all types.
    """

    in_file = settings.in_file
    atom_pot_db = settings.atom_pot_db
    spin_dn = settings.spin_dn

    # Initialize the term counter and output data list.
    # The output list is 1-indexed to match the Perl
    # original (index 0 is a placeholder).
    num_terms = 0
    out_values = [None]

    # Open the input SCF potential file.
    with open(in_file, 'r') as inpot:

        # Read past the line that lists the number of types.
        # (This was already computed from the structure.dat
        # file where the element data was also obtained.)
        # Also read past the "TOTAL__OR__SPIN_UP" label.
        inpot.readline()
        inpot.readline()

        # If we need to do the spin-down data, then read
        # past the first (spin-up/total) section entirely.
        if spin_dn:

            # Iterate over all the types in the system. We
            # will read past each of the terms for each of
            # the types.
            for type_idx in range(1, num_types + 1):

                # Read the number of terms for this type
                # and then iterate past them.
                values = prep_line(inpot)
                type_num_terms = int(values[0])
                for _ in range(type_num_terms):
                    inpot.readline()

            # Read past the "SPIN_DN" label line.
            inpot.readline()

        # Loop over each type and perform a difference
        # calculation between the SCF data and the
        # non-interacting (isolated atom) data.
        for type_idx in range(1, num_types + 1):

            # Define the isolated potential and charge
            # coefficient file path for this element type.
            element = type_elements[type_idx]
            iso_coeff_file = os.path.join(
                atom_pot_db, element, "coeff.isolated"
            )

            # Open the non-interacting (isolated atom)
            # potential file associated with this type.
            with open(iso_coeff_file, 'r') as nonpot:

                # Read the number of terms from the
                # isolated coefficient file and from the
                # SCF potential file. Both files should
                # have the same number of terms, but we
                # read both to advance the file positions.
                prep_line(nonpot)
                prep_line(inpot)

                # Read the terms of this non-interacting
                # data and compute the difference compared
                # to the SCF terms.
                while True:
                    line_np = nonpot.readline()
                    if not line_np:
                        break

                    # Read a line from the input potential
                    # file and increment the count of the
                    # number of terms read.
                    line_ip = inpot.readline()
                    num_terms += 1

                    # Prepare the data for comparison by
                    # splitting on whitespace.
                    values_np = line_np.strip().split()
                    values_ip = line_ip.strip().split()

                    # Convert to floats for arithmetic.
                    vals_np = [
                        float(v) for v in values_np
                    ]
                    vals_ip = [
                        float(v) for v in values_ip
                    ]

                    # Perform a check to make sure that the
                    # comparison is valid. The exponential
                    # alphas (column 1) must match between
                    # the SCF and isolated-atom files.
                    if vals_ip[1] != vals_np[1]:
                        sys.exit(
                            f"Compare failure: type "
                            f"{type_idx}; element "
                            f"{element}; in term = "
                            f"{vals_ip[1]}; iso term = "
                            f"{vals_np[1]}."
                        )

                    # Compute the differences:
                    #   pot_diff = SCF_pot - iso_pot
                    #   val_chg_diff = SCF_val - iso_val
                    pot_diff = vals_ip[0] - vals_np[0]
                    val_chg_diff = vals_ip[3] - vals_np[3]

                    # Create the contents of the new file.
                    # The 6 output columns are:
                    #   [0] Exponential alphas
                    #   [1] SCF Potential Coefficients
                    #   [2] Potential Difference
                    #   [3] Valence Charge Coefficients
                    #   [4] Valence Charge Difference
                    #   [5] Total Charge Coefficients
                    out_values.append([
                        vals_ip[1],  # Exponential alphas
                        vals_ip[0],  # SCF Pot Coeffs
                        pot_diff,    # Pot Difference
                        vals_ip[3],  # Val Charge Coeffs
                        val_chg_diff,  # Val Chg Diff
                        vals_ip[2],  # Total Charge Coeffs
                    ])

    return out_values, num_terms


# ============================================================
# print_temp_data: Write the temporary data file
# ============================================================

def print_temp_data(settings, out_values, num_terms):
    """Print an input file for the OLCAOrhoV Fortran program.

    This writes the temporary data file that the OLCAOrhoV
    Fortran90 program will read. Each line contains 6 columns
    of data in scientific notation format (%16.8e), matching
    the Perl original's output format.

    The columns are:
      Column 1: Exponential alphas for each Gaussian.
      Column 2: Potential coefficients for the SCF
                calculation.
      Column 3: Difference between potential coefficients
                of the SCF and isolated atom calculations.
      Column 4: Valence charge coefficients for the SCF
                calculation.
      Column 5: Difference between valence charge
                coefficients of the SCF and isolated atom
                calculations.
      Column 6: Total charge coefficients for the SCF
                calculation.

    Parameters
    ----------
    settings : ScriptSettings
        Program settings (provides out_file path).
    out_values : list of list of float
        The computed output data (1-indexed; index 0 is
        None).
    num_terms : int
        Total number of terms to write.
    """

    out_file = settings.out_file

    with open(out_file, 'w') as f:
        for term in range(1, num_terms + 1):
            line = ""
            for value in out_values[term]:
                line += f"{value:16.8e}"
            f.write(line + "\n")


# ============================================================
# run_olcao_rho_v: Execute the OLCAOrhoV Fortran program
# ============================================================

def run_olcao_rho_v(settings):
    """Run the OLCAOrhoV Fortran90 program.

    This calls the OLCAOrhoV executable to actually compute
    the charge density and potential functions on the
    requested 3D mesh. The Fortran program reads the
    temporary data file and the structure file, evaluates
    the Gaussian-fitted functions at each mesh point, and
    produces profile data and (optionally) openDX-format
    output.

    Command line parameters passed to OLCAOrhoV:
      1-3: numMeshPoints (a, b, c) — number of mesh points
           along each axis to evaluate the functions on.
      4:   num_cells — number of levels of neighboring
           replicated cells to consider. 0 = original cell
           only (1 cell), 1 = first neighbors (27 total
           cells), 2 = two levels (125 total cells).
      5:   noDX flag — when set to 1, prevents the openDX
           data from being produced. 0 means produce it.
           (The profile data is always produced.)

    Parameters
    ----------
    settings : ScriptSettings
        Program settings (provides out_file, num_mesh_points,
        num_cells, and do_open_dx).
    """

    out_file = settings.out_file
    mesh = settings.num_mesh_points
    num_cells = settings.num_cells
    no_dx_flag = 0 if settings.do_open_dx else 1

    # Build the command. The Perl original hard-coded
    # num_cells=1 and noDX=0.
    cmd = [
        "OLCAOrhoV",
        out_file,
        str(mesh[0]),
        str(mesh[1]),
        str(mesh[2]),
        str(num_cells),
        str(no_dx_flag),
    ]

    print(f"OLCAOrhoV {out_file} "
          f"{mesh[0]} {mesh[1]} {mesh[2]} "
          f"{num_cells} {no_dx_flag}")

    result = subprocess.run(
        cmd, capture_output=True, text=True
    )

    print("Program Done")

    if result.stdout:
        print(result.stdout)

    if result.returncode != 0:
        print(
            f"Error: OLCAOrhoV exited with return code "
            f"{result.returncode}.",
            file=sys.stderr,
        )
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)


# ============================================================
# main: Top-level program flow
# ============================================================

def main():
    """Main program entry point.

    The program flow follows the Perl original:
      1. Initialize settings from the rc file and command
         line (handled by ScriptSettings.__init__).
      2. Read the structure file to get the elements in the
         system and associate elements with types.
      3. Compute the difference between SCF data and
         non-interacting (isolated atom) data.
      4. Print an input file for the OLCAOrhoV Fortran
         program.
      5. Run the OLCAOrhoV Fortran program.
    """

    # Get script settings from a combination of the resource
    # control file and parameters given by the user on the
    # command line.
    settings = ScriptSettings()

    # Read the structure file to get the elements in the
    # system and associate elements with types.
    num_types, num_sites, type_elements = read_struct(
        settings
    )

    # Compute difference between SCF data and
    # non-interacting (isolated atom) data.
    out_values, num_terms = compute_diff(
        settings, num_types, type_elements
    )

    # Print an input file for the OLCAOrhoV Fortran program.
    print_temp_data(settings, out_values, num_terms)

    # Run the OLCAOrhoV Fortran program.
    run_olcao_rho_v(settings)


if __name__ == '__main__':
    # Everything before this point was a subroutine
    # definition or a request to import information from
    # external modules. Only now do we actually start running
    # the program. The purpose of this is to allow another
    # python program to import *this* script and call its
    # functions internally.
    main()
