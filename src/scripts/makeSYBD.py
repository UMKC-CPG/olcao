#!/usr/bin/env python3

# PROGRAM: makeSYBD.py
# PURPOSE: To write a data file that is usable by plotting software
#          (e.g. Origin, Veusz, gnuplot) to produce a symmetric band
#          (SYBD) diagram. The script reads three files produced by
#          an OLCAO SYBD calculation and produces a single merged
#          output file with eigenvalue data annotated with high-
#          symmetry k-point labels and positions.
#
# AUTHOR:  Paul Rulis
# LAST MODIFIED: Dec. 11, 2019 (Perl); Mar. 2026 (Python port)
#
# USAGE:
#   makeSYBD.py [-dat DATFILE] [-out OUTFILE] [-raw RAWFILE]
#               [-plot PLOTFILE] [-basis BASIS] [-edge EDGE]
#               [-help]
#
# If no file parameters are given then it is assumed that the files
#   to be used are "setup-fb.dat", "gs_sybd-fb.out", and
#   "gs_sybd-fb.raw". The default basis value is "fb", and the
#   default edge is "gs". The default output plot file will be
#   "gs_sybd-fb.plot".
#
# DETAILED DESCRIPTION:
#   The OLCAO SYBD (symmetric band) calculation produces eigenvalue
#   data along high-symmetry paths in the Brillouin zone. This
#   script post-processes that output into a format suitable for
#   direct plotting of band structure diagrams.
#
#   The script performs the following steps:
#
#   1) Reads the setup .dat input file to find the SYBD_INPUT_DATA
#      section. From that section it extracts the number of high-
#      symmetry k-points and their names. The special name "GAMMA"
#      is shortened to "G" for display purposes.
#
#   2) Reads the SYBD .out file to find the index numbers of the
#      high-symmetry k-points (lines containing "HIGH SYMMETRY
#      K POINT") and then reads further to find the position values
#      associated with those index numbers.
#
#   3) Detects path breaks: when two consecutive high-symmetry
#      k-points share the same position value, it means the band
#      path jumps from one segment to another. In that case the
#      two names are merged into a "A,B" comma-separated label
#      and the duplicate position entry is removed.
#
#   4) Reads the SYBD .raw file which contains the eigenvalue data.
#      Each k-point's data may span multiple lines (groups of 11
#      eigenvalues per line). The script concatenates these lines
#      and appends the high-symmetry k-point position and label
#      at the end of each annotated k-point line. The result is
#      written to the .plot output file.

import argparse as ap
import os
import sys
from datetime import datetime


# ----------------------------------------------------------------- #
#                       Script Settings                              #
# ----------------------------------------------------------------- #

class ScriptSettings():
    """The instance variables of this object are the user settings
       that control the program. The variable values are pulled from
       a list that is created within a resource control file and that
       are then reconciled with command line parameters."""


    def __init__(self):
        """Define default values for the parameters by pulling
        them from the resource control file in the default
        location: $OLCAO_RC/makeSYBDrc.py or from the current
        working directory if a local copy of makeSYBDrc.py is
        present."""

        # Read default variables from the resource control file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit(
                "Error: $OLCAO_RC is not set. See instructions."
            )
        sys.path.insert(1, rc_dir)
        from makeSYBDrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults file.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc file.
        self.reconcile(args)

        # At this point, the command line parameters are set and
        #   accepted. When this initialization subroutine returns
        #   the script will start running. So, we use this as a
        #   good spot to record the command line parameters that
        #   were used.
        self.recordCLP()


    def assign_rc_defaults(self, default_rc):

        # File name parameters. When empty, the default names
        #   are constructed from the basis and edge parameters.
        self.dat_file = default_rc["dat_file"]
        self.out_file = default_rc["out_file"]
        self.raw_file = default_rc["raw_file"]
        self.plot_file = default_rc["plot_file"]

        # Naming component parameters.
        self.basis = default_rc["basis"]
        self.edge = default_rc["edge"]


    def parse_command_line(self):

        # Create the parser tool.
        prog_name = "makeSYBD.py"

        description_text = """
Post-process OLCAO SYBD (symmetric band) calculation output
into a plot-ready file for band structure diagrams.

Reads three files from an OLCAO SYBD calculation:
  - A setup .dat file (k-point names and path definition)
  - A SYBD .out file (k-point index numbers and positions)
  - A SYBD .raw file (eigenvalue data for each k-point)

Produces a single .plot file with eigenvalue data annotated
with high-symmetry k-point labels and positions, suitable
for direct use with Origin, Veusz, gnuplot, or similar
plotting software.

When two consecutive high-symmetry k-points share the same
position (indicating a path break in the Brillouin zone),
their labels are merged into a comma-separated "A,B" form.
"""

        epilog_text = """
Originally by Paul Rulis.
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults given in ./makeSYBDrc.py or $OLCAO_RC/makeSYBDrc.py.
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

        # Define the dat_file argument.
        parser.add_argument(
            '-dat', dest='dat_file', type=str,
            default=self.dat_file,
            help=(
                'The setup .dat input file that contains the '
                'SYBD_INPUT_DATA section with high-symmetry '
                'k-point definitions. '
                f'Default: {self.dat_file or "setup-<basis>.dat"}'
            ),
        )

        # Define the out_file argument.
        parser.add_argument(
            '-out', dest='out_file', type=str,
            default=self.out_file,
            help=(
                'The SYBD .out output file containing high-'
                'symmetry k-point index numbers and position '
                'values. '
                f'Default: {self.out_file or "<edge>_sybd-<basis>.out"}'
            ),
        )

        # Define the raw_file argument.
        parser.add_argument(
            '-raw', dest='raw_file', type=str,
            default=self.raw_file,
            help=(
                'The SYBD .raw file containing eigenvalue data '
                'for every k-point along the band path. '
                f'Default: {self.raw_file or "<edge>_sybd-<basis>.raw"}'
            ),
        )

        # Define the plot_file argument.
        parser.add_argument(
            '-plot', dest='plot_file', type=str,
            default=self.plot_file,
            help=(
                'The output .plot file to write. This file '
                'contains the merged eigenvalue data with '
                'high-symmetry k-point annotations. '
                f'Default: {self.plot_file or "<edge>_sybd-<basis>.plot"}'
            ),
        )

        # Define the basis argument.
        parser.add_argument(
            '-basis', dest='basis', type=str,
            default=self.basis,
            help=(
                'The basis set tag used to construct default '
                'file names (e.g. fb, mb, eb). '
                f'Default: {self.basis}'
            ),
        )

        # Define the edge argument.
        parser.add_argument(
            '-edge', dest='edge', type=str,
            default=self.edge,
            help=(
                'The edge tag used to construct default file '
                'names (e.g. gs for ground state, es for '
                'excited state). '
                f'Default: {self.edge}'
            ),
        )


    def reconcile(self, args):
        self.dat_file = args.dat_file
        self.out_file = args.out_file
        self.raw_file = args.raw_file
        self.plot_file = args.plot_file
        self.basis = args.basis
        self.edge = args.edge

        # Construct default file names from the basis and edge
        #   parameters for any file names that were not explicitly
        #   provided. The basis tag is prefixed with a hyphen
        #   (e.g., "-fb") to match the OLCAO file naming
        #   convention.
        basis_tag = f"-{self.basis}"
        sybd_name = "_sybd"

        if not self.dat_file:
            self.dat_file = f"setup{basis_tag}.dat"
        if not self.out_file:
            self.out_file = (
                f"{self.edge}{sybd_name}{basis_tag}.out"
            )
        if not self.raw_file:
            self.raw_file = (
                f"{self.edge}{sybd_name}{basis_tag}.raw"
            )
        if not self.plot_file:
            self.plot_file = (
                f"{self.edge}{sybd_name}{basis_tag}.plot"
            )


    def recordCLP(self):
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


# ----------------------------------------------------------------- #
#             Read High-Symmetry K-Point Names from .dat             #
# ----------------------------------------------------------------- #

def read_kpoint_names(dat_file):
    """Read the setup .dat file and extract high-symmetry k-point
    names from the SYBD_INPUT_DATA section.

    The SYBD_INPUT_DATA section of the setup .dat file contains:
      - A header line (skipped)
      - A line with space-separated counts whose sum gives the
        total number of high-symmetry k-points
      - One line per k-point, where the last token on each line
        is the k-point name (or "GAMMA" which is shortened to
        "G" for display)

    Args:
        dat_file: Path to the setup .dat file.

    Returns:
        A tuple of (num_high_sym_kp, k_names) where:
          - num_high_sym_kp is the total number of high-symmetry
            k-points
          - k_names is a list of their display names
    """

    with open(dat_file, "r") as f:
        # Read forward until we find the SYBD_INPUT_DATA marker.
        for line in f:
            if "SYBD_INPUT_DATA" in line:
                break

        # Extract the number of high symmetry kpoints from the
        #   next two lines. The first line after the marker is
        #   skipped (it is a header). The second line contains
        #   space-separated integers whose sum is the total
        #   number of high-symmetry k-points.
        next(f)  # Skip one line.
        count_line = next(f)
        num_high_sym_kp = sum(
            int(x) for x in count_line.split()
        )

        # Read one line per high-symmetry k-point to get the
        #   name. The special name "GAMMA" is replaced with "G"
        #   for brevity in plot labels.
        k_names = []
        for _ in range(num_high_sym_kp):
            line = next(f)
            if "GAMMA" in line:
                k_names.append("G")
            else:
                tokens = line.strip().split()
                k_names.append(tokens[-1])

    return num_high_sym_kp, k_names


# ----------------------------------------------------------------- #
#          Read High-Symmetry K-Point Positions from .out            #
# ----------------------------------------------------------------- #

def read_kpoint_positions(out_file, k_names):
    """Read the SYBD .out file to find high-symmetry k-point index
    numbers and then look up their position values.

    The SYBD output file contains lines like:
      HIGH SYMMETRY K POINT   <index>
    followed later by a table of k-point data where each line has:
      <col0>  <kpoint_number>  <position_value>  ...

    This function first collects the index numbers from the
    "HIGH SYMMETRY K POINT" lines, then reads the k-point table
    to find the position value associated with each index.

    Args:
        out_file: Path to the SYBD .out file.
        k_names: List of k-point names (used to know how many
                 indices to collect).

    Returns:
        A list of position values (as strings) corresponding to
        each high-symmetry k-point.
    """

    with open(out_file, "r") as f:

        # Start with an empty list of values for the special
        #   K point positions. First we fill it with the index
        #   numbers corresponding to the special k points.
        k_values = []
        for line in f:
            # This block collects the special k point index
            #   numbers from lines containing
            #   "HIGH SYMMETRY K POINT".
            if "HIGH SYMMETRY K POINT" in line:
                tokens = line.strip().split()
                k_values.append(tokens[-1])

                # We've already discovered all the high symmetry
                #   kpoint names from the .dat input file. Now
                #   we check if we have discovered an equal
                #   number of index numbers from the first part
                #   of the kpoint description in the SYBD output.
                #   Once we have discovered all the associated
                #   index values, then we quit.
                if len(k_values) == len(k_names):
                    break

        # Now, we read the positions of the special high symmetry
        #   kpoints using the index numbers found in the above
        #   block. (We are continuing to read from the SYBD
        #   output file.)
        current = 0  # Which high symmetry kpoint we seek.
        for line in f:
            tokens = line.strip().split()
            if len(tokens) >= 3:
                # If we found the matching index number, then
                #   replace the index with the position value.
                if tokens[1] == k_values[current]:
                    k_values[current] = tokens[2]
                    current += 1
                    if current == len(k_values):
                        break

    return k_values


# ----------------------------------------------------------------- #
#              Merge Consecutive Identical Positions                  #
# ----------------------------------------------------------------- #

def merge_path_breaks(k_values, k_names):
    """Detect and merge path breaks in the k-point list.

    Go through the list and find any cases with two consecutive
    equal position values. That implies that we have found a
    change from one path segment to another in the Brillouin
    zone. In that case we merge the two names into a comma-
    separated pair (e.g. "X,Y") and remove the duplicate
    position entry.

    As we do this, the length of the k_values and k_names lists
    will shrink. Therefore, we iterate carefully to avoid going
    beyond the end of the lists.

    Args:
        k_values: List of position values (modified in place).
        k_names: List of k-point names (modified in place).

    Returns:
        The (modified) k_values and k_names lists.
    """

    # We iterate with an explicit index because the lists shrink
    #   as we merge entries. Start from index 1 (the second
    #   element) so we can compare with the previous element.
    kp = 1
    while kp < len(k_values):
        # If the current and previous position values are equal
        #   then we have found a path break.
        if k_values[kp] == k_values[kp - 1]:
            # Remove one of the identical position values.
            k_values.pop(kp)

            # Remove both names and replace the first with a
            #   merged "A,B" form of name.
            merged_name = f"{k_names[kp-1]},{k_names[kp]}"
            k_names[kp - 1] = merged_name
            k_names.pop(kp)
            # Do not increment kp because the list shrank and
            #   the next element is now at the current index.
        else:
            kp += 1

    return k_values, k_names


# ----------------------------------------------------------------- #
#           Read Raw Data and Write Annotated Plot File               #
# ----------------------------------------------------------------- #

def write_plot_file(raw_file, plot_file, k_values, k_names):
    """Read the SYBD .raw file and write the annotated .plot file.

    The .raw file begins with a header line containing:
      <something>  <num_kpoints>  <num_states>  ...

    The eigenvalue data for each k-point may span multiple lines.
    Each k-point starts with one line, followed by additional
    continuation lines. The number of continuation lines is
    determined by the total number of states: eigenvalues are
    written 11 per line, so the total number of lines per k-point
    is ceil(num_states / 11).

    For each k-point, all its lines are concatenated (without
    newlines) into a single output line. If the k-point
    corresponds to a high-symmetry point, the position value and
    name are appended at the end of that line.

    Args:
        raw_file: Path to the SYBD .raw file.
        plot_file: Path to the output .plot file to write.
        k_values: List of high-symmetry k-point position values.
        k_names: List of high-symmetry k-point display names.
    """

    with open(raw_file, "r") as raw_f:
        # Read the header line to get the number of k-points
        #   and the number of states.
        header = raw_f.readline()
        tokens = header.split()
        num_kpoints = int(tokens[1])
        num_states_orig = int(tokens[2])

        # The eigenvalues are written 11 per line. If the total
        #   number of states is not a multiple of 11, we round
        #   up to determine the effective number of states for
        #   line-counting purposes.
        if num_states_orig % 11 != 0:
            num_states = ((num_states_orig // 11) + 1) * 11
        else:
            num_states = num_states_orig

        # Number of continuation lines per k-point (lines
        #   after the first line for each k-point).
        num_cont_lines = num_states // 11

        # Index tracking which high-symmetry k-point annotation
        #   to append next.
        current = 0

        with open(plot_file, "w") as plot_f:
            for i in range(num_kpoints):
                # Read the first line for this k-point.
                line = raw_f.readline().rstrip("\n")
                out_line = line

                # Read and concatenate continuation lines.
                for _ in range(num_cont_lines):
                    line = raw_f.readline().rstrip("\n")
                    out_line += line

                # If there are remaining high-symmetry k-point
                #   annotations to append, add the position
                #   value and name.
                if current <= len(k_values) - 1:
                    out_line += (
                        f" {k_values[current]}"
                        f" {k_names[current]}"
                    )
                    current += 1

                plot_f.write(out_line + "\n")


# ----------------------------------------------------------------- #
#                          Main Program                              #
# ----------------------------------------------------------------- #

def main():

    # Get script settings from a combination of the resource
    #   control file and parameters given by the user on the
    #   command line.
    settings = ScriptSettings()

    # Step 1: Read the high-symmetry k-point names from the
    #   setup .dat input file.
    num_high_sym_kp, k_names = read_kpoint_names(
        settings.dat_file
    )

    # Step 2: Read the high-symmetry k-point index numbers
    #   and position values from the SYBD .out file.
    k_values = read_kpoint_positions(
        settings.out_file, k_names
    )

    # Step 3: Detect path breaks (consecutive identical
    #   positions) and merge their labels.
    k_values, k_names = merge_path_breaks(k_values, k_names)

    # Step 4: Read the .raw eigenvalue data and write the
    #   annotated .plot file.
    write_plot_file(
        settings.raw_file, settings.plot_file,
        k_values, k_names
    )


if __name__ == '__main__':
    # Everything before this point was a subroutine definition
    #   or a request to import information from external modules.
    #   Only now do we actually start running the program. The
    #   purpose of this is to allow another python program to
    #   import *this* script and call its functions internally.
    main()
