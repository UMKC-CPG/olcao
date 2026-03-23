#!/usr/bin/env python3

"""vasp2skl.py -- Convert a VASP position file to an OLCAO
skeleton file.

PROGRAM:  vasp2skl
PURPOSE:  This program will read in a VASP output data file
          and produce an olcao.skl skeleton input file.

          The VASP position file (CONTCAR or POSCAR) contains
          lattice vectors in Cartesian x,y,z form and atomic
          positions in fractional (direct) coordinates.  This
          script extracts those data, computes the conventional
          cell parameters (a, b, c, alpha, beta, gamma), and
          writes them in OLCAO skeleton format.

          The element names are obtained from the POTCAR file
          by scanning for "TITEL" lines, which contain the
          element symbol used in the pseudopotential.

          Both VASP 4.6 and VASP 5.x CONTCAR/POSCAR formats
          are supported.  The VASP 5.x format includes an
          extra line with element symbols between the lattice
          vectors and the atom counts; this script detects and
          handles that automatically.

          Selective dynamics lines are also handled: if the
          line after the atom counts is not "direct" (case-
          insensitive), it is skipped and the next line is
          expected to be "direct".

USAGE:    vasp2skl.py [-i INPUT] [-a ATOMFILE] [-o OUTPUT]
                      [--help]

OPTIONS:

  -i/--input    Name of the VASP atomic position file after
                relaxation.  By default, if this option is not
                given then the value is "CONTCAR".  If CONTCAR
                is not found, the script falls back to
                "POSCAR".

  -a/--atoms    Name of the file that contains the names of
                the elements as they appear in the VASP
                position file (e.g. CONTCAR).  By default,
                if this option is not given then the value
                is "POTCAR".

  -o/--output   Name of the output OLCAO skeleton file.  By
                default, if this option is not given then the
                value is "olcao.skl.relaxed".

  --help        Print usage information and exit.
"""

import argparse as ap
import math
import os
import sys
from datetime import datetime


# ----------------------------------------------------------------
#  ScriptSettings
# ----------------------------------------------------------------

class ScriptSettings():
    """User settings that control the program.

    Variable values are pulled from a resource control file
    (vasp2sklrc.py) and reconciled with command line
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
        from vasp2sklrc import parameters_and_defaults
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

        # Name of the VASP atomic position input file.
        #   If the default "CONTCAR" does not exist in the
        #   current directory, fall back to "POSCAR".
        self.input_file = default_rc["input_file"]
        if not os.path.exists(self.input_file):
            self.input_file = "POSCAR"

        # Name of the file containing element names
        #   (typically POTCAR with its TITEL lines).
        self.atom_file = default_rc["atom_file"]

        # Name of the output OLCAO skeleton file.
        self.output_file = default_rc["output_file"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv."""

        prog_name = "vasp2skl.py"

        description_text = """\
PROGRAM:  vasp2skl
PURPOSE:  Read a VASP output data file (CONTCAR or POSCAR)
          and produce an OLCAO skeleton (.skl) input file.
          Supports both VASP 4.6 and VASP 5.x formats, as
          well as selective dynamics.
"""

        epilog_text = """\
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults: ./vasp2sklrc.py or $OLCAO_RC/vasp2sklrc.py.
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

        # The -i option specifies the VASP position file.
        #   This is the VASP atomic position file after
        #   relaxation.  By default, if this option is not
        #   given then the value is "CONTCAR".  If CONTCAR
        #   does not exist, it defaults to "POSCAR".
        parser.add_argument(
            '-i', '--input',
            dest='input_file',
            type=str,
            default=self.input_file,
            help=(
                'VASP position file (CONTCAR or POSCAR). '
                f'Default: {self.input_file}'
            ),
        )

        # The -a option specifies the file containing the
        #   element names.  This is typically the POTCAR
        #   file, which contains TITEL lines that identify
        #   each element used in the calculation.
        parser.add_argument(
            '-a', '--atoms',
            dest='atom_file',
            type=str,
            default=self.atom_file,
            help=(
                'File with element names (POTCAR). '
                f'Default: {self.atom_file}'
            ),
        )

        # The -o option specifies the output skeleton file.
        parser.add_argument(
            '-o', '--output',
            dest='output_file',
            type=str,
            default=self.output_file,
            help=(
                'Output OLCAO skeleton file. '
                f'Default: {self.output_file}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command line values with rc defaults."""
        self.input_file = args.input_file
        self.atom_file = args.atom_file
        self.output_file = args.output_file

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
#  POTCAR element extraction
# ----------------------------------------------------------------

def read_element_list(atom_file):
    """Read element names from a VASP POTCAR file.

    The POTCAR file contains one or more pseudopotential
    blocks.  Each block has a "TITEL" line of the form:

        TITEL  = PAW_PBE Si 05Jan2001

    The element symbol is extracted from the fourth
    whitespace-delimited token (index 3 after splitting on
    whitespace, then index 0 after splitting on underscore).
    This handles cases like "Si_pv" where the element name
    is followed by a sub-type suffix separated by an
    underscore.

    The element names are lowercased to match OLCAO
    convention.

    Args:
        atom_file (str): Path to the POTCAR file.

    Returns:
        list of str: 1-indexed list of lowercase element
            names.  Index 0 is None.
    """
    element_list = [None]  # 1-indexed

    with open(atom_file, 'r') as f:
        for line in f:
            if 'TITEL' in line:
                tokens = line.split()
                # The element+subtype is typically the
                #   4th token (index 3), e.g. "Si_pv".
                #   Split on underscore to get just the
                #   element symbol.
                raw_name = tokens[3] if len(tokens) > 3 \
                    else tokens[-1]
                elem = raw_name.split('_')[0].lower()
                element_list.append(elem)

    return element_list


# ----------------------------------------------------------------
#  VASP position file reader
# ----------------------------------------------------------------

def read_vasp(input_file, element_list):
    """Read a VASP CONTCAR/POSCAR file and extract all
    important structural data.

    The VASP position file has the following format:

        Line 1:  Comment / system title
        Line 2:  Universal scaling factor
        Lines 3-5:  Lattice vectors (Cartesian x,y,z)
        Line 6:  [VASP 5.x only] Element symbols
        Line 6/7:  Number of atoms of each element
        Next line:  "Selective dynamics" (optional)
        Next line:  "Direct" or "Cartesian"
        Remaining:  Atom positions (fractional or Cartesian)

    This function handles both VASP 4.6 and VASP 5.x
    formats.  The VASP 5.x format includes an extra line
    with element symbols between the lattice vectors and
    the atom counts.  This is detected by checking whether
    the first token on that line contains alphabetic
    characters.

    Selective dynamics lines are also handled: if the line
    after the atom counts is not "direct" (case-insensitive),
    it is skipped and the next line is expected to be
    "direct".

    The lattice vectors are scaled by the universal scaling
    factor, and the cell parameters (a, b, c, alpha, beta,
    gamma) are computed from the scaled vectors using the
    standard dot-product formulas.

    Args:
        input_file (str): Path to the VASP position file.
        element_list (list of str): 1-indexed list of
            element names from the POTCAR file.

    Returns:
        dict: A dictionary containing:
            - system_title (str): Comment line from the
              VASP file.
            - a, b, c (float): Lattice magnitudes in
              Angstroms.
            - alpha, beta, gamma (float): Lattice angles
              in degrees.
            - num_atoms (int): Total number of atoms.
            - atom_element (list of str): 1-indexed list
              of element names for each atom.
            - atom_a, atom_b, atom_c (list of float):
              1-indexed lists of fractional coordinates.
    """

    with open(input_file, 'r') as f:
        # Line 1: Comment / system title.
        system_title = f.readline().strip()

        # Line 2: Universal scaling factor.
        scale = float(f.readline().split()[0])

        # Lines 3-5: Lattice vectors in Cartesian x,y,z.
        #   Each vector is scaled by the universal scaling
        #   factor.
        ax, ay, az = _read_scaled_vector(f, scale)
        bx, by, bz = _read_scaled_vector(f, scale)
        cx, cy, cz = _read_scaled_vector(f, scale)

        # Next line: either atom counts (VASP 4.6) or
        #   element symbols (VASP 5.x).  We detect the
        #   format by checking if the first token contains
        #   alphabetic characters.
        # ============= 10/28/2010 3:53 pm =============
        # These extra lines here are to make it compatible
        # with VASP 5.2 CONTCAR.
        values = f.readline().split()
        if any(c.isalpha() for c in values[0]):
            print("CONTCAR file from VASP 5.2")
            values = f.readline().split()
        else:
            print("CONTCAR file from VASP 4.6")
        # ===============================================

        # Parse the atom counts for each element.
        atom_counts = [int(v) for v in values]

        # Compute the total number of atoms.
        num_atoms = sum(atom_counts)

        # Determine the element for each atom.
        #   The number of atoms of each element is in
        #   atom_counts.  We build a 1-indexed list
        #   mapping each atom to its element name.
        #   element_list is 1-indexed (from POTCAR),
        #   while atom_counts is 0-indexed.
        atom_element = [None]  # 1-indexed
        num_elements = len(element_list) - 1  # skip [0]
        for i in range(1, num_elements + 1):
            count = atom_counts[i - 1]  # Note the -1
            for _ in range(count):
                atom_element.append(element_list[i])

        # Compute the cell parameters from the scaled
        #   lattice vectors using dot-product formulas.
        a = math.sqrt(ax*ax + ay*ay + az*az)
        b = math.sqrt(bx*bx + by*by + bz*bz)
        c = math.sqrt(cx*cx + cy*cy + cz*cz)

        alpha = math.acos(
            (bx*cx + by*cy + bz*cz) / (b * c)
        ) * 180.0 / math.pi
        beta = math.acos(
            (ax*cx + ay*cy + az*cz) / (a * c)
        ) * 180.0 / math.pi
        gamma = math.acos(
            (ax*bx + ay*by + az*bz) / (a * b)
        ) * 180.0 / math.pi

        # Read the label for the atom positions.
        # ======= Extra feature for VASP 5.2 =======
        # This also works if selective dynamics is used.
        # Split the line and check if it has "Direct".
        # If not, it is a "selective dynamics" line and
        # we need to read one more line to get "direct".
        line = f.readline()
        line_values = line.split()
        if line_values[0].lower() != "direct":
            line = f.readline()
        # ============================================

        # Read the fractional atom positions.
        atom_a = [None]  # 1-indexed
        atom_b = [None]
        atom_c = [None]
        for _ in range(num_atoms):
            vals = f.readline().split()
            atom_a.append(float(vals[0]))
            atom_b.append(float(vals[1]))
            atom_c.append(float(vals[2]))

    return {
        'system_title': system_title,
        'a': a, 'b': b, 'c': c,
        'alpha': alpha, 'beta': beta, 'gamma': gamma,
        'num_atoms': num_atoms,
        'atom_element': atom_element,
        'atom_a': atom_a,
        'atom_b': atom_b,
        'atom_c': atom_c,
    }


def _read_scaled_vector(f, scale):
    """Read one lattice vector line and scale it.

    Each lattice vector line in the VASP position file
    contains three Cartesian components (x, y, z).  These
    are multiplied by the universal scaling factor to give
    the true lattice vector in Angstroms.

    Args:
        f (file): Open file handle positioned at the line.
        scale (float): Universal scaling factor.

    Returns:
        tuple of float: (x, y, z) scaled vector components.
    """
    vals = f.readline().split()
    x = scale * float(vals[0])
    y = scale * float(vals[1])
    z = scale * float(vals[2])
    return x, y, z


# ----------------------------------------------------------------
#  Skeleton file writer
# ----------------------------------------------------------------

def print_skl(output_file, data):
    """Write the OLCAO skeleton (.skl) file.

    The skeleton file format is:

        title
        <system title>
        end
        cell
        <a> <b> <c> <alpha> <beta> <gamma>
        fractional <num_atoms>
        <element> <frac_a> <frac_b> <frac_c>
        ...
        space 1_a
        supercell 1 1 1
        prim

    The "space 1_a" line specifies space group P1 (no
    symmetry), which is appropriate because the VASP
    structure has already been fully expanded.

    The "supercell 1 1 1" line indicates no supercell
    expansion.

    The "prim" line requests that the primitive cell be
    used.

    Args:
        output_file (str): Path to the output skeleton
            file.
        data (dict): Structural data dictionary as
            returned by read_vasp().
    """
    with open(output_file, 'w') as f:
        # Print the header for the skl file.
        f.write("title\n")
        f.write(f"{data['system_title']}\n")
        f.write("end\n")
        f.write("cell\n")

        # Print the cell parameters.
        f.write(
            f"{data['a']:10.6f} "
            f"{data['b']:10.6f} "
            f"{data['c']:10.6f} "
            f"{data['alpha']:10.6f} "
            f"{data['beta']:10.6f} "
            f"{data['gamma']:10.6f}\n"
        )

        # Print out the atom positions.
        f.write(f"fractional {data['num_atoms']}\n")
        for i in range(1, data['num_atoms'] + 1):
            f.write(
                f"{data['atom_element'][i]:<4s} "
                f"{data['atom_a'][i]:15.9f} "
                f"{data['atom_b'][i]:15.9f} "
                f"{data['atom_c'][i]:15.9f}\n"
            )

        # Print out the trailing information.
        f.write("space 1_a\n")
        f.write("supercell 1 1 1\n")
        f.write("prim\n")


# ----------------------------------------------------------------
#  Main
# ----------------------------------------------------------------

def main():
    """Main entry point for vasp2skl.

    Reads a VASP CONTCAR/POSCAR file, extracts element
    names from the POTCAR file, computes cell parameters,
    and writes the result as an OLCAO skeleton file.

    The workflow mirrors the original Perl script:
      1. Parse command line and set defaults.
      2. Read element names from POTCAR (TITEL lines).
      3. Read the VASP position file (CONTCAR/POSCAR).
      4. Write the OLCAO skeleton file.
    """

    # Get script settings from a combination of the
    #   resource control file and parameters given by the
    #   user on the command line.
    settings = ScriptSettings()

    # Read the element list from the POTCAR file.
    #   Each TITEL line in the POTCAR identifies one
    #   element used in the VASP calculation.
    element_list = read_element_list(settings.atom_file)

    # Read the VASP position file and extract all
    #   important structural data.
    data = read_vasp(settings.input_file, element_list)

    # Write the OLCAO skeleton file.
    print_skl(settings.output_file, data)


if __name__ == '__main__':
    # Everything before this point was a subroutine
    #   definition or a request to import information from
    #   external modules.  Only now do we actually start
    #   running the program.  The purpose of this is to
    #   allow another python program to import *this* script
    #   and call its functions internally.
    main()
