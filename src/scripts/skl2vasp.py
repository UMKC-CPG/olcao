#!/usr/bin/env python3

"""skl2vasp.py -- Convert an OLCAO skeleton file to VASP input.

PROGRAM:  skl2vasp
PURPOSE:  This program will take an OLCAO skeleton file and
          will produce input files for executing a VASP
          calculation with a few intelligent options.

          Currently, only the POSCAR file is generated.  The
          POTCAR, KPOINTS, and INCAR files were produced by
          the original Perl script but are not yet ported to
          the Python structure_control module.  The command
          line options for potential type, sub-type, GW flag,
          Gamma k-point, job type, and cell type are parsed
          and stored for future use when those files are
          implemented.

USAGE:    skl2vasp.py [-j JOBTYPE] [-p POTTYPE [-s SUBTYPE]]
                      [--gw] [--gamma] [-t CELLTYPE]
                      [-i INPUT] [--help]

OPTIONS:

  -j/--job      Type of calculation to set up.  This will
                choose optimal parameters for some of the
                input files.  Possible options are:
                "relaxfull", "relaxion1", "relaxion2",
                "relaxvol", "cij", "phonon", "gpt".
                Default: "relaxfull".

  -p/--pot      Base type of pseudopotential.  Options:
                "potLDA"      = Ultrasoft LDA
                "potGGA"      = Ultrasoft GGA
                "potpawLDA"   = PAW LDA (VASP 4.x)
                "potpawGGA"   = PAW GGA (VASP 4.x)
                "potpawPBE"   = PAW PBE (VASP 4.x)
                "potpawLDA5x" = PAW LDA (VASP 5.x)
                "potpawPBE5x" = PAW PBE (VASP 5.x)
                NOTE: Do not replace "x" with a number.
                Default: "potpawPBE".

  -s/--subpot   Sub-type of the potential.  Many options are
                valid only for particular elements or base
                types.  User responsibility to know what they
                are asking for.  Options:
                "s"    = Soft
                "h"    = Hard
                "pv"   = Semi-core p states as valence
                "sv"   = Semi-core s and p as valence
                "d"    = Semi-core d states as valence
                Default: "" (standard potential).

  --gw          Append the GW sub-sub-type to the potential
                sub-type.

  --gamma       Use a Gamma-centered k-point instead of a
                general Monkhorst-Pack mesh.

  -t/--celltype Override the Bravais lattice cell type that
                would be auto-detected from the space group.
                Options: "tric", "mono", "ortho", "tet",
                "trig", "hex", "sc", "bcc", "fcc".
                Default: auto-detect from skeleton file.

  -i/--input    Path to the OLCAO skeleton input file.
                Default: "olcao.skl".

  --help        Print usage information and exit.
"""

import argparse as ap
import os
import sys
from datetime import datetime

from structure_control import StructureControl


# ----------------------------------------------------------------
#  Potential sub-type string mapping.
#
#  This dictionary maps user-supplied sub-type abbreviations
#  to the directory-name suffixes used in the VASP potential
#  database.  Some sub-types gain an underscore prefix while
#  others (numeric fractions, etc.) do not.
# ----------------------------------------------------------------
SUB_POT_MAP = {
    "s": "_s",
    "h": "_h",
    "sv": "_sv",
    "pv": "_pv",
    "d": "_d",
    "soft": "_soft",
    "200ev": "_200ev",
    ".75": ".75",
    "1.25": "1.25",
    "1.5": "1.5",
    ".5": ".5",
    ".33": ".33",
    ".66": ".66",
    "3": "_3",
    "2": "_2",
    "new": "_new",
}


# ----------------------------------------------------------------
#  ScriptSettings
# ----------------------------------------------------------------

class ScriptSettings():
    """User settings that control the program.

    Variable values are pulled from a resource control file
    (skl2vasprc.py) and reconciled with command line parameters.
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
        from skl2vasprc import parameters_and_defaults
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

        # The type of VASP job the input files are for.
        self.job_type = default_rc["job_type"]

        # The base type of potential to use.
        self.pot_type = default_rc["pot_type"]

        # The sub-type of the potential to use.
        self.sub_pot_type = default_rc["sub_pot_type"]

        # Flag to include the GW sub-sub-type.
        self.gw = default_rc["gw"]

        # Flag to use a Gamma-centered k-point.
        self.gamma_kpoint = default_rc["gamma_kpoint"]

        # Bravais lattice cell type override.
        self.cell_type = default_rc["cell_type"]

        # Path to the OLCAO skeleton input file.
        self.input_file = default_rc["input_file"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv."""

        prog_name = "skl2vasp.py"

        description_text = """\
PROGRAM:  skl2vasp
PURPOSE:  Take an OLCAO skeleton file and produce VASP
          input files (currently POSCAR only).  Options
          for potential type, job type, and k-points are
          parsed and stored for future expansion.
"""

        epilog_text = """\
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults: ./skl2vasprc.py or $OLCAO_RC/skl2vasprc.py.
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

        # The -j option defines the type of VASP
        #   calculation.  Options: "relaxfull",
        #   "relaxion1", "relaxion2", "relaxvol", "cij",
        #   "phonon", "gpt".
        parser.add_argument(
            '-j', '--job',
            dest='job_type',
            type=str,
            default=self.job_type,
            help=(
                'VASP job type. '
                f'Default: {self.job_type}'
            ),
        )

        # The -p option defines the base potential type.
        #   Options: "potLDA", "potGGA", "potpawLDA",
        #   "potpawGGA", "potpawPBE", "potpawLDA5x",
        #   "potpawPBE5x".
        parser.add_argument(
            '-p', '--pot',
            dest='pot_type',
            type=str,
            default=self.pot_type,
            help=(
                'Base potential type. '
                f'Default: {self.pot_type}'
            ),
        )

        # The -s option defines the potential sub-type.
        #   Options: "s" (soft), "h" (hard), "pv"
        #   (semi-core p), "sv" (semi-core s+p), "d"
        #   (semi-core d), "soft", "200ev", etc.
        #   In the event that a sub type is not available
        #   for a particular element, the standard type
        #   will be used.
        parser.add_argument(
            '-s', '--subpot',
            dest='sub_pot_type',
            type=str,
            default=self.sub_pot_type,
            help=(
                'Potential sub-type. '
                f'Default: "{self.sub_pot_type}"'
            ),
        )

        # The --gw option appends the GW sub-sub-type
        #   to the potential sub-type string.
        parser.add_argument(
            '--gw',
            dest='gw',
            action='store_true',
            default=self.gw,
            help=(
                'Append GW sub-type. '
                f'Default: {self.gw}'
            ),
        )

        # The --gamma option requests that the calculation
        #   use a Gamma-centered k-point instead of a
        #   general Monkhorst-Pack mesh.
        parser.add_argument(
            '--gamma',
            dest='gamma_kpoint',
            action='store_true',
            default=self.gamma_kpoint,
            help=(
                'Use Gamma k-point. '
                f'Default: {self.gamma_kpoint}'
            ),
        )

        # The -t option allows the user to override the
        #   cell type that would otherwise be determined
        #   from the space group label.  This may be
        #   useful when a calculation is being done on a
        #   crystal system (say hexagonal) but in which
        #   some symmetry has been broken (e.g. now
        #   triclinic) and the user still wants to proceed
        #   as if the cell was hexagonal.
        #   Options: "tric", "mono", "ortho", "tet",
        #   "trig", "hex", "sc", "bcc", "fcc".
        parser.add_argument(
            '-t', '--celltype',
            dest='cell_type',
            type=str,
            default=self.cell_type,
            help=(
                'Override Bravais lattice cell type. '
                f'Default: "{self.cell_type}"'
            ),
        )

        # The -i option specifies the skeleton input file.
        parser.add_argument(
            '-i', '--input',
            dest='input_file',
            type=str,
            default=self.input_file,
            help=(
                'OLCAO skeleton input file. '
                f'Default: {self.input_file}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command line values with rc defaults."""
        self.job_type = args.job_type
        self.pot_type = args.pot_type
        self.sub_pot_type = args.sub_pot_type
        self.gw = args.gw
        self.gamma_kpoint = args.gamma_kpoint
        self.cell_type = args.cell_type
        self.input_file = args.input_file

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
#  Helper functions
# ----------------------------------------------------------------

def resolve_sub_pot_type(sub_pot_type, gw):
    """Convert the user-supplied potential sub-type
    abbreviation into the directory-name suffix used in the
    VASP potential database.

    Some sub-types gain an underscore prefix (e.g. "s" becomes
    "_s") while others (numeric fractions like ".75") are used
    as-is.  If the GW flag is set, "_GW" is appended to the
    result.

    Args:
        sub_pot_type (str): User-supplied sub-type
            abbreviation (e.g. "s", "pv", ".75").
        gw (bool): If True, append "_GW" to the result.

    Returns:
        str: The resolved sub-type string for use in VASP
            potential directory lookups.
    """
    # Map the abbreviation to the directory suffix.
    resolved = SUB_POT_MAP.get(
        sub_pot_type, sub_pot_type
    )

    # Append the GW sub-sub-type if requested.
    if gw:
        resolved = resolved + "_GW"

    return resolved


def get_cell_type_from_sg(sg_name):
    """Determine the crystal cell type from a space group
    designation by looking it up in the OLCAO space group
    database.

    The approach:
      1. Open the space group definition file at
         $OLCAO_DATA/spaceDB/<sg_name>.
      2. Read the space group number from the second line.
      3. Map the number to a cell type using the standard
         crystallographic ranges.

    Args:
        sg_name (str): Space group designation as it appears
            in the skeleton file (e.g. "1_a", "225_a").

    Returns:
        str: Cell type string.  One of "tric", "mono",
            "ortho", "tet", "trig", "hex", "cubic".
    """
    data_dir = os.getenv('OLCAO_DATA', '')
    sg_path = os.path.join(data_dir, 'spaceDB', sg_name)

    with open(sg_path, 'r') as f:
        f.readline()                        # skip line 1
        sg_number = int(f.readline().split()[0])

    # Determine the cell type from the space group number.
    if sg_number <= 2:
        return "tric"
    elif sg_number <= 15:
        return "mono"
    elif sg_number <= 74:
        return "ortho"
    elif sg_number <= 142:
        return "tet"
    elif sg_number <= 167:
        return "trig"
    elif sg_number <= 194:
        return "hex"
    elif sg_number <= 230:
        return "cubic"
    else:
        sys.exit(
            "Error: space group number > 230."
        )


def pre_sort_skl(input_file, output_file, cell_type):
    """Create a copy of the OLCAO skeleton file with atoms
    sorted alphabetically by element name.

    VASP requires that atoms of the same element be grouped
    together in the POSCAR file.  This function creates a
    sorted copy of the skeleton file so that when
    StructureControl reads it, the atoms are already in
    element-sorted order.

    The approach:
      1. Read and immediately output the header information
         (everything up to and including the frac/cart line).
      2. Read the atom list, sort it case-insensitively, and
         write the sorted lines.
      3. Read and immediately output the tail information,
         extracting the cell type from the space group
         designation if the cell type was not overridden on
         the command line.

    Args:
        input_file (str): Path to the OLCAO skeleton file.
        output_file (str): Path to write the sorted copy.
        cell_type (str): User-specified cell type override.
            If empty, the cell type will be determined from
            the space group designation in the skeleton file.

    Returns:
        str: The resolved cell type string.
    """
    resolved_cell_type = cell_type

    with open(input_file, 'r') as fin, \
         open(output_file, 'w') as fout:

        # Read and regurgitate the header information.
        #   Once we reach the part with the atomic
        #   coordinates (a line containing "frac" or
        #   "cart"), we shift to collecting the atom
        #   lines for sorting.  We also grab the number
        #   of atom lines in the file from this line.
        num_atoms = 0
        for line in fin:
            fout.write(line)
            if 'frac' in line or 'cart' in line:
                values = line.split()
                num_atoms = int(values[1])
                break

        # Collect the atom lines for sorting.
        atom_lines = []
        for _ in range(num_atoms):
            atom_lines.append(fin.readline())

        # Sort the atom lines case-insensitively so that
        #   atoms of the same element are grouped together.
        sorted_lines = sorted(
            atom_lines, key=str.casefold
        )

        # Print the sorted atom lines to the output file.
        for line in sorted_lines:
            fout.write(line)

        # Read and regurgitate the tail information,
        #   keeping track of the space group request so
        #   we can determine the cell type if it was not
        #   specified on the command line.
        for line in fin:
            fout.write(line)

            # Get the cell type from the space group
            #   designation if it was not defined on the
            #   command line.
            if ('space' in line
                    and not resolved_cell_type):
                values = line.split()
                if len(values) >= 2:
                    resolved_cell_type = (
                        get_cell_type_from_sg(values[1])
                    )

    return resolved_cell_type


# ----------------------------------------------------------------
#  Main
# ----------------------------------------------------------------

def main():
    """Main entry point for skl2vasp.

    Reads an OLCAO skeleton file, pre-sorts atoms by element,
    determines the cell type, and writes a VASP POSCAR file.

    NOTE: The original Perl skl2vasp also generated POTCAR,
    KPOINTS, and INCAR files.  Those are not yet ported in
    structure_control.py.  The potential type, sub-type, GW,
    Gamma k-point, and job type options are parsed and stored
    for future use when those files are implemented.
    """

    # Get script settings from a combination of the resource
    #   control file and parameters given by the user on the
    #   command line.
    settings = ScriptSettings()

    # Resolve the potential sub-type string for directory
    #   lookup (will be used when POTCAR generation is
    #   implemented).
    settings.sub_pot_type = resolve_sub_pot_type(
        settings.sub_pot_type, settings.gw
    )

    # Pre-sort the element list in the skeleton file.
    #   VASP requires atoms of the same element to be
    #   grouped together.  This creates a sorted copy of
    #   the skeleton file and determines the cell type
    #   from the space group if not overridden.
    sorted_skl = settings.input_file + ".sorted"
    settings.cell_type = pre_sort_skl(
        settings.input_file,
        sorted_skl,
        settings.cell_type,
    )

    # Read the sorted skeleton file.
    sc = StructureControl()
    sc.read_input_file(
        sorted_skl, use_file_species=False
    )

    # Print the VASP POSCAR file.
    sc.print_vasp(filename='POSCAR')

    # Clean up the temporary sorted skeleton file.
    if os.path.exists(sorted_skl):
        os.remove(sorted_skl)


if __name__ == '__main__':
    # Everything before this point was a subroutine
    #   definition or a request to import information from
    #   external modules.  Only now do we actually start
    #   running the program.  The purpose of this is to
    #   allow another python program to import *this* script
    #   and call its functions internally.
    main()
