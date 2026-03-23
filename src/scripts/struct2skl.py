#!/usr/bin/env python3

"""struct2skl.py -- Convert a structure.dat file to OLCAO skeleton.

PROGRAM:  struct2skl
PURPOSE:  This program will convert a file containing a structure
          in the format of an OLCAO structure.dat file into an
          OLCAO skeleton file with the atoms all rearranged such
          that they are sorted according to their elements.  This
          will also produce a mapping between the structure file
          and the skl atoms.

USAGE:    struct2skl.py [-i INPUT] [-o OUTPUT] [-m MAPFILE]
                        [--structtypes] [--makemol] [--help]

OPTIONS:

  -i/--input      Path to the structure.dat input file.  If
                  this option is not given, then the default
                  value for the input file is structure.dat.
  -o/--output     Path to the OLCAO skeleton output file.  If
                  this option is not given, then the default
                  value for the output file is olcao.skl.
  -m/--map        Path to the mapping file that records the
                  correspondence between structure.dat atoms
                  and SKL atoms.
  --structtypes   Use the species types from the structure.dat
                  file in the skeleton output.  The default is
                  to have no types assigned (i.e., all atoms
                  have type # = 1).
  --makemol       Export the structure as a cluster (molecule)
                  rather than inside a periodic cell.
                  NOTE: This option is documented in the
                  original Perl script's help text but was
                  not implemented in the Perl code.  It is
                  included here for future use.
  --help          Print usage information and exit.
"""

import argparse as ap
import os
import sys
from datetime import datetime

from structure_control import StructureControl


# ----------------------------------------------------------------
#  ScriptSettings
# ----------------------------------------------------------------

class ScriptSettings():
    """User settings that control the program.

    Variable values are pulled from a resource control file
    (struct2sklrc.py) and reconciled with command line
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
        from struct2sklrc import parameters_and_defaults
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

        # Path to the structure.dat input file.
        self.input_file = default_rc["input_file"]

        # Path to the OLCAO skeleton output file.
        self.output_file = default_rc["output_file"]

        # Path to the struct-to-SKL atom map file.
        self.map_file = default_rc["map_file"]

        # Flag to use structure.dat species types.
        self.struct_types = default_rc["struct_types"]

        # Flag to export as cluster (molecule).
        self.make_mol = default_rc["make_mol"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv."""

        prog_name = "struct2skl.py"

        description_text = """\
PROGRAM:  struct2skl
PURPOSE:  Convert an OLCAO structure.dat file into an
          OLCAO skeleton (.skl) file with atoms sorted
          by element.  Also produces a mapping file
          between structure.dat and SKL atom indices.
"""

        epilog_text = """\
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults: ./struct2sklrc.py or $OLCAO_RC/struct2sklrc.py.
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

        # The -i option specifies the structure.dat input.
        parser.add_argument(
            '-i', '--input',
            dest='input_file',
            type=str,
            default=self.input_file,
            help=(
                'structure.dat input file. '
                f'Default: {self.input_file}'
            ),
        )

        # The -o option specifies the skeleton output file.
        parser.add_argument(
            '-o', '--output',
            dest='output_file',
            type=str,
            default=self.output_file,
            help=(
                'OLCAO skeleton output file. '
                f'Default: {self.output_file}'
            ),
        )

        # The -m option specifies the atom map file.
        parser.add_argument(
            '-m', '--map',
            dest='map_file',
            type=str,
            default=self.map_file,
            help=(
                'Struct-to-SKL atom map file. '
                f'Default: {self.map_file}'
            ),
        )

        # The --structtypes option requests that the
        #   skeleton file use the species types defined
        #   in the structure.dat file.  The default is
        #   to have no types assigned (all atoms of the
        #   same element share a single species, i.e.
        #   type # = 1).
        parser.add_argument(
            '--structtypes',
            dest='struct_types',
            action='store_true',
            default=self.struct_types,
            help=(
                'Use structure.dat species types. '
                f'Default: {self.struct_types}'
            ),
        )

        # The --makemol option requests that the given
        #   structure file be exported such that the
        #   atoms are not inside of a periodic cell.
        #   Instead, they will be in a cluster
        #   configuration.
        #   NOTE: This option was documented in the
        #   original Perl script's help text but was not
        #   implemented.  It is included here as a
        #   placeholder for future development.
        parser.add_argument(
            '--makemol',
            dest='make_mol',
            action='store_true',
            default=self.make_mol,
            help=(
                'Export as cluster (not implemented). '
                f'Default: {self.make_mol}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command line values with rc defaults."""
        self.input_file = args.input_file
        self.output_file = args.output_file
        self.map_file = args.map_file
        self.struct_types = args.struct_types
        self.make_mol = args.make_mol

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
#  Main
# ----------------------------------------------------------------

def main():
    """Main entry point for struct2skl.

    Reads an OLCAO structure.dat file, converts it to an
    OLCAO skeleton (.skl) file, and writes the atom index
    mapping between the two formats.
    """

    # Get script settings from a combination of the resource
    #   control file and parameters given by the user on the
    #   command line.
    settings = ScriptSettings()

    # Create the StructureControl object and read the
    #   structure.dat style input file.  The
    #   use_file_species flag controls whether the species
    #   types from the structure.dat file are preserved in
    #   the skeleton output.
    sc = StructureControl()
    sc.read_input_file(
        settings.input_file,
        use_file_species=settings.struct_types,
    )

    # Get the title for the system and prepend a note about
    #   the source file.
    title = (
        f"Generated from {settings.input_file} "
        f"structure file:\n{sc.title}\n"
    )

    # Print the OLCAO skeleton (.skl) output file.
    sc.print_olcao(
        filename=settings.output_file,
        title=title,
        style="cartesian",
    )

    # Print the data file that maps the structure.dat
    #   nature of each atom to the SKL atoms.
    sc.print_olcao_map(filename=settings.map_file)


if __name__ == '__main__':
    # Everything before this point was a subroutine
    #   definition or a request to import information from
    #   external modules.  Only now do we actually start
    #   running the program.  The purpose of this is to
    #   allow another python program to import *this* script
    #   and call its functions internally.
    main()
