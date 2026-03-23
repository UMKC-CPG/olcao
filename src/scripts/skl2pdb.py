#!/usr/bin/env python3

"""skl2pdb.py -- Convert an OLCAO skeleton file to PDB format.

PROGRAM:  skl2pdb
PURPOSE:  This program will convert an OLCAO skeleton .skl file
          into a Protein Data Bank .pdb input file with the
          atoms all rearranged such that they are sorted
          according to their elements.  This will also produce
          a mapping between the pdb file atoms and the skl
          atoms.

USAGE:    skl2pdb.py [-i INPUT] [-o OUTPUT] [-m MAPFILE]
                     [--skltypes] [--help]

OPTIONS:

  -i/--input    Path to the OLCAO skeleton input file.  If
                this option is not given, then the default
                value for the input file is olcao.skl.
  -o/--output   Path to the PDB output file.  If this option
                is not given, then the default value for the
                output file is olcao.pdb.
  -m/--map      Path to the mapping file that records the
                correspondence between SKL and PDB atom
                indices.
  --skltypes    Use the species types defined in the skeleton
                file as atom types in the PDB output.  The
                default is to have no types assigned.
  --help        Print usage information and exit.
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
    (skl2pdbrc.py) and reconciled with command line parameters.
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
        from skl2pdbrc import parameters_and_defaults
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

        # Path to the OLCAO skeleton input file.
        self.input_file = default_rc["input_file"]

        # Path to the PDB output file.
        self.output_file = default_rc["output_file"]

        # Path to the SKL-to-PDB atom map file.
        self.map_file = default_rc["map_file"]

        # Flag to use skeleton species types in PDB output.
        self.skl_types = default_rc["skl_types"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv."""

        prog_name = "skl2pdb.py"

        description_text = """\
PROGRAM:  skl2pdb
PURPOSE:  Convert an OLCAO skeleton (.skl) file into a
          Protein Data Bank (.pdb) file with the atoms
          sorted by element.
"""

        epilog_text = """\
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults: ./skl2pdbrc.py or $OLCAO_RC/skl2pdbrc.py.
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

        # The -o option specifies the PDB output file.
        parser.add_argument(
            '-o', '--output',
            dest='output_file',
            type=str,
            default=self.output_file,
            help=(
                'PDB output file. '
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
                'SKL-to-PDB atom map file. '
                f'Default: {self.map_file}'
            ),
        )

        # The --skltypes option requests that the PDB file
        #   use the species types defined in the skeleton
        #   file.  The default is to have no types assigned.
        parser.add_argument(
            '--skltypes',
            dest='skl_types',
            action='store_true',
            default=self.skl_types,
            help=(
                'Use skeleton species types in PDB. '
                f'Default: {self.skl_types}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command line values with rc defaults."""
        self.input_file = args.input_file
        self.output_file = args.output_file
        self.map_file = args.map_file
        self.skl_types = args.skl_types

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
    """Main entry point for skl2pdb.

    Reads an OLCAO skeleton file and converts it to PDB format.
    """

    # Get script settings from a combination of the resource
    #   control file and parameters given by the user on the
    #   command line.
    settings = ScriptSettings()

    # Create the StructureControl object and read the OLCAO
    #   skeleton input file.  The use_file_species flag
    #   controls whether the skeleton species types are
    #   preserved in the PDB output.
    sc = StructureControl()
    sc.read_input_file(
        settings.input_file,
        use_file_species=settings.skl_types,
    )

    # Print the PDB output file.
    sc.print_pdb(filename=settings.output_file)


if __name__ == '__main__':
    # Everything before this point was a subroutine
    #   definition or a request to import information from
    #   external modules.  Only now do we actually start
    #   running the program.  The purpose of this is to
    #   allow another python program to import *this* script
    #   and call its functions internally.
    main()
