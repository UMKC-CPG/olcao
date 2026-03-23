#!/usr/bin/env python3

"""pdb2skl.py -- Convert a PDB file to an OLCAO skeleton file.

PROGRAM:  pdb2skl
PURPOSE:  This program will convert a Protein Data Bank .pdb
          file into an olcao.skl input file with the atoms all
          rearranged such that they are sorted according to
          their elements.  This will also produce a mapping
          between the pdb file atoms and the skl atoms.

USAGE:    pdb2skl.py [-i INPUT] [-o OUTPUT] [-m MAPFILE]
                     [-b BUFFER] [--pdbtypes] [--help]

OPTIONS:

  -i/--input    Path to the PDB input file.  If this option is
                not given, then the default value for the input
                file is model.pdb.
  -o/--output   Path to the OLCAO skeleton output file.  If
                this option is not given, then the default value
                for the output file is olcao.skl.
  -m/--map      Path to the mapping file that records the
                correspondence between PDB atoms and SKL atoms.
  -b/--buffer   Size of the buffer (Angstroms) to add around
                the molecule if the PDB file does not define a
                cell size.  Default = 0 A.
  --pdbtypes    Use the atom types defined in the PDB file as
                species tags in the skeleton file.  The default
                is to have no types assigned (all atoms of the
                same element share a single species).
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
    (pdb2sklrc.py) and reconciled with command line parameters.
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
        from pdb2sklrc import parameters_and_defaults
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

        # Path to the PDB input file.
        self.input_file = default_rc["input_file"]

        # Path to the OLCAO skeleton output file.
        self.output_file = default_rc["output_file"]

        # Path to the PDB-to-SKL atom map file.
        self.map_file = default_rc["map_file"]

        # Buffer size (Angstroms) around the molecule if
        #   the PDB file does not define a cell.
        self.buffer = default_rc["buffer"]

        # Flag to use PDB atom types as species tags.
        self.pdb_types = default_rc["pdb_types"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv."""

        prog_name = "pdb2skl.py"

        description_text = """\
PROGRAM:  pdb2skl
PURPOSE:  Convert a Protein Data Bank .pdb file into an
          OLCAO skeleton (.skl) input file with the atoms
          sorted by element.  Also produces a mapping file
          between PDB and SKL atom indices.
"""

        epilog_text = """\
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults: ./pdb2sklrc.py or $OLCAO_RC/pdb2sklrc.py.
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

        # The -i option specifies the PDB input file.
        parser.add_argument(
            '-i', '--input',
            dest='input_file',
            type=str,
            default=self.input_file,
            help=(
                'PDB input file. '
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
                'PDB-to-SKL atom map file. '
                f'Default: {self.map_file}'
            ),
        )

        # The -b option defines the buffer size
        #   (Angstroms) if the PDB file does not define
        #   a cell size.
        parser.add_argument(
            '-b', '--buffer',
            dest='buffer',
            type=float,
            default=self.buffer,
            help=(
                'Buffer size in Angstroms around the '
                'molecule when no cell is defined. '
                f'Default: {self.buffer}'
            ),
        )

        # The --pdbtypes option requests that the skeleton
        #   file use the atom types defined in the PDB file
        #   as species tags.  The default is to have no
        #   types assigned.
        parser.add_argument(
            '--pdbtypes',
            dest='pdb_types',
            action='store_true',
            default=self.pdb_types,
            help=(
                'Use PDB atom types as species tags. '
                f'Default: {self.pdb_types}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command line values with rc defaults."""
        self.input_file = args.input_file
        self.output_file = args.output_file
        self.map_file = args.map_file
        self.buffer = args.buffer
        self.pdb_types = args.pdb_types

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
    """Main entry point for pdb2skl.

    Reads a PDB file, converts it to an OLCAO skeleton file,
    and writes the atom index mapping between the two formats.
    """

    # Get script settings from a combination of the resource
    #   control file and parameters given by the user on the
    #   command line.
    settings = ScriptSettings()

    # Create the StructureControl object and read the PDB
    #   input file.  The use_file_species flag controls
    #   whether PDB atom types are preserved as species tags
    #   in the skeleton output.
    sc = StructureControl()
    sc.read_input_file(
        settings.input_file,
        use_file_species=settings.pdb_types,
    )

    # Print the OLCAO skeleton (.skl) output file.
    sc.print_olcao(
        filename=settings.output_file,
        title="Generated from PDB file.",
        style="cartesian",
    )

    # Print the data file that maps the PDB nature of each
    #   atom to the SKL atoms.
    sc.print_olcao_map(filename=settings.map_file)


if __name__ == '__main__':
    # Everything before this point was a subroutine
    #   definition or a request to import information from
    #   external modules.  Only now do we actually start
    #   running the program.  The purpose of this is to
    #   allow another python program to import *this* script
    #   and call its functions internally.
    main()
