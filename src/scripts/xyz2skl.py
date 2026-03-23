#!/usr/bin/env python3

"""xyz2skl.py -- Convert an XYZ coordinate file to an OLCAO
skeleton file.

PROGRAM:  xyz2skl
PURPOSE:  This program will convert a file containing a
          structure in pure XYZ direct space coordinates into
          an OLCAO skeleton file with the atoms all rearranged
          such that they are sorted according to their
          elements.  This will also produce a mapping between
          the XYZ file and the SKL atoms.

          The conversion uses StructureControl to read the XYZ
          file (which handles the coordinate parsing, element
          identification, and cell parameter computation) and
          then writes the OLCAO skeleton format and a map file
          that records the correspondence between atom indices
          in the original XYZ file and the output skeleton.

USAGE:    xyz2skl.py [-i INPUT] [-o OUTPUT] [-m MAPFILE]
                     [--xyz-types] [--help]

OPTIONS:

  -i/--input    Name of the input file containing the
                structure in pure XYZ direct space
                coordinates.  If this option is not given,
                the default value is "model.xyz".

  -o/--output   Name of the output OLCAO skeleton file.
                If this option is not given, the default
                value is "olcao.skl".

  -m/--map      Name of the mapping file that stores the
                correspondence between the XYZ file atoms
                and the SKL atoms.  If this option is not
                given, the default value is "sklXYZ.map".

  --xyz-types   When given, the skeleton file will use the
                atom types as defined in the XYZ file.  The
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
    (xyz2sklrc.py) and reconciled with command line
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
        from xyz2sklrc import parameters_and_defaults
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

        # Name of the XYZ input file.
        self.input_file = default_rc["input_file"]

        # Name of the output OLCAO skeleton file.
        self.output_file = default_rc["output_file"]

        # Name of the XYZ-to-SKL mapping file.
        self.map_file = default_rc["map_file"]

        # Whether to use XYZ file types in the skeleton.
        self.xyz_types = default_rc["xyz_types"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv."""

        prog_name = "xyz2skl.py"

        description_text = """\
PROGRAM:  xyz2skl
PURPOSE:  Convert a file containing a structure in pure XYZ
          direct space coordinates into an OLCAO skeleton
          file with atoms sorted by element.  Also produces
          a mapping between the XYZ and SKL atom indices.
"""

        epilog_text = """\
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults: ./xyz2sklrc.py or $OLCAO_RC/xyz2sklrc.py.
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

        # The -i option specifies the XYZ input file.
        #   This is the file containing the structure in
        #   pure XYZ direct space coordinates.
        parser.add_argument(
            '-i', '--input',
            dest='input_file',
            type=str,
            default=self.input_file,
            help=(
                'XYZ input file. '
                f'Default: {self.input_file}'
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

        # The -m option specifies the mapping file that
        #   records the correspondence between the XYZ file
        #   atoms and the SKL atoms.
        parser.add_argument(
            '-m', '--map',
            dest='map_file',
            type=str,
            default=self.map_file,
            help=(
                'XYZ-to-SKL atom mapping file. '
                f'Default: {self.map_file}'
            ),
        )

        # The --xyz-types option requests that the skeleton
        #   file use the atom types as defined in the XYZ
        #   file.  The default is to have no types assigned.
        parser.add_argument(
            '--xyz-types',
            dest='xyz_types',
            action='store_true',
            default=self.xyz_types,
            help=(
                'Use atom types from XYZ file. '
                f'Default: {self.xyz_types}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command line values with rc defaults."""
        self.input_file = args.input_file
        self.output_file = args.output_file
        self.map_file = args.map_file
        self.xyz_types = args.xyz_types

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
    """Main entry point for xyz2skl.

    Reads an XYZ coordinate file via StructureControl,
    prepends a provenance note to the title, writes the
    OLCAO skeleton file, and writes the XYZ-to-SKL atom
    mapping file.

    The workflow mirrors the original Perl script:
      1. Parse command line and set defaults.
      2. Read the XYZ input file via StructureControl.
      3. Prepend "Generated from XYZ file:" to the title.
      4. Write the OLCAO skeleton file (Cartesian coords).
      5. Write the atom index mapping file.
    """

    # Get script settings from a combination of the
    #   resource control file and parameters given by the
    #   user on the command line.
    settings = ScriptSettings()

    # Read the XYZ input file and store/calculate all the
    #   important information via StructureControl.
    sc = StructureControl()
    sc.read_input_file(
        settings.input_file,
        use_file_species=settings.xyz_types,
    )

    # Get the title for the system and prepend a
    #   provenance note indicating it came from an XYZ
    #   file.
    title = ""
    if sc.system_title and len(sc.system_title) > 1:
        title = sc.system_title[1].strip()
    title = f"Generated from XYZ file: {title}"

    # Print the OLCAO skeleton file using Cartesian
    #   coordinates.  The XYZ file contains Cartesian
    #   coordinates by definition, so "cartesian" is the
    #   natural coordinate style for the output.
    sc.print_olcao(
        filename=settings.output_file,
        title=title,
        style="cartesian",
    )

    # Print the data file that maps the XYZ nature of
    #   each atom to the SKL atoms.  This mapping is
    #   useful for post-processing when the user needs to
    #   correlate results back to the original XYZ atom
    #   ordering.
    sc.print_olcao_map(filename=settings.map_file)


if __name__ == '__main__':
    # Everything before this point was a subroutine
    #   definition or a request to import information from
    #   external modules.  Only now do we actually start
    #   running the program.  The purpose of this is to
    #   allow another python program to import *this* script
    #   and call its functions internally.
    main()
