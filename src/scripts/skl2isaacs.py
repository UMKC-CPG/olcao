#!/usr/bin/env python3

"""skl2isaacs.py -- Convert an OLCAO skeleton file to ISAACS
input files.

PROGRAM:  skl2isaacs
PURPOSE:  This program will convert an OLCAO skeleton .skl
          file into a set of ISAACS input files: an ISAACS
          XML project file (.ipf) and a Chem3D coordinate
          file (.chem3d).

          ISAACS (Interactive Structure Analysis of Amorphous
          and Crystalline Structures) is a tool for analysing
          structural properties of both crystalline and
          disordered materials.  It uses a two-file scheme:
          the .ipf file contains chemistry, box, PBC, and
          cutoff information in XML format, while the .c3d
          file contains Cartesian atom coordinates.

          The conversion uses StructureControl to read the
          skeleton file (which handles coordinate parsing,
          element identification, and cell parameter
          computation) and then delegates the ISAACS output
          to StructureControl.print_isaacs().

USAGE:    skl2isaacs.py [-i INPUT] [-o OUTPUTROOT] [--help]

OPTIONS:

  -i/--input    Name of the input OLCAO skeleton file.
                If this option is not given, the default
                value is "olcao.skl".

  -o/--output   Root name for the output ISAACS files.
                Two files are produced:
                  <root>.ipf    — ISAACS XML project file
                  <root>.chem3d — Chem3D coordinate file
                If this option is not given, the default
                root is "olcao".

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
    (skl2isaacsrc.py) and reconciled with command line
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
        from skl2isaacsrc import parameters_and_defaults
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

        # Name of the input OLCAO skeleton file.
        self.input_file = default_rc["input_file"]

        # Root name for the ISAACS output files.
        self.output_root = default_rc["output_root"]

    def parse_command_line(self):
        """Build the argument parser and parse sys.argv."""

        prog_name = "skl2isaacs.py"

        description_text = """\
PROGRAM:  skl2isaacs
PURPOSE:  Convert an OLCAO skeleton (.skl) file into a set
          of ISAACS input files: an ISAACS XML project file
          (.ipf) and a Chem3D coordinate file (.chem3d).
"""

        epilog_text = """\
Please contact Paul Rulis (rulisp@umkc.edu) with questions.
Defaults: ./skl2isaacsrc.py or $OLCAO_RC/skl2isaacsrc.py.
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

        # The -i option specifies the input skeleton file.
        parser.add_argument(
            '-i', '--input',
            dest='input_file',
            type=str,
            default=self.input_file,
            help=(
                'Input OLCAO skeleton file. '
                f'Default: {self.input_file}'
            ),
        )

        # The -o option specifies the root name for the
        #   ISAACS output files.  Two files are produced:
        #   <root>.ipf (ISAACS XML project file) and
        #   <root>.chem3d (Chem3D coordinate file).
        parser.add_argument(
            '-o', '--output',
            dest='output_root',
            type=str,
            default=self.output_root,
            help=(
                'Root name for ISAACS output files '
                '(<root>.ipf and <root>.chem3d). '
                f'Default: {self.output_root}'
            ),
        )

    def reconcile(self, args):
        """Reconcile command line values with rc defaults."""
        self.input_file = args.input_file
        self.output_root = args.output_root

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
    """Main entry point for skl2isaacs.

    Reads an OLCAO skeleton file via StructureControl and
    writes the ISAACS XML project file (.ipf) and the
    Chem3D coordinate file (.chem3d).

    The workflow mirrors the original Perl script:
      1. Parse command line and set defaults.
      2. Read the skeleton input file via StructureControl.
      3. Write the ISAACS output files.
    """

    # Get script settings from a combination of the
    #   resource control file and parameters given by the
    #   user on the command line.
    settings = ScriptSettings()

    # Read the skeleton input file and store/calculate
    #   all important information via StructureControl.
    sc = StructureControl()
    sc.read_input_file(
        settings.input_file,
        use_file_species=False,
    )

    # Print the ISAACS input files.
    #   StructureControl.print_isaacs() takes a base
    #   filename and produces two files:
    #   <base>.ipf and <base>.c3d.
    sc.print_isaacs(filename=settings.output_root)


if __name__ == '__main__':
    # Everything before this point was a subroutine
    #   definition or a request to import information from
    #   external modules.  Only now do we actually start
    #   running the program.  The purpose of this is to
    #   allow another python program to import *this* script
    #   and call its functions internally.
    main()
