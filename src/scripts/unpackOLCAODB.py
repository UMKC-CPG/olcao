#!/usr/bin/env python3

"""PROGRAM: unpackOLCAODB.py

PURPOSE: This program unpacks the database files that come with the
   OLCAO package. The databases include atomic basis sets, atomic
   potentials, precursor data, space group definitions, and
   symmetry band (SYBD) k-point path data. These databases are
   stored as compressed tarballs (.tgz) in the directory pointed
   to by the $OLCAO_DATA environment variable and must be
   unpacked before the OLCAO package can be used.

   After unpacking the space group database (spaceDB.tgz), this
   script also runs the "remake" script inside the spaceDB/
   directory to regenerate the individual space group files from
   the packed representation.

USAGE: unpackOLCAODB.py [-h]

   This script takes no arguments beyond the optional -h for
   help. It relies entirely on the $OLCAO_DATA environment
   variable to locate the database directory.

DATABASES UNPACKED:
   atomicBDB.tgz   - Atomic basis set database
   atomicPDB.tgz   - Atomic potential database
   precursorDB.tgz  - Precursor data for atomic SCF
   spaceDB.tgz     - Space group symmetry operations
   sybdDB.tgz      - Symmetry band k-point path definitions
"""

import argparse as ap
import os
import subprocess
import sys
from datetime import datetime


# Define the main class that holds script data structures and
#   settings.
class ScriptSettings():
    """The instance variables of this object are the user settings
       that control the program. For this script there are no
       user-configurable parameters; the class exists to maintain
       a consistent structure with other OLCAO Python scripts and
       to handle command line parsing (help text) and command
       logging."""


    def __init__(self):
        """Initialize settings by parsing the command line and
        recording the invocation."""

        # Parse the command line (help text only for this
        #   script).
        self.parse_command_line()

        # Record the command line parameters that were used.
        self.recordCLP()


    def parse_command_line(self):
        """Set up and run the argument parser. This script has no
        arguments beyond -h/--help, but the parser provides
        consistent help text and usage information."""

        # Create the parser tool.
        prog_name = "unpackOLCAODB.py"

        description_text = """
Version: 1.0 (March 2026)
Requires: $OLCAO_DATA environment variable set to the OLCAO
          data directory containing the .tgz database files.

Unpack the OLCAO database files (.tgz archives) located in the
$OLCAO_DATA directory. This must be done once after installing
or updating the OLCAO package before any calculations can be
run. The databases provide atomic basis sets, potentials,
precursor data, space group definitions, and symmetry band
k-point path information.

After unpacking, the space group database is rebuilt by running
the 'remake' script inside the spaceDB/ directory.
"""

        epilog_text = """
Please contact rulisp@umkc.edu regarding questions.
"""

        parser = ap.ArgumentParser(
                prog=prog_name,
                formatter_class=ap.RawDescriptionHelpFormatter,
                description=description_text,
                epilog=epilog_text)

        # No arguments to add for this script.

        # Parse the arguments (validates that no unexpected
        #   arguments were given).
        parser.parse_args()


    def recordCLP(self):
        """Record the command line invocation to the 'command'
        file for reproducibility tracking."""
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime(
                    "%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write("Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


def unpack_databases():
    """Unpack all OLCAO database tarballs in the $OLCAO_DATA
    directory and rebuild the space group database.

    The function changes to the $OLCAO_DATA directory, extracts
    each .tgz archive, then enters the spaceDB/ subdirectory to
    run the 'remake' script that regenerates individual space
    group files.

    Each tarball is extracted with 'tar -xzf'. If any extraction
    or the remake step fails, the program exits with an error
    message indicating which step failed."""

    # Determine the OLCAO data directory from the environment.
    olcao_data = os.getenv('OLCAO_DATA')
    if not olcao_data:
        sys.exit(
            "Error: $OLCAO_DATA is not set. "
            "Please set it to the OLCAO data directory."
        )

    # Go to the $OLCAO_DATA directory (which is where the
    #   databases are).
    try:
        os.chdir(olcao_data)
    except OSError as e:
        sys.exit(
            f"Error: Cannot change to $OLCAO_DATA "
            f"directory '{olcao_data}': {e}"
        )

    # Define the list of database tarballs to unpack.
    databases = [
        "atomicBDB.tgz",
        "atomicPDB.tgz",
        "precursorDB.tgz",
        "spaceDB.tgz",
        "sybdDB.tgz",
    ]

    # Issue the necessary unpacking commands.
    for db in databases:
        result = subprocess.run(
                ["tar", "-xzf", db],
                capture_output=True, text=True)
        if result.returncode != 0:
            sys.exit(f"Failed to unpack {db}: "
                     f"{result.stderr.strip()}")

    # Enter into the space group database directory.
    try:
        os.chdir("spaceDB")
    except OSError as e:
        sys.exit(
            f"Error: Cannot change to spaceDB/ "
            f"directory: {e}"
        )

    # Recreate the space groups.
    result = subprocess.run(
            ["./remake"],
            capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(
            "Failed to remake the space group "
            f"database: {result.stderr.strip()}"
        )

    # Return to the original working directory (the OLCAO
    #   data directory).
    os.chdir("..")


def main():

    # Get script settings (parse help text, record invocation).
    settings = ScriptSettings()

    # Unpack the OLCAO databases.
    unpack_databases()

    # Report success.
    print("OLCAO databases unpacked successfully.")


if __name__ == '__main__':
    # Everything before this point was a subroutine definition
    #   or a request to import information from external
    #   modules. Only now do we actually start running the
    #   program. The purpose of this is to allow another python
    #   program to import *this* script and call its functions
    #   internally.
    main()
