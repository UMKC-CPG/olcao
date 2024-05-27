#!/usr/bin/env python3

import argparse as ap
import os
import sys
import pandas as pd
from datetime import datetime


# Define the main class that holds script data structures and settings.
class ScriptSettings():
    """The instance variables of this object are the user settings that
       control the program. The variable values are pulled from a list
       that is created within a resource control file and that are then
       reconciled with command line parameters."""


    def __init__(self):
        """Define default values for the graph parameters by pulling them
        from the resource control file in the default location:
        $OLCAO_RC/collectSYBDrc.py or from the current working directory if a
        local copy of collectSYBDrc.py is present."""

        # Read default variables from the resource control file.
        sys.path.insert(1, os.getenv('OLCAO_RC'))
        from collectSYBDrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults file.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc file.
        self.reconcile(args)

        # At this point, the command line parameters are set and accepted.
        #   When this initialization subroutine returns the script will
        #   start running. So, we use this as a good spot to record the
        #   command line parameters that were used.
        self.recordCLP()


    def assign_rc_defaults(self, default_rc):

        # First group of default variables.
        self.hskp = default_rc["hskp"]


    def parse_command_line(self):
    
        # Create the parser tool.
        prog_name = "collectSYBD.py"

        description_text = """
Version: 14d-02m-2024y

The purpose of this program is to gather a set of symmetric band structure
(SYBD) plots and then parse them for specific information that will be used
to create a new specific type of plot.
"""

        epilog_text = """
Please contact Paul Rulis (rulisp@umkc.edu) regarding questions.
Defaults are given in ./collectSYBDrc.py or $OLCAO_RC/collectSYBDrc.py.
"""

        parser = ap.ArgumentParser(prog = prog_name,
                formatter_class=ap.RawDescriptionHelpFormatter,
                description = description_text,
                epilog = epilog_text)
    
        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()


    def add_parser_arguments(self, parser):
    
        # Define the XYZa_list argument.
        parser.add_argument('-d', '--delta', nargs=1, dest='hskp',
                            type=str, default=self.hskp,
                            help='Argument hskp: high symmetry k-point. ' +
                            f'Default: {self.hskp}')
    

    def reconcile(self, args):
        self.hskp = args.hskp


    def recordCLP(self):
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write(f"Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


# Traverse the direcetory tree from the current working directory and
#   look for any file that matches the target. Then, save that file and
#   its directory location in a sorted list.
def make_SYBD_file_list():
    target_file = "gs_sybd-fb.plot"

    # Initialize the list of sybd files.
    sybd_file_list = []

    # Traverse the directory tree from the current working directory.
    for (root, dirs, files) in os.walk(os.getcwd()):

        # Add any matching file to a list.
        if (target_file in files):
            sybd_file_list.append(f"{root}/{target_file}")

    # Return the sorted list.
    return sybd_file_list.sort()


def get_hskp_data(hskp, sybd_file_list):
    # Consider each sybd file.
    for sybd_file in sybd_file_list:

        line = []
        hskp_pos = 0.0
        with open(sybd_file, 'r') as f:
            curr_hskp = ""
            while (curr_hskp != hskp):

                line = f.readline().split()
                curr_hskp = line[len(line)-1]
                print(curr_hskp, hskp)

            hskp_pos = line[len(line)-2]
            if (hskp_pos == 0.0):


        with open(sybd_file, 'r') as f:
            for


        # Import the sybd file into a Pandas data frame.
        sybd = pd.read_table(sybd_file, header=None, delim_whitespace=True)
        print(sybd)

        # Get the rows of data that contain the target hskp name.
        #hskp_data = sybd.loc()


def main():

    # Get script settings from a combination of the resource control file
    #   and parameters given by the user on the command line.
    settings = ScriptSettings()

    # Create a list of all the SYBD plot files in subdirectories of the pwd.
    sybd_file_list = make_SYBD_file_list()

    # Open each file in the order obtained and extract the row associated
    #   with the targeted high symmetry k-point.
    get_hskp_data(settings.hskp, sybd_file_list)

    # Print a new "SYBD plot" file that shows a sequence of energy eigenvalue
    #   changes for the targeted hskp.
    #print_delta_hskp(settings.hskp)


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. The purpose of this is to allow another
    #   python program to import *this* script and call its functions
    #   internally.
    main()
