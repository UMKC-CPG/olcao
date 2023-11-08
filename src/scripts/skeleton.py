#!/usr/bin/env python3

import argparse as ap
import os
import sys
from datetime import datetime

def print_help():
    output = """
#Script to XYZ
#Last updated Month, DD, YYYY
########################################################################
#This python script will execute the necessary steps to ...
# Note that default parameters are define by the file "$OLCAO_RC/XYZrc".
#
#-----------------------------------------------------------------------
#
#USAGE:  XYZ [-x a1 a2] [-y b1 b2 b3] [-z c] filename
#OR:     XYZ --help | -h
#
#If the "-x" option is given then a1 and a2 are ...
#If the "-y" option is given then b1, b2, b3 are ...
#If the "-z" option is given then c is ...
#The "filename" is required and has no default value here.
#
#
#DEFAULTS: Check the XYZrc file in either the local directory or in
#   the $OLCAO_RC directory for total confirmation of default values.
#For -x: The default is probably ...
#For -y: The default is probably ...
#For -z: The default is probably ...
#
#-----------------------------------------------------------------------
#
#REQUIREMENTS
#
#You need to have ...
#You need to be in ...
#
"""
    print (output)


# Define the main class that holds script data structures and settings.
class ScriptSettings():
    """The instance variables of this object are the user settings that
       control the program. The variable values are pulled from a list
       that is created within a resource control file and that are then
       reconciled with command line parameters."""


    def __init__(self):
        """Define default values for the graph parameters by pulling them
        from the resource control file in the default location:
        $OLCAO_RC/XYZrc.py or from the current working directory if a local
        copy of XYZrc.py is present."""

        # Read default variables from the resource control file.
        sys.path.insert(1, os.getenv('OLCAO_RC'))
        from XYZrc import parameters_and_defaults
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
        self.a_list = [default_rc[0], default_rc[1]]

        # Second group of default variables.
        self.b_list = [default_rc[2], default_rc[3], default_rc[4]]

        # Third group of default variables.
        self.c = default_rc[5]


    def parse_command_line(self):
    
        # Create the parser tool.
        prog_name = "<Program Name>"
        description_text = """
<Description Text>
"""
        epilog_text = """
<Epilog Text>
"""
        parser = ap.ArgumentParser(prog = prog_name,
                                   description = description_text,
                                   epilog = epilog_text)
    
        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()


    def add_parser_arguments(self, parser):
    
        # Define the XYZa_list argument.
        parser.add_argument('-x', '--XYZx', nargs=2, dest='a_list',
                            type=str, default=self.a_list,
                            help='Arguments a_list. Default: ' +
                            f'{self.a_list}')
    
        # Define the XYZb_list argument.
        parser.add_argument('-y', '--XYZb', nargs=3, dest='b_list',
                            type=float, default=self.b_list,
                            help='Arguments b_list. ' +
                            f'Default: {self.b_list}')
    
        # Define the XYZc argument.
        parser.add_argument('-c', '--XYZc', dest='c', type=float,
                            default=self.c, help='Argument c. ' +
                            f'Default: {self.c}')
    
        # Define the positional arguments.
        parser.add_argument('filename', 
                            help='Name of the file to operate on.')


    def reconcile(self, args):
        self.a_list = args.a_list.copy()
        self.b_list = args.b_list.copy()
        self.c = args.c


    def recordCLP(self):
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write(f"Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


def main():

    # Get script settings from a combination of the resource control file
    #   and parameters given by the user on the command line.
    settings = ScriptSettings()

    # Start executing the main activities of the program.

    # Finalize the program activities and quit.


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. The purpose of this is to allow another
    #   python program to import *this* script and call its functions
    #   internally.
    main()
