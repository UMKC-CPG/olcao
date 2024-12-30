#!/usr/bin/env python3

import argparse as ap
import os
import sys
from datetime import datetime


# Define the main class that holds script data structures and settings.
class ScriptSettings():
    """The instance variables of this object are the user settings that
       control the program. The variable values are pulled from a list
       that is created within a resource control file and that are then
       reconciled with command line parameters."""

       # Declare class variables.
       a_list = 0
       b_list = 0
       c = 0


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
        a_list = [default_rc["some"], default_rc["thing"]]

        # Second group of default variables.
        b_list = [default_rc["more"], default_rc["other"],
                default_rc["things"]]

        # Third group of default variables.
        c = default_rc["final_thing"]


    def parse_command_line(self):
    
        # Create the parser tool.
        prog_name = "<Program Name>"

        description_text = """
<Meta data: Version, date of last edit, relevant URL, requirements.>
<Description Text: Purpose, capabilities, limitations.>
"""

        epilog_text = """
Please contact <name> (<email>) regarding questions.
Defaults are given in ./XYZrc.py or $OLCAO_RC/XYZrc.py.
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
        parser.add_argument('-x', '--XYZx', nargs=2, dest='a_list',
                            type=str, default=a_list,
                            help=f'Argument a_list. Default: {a_list}')
    
        # Define the XYZb_list argument.
        parser.add_argument('-y', '--XYZb', nargs=3, dest='b_list',
                            type=float, default=b_list,
                            help=f'Argument b_list. Default: {b_list}')
    
        # Define the XYZc argument.
        parser.add_argument('-c', '--XYZc', dest='c', type=float,
                            default=c, help=f'Argument c. Default: {c}')
    
        # Define the positional arguments.
        parser.add_argument('filename', 
                            help='Name of the file to operate on.')


    def reconcile(self, args):
        a_list = args.a_list.deepcopy()
        b_list = args.b_list.deepcopy()
        c = args.c


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
