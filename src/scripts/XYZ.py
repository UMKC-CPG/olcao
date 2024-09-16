#!/usr/bin/env python3

import XYZSettings as sets

def main():

    # Get script settings from a combination of the resource control file
    #   and parameters given by the user on the command line.
    settings = sets.ScriptSettings()

    # Start executing the main activities of the program.

    # Finalize the program activities and quit.


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. The purpose of this is to allow another
    #   python program to import *this* script and call its functions
    #   internally.
    main()
