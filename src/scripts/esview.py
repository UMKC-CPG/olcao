!#/usr/bin/env python3

import argparse as ap
import os
import sys
import numpy as np
import vedo as v
import math as m


class ScriptSettings():
    """The instance variables of this object are the user settings that
       control the program. The variable values are pulled from a list
       that is created within a resource control file and that are then
       reconciled with command line parameters."""

    def __init__(self):
        # Read default variable values from the resource control file.
        sys.path.insert(1, os.getenv('OLCAO_DIR') + "/.olcao")
        from esviewrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults file.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line(default_rc)

        # Reconcile default settings in the resource control file
        #   with any command line arguments that were given.
        self.reconcile(args, default_rc)
