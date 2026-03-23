#!/usr/bin/env python3

"""dump2sklrc.py -- Resource control (default parameters) for dump2skl.py.

This file defines the default values for every command line parameter
accepted by dump2skl.py.  A copy is installed to $OLCAO_RC during the
CMake install step.  Users may place a local copy in their working
directory to override the installed defaults on a per-project basis.

To see what each parameter does, run:  dump2skl.py --help
"""


import os


def parameters_and_defaults():
    param_dict = {
        # Path to the LAMMPS dump file (required; no useful
        # default).
        "dump_file": "",

        # Path to the LAMMPS data file (required; no useful
        # default).
        "data_file": "",

        # Frame number to extract.  1 = first frame,
        # -1 = last frame.  Set to None when using timestep
        # instead.
        "frame": None,

        # Timestep number to extract (must match a value in
        # the dump file).  Set to None when using frame
        # instead.
        "timestep": None,

        # System name written into the skeleton file header.
        "name": "name",

        # T (F) = Add (don't add) 100 Angstrom padding on
        # every side of the cell to suppress periodic
        # interactions.
        "non_periodic": False,
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
