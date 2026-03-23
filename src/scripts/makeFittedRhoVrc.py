#!/usr/bin/env python3

"""Resource control file for makeFittedRhoV.py.

This file provides default parameter values for the
makeFittedRhoV.py script. Users can customize defaults by
placing a copy of this file in the current working directory
or in the directory pointed to by the $OLCAO_RC environment
variable.

The parameters control the input/output file names, the
database path for isolated atom potential/charge coefficient
files, the spin polarization setting, and whether openDX data
should be generated.

Parameters
----------
struct_file : str
    Path to the input structure file. Default: "structure.dat".
in_file : str
    Path to the SCF charge and potential data file produced
    by the OLCAO SCF calculation. Default: "gs_scfV-fb.dat".
out_file : str
    Path to the temporary output file that will be created
    for the OLCAOrhoV Fortran program. Default:
    "tempPotRho.dat".
atom_pot_db : str
    Path to the database directory containing isolated atom
    potential/charge coefficient files. Each element has its
    own subdirectory containing a "coeff.isolated" file.
    Default: "$OLCAO_DATA/atomicPDB".
spin_dn : bool
    If True, use the spin-down section of the potential file
    instead of the total/spin-up section. This is relevant
    only for spin-polarized calculations. Default: False.
do_open_dx : bool
    If True, generate openDX-format plottable data in
    addition to profile data. Default: True.
num_mesh_points : list of int
    Default mesh points along a, b, c axes. A value of [0,
    0, 0] indicates that mesh points *must* be provided on
    the command line. Default: [0, 0, 0].
num_cells : int
    Number of levels of neighboring replicated cells to
    consider when evaluating the charge density and potential
    functions. 0 = original cell only, 1 = first neighbors
    (27 total cells), 2 = first two levels (125 total cells).
    Default: 1.
"""


import os


def parameters_and_defaults():
    """Return a dictionary of default parameter values.

    Returns
    -------
    dict
        Keys are parameter names; values are the defaults.
        Each entry is annotated with its type and meaning.
    """

    # Resolve the database path from the environment.
    olcao_data = os.getenv('OLCAO_DATA', '')
    default_db = os.path.join(olcao_data, "atomicPDB")

    param_dict = {
        # Path to the input structure file.
        "struct_file": "structure.dat",  # String

        # Path to the SCF potential/charge data file.
        "in_file": "gs_scfV-fb.dat",  # String

        # Path to the temporary output file for OLCAOrhoV.
        "out_file": "tempPotRho.dat",  # String

        # Path to the isolated atom potential database.
        "atom_pot_db": default_db,  # String

        # T (F) = use spin-down (spin-up/total) section.
        "spin_dn": False,  # Boolean

        # T (F) = do (don't) generate openDX data.
        "do_open_dx": True,  # Boolean

        # Mesh points along a, b, c. [0,0,0] = must
        # be specified on the command line.
        "num_mesh_points": [0, 0, 0],  # List of int

        # Number of levels of neighboring replicated
        # cells. 0=original only, 1=27 cells, 2=125.
        "num_cells": 1,  # Integer
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
