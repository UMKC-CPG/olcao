#!/usr/bin/env python3

"""condenserc.py -- Default parameters for condense.py.

This resource control file provides a single function that returns a
dictionary of all default parameters used by condense.py.  Users may
place a customized copy of this file in the current working directory
or in the directory pointed to by $OLCAO_RC to override any defaults.

The condense script creates all the necessary input files for running
a LAMMPS condensation simulation on a set of molecules with defined
reaction types.  The parameters here control file paths, simulation
cell geometry, initial molecular speeds, and related options.
"""


def parameters_and_defaults():
    """Return a dictionary of default parameters for condense.py.

    Returns
    -------
    dict
        Keys are parameter names; values are the defaults.
    """
    param_dict = {

        # ---- Input file ----
        # Name of the condense input file that specifies the
        # composition, cell size, stages, and reactions.
        "input_file":       "condense.in",

        # ---- Force collision mode ----
        # When True, the lammps.dat input file will arrange
        # molecules so that they are put on a forced collision
        # course near the center of the cell.  This option only
        # works when the number of molecules is three or fewer.
        "force_collision":  False,

        # ---- Simulation cell ----
        # Default cubic cell side length in Angstroms.
        "cell_size":        100.0,

        # ---- Molecular speed ----
        # Maximum initial speed of each molecule
        # (Angstroms/time-step).
        "max_speed":        10.0,

        # ---- Directory names ----
        # Name of the directory where combined reaction template
        # files are stored.
        "rxn_template_dir": "reactionTemplates",
    }

    return param_dict


if __name__ == "__main__":
    for k, v in parameters_and_defaults().items():
        print(f"  {k:24s} = {v!r}")
