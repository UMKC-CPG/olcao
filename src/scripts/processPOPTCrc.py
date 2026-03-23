#!/usr/bin/env python3

"""Resource control file for processPOPTC.py.

This file provides default parameter values for the
processPOPTC.py script. Users can customize defaults by
placing a copy of this file in the current working directory
or in the directory pointed to by the $OLCAO_RC environment
variable.

The parameters control which spin channel is processed and
the file names used for intermediate and output data.

Parameters
----------
spin : int
    Spin channel to process. 1 = spin up (or total for
    non-spin-polarized calculations). 2 = spin down.
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

    param_dict = {
        # Spin channel. 1 = up/total, 2 = down.
        "spin": 1,  # Integer
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
