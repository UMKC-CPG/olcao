#!/usr/bin/env python3

"""Resource control file for makePDOS.py.

This file provides default parameter values for the makePDOS.py
script. Users can customize defaults by placing a copy of this
file in the current working directory or in the directory pointed
to by the $OLCAO_RC environment variable.

The parameters control file names, output behavior, and the
orbital decomposition mode used when no explicit control file
is provided.

Parameters
----------
control_file : str
    Path to the user-created control file that defines how to
    collect and organize data from the raw PDOS file. An empty
    string means no control file is used and internal defaults
    will be applied. (See makePDOS.py -h for control file
    format details.)
raw_data_file : str
    Path to the raw PDOS data file produced by the OLCAO DOS
    program. Default: "gs_dos-fb.p.raw".
output_file : str
    Path to the output file where the processed PDOS curves
    will be written. Default: "PDOS.plot".
neg_to_zero : bool
    If True, any negative values in the final output will be
    replaced with 0.0. This can be useful for cleaning up
    small numerical artifacts. Default: False.
xanes : bool
    If True, the internally generated default control settings
    will decompose the PDOS curves into s+d and p states (a
    common decomposition for XANES analysis). Only applies
    when no explicit control file is given. Default: False.
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
        # Path to the control file. Empty string means use
        # internal defaults.
        "control_file": "",        # String

        # Path to the raw PDOS data file.
        "raw_data_file": "gs_dos-fb.p.raw",  # String

        # Path to the output file.
        "output_file": "PDOS.plot",  # String

        # If True, clamp negative values to zero in output.
        "neg_to_zero": False,      # Boolean

        # If True, use XANES-style s+d / p decomposition
        # for the default (no control file) case.
        "xanes": False,            # Boolean
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
