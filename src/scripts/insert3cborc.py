#!/usr/bin/env python3

"""Resource control file for insert3cbo.py.

This file provides default parameter values for the
insert3cbo.py script. Users can customize defaults by
placing a copy of this file in the current working
directory or in the directory pointed to by the $OLCAO_RC
environment variable.

The parameters control the input and output file names
used when merging two-center and three-center bond order
data.

Parameters
----------
data_3c : str
    Path to the three-center bond order raw data file.
    Default: "gs_bond-mb.3c.raw".
data_2c : str
    Path to the two-center bond order raw data file.
    Default: "gs_bond-mb.raw".
struct_data : str
    Path to the structure data file.
    Default: "structure.dat".
raw_23c_out : str
    Path to the output file for the combined two- and
    three-center bond order data.
    Default: "gs_bond-mb.23cbo.raw".
struct_23c_out : str
    Path to the output structure file with fake centroid
    atoms added.
    Default: "structure.23cbo.dat".
"""


import os


def parameters_and_defaults():
    """Return a dictionary of default parameter values.

    Returns
    -------
    dict
        Keys are parameter names; values are the
        defaults. Each entry is annotated with its type
        and meaning.
    """

    param_dict = {
        # Path to the 3-center bond order data file.
        "data_3c": "gs_bond-mb.3c.raw",  # String

        # Path to the 2-center bond order data file.
        "data_2c": "gs_bond-mb.raw",  # String

        # Path to the structure data file.
        "struct_data": "structure.dat",  # String

        # Path to the combined 2+3 center bond order
        # output file.
        "raw_23c_out": "gs_bond-mb.23cbo.raw",  # String

        # Path to the output structure file with fake
        # centroid atoms.
        "struct_23c_out": "structure.23cbo.dat",  # String
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
