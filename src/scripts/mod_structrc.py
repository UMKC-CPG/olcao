#!/usr/bin/env python3

"""mod_structrc.py -- Default parameters for mod_struct.py.

This resource control file provides a single function that returns a
dictionary of all default parameters used by modStruct.  Users may
place a customized copy of this file in the current working directory
or in the directory pointed to by $OLCAO_RC to override any defaults.
"""

import os


def parameters_and_defaults():
    param_dict = {

        # ---- Input / Output files ----
        "in_file":  "olcao.skl",    # Default input structure file.
        "out_file": "olcao.skl.new",  # Default output structure file.

        # ---- ABC-XYZ axis assignment order ----
        # These define how the lattice vectors (a,b,c) map to
        # Cartesian axes (x,y,z).  The defaults cause a to align
        # with x, b to be in the xy plane, and c to be arbitrary.
        "abc_order": [1, 2, 3],  # a=1, b=2, c=3
        "xyz_order": [1, 2, 3],  # x=1, y=2, z=3
    }
    return param_dict


if __name__ == "__main__":
    for k, v in parameters_and_defaults().items():
        print(f"  {k:24s} = {v!r}")
