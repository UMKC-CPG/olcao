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
        # Stored in the 1-indexed layout used everywhere else in
        # the OLCAO Python scripts: slot 0 is an unused None
        # sentinel, slots 1..3 hold the three axis labels, matching
        # Perl's ``@abc[1..3]`` / ``@xyz[1..3]`` convention.
        "abc_order": [None, 1, 2, 3],  # slot 1=a, 2=b, 3=c
        "xyz_order": [None, 1, 2, 3],  # slot 1=x, 2=y, 3=z
    }
    return param_dict


if __name__ == "__main__":
    for k, v in parameters_and_defaults().items():
        print(f"  {k:24s} = {v!r}")
