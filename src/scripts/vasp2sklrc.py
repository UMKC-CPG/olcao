#!/usr/bin/env python3

# Resource control file for vasp2skl.py.
#
# This file defines the default parameter values for the
# vasp2skl script.  A copy of this file is installed to
# $OLCAO_RC during the build process.  Users may place a
# modified copy in their current working directory to
# override the installed defaults on a per-project basis.
#
# To customize: copy this file into your working directory
# and edit the values below.  Command line arguments will
# still override any values set here.


def parameters_and_defaults():
    param_dict = {
        # String: Name of the VASP atomic position input
        #   file.  Typically "CONTCAR" (after relaxation)
        #   or "POSCAR" (initial positions).  If "CONTCAR"
        #   is not found, the script falls back to "POSCAR".
        "input_file": "CONTCAR",

        # String: Name of the file containing element names
        #   (typically the VASP POTCAR file, which has
        #   "TITEL" lines identifying each element).
        "atom_file": "POTCAR",

        # String: Name of the output OLCAO skeleton file.
        "output_file": "olcao.skl.relaxed",
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
