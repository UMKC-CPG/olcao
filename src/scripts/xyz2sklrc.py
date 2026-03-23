#!/usr/bin/env python3

# Resource control file for xyz2skl.py.
#
# This file defines the default parameter values for the
# xyz2skl script.  A copy of this file is installed to
# $OLCAO_RC during the build process.  Users may place a
# modified copy in their current working directory to
# override the installed defaults on a per-project basis.
#
# To customize: copy this file into your working directory
# and edit the values below.  Command line arguments will
# still override any values set here.


def parameters_and_defaults():
    param_dict = {
        # String: Name of the input file containing the
        #   structure in pure XYZ direct space coordinates.
        "input_file": "model.xyz",

        # String: Name of the output OLCAO skeleton file.
        "output_file": "olcao.skl",

        # String: Name of the mapping file that records
        #   the correspondence between XYZ and SKL atoms.
        "map_file": "sklXYZ.map",

        # Boolean: If True, the skeleton file will use the
        #   atom types as defined in the XYZ file.  If
        #   False (default), no types are assigned.
        "xyz_types": False,
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
