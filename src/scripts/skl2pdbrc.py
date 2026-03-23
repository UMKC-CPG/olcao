#!/usr/bin/env python3

# Resource control file for skl2pdb.py.
#
# This file defines the default parameter values for the
# skl2pdb script. A copy of this file is installed to
# $OLCAO_RC during the build process. Users may place a
# modified copy in their current working directory to
# override the installed defaults on a per-project basis.
#
# To customize: copy this file into your working directory
# and edit the values below. Command line arguments will
# still override any values set here.


def parameters_and_defaults():
    param_dict = {
        # String: Path to the OLCAO skeleton input file.
        "input_file": "olcao.skl",

        # String: Path to the PDB output file.
        "output_file": "olcao.pdb",

        # String: Path to the SKL-to-PDB atom map file.
        "map_file": "sklPDB.map",

        # Boolean: If True, use the species types defined
        #   in the skeleton file as atom types in the PDB
        #   output.  If False, all atoms of the same
        #   element share a single type.
        "skl_types": False,
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
