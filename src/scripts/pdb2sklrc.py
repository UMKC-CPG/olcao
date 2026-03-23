#!/usr/bin/env python3

# Resource control file for pdb2skl.py.
#
# This file defines the default parameter values for the
# pdb2skl script. A copy of this file is installed to
# $OLCAO_RC during the build process. Users may place a
# modified copy in their current working directory to
# override the installed defaults on a per-project basis.
#
# To customize: copy this file into your working directory
# and edit the values below. Command line arguments will
# still override any values set here.


def parameters_and_defaults():
    param_dict = {
        # String: Path to the PDB input file.
        "input_file": "model.pdb",

        # String: Path to the OLCAO skeleton output file.
        "output_file": "olcao.skl",

        # String: Path to the PDB-to-SKL atom map file.
        "map_file": "skl2PDB.map",

        # Float: Buffer size (Angstroms) added around the
        #   molecule if the PDB file does not define a cell.
        #   Default: 0.0 (no buffer).
        "buffer": 0.0,

        # Boolean: If True, use the atom types defined in
        #   the PDB file as species tags in the skeleton
        #   file.  If False, all atoms of the same element
        #   share a single species (type 1).
        "pdb_types": False,
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
