#!/usr/bin/env python3

# Resource control file for struct2skl.py.
#
# This file defines the default parameter values for the
# struct2skl script. A copy of this file is installed to
# $OLCAO_RC during the build process. Users may place a
# modified copy in their current working directory to
# override the installed defaults on a per-project basis.
#
# To customize: copy this file into your working directory
# and edit the values below. Command line arguments will
# still override any values set here.


def parameters_and_defaults():
    param_dict = {
        # String: Path to the structure.dat input file.
        "input_file": "structure.dat",

        # String: Path to the OLCAO skeleton output file.
        "output_file": "olcao.skl",

        # String: Path to the struct-to-SKL atom map file.
        "map_file": "sklStruct.map",

        # Boolean: If True, use the species types from the
        #   structure.dat file in the skeleton output.
        #   If False, all atoms of the same element share
        #   a single species (type 1).
        "struct_types": False,

        # Boolean: If True, export the structure as a
        #   cluster (molecule) rather than inside a
        #   periodic cell.
        #   NOTE: This option is documented in the original
        #   Perl script's help text but was not implemented
        #   in the Perl code.  It is included here for
        #   future use.
        "make_mol": False,
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
