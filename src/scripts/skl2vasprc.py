#!/usr/bin/env python3

# Resource control file for skl2vasp.py.
#
# This file defines the default parameter values for the
# skl2vasp script. A copy of this file is installed to
# $OLCAO_RC during the build process. Users may place a
# modified copy in their current working directory to
# override the installed defaults on a per-project basis.
#
# To customize: copy this file into your working directory
# and edit the values below. Command line arguments will
# still override any values set here.


def parameters_and_defaults():
    param_dict = {
        # String: Type of VASP calculation to set up.
        #   Valid options: "relaxfull", "relaxion1",
        #   "relaxion2", "relaxvol", "cij", "phonon",
        #   "gpt".
        "job_type": "relaxfull",

        # String: Base type of pseudopotential.
        #   Valid options: "potLDA", "potGGA",
        #   "potpawLDA", "potpawGGA", "potpawPBE",
        #   "potpawLDA5x", "potpawPBE5x".
        "pot_type": "potpawPBE",

        # String: Sub-type of the potential (e.g. "s",
        #   "h", "pv", "sv", "d", "soft", "200ev", etc.)
        #   Empty string means use the standard potential.
        "sub_pot_type": "",

        # Boolean: If True, append the GW sub-sub-type
        #   to the potential sub-type.
        "gw": False,

        # Boolean: If True, use a Gamma-centered k-point
        #   mesh instead of a general Monkhorst-Pack mesh.
        "gamma_kpoint": False,

        # String: Override the Bravais lattice cell type
        #   that would otherwise be determined from the
        #   space group.  Valid options: "tric", "mono",
        #   "ortho", "tet", "trig", "hex", "sc", "bcc",
        #   "fcc".  Empty string means auto-detect.
        "cell_type": "",

        # String: Path to the OLCAO skeleton input file.
        "input_file": "olcao.skl",
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
