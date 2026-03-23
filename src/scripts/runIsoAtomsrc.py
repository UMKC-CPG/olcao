#!/usr/bin/env python3

# Resource control file for runIsoAtoms.py.
#
# This file defines the default parameter values for the
# runIsoAtoms script.  A copy of this file is installed to
# $OLCAO_RC during the build process.  Users may place a
# modified copy in their current working directory to
# override the installed defaults on a per-project basis.
#
# To customize: copy this file into your working directory
# and edit the values below.  Command line arguments will
# still override any values set here.


def parameters_and_defaults():
    param_dict = {
        # Integer: Atomic Z number of a single element to
        #   process.  When set to 0 (default), all elements
        #   in the database are processed.  When set to a
        #   specific Z number (e.g. 14 for silicon), only
        #   that element is processed.
        "target_z": 0,

        # Boolean: If True, the -nocore option is passed
        #   to makeinput so that the all-electron electronic
        #   structure is computed.  (It is somewhat counter-
        #   intuitive that an option called -nocore should
        #   yield an all-electron result.  But, the idea is
        #   that -nocore *will not* include any core orbitals
        #   in the list of orbitals that are orthogonalized
        #   away.  Therefore, all of the core orbitals *will
        #   be* included in the list of orbitals to compute.)
        "nocore": False,
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
