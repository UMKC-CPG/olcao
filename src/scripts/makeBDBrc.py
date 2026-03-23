#!/usr/bin/env python3

# Resource control file for makeBDB.py.
#
# This file defines the default parameter values for the
# makeBDB script. A copy of this file is installed to $OLCAO_RC
# during the build process. Users may place a modified copy in
# their current working directory to override the installed
# defaults on a per-project basis.
#
# To customize: copy this file into your working directory and
# edit the values below. Command line arguments will still
# override any values set here.


def parameters_and_defaults():
    param_dict = {
        # String: Which step(s) to perform.
        #   "all"  = all three steps (SCF, fit, contract).
        #   "scf"  = only the atomic SCF step.
        #   "fit"  = only the Gaussian fitting step.
        #   "cont" = only the contraction file step.
        "step": "all",

        # Integer: Target element Z number. 0 means process
        #   all elements. A non-zero value restricts
        #   processing to that single element.
        "element": 0,

        # Integer: Number of parallel processes (forks) to
        #   use when performing the Gaussian fitting of
        #   numerical potentials. 1 means serial execution.
        "fork": 1,

        # Boolean: If True, overwrite any previously obtained
        #   results for the requested operations.
        "force": False,

        # Boolean: If True, print progress reports as the
        #   program works.
        "verbose": False,
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
