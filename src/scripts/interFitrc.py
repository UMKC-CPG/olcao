#!/usr/bin/env python3


import os


def parameters_and_defaults():
    param_dict = {
            "comp" : 1,        # Int = Number of components (1, 2, or 4).
            "shrink_rc" : 0.0,  # Real = Shrinking function critical radius.
            "shrink_sig" : 0.0, # Real = Shrinking function sigma parameter.
            "tolerance" : 1e-3, # Real = RMSE tolerance for fitting acceptance.
            }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
