#!/usr/bin/env python3


import os


def parameters_and_defaults():
    param_dict = {
            "batch" : False, # T (F) = Generate SBATCH scripts (run directly).
            "num_cores" : 1, # Int = Number of cores for SBATCH grouping.
            }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
