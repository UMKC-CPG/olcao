#!/usr/bin/env python3


import os


def parameters_and_defaults():
    param_list = [
            True, # T (F) = show (do not show) reciprocal lattice points
            True, # T (F) = show (do not show) 
            False, # T (F) = use (do not use) vectorization; !production only!
            4, # highest angular momentum (0=s; 1=p; 2=d; 3=f; 4=g; etc.)
            False, # T (F) = include (exclude) two-center overlap
            False, # T (F) = include (exclude) two-center kinetic energy
            False, # T (F) = include (exclude) three-center nuclear attraction
            False, # T (F) = include (exclude) three-center overlap (with c=s)
            False  # T (F) = include (exclude) two-center momentum
            ]
    return param_list


if __name__ == '__main__':
    print(parameters_and_defaults())
