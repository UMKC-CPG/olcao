#!/usr/bin/env python3


import os


def parameters_and_defaults():
    param_list = [
            False, # T (F) = include (exclude) two-center overlap
            False, # T (F) = include (exclude) two-center kinetic energy
            False, # T (F) = include (exclude) three-center nuclear attraction
            False, # T (F) = include (exclude) three-center overlap (with c=s)
            False  # T (F) = include (exclude) two-center momentum
            ]
    return param_list


if __name__ == '__main__':
    print(parameters_and_defaults())
