#!/usr/bin/env python3


import os


def parameters_and_defaults():
    param_list = [
            False, # T (F) = include (exclude) two-center overlap
            False, # T (F) = include (exclude) two-center kinetic energy
            False, # T (F) = include (exclude) three-center nuclear attraction
            False, # T (F) = include (exclude) three-center overlap (with c=s)
            False, # T (F) = include (exclude) two-center momentum
            False, # T (F) = include (exclude) mass velocity
            False, # T (F) = include (exclude) dipole moment term
            False, # T (F) = include (exclude) kinetic energy derivative term
            False, # T (F) = include (exclude) nuclear pot. derivative term
            False  # T (F) = include (exclude) electronic pot. derivative term
            ]
    return param_list


if __name__ == '__main__':
    print(parameters_and_defaults())
