#!/usr/bin/env python3


import os


def parameters_and_defaults():
    param_dict = {
            # T (F) = Optical properties computed in serial
            #   (simultaneous) x,y,z fashion. Serial conserves
            #   memory at the cost of time.
            "serialxyz" : False,
            # T (F) = Run Fortran programs in (without) a
            #   valgrind environment with --leak-check.
            "valgrind" : False,
            }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
