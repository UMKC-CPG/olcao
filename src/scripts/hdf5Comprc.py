#!/usr/bin/env python3


import os


def parameters_and_defaults():
    param_dict = {
            "file1" : "gs_pscf-fb.hdf5", # File 1
            "file2" : "gs_pscf-fb.hdf5", # File 1
            "dir1" : ".", # Dir to find file 1 intermediate.
            "dir2" : ".", # Dir to find file 2 intermediate.
            "group1" : "00001_00001_00001/atomIntgGroup/atomOverlap", # Group1
            "group2" : "00001_00001_00001/atomIntgGroup/atomOverlap", # Group2
            "dataset1" : "0000001", # Dataset1
            "dataset2" : "0000002", # Dataset2
            "unpack" : False, # Set to True to unpack.
            "imaginary" : False, # Set True for imag. part of packed matrix.
            "magnitude" : False # Set True for absolute magnitude of complex.
            }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
