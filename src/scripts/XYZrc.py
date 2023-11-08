#!/usr/bin/env python3


import os


def parameters_and_defaults():
    param_list = [
            "a_arg1",
            "a_arg2",
            1.0,
            2.0,
            3.0,
            False
            ]
    return param_list


if __name__ == '__main__':
    print(parameters_and_defaults())
