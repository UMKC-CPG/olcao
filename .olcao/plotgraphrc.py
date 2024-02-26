#!/usr/bin/env python3


import math as m
import os


def parameters_and_defaults():
    param_dict = {
            "title" : "Plotted by " + os.getenv('USER'),
            "print_command" : True,
            "page_orientation" : "portrait",  # "portrait"; "landscape"
            "fig_format" : "letter",  # "letter"; "A4"
            "fig_height" : 80,
            "fig_width" : 40,
            "fig_type" : "general",  # general; dos, optc, sybd, loci, ...
            "link_subplots" : True,
            "multi_x_cols" : False,
            "x_col" : 1,
            "y_col" : "",
            "x_min" : m.inf,
            "x_max" : m.inf,
            "y_min" : m.inf,
            "y_max" : m.inf,
            "x_axis_inc" : 5,
            "x_axis_minor_ticks" : 4,
            "y_axis_major_ticks" : 10,
            "y_axis_minor_ticks" : 4,
            "subplots_per_fig" : "1",
            "subplot_separation" : 10,
            "curves_per_subplot" : "1",
            "curve_thickness" : "750",
            "curve_style" : "0",
            "curve_color" : "1",
            "curve_separation" : "10"
            }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
