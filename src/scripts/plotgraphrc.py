#!/usr/bin/env python3


import math as m
import os


def parameters_and_defaults():
    param_dict = {
            "display" : "mpl",  # -d mpl; veusz, plotly, ...
            "infile" : "gs_dos-fb.t.plot",  # -i
            "outfile" : "Xscript.py",  # -o; If not given, X=display
            "title" : "Plotted by " + os.getenv('USER'),  # -t
            "print_command" : False,  # -pc
            "print_legend" : False,  # -pl
            "page_orientation" : "portrait",  # -po "portrait"; "landscape"
            "fig_format" : "letter",  # -ff "letter"; "A4"
            "fig_height" : 80,  # -fh
            "fig_width" : 40,  # -fw
            "fig_type" : "general",  # -ft general; dos, optc, sybd, loci, ...
            "link_subplots" : True,  # -l
            "multi_x_cols" : False,  # -mx
            "x_col" : 1,  # -xc
            "y_col" : "",  # -yc
            "x_min" : m.inf,  # -xi; Use m.inf to autoscale.
            "x_max" : m.inf,  # -xa; Use m.inf to autoscale.
            "y_min" : m.inf,  # -yi; Use m.inf to autoscale.
            "y_max" : m.inf,  # -ya; Use m.inf to autoscale.
            "x_multiple" : 5,  # -xm
            "y_multiple" : 1,  # -ym
            "x_axis_inc" : 5,  # -xai
            "x_axis_minor_ticks" : 4,  # -xmi
            "y_axis_major_ticks" : 10,  # -yma
            "y_axis_minor_ticks" : 4,  # -ymi
            "linked_x_axes" : False,  # -lx
            "linked_y_axes" : False,  # -ly
            "subplot_separation" : 10,  # -ss
            "subplots_per_fig" : [1],  # -sp
            "curves_per_subplot" : [1],  # -cs
            "curve_separation" : [10],  # -cp
            "curve_width" : [1],  # -ct
            "curve_width_start" : 0.5,  # -cta
            "curve_width_step" : 0.25,  # -cte
            "curve_width_size" : 1,  # -cti
            "curve_style" : [],  # -cy
            "curve_style_start" : 0,  # -cya
            "curve_style_step" : 1,  # -cye
            "curve_style_size" : 1,  # -cyi
            "curve_color" : [],  # -cc
            "curve_color_start" : 18,  # -cca
            "curve_color_step" : 1,  # -cce
            "curve_color_size" : 1,  # -cci
            "curve_mark" : [],  # -cm
            "curve_mark_start" : 0,  # -cma
            "curve_mark_step" : 1,  # -cme
            "curve_mark_size" : 1  # -cmi
            }
    # Note that the colors use:  Use https://xkcd.com/color/rgb/
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
