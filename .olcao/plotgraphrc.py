#!/usr/bin/env python3


import os


def parameters_and_defaults():
    param_list = [
            "Plotted by " + os.getenv('USER'), # title
            "", # print_command
            80, # page_height
            40, # page_width
            0, # linked_plots
            0, # plot_type
            0, # multi_x_cols
            1, # x_col
            [2], # y_col
            "", # x_min
            "", # x_max
            "", # y_min
            "", # y_max
            5, # x_axis_inc
            4, # x_axis_minor_ticks
            10, # y_axis_major_ticks
            4, # y_axis_minor_ticks
            [1], # plots_per_image
            10, # plot_separation
            [1], # lines_per_plot
            [750], # line_thickness
            [0], # line_style
            [1], # line_color
            [10] # line_separation
            ]
    return param_list


if __name__ == '__main__':
    print(parameters_and_defaults())
