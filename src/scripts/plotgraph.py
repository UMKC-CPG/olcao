#!/usr/bin/env python3

import argparse as ap
import pandas as pd
import os
import sys
import copy
import matplotlib.pyplot as plt
import math as m
import subprocess
from datetime import datetime


# Define the main class that holds plot structure and information settings.
class ScriptSettings():
    """The instance variables of this object are the user settings that
       control the program. The variable values are pulled from a list
       that is created within a resource control file and that are then
       reconciled with command line parameters."""


    def __init__(self):
        """Define default values for the graph parameters by pulling them
        from the resource control file $OLCAO_RC/plotgraphrc.py"""

        # Read default variables from the plotgraph resource control file.
        sys.path.insert(1, os.getenv('OLCAO_RC'))
        from plotgraphrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults file.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line(default_rc)

        # Reconcile default settings in the resource control file with
        #   any command line arguments that were given.
        self.reconcile(args)

        # Now that the arguments (command line or default from the rc file)
        #   have been assigned, we can do some follow-up processing such as
        #   expanding some expressions into arrays or noting if axis limits
        #   were explicitly specified on either the command line or the rc
        #   file.
        self.process_settings()

        # Record the command line parameters.
        self.recordCLP()


    def assign_rc_defaults(self, default_rc):

        # Graphic display method.
        self.display = default_rc["display"]

        # Filenames
        self.infile = default_rc["infile"]
        self.outfile = default_rc["outfile"]

        # Supporting display information
        self.title = default_rc["title"]
        self.print_command = default_rc["print_command"]
        self.print_legend = default_rc["print_legend"]

        # Plot formats
        self.fig_format = default_rc["fig_format"]
        self.page_orientation = default_rc["page_orientation"]
        self.fig_height = default_rc["fig_height"]
        self.fig_width = default_rc["fig_width"]
        self.fig_type = default_rc["fig_type"]
        self.link_subplots = default_rc["link_subplots"]

        # Data columns
        self.multi_x_cols = default_rc["multi_x_cols"]
        self.x_col = default_rc["x_col"]
        self.y_col = default_rc["y_col"]

        # Ranges and ticks
        self.x_min = default_rc["x_min"]
        self.x_max = default_rc["x_max"]
        self.y_min = default_rc["y_min"]
        self.y_max = default_rc["y_max"]
        self.x_multiple = default_rc["x_multiple"]
        self.y_multiple = default_rc["y_multiple"]
        self.x_axis_inc = default_rc["x_axis_inc"]
        self.x_axis_minor_ticks = default_rc["x_axis_minor_ticks"]
        self.y_axis_major_ticks = default_rc["y_axis_major_ticks"]
        self.y_axis_minor_ticks = default_rc["y_axis_minor_ticks"]
        self.linked_x_axes = default_rc["linked_x_axes"]
        self.linked_y_axes = default_rc["linked_y_axes"]

        # Subplot and curve variables
        self.subplot_separation = default_rc["subplot_separation"]
        self.subplots_per_fig = default_rc["subplots_per_fig"]
        self.curves_per_subplot = default_rc["curves_per_subplot"]
        self.curve_separation = default_rc["curve_separation"]
        self.curve_width = default_rc["curve_width"]
        self.curve_width_start = default_rc["curve_width_start"]
        self.curve_width_step = default_rc["curve_width_step"]
        self.curve_width_size = default_rc["curve_width_size"]
        self.curve_style = default_rc["curve_style"]
        self.curve_style_start = default_rc["curve_style_start"]
        self.curve_style_step = default_rc["curve_style_step"]
        self.curve_style_size = default_rc["curve_style_size"]
        self.curve_color = default_rc["curve_color"]
        self.curve_color_start = default_rc["curve_color_start"]
        self.curve_color_step = default_rc["curve_color_step"]
        self.curve_color_size = default_rc["curve_color_size"]
        self.curve_mark = default_rc["curve_mark"]
        self.curve_mark_start = default_rc["curve_mark_start"]
        self.curve_mark_step = default_rc["curve_mark_step"]
        self.curve_mark_size = default_rc["curve_mark_size"]

        #    global dos = default_rc[24]
        #    global points = default_rc[25]


    def parse_command_line(self, default_rc):
    
        # Create the parser tool.
        parser = ap.ArgumentParser(description='Control parameters')
    
        # Add arguments to the parser.
        self.add_parser_arguments(parser, default_rc)

        # Parse the arguments.
        args = parser.parse_args()

        # Return the result.
        return args


    def list_of_strings(self, argument):
        return argument.split(',')


    def add_parser_arguments(self, parser, default_rc):

        # Define the selected display library to use.
        parser.add_argument('-d', '--display', dest='display',
                            choices=['mpl', 'veusz'],
                            default=self.display, help='Display library to ' +
                            f'use. Default: {self.display}.')

        # Define the input file and output file arguments.
        parser.add_argument('-i', '--infile', dest='infile', type=ascii,
                            default=self.infile, help='Input file name. ' +
                            f'Default: {self.infile}.')
        parser.add_argument('-o', '--outfile', dest='outfile', type=ascii,
                            default=self.outfile, help='Output file name. ' +
                            f'Default: {self.outfile}.')
    
        # Define the title argument.
        parser.add_argument('-t', '--title', dest='title', type=ascii,
                            default=self.title,
                            help='Title for the figure set. Default: ' +
                            f'{self.title}.')
    
        # Define the flag to print the command used to create the figure.
        if (self.print_command == True):
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-pc', '--printcmd', action=store_action,
                            dest='print_command', default=self.print_command,
                            help='Flag to print the command line. Default: ' +
                            f'{self.print_command}.')
    
        # Define the flag to print a legend into the figure.
        if (self.print_legend == True):
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-pl', '--printlgd', action=store_action,
                            dest='print_legend', default=self.print_legend,
                            help='Flag to print a legend. Default: ' +
                            f'{self.print_legend}.')

        # Define the page orientation.
        parser.add_argument('-po', '--pageo', dest='page_orientation',
                            choices=['portrait', 'landscape'],
                            default=self.page_orientation,
                            help='Page orientation. Default: ' +
                            f'{self.page_orientation}.')

        # Define the fig format.
        parser.add_argument('-ff', '--figformat', dest='fig_format',
                            choices=['letter', 'A4'], default=self.fig_format,
                            help=f'Fig format. Default: {self.fig_format}.')

        # Define the fig height.
        parser.add_argument('-fh', '--figh', dest='fig_height', type=float,
                            default=self.fig_height, help='Fig height. ' +
                            f'Default: {self.fig_height} percent.')

        # Define the fig width.
        parser.add_argument('-fw', '--figw', dest='fig_width', type=float,
                            default=self.fig_width, help='Fig width. ' +
                            f'Default: {self.fig_width} percent.')

        # Define the figure type.
        parser.add_argument('-ft', '--figt', dest='fig_type',
                            choices=['general','dos','optc','sybd'],
                            default=self.fig_type,
                            help='Flag to request a specific figure type. ' +
                            f'Default: {self.fig_type}.')

        # Define the flag to request linked subplots.
        if (self.link_subplots == True):
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-l', '--link_subplots', action=store_action,
                            dest='link_subplots', default=self.link_subplots,
                            help='Flag to request linked subplots. Default: '+
                            f'{self.link_subplots}.')

        # Define the flag to indicate that the input has multiple x columns.
        if (self.multi_x_cols == True):
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-mx', '--multix', action=store_action,
                            dest='multi_x_cols', default=self.multi_x_cols,
                            help='Flag to indicate that input has multiple ' +
                            f'x columns. Default: {self.multi_x_cols}.')

        # Select the x_column.
        parser.add_argument('-xc', '--xcol', dest='x_col', type=int,
                            default=self.x_col, help='X column. Default: ' +
                            f'{self.x_col}.')

        # Select the y_column(s).
        parser.add_argument('-yc', '--ycol', dest='y_col', type=ascii,
                            default=self.y_col, help='Y column(s) in array ' +
                            f'form with .. options. Default: {self.y_col}.')

        # Set the x_min
        if (self.x_min == m.inf):
            x_min_default = "Data min"
        else:
            x_min_default = self.x_min
        parser.add_argument('-xi', '--xmin', dest='x_min', type=float,
                            default=self.x_min, help='X min value. ' +
                            f'Default: {x_min_default}.')

        # Set the x_max
        if (self.x_max == m.inf):
            x_max_default = "Data max"
        else:
            x_max_default = self.x_max
        parser.add_argument('-xa', '--xmax', dest='x_max', type=float,
                            default=self.x_max, help='X max value. ' +
                            f'Default: {x_max_default}.')

        # Set the y_min
        if (self.y_min == m.inf):
            y_min_default = "Data min"
        else:
            y_min_default = self.y_min
        parser.add_argument('-yi', '--ymin', dest='y_min', type=float,
                            default=self.y_min, help='Y min value. ' +
                            f'Default: {y_min_default}.')

        # Set the y_max
        if (self.y_max == m.inf):
            y_max_default = "Data max"
        else:
            y_max_default = self.y_max
        parser.add_argument('-ya', '--ymax', dest='y_max', type=float,
                            default=self.y_max, help='Y max value. ' +
                            f'Default: {y_max_default}.')

        # Set the default x_multiple for axis rounding.
        parser.add_argument('-xm', '--xmult', dest='x_multiple',
                            type=float, default=self.x_multiple, help='X ' +
                            f'axis multiple to round limits to. Default: ' +
                            f'{self.x_multiple}.')

        # Set the default y_multiple for axis rounding.
        parser.add_argument('-ym', '--ymult', dest='y_multiple',
                            type=float, default=self.y_multiple, help='Y ' +
                            f'axis multiple to round limits to. Default: ' +
                            f'{self.y_multiple}.')

        # Set the default x axis increment.
        parser.add_argument('-xai', '--xaxisinc', dest='x_axis_inc',
                            type=float, default=self.x_axis_inc, help='X ' +
                            f'axis increment. Default: {self.x_axis_inc}.')

        # Set the default x axis minor ticks.
        parser.add_argument('-xmi', '--xminor', dest='x_axis_minor_ticks',
                            type=int, default=self.x_axis_minor_ticks,
                            help='X axis increment. Default: ' +
                            f'{self.x_axis_minor_ticks}.')

        # Set the default y axis major ticks.
        parser.add_argument('-yma', '--ymajor', dest='y_axis_major_ticks',
                            type=int, default=self.y_axis_major_ticks,
                            help='Y axis major ticks. Default: ' +
                            f'{self.y_axis_major_ticks}.')

        # Set the default y axis minor ticks.
        parser.add_argument('-ymi', '--yminor', dest='y_axis_minor_ticks',
                            type=int, default=self.y_axis_minor_ticks,
                            help='Y axis minor ticks. Default: ' +
                            f'{self.y_axis_minor_ticks}.')

        # Set the default behavior for whether to link x axes or not.
        if (self.linked_x_axes == True):
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-lx', '--linkx', action=store_action,
                            dest='linked_x_axes', default=self.linked_x_axes,
                            help='Flag to indicate that the x axes should ' +
                            f'be linked. Default: {self.linked_x_axes}.')

        # Set the default behavior for whether to link y axes or not.
        if (self.linked_y_axes == True):
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-ly', '--linky', action=store_action,
                            dest='linked_y_axes', default=self.linked_y_axes,
                            help='Flag to indicate that the y axes should ' +
                            f'be linked. Default: {self.linked_y_axes}.')

        # Define the plot separation.
        parser.add_argument('-ss', '--subplotsep', dest='subplot_separation',
                            type=int, default=self.subplot_separation,
                            help='Subplot separation. ' +
                            f'Default: {self.subplot_separation}.')

        # Set the number of subplots per figure.
        parser.add_argument('-sp', '--subplotsperpage',
                            dest='subplots_per_fig', type=int,
                            default=self.subplots_per_fig,
                            nargs='+',help='Subplots per figure. ' +
                            f'Default: {self.subplots_per_fig}.')

        # Define the number of curves per subplot.
        parser.add_argument('-cs', '--curves', dest='curves_per_subplot',
                            type=int, default=self.curves_per_subplot,
                            nargs='+',help='Curves per subplot. ' +
                            f'Default: {self.curves_per_subplot}.')

        # Define the curve width.
        parser.add_argument('-ct', '--curvewidth', dest='curve_width',
                            type=float, default=self.curve_width,
                            nargs='+', help='Curve width. ' +
                            f'Default: {self.curve_width}.')
        parser.add_argument('-cta', '--curvewidthstart',
                            dest='curve_width_start', type=float,
                            default=self.curve_width_start,
                            help='Starting curve width. Default: ' +
                            f'{self.curve_width_start}.')
        parser.add_argument('-cte', '--curvewidthstep',
                            dest='curve_width_step', type=float,
                            default=self.curve_width_step,
                            help='Curve width step size. Default: ' +
                            f'{self.curve_width_step}.')
        parser.add_argument('-cti', '--curvewidthsize',
                            dest='curve_width_size', type=int,
                            default=self.curve_width_size,
                            help='Length of size variation list. Default: ' +
                            f'{self.curve_width_size}.')

        # Define the curve style.
        parser.add_argument('-cy', '--curvestyle', dest='curve_style',
                            type=self.list_of_strings,
                            default=self.curve_style,
                            help='Curve style. Note that -cy cannot take ' +
                            'shorthand input (e.g., -. or -). You need to' +
                            'use the long option or edit the rc file ' +
                            f'for that instead. Default: {self.curve_style}.')
        parser.add_argument('-cya', '--curvestylestart',
                            dest='curve_style_start', type=int,
                            default=self.curve_style_start, help='Starting ' +
                            'curve style index in curve_style.dat. ' +
                            f'Default: {self.curve_style_start}')
        parser.add_argument('-cye', '--curvestylestep',
                            dest='curve_style_step', type=int,
                            default=self.curve_style_step, help='Indices ' +
                            ' to step in curve_style.dat. ' +
                            f'Default: {self.curve_style_step}')
        parser.add_argument('-cyi', '--curvestylesize',
                            dest='curve_style_size', type=int,
                            default=self.curve_style_size, help='Length ' +
                            ' of the list to build from curve_styles.dat. ' +
                            f'Default: {self.curve_style_size}')

        # Define the curve color.
        parser.add_argument('-cc', '--curvecolor', dest='curve_color',
                            type=ascii, default=self.curve_color,
                            nargs='+',help='Curve color. ' +
                            f'Default: {self.curve_color}.')
        parser.add_argument('-cca', '--curvecolorstart',
                            dest='curve_color_start', type=int,
                            default=self.curve_color_start, help='Starting ' +
                            'curve color index in curve_color.dat. ' +
                            f'Default: {self.curve_color_start}')
        parser.add_argument('-cce', '--curvecolorstep',
                            dest='curve_color_step', type=int,
                            default=self.curve_color_step, help='Indices ' +
                            ' to step in curve_color.dat. ' +
                            f'Default: {self.curve_color_step}')
        parser.add_argument('-cci', '--curvecolorsize',
                            dest='curve_color_size', type=int,
                            default=self.curve_color_size, help='Length ' +
                            ' of the list to build from curve_colors.dat. ' +
                            f'Default: {self.curve_color_size}')

        # Define the curve mark.
        parser.add_argument('-cm', '--curvemark', dest='curve_mark',
                            type=ascii, default=self.curve_mark,
                            nargs='+',help='Curve mark. ' +
                            f'Default: {self.curve_mark}.')
        parser.add_argument('-cma', '--curvemarkstart',
                            dest='curve_mark_start', type=int,
                            default=self.curve_mark_start, help='Starting ' +
                            'curve mark index in curve_mark.dat. ' +
                            f'Default: {self.curve_mark_start}')
        parser.add_argument('-cme', '--curvemarkstep',
                            dest='curve_mark_step', type=int,
                            default=self.curve_mark_step, help='Indices ' +
                            ' to step in curve_mark.dat. ' +
                            f'Default: {self.curve_mark_step}')
        parser.add_argument('-cmi', '--curvemarksize',
                            dest='curve_mark_size', type=int,
                            default=self.curve_mark_size, help='Length ' +
                            ' of the list to build from curve_marks.dat. ' +
                            f'Default: {self.curve_mark_size}')

        # Define the curve separation.
        parser.add_argument('-cp', '--curvesep', dest='curve_separation',
                            type=float, default=self.curve_separation,
                            nargs='+',help='Curve separation. ' +
                            f'Default: {self.curve_separation}.')

        # Execute the argument parsing.
        args = parser.parse_args()
    
        # Return the results.
        return args


    def reconcile(self, args):

        # The args now holds either the parameter as given on the command line
        #   or the default value pulled from the RC file. Thus, we just copy
        #   the args into the self variables.
        self.display = args.display.strip("'")
        self.infile = args.infile.strip("'")
        self.outfile = args.outfile.strip("'")
        self.title = args.title.strip("'")
        self.print_command = args.print_command
        self.print_legend = args.print_legend
        self.page_orientation = args.page_orientation.strip("'")
        self.fig_format = args.fig_format.strip("'")
        self.fig_height = args.fig_height
        self.fig_width = args.fig_width
        self.fig_type = args.fig_type
        self.link_subplots = args.link_subplots
        self.multi_x_cols = args.multi_x_cols
        self.x_col = args.x_col
        self.y_col = args.y_col.strip("'")
        self.x_min = args.x_min
        self.x_max = args.x_max
        self.y_min = args.y_min
        self.y_max = args.y_max
        self.x_multiple = args.x_multiple
        self.y_multiple = args.y_multiple
        self.x_axis_inc = args.x_axis_inc
        self.x_axis_minor_ticks = args.x_axis_minor_ticks
        self.y_axis_major_ticks = args.y_axis_major_ticks
        self.y_axis_minor_ticks = args.y_axis_minor_ticks
        self.linked_x_axes = args.linked_x_axes
        self.linked_y_axes = args.linked_y_axes
        self.subplot_separation = args.subplot_separation
        self.subplots_per_fig = copy.deepcopy(args.subplots_per_fig)
        self.curves_per_subplot = copy.deepcopy(args.curves_per_subplot)
        self.curve_separation = copy.deepcopy(args.curve_separation)
        self.curve_width = copy.deepcopy(args.curve_width)
        self.curve_width_start = args.curve_width_start
        self.curve_width_step = args.curve_width_step
        self.curve_width_size = args.curve_width_size
        self.curve_style = copy.deepcopy(args.curve_style)
        self.curve_style_start = args.curve_style_start
        self.curve_style_step = args.curve_style_step
        self.curve_style_size = args.curve_style_size
        self.curve_color = copy.deepcopy(args.curve_color)
        self.curve_color_start = args.curve_color_start
        self.curve_color_step = args.curve_color_step
        self.curve_color_size = args.curve_color_size
        self.curve_mark = copy.deepcopy(args.curve_mark)
        self.curve_mark_start = args.curve_mark_start
        self.curve_mark_step = args.curve_mark_step
        self.curve_mark_size = args.curve_mark_size


    def process_settings(self):

        # The options defining the visuals for the lines may need to be
        #   processed. Example visual aspects include the line thickness,
        #   the line color, and the solid/dotted/dashed nature of the line.
        #   For each aspect, the command line could include an explicit
        #   list of the options, or it could include a "start", "step", and
        #   "size" that are used to construct the list of options.
        # So, if the explicit list is *not* given and thus empty by default,
        #   then we need to construct it using the start/step/size.

        # Do the curve width first because it is the easiest.
        if (len(self.curve_width) == 0):
            for i in range(self.curve_width_size):
                self.curve_width.append(self.curve_width_start +
                        self.curve_width_step * i)

        # The remaining visuals have data files that define the possible
        #   values. So, we need to be able to open those files from the
        #   environment variable OLCAO_DATA directory.
        OLCAO_DATA = os.getenv('OLCAO_DATA')

        # The curve styles must pull data from curve_styles.dat.
        with open(f"{OLCAO_DATA}/curve_styles.dat", "r") as dat:
            styles = dat.read().splitlines()
        styles = [eval(x) for x in styles]

        # Append styles from curve_styles.dat starting at the requested
        #   index and stepping by the given step while remaining within
        #   the array that contains all styles from curve_styles.dat.
        for x in range(self.curve_style_size):
            self.curve_style.append(styles[(self.curve_style_start +
                    self.curve_style_step * x) % len(styles)])

        # The curve colors must pull data from curve_colors.dat. Also,
        #   The data file contains additional #hex representations of the
        #   colors that we will ignore. We just want the fun color name.
        #   Fortunately, after the color name, there is always a tab
        #   character so we just take the slice of the line string from
        #   the first character to (but not including) the first tab.
        with open(f"{OLCAO_DATA}/curve_colors.dat", "r") as dat:
            colors = dat.read().splitlines()
        colors = [x[:x.find('\t')] for x in colors]

        # Append colors from curve_colors.dat starting at the requested
        #   index and stepping by the given step while remaining within
        #   the array that contains all colors from curve_colors.dat.
        for x in range(self.curve_color_size):
            self.curve_color.append('xkcd:' + colors[
                    (self.curve_color_start + self.curve_color_step * x)
                    % len(colors)])

        # The curve marks must pull data from curve_marks.dat.
        with open(f"{OLCAO_DATA}/curve_marks.dat", "r") as dat:
            marks = dat.read().splitlines()

        # Append marks from curve_marks.dat starting at the requested
        #   index and stepping by the given step while remaining within
        #   the array that contains all marks from curve_marks.dat.
        for x in range(self.curve_mark_size):
            self.curve_mark.append(marks[(self.curve_mark_start +
                    self.curve_mark_step * x) % len(marks)])


        # The color command line string may contain spaces and hence was
        #   enclosed in quotes that now need to be stripped out. Do that.
        self.curve_color = [e.strip("'") for e in self.curve_color]

        # Once the command line is parsed, if any axis limit (x,y ; min,max)
        #   is recognized as m.inf then we need to set a flag indicating that
        #   this axis limit should be determined automatically. We start by
        #   assuming that none of the axis limits were specified.
        self.x_min_exists = False
        self.x_max_exists = False
        self.y_min_exists = False
        self.y_max_exists = False
        if (self.x_min != m.inf):
            self.x_min_exists = True
        if (self.x_max != m.inf):
            self.x_max_exists = True
        if (self.y_min != m.inf):
            self.y_min_exists = True
        if (self.y_max != m.inf):
            self.y_max_exists = True

        # Internally, it is easier for multi_x_cols to be an integer 1 for
        #   False, and integer 2 for True.
        if (self.multi_x_cols == False):
            self.multi_x_cols = 1
        else:
            self.multi_x_cols = 2

        # If the outfile contains an "X" as the first character, it will be
        #   replaced with the display mode to clarify the purpose of the
        #   script.
        if (self.outfile == "Xscript.py"):
            self.outfile = self.outfile.replace("X", self.display, 1)


    def recordCLP(self):
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write(f"Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


# Convert a string of space separated integers mixed with #..# style
#   expressions into a list of integers where the #..# style expressions
#   "fill in" the numbers between the first and last. E.g., 3..5 becomes
#   a string of the form "3 4 5".
def string_2_list(in_string):

    # Get each space separated element in the given in_string.
    string_terms = in_string.strip("'").split()

    # Initialize an empty string that will hold the final list
    #   of all integers.
    expanded_string_terms = ""

    # Consider each term.
    for term in (string_terms):

        # If the term contains a "..", then we need to expand the
        #   expression into a sequence of integers.
        if (".." in term):

            # Find the beginning and ending limits of the range and
            #   store as integers.
            range_limits = term.split("..")

            # Create a range of integers and then convert to a list of
            #   strings (list comprehension), and finally, join them.
            expanded_string_terms += " " + " ".join([str(x) for x in
                    range(int(range_limits[0]), int(range_limits[1])+1)])
        else:
            # Otherwise, just add the term to the list.
            expanded_string_terms += " " + term

    # Remove the leading space.
    expanded_string_terms = expanded_string_terms[1:]

    # Convert the expanded string into a list of integers.
    out_list = [int(e) for e in expanded_string_terms.split()]

    # Return the list we created.
    return out_list


def read_data_headers(settings):

    # The columns headers are usually stored in the first line of the data
    #   file. The exception is symmteric band structure plot (sybd). In that
    #   case the first line is just data.
    # In any case, we need to know the number of columns so we get the first
    #   line as if we were going to use the whole thing as headers.
    with open(settings.infile, "r") as f:
        all_headers = f.readline().split()

    # If we are doing sybd, then we redefine the headers with index numbers.
    #   Note that the last two columns of the input for sybd contain special
    #   information for the high symmetry k-points. We will also prepend
    #   "Energy" as the header for the first column.
    if (settings.fig_type == "sybd"):
        all_headers = [str(x+1) for x in range(len(all_headers)-3)]
        all_headers.insert(0, "Energy")

    # Now that we know the number of columns, we can expand the y_col
    #   variable into a list of integers. In the special case where the
    #   y_col was not provided (so that the default of "all columns" is
    #   used, we simply store the relevant columns.
    if (not settings.y_col):
        if (settings.multi_x_cols == 2):  # True
            # Get all even numbered columns from two onward.
            settings.y_col = [x for x in range(2, len(all_headers)+1, 2)]
        else:
            # Get all columns from two onward.
            settings.y_col = [x for x in range(2, len(all_headers)+1)]
    else:
        settings.y_col = string_2_list(settings.y_col)

    # Extract only the headers from the requested y_cols. Note that the
    #   index number used to refer to each y_col starts counting from 1
    #   as the first column in the input file. Thus, for a simple two-
    #   column x-y data set, the y_col value is 2 and consequently we
    #   need to subtract one for the python array index below.
    y_col_headers = [all_headers[col-1] for col in settings.y_col]

    # For each y_col, extract the associated x_col header. Note, for the
    #   same reason as described above, we will subtract 2 when getting
    #   the x_col values.
    if (settings.multi_x_cols == 2):  # True
        x_col_headers = [all_headers[col-2] for col in settings.y_col]
    else:
        x_col_headers = [all_headers[0] for col in settings.y_col]


    return (x_col_headers, y_col_headers)


def print_veusz_header(settings):
    with open(settings.outfile, "w") as s:
        s.write("""#!/usr/bin/env python3
""")


def print_mpl_header(settings):
    with open(settings.outfile, "w") as s:
        s.write("""#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
"""
        )

        if (settings.fig_type == "sybd"):
            # Define the column names for all the SYBD plot data.
            last_idx = settings.curves_per_subplot[0]
            s.write(f"col_names = [str(x+1) for x in range({last_idx})]\n")
            s.write("col_names.insert(0, 'Energy')\n")
            s.write("col_names.append('sykp_value')\n")
            s.write("col_names.append('sykp_name')\n")

            # Read in the SYBD data and meta data in one shot.
            s.write("""

data = pd.read_csv(""" + f"'{settings.infile}', header = None,\n" +
"""                   names = col_names, usecols = col_names, sep='\\\\s+')\n\n"""
            )
            
            # Extract out the SYBD meta data.
            s.write(f"meta = data.iloc[:, {last_idx+1}:{last_idx+3}].dropna()\n")

            # Extract out the SYBD data.
            s.write(f"data = data.iloc[:, 0:{last_idx+1}]\n")

        else:
            s.write("""

data = pd.read_csv(""" + f"'{settings.infile}'" + """, sep='\\\\s+')

"""
            )

        # Create the arrays that hold the figures, subplots, axes, and curves.
        s.write("""# Arrays to hold figures, subplots, axes, and curves.
figs = []
subplots = []
curves = []
"""
                )



def adjust_for_fig_type(settings):
    if (settings.fig_type == "sybd"):
        settings.curves_per_subplot = [len(settings.y_col)]


def compute_fig_size(settings):

    # Make short-hand copies of the page orientation and fig format.
    page_orientation = settings.page_orientation
    fig_format = settings.fig_format
    fig_height = settings.fig_height
    fig_width = settings.fig_width

    # Define an array to hold the page_size.
    page_size = []

    # Compute the page size from the supplied format.
    if (fig_format == "letter"):
        page_size = [8.5, 11.0]  # Portrait by default.
    elif (fig_format == "A4"):
        page_size = [8.27, 11.69]  # Portrait by default.
    else:
        print(f"Unknown fig_format: {fig_format}.")
        quit()

    # Change the orientation of the page as desired.
    if (page_orientation == "landscape"):
        page_size = [page_size[1], page_size[0]]
    elif (page_orientation != "portrait"):
        print(f"Unknown page_orientation: {page_orientation}.")
        quit()

    # Compute the figure size from the page size and figure use %.
    return [page_size[0] * fig_height / 100.0,
            page_size[1] * fig_width / 100.0]


def draw_sybd_verticals(s):
    s.write(f"plt.vlines(meta['sykp_value'], y_min, y_max)\n")
    #s.write(f"for idx, row in meta.iterrows():")


def modify_sybd_horz_axis(s):
    s.write("plt.tick_params(axis='x', which='both', bottom=False, " +
            "labelbottom=False)\n")
    s.write("for (val, name) in zip(meta['sykp_value'], meta['sykp_name']):\n")
    s.write("    plt.annotate(name, (val, y_min), xytext=(-0.5, -1.1), " +
            "textcoords='offset fontsize')\n")


def apply_x_limits(settings, s, wrap_subplot, curve):
    if (settings.linked_x_axes == False):
        if (settings.x_min_exists == False or settings.x_max_exists == False):
            # See below for the y_min, y_max description. This one is almost
            #   the same except that we need to account for the multi_x_cols
            #   flag. For only one x_col, we always use the first column. For
            #   multiple x columns, we need to build a list of them.
            if (settings.multi_x_cols == 1):
                curr_cols = [0]
            else:
                curr_cols = []
                for idx in range(settings.curves_per_subplot[wrap_subplot]):
                    curr_cols.append(settings.x_col[curve+idx]-1)
            s.write(f"data_temp = data.iloc[:, {curr_cols}]\n")


        if (settings.x_min_exists == False):
            s.write(f"x_min = data_temp.min()\n")
            s.write(f"if (type(x_min) is not float and type(x_min) is not int):\n")
            s.write(f"    x_min = x_min.min()\n")
            s.write(f"if (x_min != 0):\n")
            s.write(f"    x_min = (int(x_min/{settings.x_multiple}) - 1) " +
                    f" * {settings.x_multiple}\n")
        else:
            s.write(f"x_min = {settings.x_min}\n")

        if (settings.x_max_exists == False):
            s.write(f"x_max = data_temp.max()\n")
            s.write(f"if (type(x_max) is not float and type(x_max) is not int):\n")
            s.write(f"    x_max = x_max.max()\n")
            if (settings.fig_type != "sybd"):
                s.write(f"if (x_max != 0):\n")
                s.write(f"    x_max = (int(x_max/{settings.x_multiple}) + 1) " +
                        f"* {settings.x_multiple}\n")
        else:
            s.write(f"x_max = {settings.x_max}\n")

        # Use the computed or assigned x_min and x_max.
        s.write(f"plt.xlim((x_min, x_max))\n")


def apply_y_limits(settings, s, wrap_subplot, curve):
    if (settings.linked_y_axes == False):
        if (settings.y_min_exists == False or settings.y_max_exists == False):
            # We have to find the maximum and/or the minimum y value in the
            #   data set for this subplot.

            # Construct a dataframe of the columns that will be included in
            #   this subplot. Start with the current curve and end with the
            #   last curve that will be added to this subplot. The last curve
            #   is determined by either being the last curve in this subplot's
            #   list of curves_per_subplot or by being the last curve in the
            #   whole data set (i.e., we request curve_per_subplot number of
            #   curves, but the data set only has a few remaining columns
            #   before we reach the end of the data set.)
            curr_cols = []
            for idx in range(settings.curves_per_subplot[wrap_subplot]):
                curr_cols.append(settings.y_col[curve+idx]-1)
            s.write(f"data_temp = data.iloc[:, {curr_cols}]\n")

        # For each of the y_min and y_max, if no value was specified, then
        #   get the min/max of the above constructed dataset and then extend
        #   it to align with the desired y_multiple. That will make the
        #   subplot on a "round" number, i.e., the y_multiple.
        # If a y_min or y_max was specified, then just assign it. 
        if (settings.y_min_exists == False):
            s.write(f"y_min = data_temp.min()\n")
            s.write(f"if (type(y_min) is not float and type(y_min) is not int):\n")
            s.write(f"    y_min = y_min.min()\n")
            s.write(f"if(y_min != 0):\n")
            s.write(f"    y_min = (int(y_min/{settings.y_multiple}) - 1) " +
                    f" * {settings.y_multiple}\n")
        else:
            s.write(f"y_min = {settings.y_min}\n")

        if (settings.y_max_exists == False):
            s.write(f"y_max = data_temp.max()\n")
            s.write(f"if (type(y_max) is not float and type(y_max) is not int):\n")
            s.write(f"    y_max = y_max.max()\n")
            s.write(f"if(y_max != 0):\n")
            s.write(f"    y_max = (int(y_max/{settings.y_multiple}) + 1) " +
                    f"* {settings.y_multiple}\n")
        else:
            s.write(f"y_max = {settings.y_max}\n")

        # Use the computed or assigned y_min and y_max.
        s.write(f"plt.ylim((y_min, y_max))\n")


def print_figs_subplots_curves(settings, x_col_headers, y_col_headers):

    # To create our set of figures, their subplots, and their curves, we need
    #   to consider a few different things. We assume that by this time, the
    #   list of columns to use has been prepared. Now, we just go through the
    #   columns and add curves, subplots, and figures as needed.

    # Key variables to understand are:
    make_new_fig = False
    make_new_subplot = False
    curr_subplot_curve = 1
    curr_subplot = 1
    wrap_subplot = 0
    wrap_fig = 0

    # Initialize counts of the total # of figures, subplots, and curves.
    num_total_figs = 0
    num_total_subplots = 0
    num_total_curves = 0

    # curves_per_subplot: This list defines the number of curves to allow
    #   for a series of subplots. If the list has a length of three, with
    #   values of 4, 1, 2, then the first subplot will have four curves, the
    #   second subplot will have one curve, and the third subplot will have
    #   two curves. If there are more columns of data, then the fourth
    #   subplot will wrap back around and expect four curves again.
    # wrap_subplot: This is an integer that equals an index number for the
    #   curves_per_subplot list. The value in the list at the wrap_subplot
    #   index says how many curves our current subplot can hold. But, the main
    #   purpose of the wrap_subplot variable is to keep track of which subplot
    #   within the list of subplots we are currently working through (i.e.,
    #   adding curves to) so that when we go to the next subplot we can figure
    #   out how many curves the subplot will have and if we should wrap back
    #   to the beginning of the curves_per_subplot list.
    # num_subplots_reset: This is the length of the curves_per_subplot list.
    #   When the wrap_subplot index number is equal to this number, then we
    #   should make wrap_subplot go back to list index zero.
    # curr_subplot_curve: This integer counts the number of curves that we have
    #   added to whatever current subplot we are working on. (I.e., it counts
    #   the number of curves that have been added to the subplot that is
    #   currently indexed by wrap_subplot.) This number goes from one to the
    #   number stored in curves_per_subplot[wrap_subplot].

    # subplots_per_fig: This is just like curves_per_subplot, except that it
    #   applies to subplots and figures respectively. This list says how many
    #   subplots should appear in a sequence of figures. I.e., a list of
    #   4, 1, 2 says that the first figure has four subplots, the second
    #   figure has one subplot, and the third figure as two subplots. If
    #   another figure is added beyong the third, then we wrap back around and
    #   this figure will have four subplots.
    # wrap_fig: This is an integer that equals an index number for the
    #   subplots_per_fig list. The value in the list at the wrap_fig index
    #   says how many subplots this figure will hold. But, the main purpose
    #   of the wrap_fig variable is to keep track of which figure within the
    #   list we are currently working through (i.e., adding subplots to) so
    #   that when we go to the next figure, we can determine how many subplots
    #   the figure will have and if we should wrap back to the beginning of
    #   the subplots_per_fig list.
    # num_figs_reset: This is the length of the subplots_per_fig list. When
    #   the wrap_fig index number is equal to this number, then we should make
    #   wrap_fig go back to the list index zero.
    # curr_subplot: This integer counts the number of subplots that we have
    #   added to whatever current figure we are working on. (I.e., it counts
    #   the number of subplots that have been added to the figure that is
    #   currently indexed by wrap_fig.) This number goes from one to the
    #   number stored in subplots_per_fig[wrap_fig].

    # Compute the figure size compared to the page.
    fig_size = compute_fig_size(settings)

    # Define convenient shorthand variables for the number of curves to plot,
    #   and the track when the wrap variables should be reset back to the zero
    #   index.
    num_curves = len(settings.y_col)
    num_subplots_reset = len(settings.curves_per_subplot)
    num_figs_reset = len(settings.subplots_per_fig)

    # Reopen the file that the script will be written to using the append
    #   option.
    s = open(settings.outfile, "a")

    # Consider each requested curve and plot it in the appropriate figure
    #   and subplot. (Make new figures and new subplots as needed.)
    for curve in range(num_curves):

        # Manage the special case of the first figure, subplot, and curve.
        if (curve == 0):
            make_new_fig = True  # Trigger the creation of a figure
            make_new_subplot = True  # Trigger the creation of a subplot
            curr_subplot_curve = 1  # Start counting curves within this subplot
            curr_subplot = 1  # ID the newly created subplot
            wrap_subplot = 0  # List index for num. curves per subplot
            wrap_fig = 0  # List index for num. subplots per figure
        else:
            # Assume that we do not need to make a new figure or subplot
            make_new_fig = False
            make_new_subplot = False

            # If we have already added all the curves to this subplot that
            #   should be added, then we'll need to start making a new subplot.
            if (curr_subplot_curve ==
                    settings.curves_per_subplot[wrap_subplot]):

                # Add a legend to the previous subplot if needed.
                if (settings.print_legend == True):
                    s.write(f"plt.legend()\n")

                # Update indices for tracking progress.
                make_new_subplot = True  # Trigger creation of a new subplot
                curr_subplot_curve = 1  # Reset count of curves for new subplot
                curr_subplot += 1  # Increment subplot ID in this figure
                wrap_subplot += 1  # Inc. list index for num. curves per subp.

                # In the event that we have reached the end of the list of the
                #   number of subplots per figure, then we need to go back to
                #   the beginning of the list.
                if (wrap_subplot == num_subplots_reset):
                    wrap_subplot = 0

                # If we have already added all the subplots to this figure
                #   that should be added, then we'll need to start making
                #   a new figure.
                if (curr_subplot > settings.subplots_per_fig[wrap_fig]):
                    make_new_fig = True
                    curr_subplot = 1
                    wrap_fig += 1
                    if (wrap_fig == num_figs_reset):
                        wrap_fig = 0
            else:
                # Increment the count of the number of curves in this subplot.
                curr_subplot_curve += 1

        if (make_new_subplot == True):

            # Before we create a new figure and/or subplot, we need to refine
            #   the axes according to the command line input and/or default
            #   values.


            # Increment count of the total number of subplots in all figs.
            num_total_subplots += 1

            # If this subplot needs to be placed in a new figure, then do it.
            if (make_new_fig == True):

                # Increment the count of the total number of figures.
                num_total_figs += 1

                # Create the figure and append it to a list of figures.
                s.write(f"figs.append(plt.figure({num_total_figs}, " +
                        f"figsize={fig_size}))\n")

                # Create the subplot entity (that actually contains all
                #   subplots) in the figure. Note, this also "selects" the
                #   index #1 position of the subplot as the next subplot
                #   that will be filled with a "plot" method in pyplot.
                s.write(f"subplots.append(plt.subplot(" +
                        f"{settings.subplots_per_fig[wrap_subplot]}," +
                        " 1, 1))\n")

            else:
                # Select the next subplot index of the subplot.
                s.write(f"subplots.append(plt.subplot(" +
                        f"{settings.subplots_per_fig[wrap_subplot]}, 1, " +
                        f"{curr_subplot}), )\n")

            # Now that we have created the subplot, we can apply the limits
            #   to the axes.
            apply_x_limits(settings, s, wrap_subplot, curve)
            apply_y_limits(settings, s, wrap_subplot, curve)

        # Now, we just add the next curve to whatever subplot we selected.
        #   We also pre-compute the indices at which the color, style, and
        #   line width values should be found in their respective lists.
        color_index = (curr_subplot_curve - 1)%len(settings.curve_color)
        mark_index = (curr_subplot_curve - 1)%len(settings.curve_mark)
        style_index = (curr_subplot_curve - 1)%len(settings.curve_style)
        width_index = (curr_subplot_curve -1)%len(settings.curve_width)
        s.write("curves.append(plt.plot(" +
                f"data['{x_col_headers[curve]}'], " +
                f"data['{y_col_headers[curve]}'], " +
                f"label='{y_col_headers[curve]}', " +
                f"color='{settings.curve_color[color_index]}', " +
                f"marker='{settings.curve_mark[mark_index]}', " +
                f"linestyle={settings.curve_style[style_index]}, " +
                f"linewidth={settings.curve_width[width_index]}))\n")

    # In the event that a SYBD plot was created, we now need to do a few
    #   things to improve the visuals. (1) Add the vertical lines using the
    #   meta data. (2) Modify the labels.
    if (settings.fig_type == "sybd"):
        draw_sybd_verticals(s)
        modify_sybd_horz_axis(s)

    # Add a legend to the previous subplot if needed.
    if (settings.print_legend == True):
        s.write(f"plt.legend()\n")

    # Show the plot(s)
    s.write("plt.show()\n")

    # Close the file.
    s.close()

    # Final script activities.
    if (settings.display == "mpl"):
        # Set permissions to rwx for the user and rx for the group.
        os.chmod(settings.outfile, 0o750)

        # Execute the script.
        subprocess.run([f"./{settings.outfile}"])

    
def main():

    # Define and initialize all the global variables.
    num_lines = 0

    # Set default values from the user resource control file.
    settings = ScriptSettings()

    # Get the header for each requested column.
    (x_col_headers, y_col_headers) = read_data_headers(settings)

    # Depending on the type of plot some automatic adjustments should be made
    #   to the figure.
    adjust_for_fig_type(settings)

    # If the user wants a veusz visual then make the header for it.
    #   Otherwise, make a matplotlib based header.
    if (settings.display == "mpl"):
        print_mpl_header(settings)
    elif (settings.display == "veusz"):
        print_veusz_header(settings)
    else:
        print ("Unknown display style. Should not get here.")
        exit()

    # Compute the number of lines and plots to make.
    print_figs_subplots_curves(settings, x_col_headers, y_col_headers)


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program.
    main()

    #from cProfile import Profile
    #from pstats import SortKey, Stats

    #with Profile() as profile:
    #    print(f"{main() = }")
    #    (
    #        Stats(profile)
    #        .strip_dirs()
    #        .sort_stats(SortKey.CALLS)
    #        .print_stats()
    #    )
