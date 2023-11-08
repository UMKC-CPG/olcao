#!/usr/bin/env python3

import argparse as ap
import pandas as pd
import os
import sys

def print_help():
    output = """
#Script to plot data v0.5  Paul Rulis
#Last update Mar. 04, 2022
######################################################################
#This python script will execute the necessary steps to plot a graph to
#	the screen using either matplotlib or veusz. Note that plot default
#   parameters are define by the file "$OLCAO_DIR/.olcao/plotgraphrc".
#
#---------------------------------------------------------------------
#
#USAGE:  plotgraph [-x n] [-y m1 [m2..m3] ... ] [-xrange xmin xmax]
#		   [-yrange ymin ymax] [-c] [-noef] [-eos] [-sybd] [-li]
#		   [-iter] [-pdos] [-title "TITLE"] [-mxtics x] [-mytics y]
#          [-points] [-shift shift] [-log] [-v]
#          filename1 [filename2]
#OR:     plotgraph -help
#
#filename1 is the name of the file in the current directory to be
#    plotted. If filename2 is given then it is assumed you are
#    plotting a spin-polarized PDOS. In this case filename1 is
#    taken to be the UP-SPIN and filename2 is the DOWN-SPIN. They
#    will be plotted on the top/bottom respectivly. Whatever
#    settings you provide in terms of range and data columns to
#    plot will be applied to both files.
#If the "-x" option is given then n is the column number to be used as
#    the x-axis in the graph. If the "-y" option is given then m
#    is the column number to be used as the y-axis in the graph.
#    You can optionally give an "m2" or an "m2..m3" which will cause
#    the program to also plot these columns or range of columns. It
#    is assumed that the data is formatted in space seperated columns
#    of data. If nether "-x" or "-y" are given then it is assumed
#    that the first column will be the x-axis and that the second
#    column will be the y-axis. If a "0" is given for the -y option then
#    every single column from 1 to the number of columns will be plotted.
#The -xrange parameter is pretty self explanitory.  If provided it will
#    override the default behavior for deciding the range of the x axis to
#    plot. In the case of pdos and localization index (li) plots the default
#    is -25 to 25. The default behavior otherwise is to use autoscaling.
#The -yrange is also pretty self explanitory.
#    If it isn't given then the data files will be searched for
#    the highest y value and the y-range will be 10% greater than
#    that. The exception is if the -eos option is given. In that
#    case the y range is computed from the 2nd input file. The other
#    exception is for sybd plots. In that case the default is -20 to 20.
#The -noef option will turn off the drawing of a straight line at the fermi
#    energy for pdos and sybd plots. Otherwise, the plotting program will
#    assume that you wish to have a straight line at the 0 mark. This
#    represents the Fermi Energy.
#The -eos command tells the script that you want to plot energy
#    vs. volume.  This requires two filenames.  One contains
#    the data points for the curve made by the eos program, and
#    the other contains the data points given to the eos program.
#    They must be given in that order.  If no filenames are given
#    then the following files are assumed:  eos.plot eos.points.
#The -title option followed by a quoted string will be printed as a graph
#    title in place of the usual "Printed by $USER"
#The -c option will prevent the command used to call plotgraph from being
#    printed as a part of the title.
#The -mxtics switch followed by an integer number will produce minor tic
#    marks on the x axis.  The number integer number is the number of
#    divisions, NOT the number of tics.  The default value is 10.  To
#    not have any minor tic marks just give a value of 1.
#The -mytics behaves exactly the same as -mxtics except it applies to the
#    y axis instead.
#The -points option will make the printout use little dots for each data
#    point instead of connecting the data points with lines.
#The -shift option takes one argument which is the amount to shift each line
#    by in a multi-line plot.
#The -log option will turn on a log plot for the y axis.
#The -v option will construct a script that can be used in veusz. Otherwise
#    the program will create and immediately show a plot using matplotlib.
#If the "-help" switch is given then the program will print some (hopefully)
#    useful information about how to use it and then exit.
#
#
#---------------------------------------------------------------------
#
#REQUIREMENTS
#
#You need to be in the same directory as the file you wish to graph.
#	Maybe in the future it will also accept the data to be stored
#	in a subdirectory such as DOS for pdos plots.
#
"""
    print (output)


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
        self.reconcile(args, default_rc)

    def assign_rc_defaults(self, default_rc):

        # Supporting display information
        self.title = default_rc[0]
        self.print_command = default_rc[1]

        # Plot formats
        self.page_height = default_rc[2]
        self.page_width = default_rc[3]
        self.link_plots = default_rc[4]
        self.plot_type = default_rc[5]

        # Data columns
        self.multi_x_cols = default_rc[6]
        self.x_col = default_rc[7]
        self.y_col = default_rc[8]

        # Ranges and ticks
        self.x_min = default_rc[9]
        self.x_max = default_rc[10]
        self.y_min = default_rc[11]
        self.y_max = default_rc[12]
        self.x_axis_inc = default_rc[13]
        self.x_axis_minor_ticks = default_rc[14]
        self.y_axis_major_ticks = default_rc[15]
        self.y_axis_minor_ticks = default_rc[16]

        # Subplots and lines variables
        self.plots_per_page = default_rc[17]
        self.plot_separation = default_rc[18]
        self.lines_per_plot = default_rc[19]
        self.line_thickness = default_rc[20]
        self.line_style = default_rc[21]
        self.line_color = default_rc[22]
        self.line_separation = default_rc[23]

        #    global dos = default_rc[24]
        #    global points = default_rc[25]


    def parse_command_line(self, default_rc):
    
        # Create the parser tool.
        parser = ap.ArgumentParser(description='Control parameters')
    
        # Add arguments to the parser.
        self.add_parser_arguments(parser, default_rc)

        # Parse the arguments and return the results.
        return parser.parse_args()

    def add_parser_arguments(self, parser, default_rc):
    
        # Define the title argument.
        parser.add_argument('-t', '--title', dest='title', type=ascii,
                            default=self.title,
                            help='Title for the plot. Default: ' +
                            f'{self.title}')
    
        # Define the flag to print the command used to create the figure.
        if (self.print_command == True):
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-pc', '--printcmd', action=store_action,
                            dest='print_command', default=self.print_command,
                            help='Flag to print the command line. Default: ' +
                            f'{self.print_command}')
    
        # Define the page height.
        parser.add_argument('-ph', '--pageh', dest='page_height', type=float,
                            default=self.page_height, help='Page height. ' +
                            f'Default: {self.page_height}')
    
        # Define the page width.
        parser.add_argument('-pw', '--pagew', dest='page_width', type=float,
                            default=self.page_width, help='Page width. ' +
                            f'Default: {self.page_width}')
    
        # Define the flag to request linked plots.
        if (self.link_plots == True):
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-l', '--link', action=store_action,
                            dest='link_plots', default=self.link_plots,
                            help='Flag to request linked plots. Default: ' +
                            f'{self.link_plots}')
    
        # Define the plot type.
        parser.add_argument('-pt', '--plott', dest='plot_type', type=int,
                            choices=[0,1,2], default=self.plot_type,
                            help='Flag to request a specific plot type. ' +
                            'Options: 0=generic; 1=pdos; 2=optc. Default: ' +
                            f'{self.plot_type}')
    
        # Define the flag to indicate that the input has multiple x columns.
        if (self.multi_x_cols == True):
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-mx', '--multix', action=store_action,
                            dest='multi_x_cols', default=self.multi_x_cols,
                            help='Flag to indicate that input has multiple ' +
                            f'x columns. Default: {self.multi_x_cols}')
    
        # Select the x_column.
        parser.add_argument('-xc', '--xcol', dest='x_col', type=int,
                            default=self.x_col, help='X column. Default: ' +
                            f'{self.x_col}')
    
        # Select the y_column(s).
        parser.add_argument('-yc', '--ycol', dest='y_col_string', type=ascii,
                            default=self.y_col, help='Y column(s) in array ' +
                            f'form with .. options. Default: {self.y_col}')
    
        # Set the x_min
        if (self.x_min == ""):
            x_min_default = "Data min"
        else:
            x_min_default = self.x_min
        parser.add_argument('-xi', '--xmin', dest='x_min', type=float,
                            default=self.x_min, help='X min value. ' +
                            f'Default: {x_min_default}')

        # Set the x_max
        if (self.x_max == ""):
            x_max_default = "Data max"
        else:
            x_max_default = self.x_max
        parser.add_argument('-xa', '--xmax', dest='x_max', type=float,
                            default=self.x_max, help='X max value. ' +
                            f'Default: {x_max_default}')

        # Set the y_min
        if (self.y_min == ""):
            y_min_default = "Data min"
        else:
            y_min_default = self.y_min
        parser.add_argument('-yi', '--ymin', dest='y_min', type=float,
                            default=self.y_min, help='Y min value. ' +
                            f'Default: {y_min_default}')

        # Set the y_max
        if (self.y_max == ""):
            y_max_default = "Data max"
        else:
            y_max_default = self.y_max
        parser.add_argument('-ya', '--ymax', dest='y_max', type=float,
                            default=self.y_max, help='Y max value. ' +
                            f'Default: {y_max_default}')

        # Set the default x axis increment.
        parser.add_argument('-xai', '--xaxisinc', dest='x_axis_inc',
                            type=int, default=self.x_axis_inc, help='X ' +
                            f'axis increment. Default: {self.x_axis_inc}')

        # Set the default x axis minor ticks.
        parser.add_argument('-xm', '--xminor', dest='x_axis_minor_ticks',
                            type=int, default=self.x_axis_minor_ticks,
                            help='X axis increment. Default: ' +
                            f'{self.x_axis_minor_ticks}')

        # Set the default y axis major ticks.
        parser.add_argument('-yma', '--ymajor', dest='y_axis_major_ticks',
                            type=int, default=self.y_axis_major_ticks,
                            help='Y axis major ticks. Default: ' +
                            f'{self.y_axis_major_ticks}')

        # Set the default y axis minor ticks.
        parser.add_argument('-ymi', '--yminor', dest='y_axis_minor_ticks',
                            type=int, default=self.y_axis_minor_ticks,
                            help='Y axis minor ticks. Default: ' +
                            f'{self.y_axis_minor_ticks}')

        # Set the number of plots per page.
        parser.add_argument('-pp', '--plotsperpage', dest='plots_per_page',
                            type=ascii, default=self.plots_per_page,
                            help='Plots per page. ' +
                            f'Default: {self.plots_per_page}')

        # Define the plot separation.
        parser.add_argument('-ps', '--plotsep', dest='plot_separation',
                            type=ascii, default=self.plot_separation,
                            help='Plot separation. ' +
                            f'Default: {self.plot_separation}')

        # Define the number of lines per plot.
        parser.add_argument('-lpp', '--lines', dest='lines_per_plot',
                            type=ascii, default=self.lines_per_plot,
                            help='Lines per plot. ' +
                            f'Default: {self.lines_per_plot}')

        # Define the line thickness.
        parser.add_argument('-lt', '--linethick', dest='line_thickness',
                            type=ascii, default=self.line_thickness,
                            help='Line thickness. ' +
                            f'Default: {self.line_thickness}')

        # Define the line style.
        parser.add_argument('-ls', '--linestyle', dest='line_style',
                            type=ascii, default=self.line_style,
                            help='Line style. ' +
                            f'Default: {self.line_style}')

        # Define the line color.
        parser.add_argument('-lc', '--linecolor', dest='line_color',
                            type=ascii, default=self.line_color,
                            help='Line color. ' +
                            f'Default: {self.line_color}')

        # Define the line separation.
        parser.add_argument('-lp', '--linesep', dest='line_separation',
                            type=ascii, default=self.line_separation,
                            help='Line separation. ' +
                            f'Default: {self.line_separation}')

        # Define the positional arguments.
        parser.add_argument('filename', required=True, 
                            help='Name of the file to plot.')

        # Execute the argument parsing.
        args = parser.parse_args()
    
        # Return the results.
        return(args)

#    def reconcile(self, args, default_rc):
#        # If the argument *was* given on the command line, then we use it.
#        #   We determine if it was given on the command line by the fact that
#        #   the *args* 
#    
#        # Return the settings that incorporate the command line arguments.
#        return settings


def read_data():
    return pd.read_csv(args.filename1, sep='\s+')



def compute_num_lines_and_plots():

    # We need to produce a number (num_lines) that is equal to the total number
    #   of lines that will be plotted based on the given data. We also must
    #   produce a number (num_plots) that is equal to the total number of plots
    #   that will be created based on the given data or the command line. Note
    #   that a "plot" is an x-y axis with some number of lines. There may be
    #   multiple plots per page and multiple lines per plot.

    if (args.multi_x_cols):
        num_cols = 0  # The x axis for each plot will be counted separately.
    else:
        num_cols = 1  # One x axis for all plots so we count it at the outset.

    # Initialize the count of the number of lines and plots.
    num_lines = 0
    num_plots = 0

    wrap_plot = 1  # Plot number from the set {1..num_plots_reset}
    wrap_page = 1  # Page number from the set {1..num_pages_reset}
    curr_line = 1  # Line number from the set {1..lines_per_plot[wrap_plot]}
    curr_plot = 1  # Plot number from the set {1..plots_per_page[wrap_page]}

    num_pages_reset = size(plots_per_page)
    num_plots_reset = size(lines_per_plot)

    while(True):

        # Add a new line to the counter for this plot and the total counter.
        curr_line += 1
        num_lines += 1

        # Count that we are at the next column.
        num_cols += 1

        # Check if we passed the number of lines for this plot. If so, then we
        #   need to go to the next plot in the lines_per_plot array
        #   (wrap_plot) and we need to go to the next plot in the
        #   plots_per_page for the current page (curr_plot), and we will also
        #   need to increment the total number of plots (num_plots). Of course
        #   we will also have to reset the line counter for the next plot
        #   (curr_line).
        # We also need to check if a "x" column for this plot existed so that
        #   we can track the number of colums looked at so far so as to make
        #   sure that we don't go over the total number of columns in the data
        #   set (or that if we do, we can count back the right number of lines
        #   etc.) Also, we must check if we have gone through all the plots in
        #   the set {1..num_plots_reset}.
        if (curr_line > lines_per_plot[wrap_plot]):
            wrap_plot += 1
            curr_plot += 1
            num_plots += 1
            curr_line = 1
            if (args.multi_x_cols):
                num_cols += 1
            if (wrap_plot > num_plots_reset):
                wrap_plot = 1

            # Check if we passed the number of plots for this page. If so,
            #   then we need to go to the next page in the plots_per_page
            #   array (wrap_page). We will also need to reset the plot counter
            #   for the next page (curr_plot).
            if (curr_plot > plots_per_page[wrap_page]):
                wrap_page += 1
                curr_plot = 1
                if (wrap_page > num_pages_reset):
                    wrap_page = 1

        # Check if we are done by comparing the current count of the number
        #   of columns to the number of columns in the data frame (obtained
        #   by looking at the data frame shape [rows,cols]).
        if (num_cols >= data.shape[1]):

            # In the special case where someone "accidentally" asks for one
            #   plot with lots of lines in it but where the data file only
            #   has a few columns worth of data in it (i.e., less data than
            #   the number of requested lines), then the num_plots variable
            #   will not be set correctly and so that needs to be fixed now.
            if (num_plots == 0):
                num_plots = 1

                # Additionally, because the num_plots variable was never
                #   incremented, the multi_x_cols case will not have the
                #   correct column count in the same num_plot == 0 situation.
                if (args.multi_x_cols):
                    num_cols += 1

            # Now we can break because we are done and all columns have been
            #   counted.
            break


def main():

    # Set default values from the user resource control file.
    settings = ScriptSettings()

    # If the user asked for help then print it and exit.
    if (args.help_flag):
        print_help()
        exit()

    # Read the data into memory. (Presently, the script reads everything into
    #   memory all at once. In the future it may read only the necessary
    #   data columns.)
    data = read_data()

    # Compute the number of lines and plots to make.
    compute_num_lines_and_plots(settings)

    # If the user wants a veusz visual then make it. Otherwise, just show a
    #   matplotlib based graph.
    if (args.veusz_flag):
        print_veusz()
    else:
        show_mpl()


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program.
    main()
