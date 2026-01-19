#!/usr/bin/env python3

import argparse as ap
import os
import sys
from datetime import datetime
import numpy as np
import math as m
import vedo as v


# Define the main class that holds script data structures and settings.
class ScriptSettings():
    """The instance variables of this object are the user settings that
       control the program. The variable values are pulled from a list
       that is created within a resource control file and that are then
       reconciled with command line parameters."""

    def __init__(self):
        """Define default values for the graph parameters by pulling them
        from the resource control file in the default location:
        $OLCAO_RC/viewCellrc.py or from the current working directory if a local
        copy of viewCellrc.py is present."""

        # Read default variables from the resource control file.
        sys.path.insert(1, os.getenv('OLCAO_RC'))
        from viewCellrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults file.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc file.
        self.reconcile(args)

        # At this point, the command line parameters are set and accepted.
        #   When this initialization subroutine returns the script will
        #   start running. So, we use this as a good spot to record the
        #   command line parameters that were used.
        self.recordCLP()


    def assign_rc_defaults(self, default_rc):

        # Set default values.
        self.bz_to_show = default_rc["bz"]
        self.real_alpha = default_rc["real_alpha"]
        self.recip_alpha = default_rc["recip_alpha"]
        self.bz_alpha = default_rc["bz_alpha"]
        self.mesh_kp_alpha = default_rc["mesh_kp_alpha"]
        self.fold_kp_alpha = default_rc["fold_kp_alpha"]
        self.path_kp_alpha = default_rc["path_kp_alpha"]


    def parse_command_line(self):
    
        # Create the parser tool.
        prog_name = "viewCell"

        description_text = """
Version 0.1
Program to view the cell and kpoints used in OLCAO.
"""

        epilog_text = """
Please contact Paul Rulis (rulisp@umkc.edu) regarding questions.
Defaults are given in ./viewCellrc.py or $OLCAO_RC/viewCellrc.py.
"""

        parser = ap.ArgumentParser(prog = prog_name,
                formatter_class=ap.RawDescriptionHelpFormatter,
                description = description_text,
                epilog = epilog_text)
    
        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()


    def add_parser_arguments(self, parser):

        # Define the bz argument.
        parser.add_argument('-bz', '--brillouinzone',
                            dest='bz_to_show', type=int,
                            default=self.bz_to_show,
                            help=f'Argument bz_to_show. ' +
                            f'Default: {self.bz_to_show}')

        # Define the real_alpha argument.
        parser.add_argument('-rla', '--real_alpha',
                            dest='real_alpha', type=float,
                            default=self.real_alpha,
                            help=f'Argument real_alpha. ' +
                            f'Default: {self.real_alpha}')

        # Define the recip_alpha argument.
        parser.add_argument('-rca', '--recip_alpha',
                            dest='recip_alpha', type=float,
                            default=self.recip_alpha,
                            help=f'Argument recip_alpha. ' +
                            f'Default: {self.recip_alpha}')

        # Define the bz_alpha argument.
        parser.add_argument('-bza', '--bz_alpha',
                            dest='bz_alpha', type=float,
                            default=self.bz_alpha,
                            help=f'Argument bz_alpha. ' +
                            f'Default: {self.bz_alpha}')

        # Define the mesh_kp_alpha argument.
        parser.add_argument('-mkpa', '--mesh_kp_alpha',
                            dest='mesh_kp_alpha', type=float,
                            default=self.mesh_kp_alpha,
                            help=f'Argument mesh_kp_alpha. ' +
                            f'Default: {self.mesh_kp_alpha}')

        # Define the fold_kp_alpha argument.
        parser.add_argument('-fkpa', '--fold_kp_alpha',
                            dest='fold_kp_alpha', type=float,
                            default=self.fold_kp_alpha,
                            help=f'Argument fold_kp_alpha. ' +
                            f'Default: {self.fold_kp_alpha}')

        # Define the path_kp_alpha argument.
        parser.add_argument('-pkpa', '--path_kp_alpha',
                            dest='path_kp_alpha', type=float,
                            default=self.path_kp_alpha,
                            help=f'Argument path_kp_alpha. ' +
                            f'Default: {self.path_kp_alpha}')


    def reconcile(self, args):
        bz_to_show = args.bz_to_show
        real_alpha = args.real_alpha
        recip_alpha = args.recip_alpha
        bz_alpha = args.bz_alpha
        mesh_kp_alpha = args.mesh_kp_alpha
        fold_kp_alpha = args.fold_kp_alpha
        path_kp_alpha = args.path_kp_alpha


    def recordCLP(self):
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write(f"Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


class UserInterface():

    def __init__(self, settings, data, plt):

        slider_tab_length = 0.0025

        # Create a slider to control transparency of the real space cell.
        self.rl_slider = plt.add_slider(self.rl_slider_func,
                pos=[(0.025, 0.05), (0.175, 0.05)], title="Real Cell",
                xmin=0.0, xmax=1.0, value=data.real_cell.alpha(),
                slider_length=slider_tab_length)

        # Create a slider to control transparency of the recip space cell.
        self.rc_slider = plt.add_slider(self.rc_slider_func,
                pos=[(0.225, 0.05), (0.375, 0.05)], title="Recip Cell",
                xmin=0.0, xmax=1.0, value=data.recip_cell.alpha(),
                slider_length=slider_tab_length)

        # Create a slider to control transparency of the BZ cell.
        self.bz_slider = plt.add_slider(self.bz_slider_func,
                pos=[(0.425, 0.05), (0.575, 0.05)], title="Brillouin",
                xmin=0.0, xmax=1.0, value=data.bz_cell.alpha(),
                slider_length=slider_tab_length)

        # Create a slider to control transparency of the mesh kpoints.
        self.mesh_kp_slider = plt.add_slider(self.mesh_kp_slider_func,
                pos=[(0.025, 0.15), (0.175, 0.15)], title="Mesh KP",
                xmin=0.0, xmax=1.0, value=data.mesh_kp.alpha(),
                slider_length=slider_tab_length)

        # Create a slider to control transparency of the folded kpoints.
        self.fold_kp_slider = plt.add_slider(self.fold_kp_slider_func,
                pos=[(0.225, 0.15), (0.375, 0.15)], title="Folded KP",
                xmin=0.0, xmax=1.0, value=data.fold_kp.alpha(),
                slider_length=slider_tab_length)

        # Create a slider to control transparency of the path kpoints.
        self.path_kp_slider = plt.add_slider(self.path_kp_slider_func,
                pos=[(0.425, 0.15), (0.575, 0.15)], title="Path KP",
                xmin=0.0, xmax=1.0, value=data.path_kp.alpha(),
                slider_length=slider_tab_length)

    def rl_slider_func(self, widget, event):
        data.real_cell.alpha(widget.value)
        if (data.real_cell.alpha() == 0.0):
            data.real_cell.toggle()

    def rc_slider_func(self, widget, event):
        data.recip_cell.alpha(widget.value)
        if (data.recip_cell.alpha() == 0.0):
            data.recip_cell.toggle()

    def bz_slider_func(self, widget, event):
        data.bz_cell.alpha(widget.value)
        if (data.bz_cell.alpha() == 0.0):
            data.bz_cell.toggle()

    def mesh_kp_slider_func(self, widget, event):
        data.mesh_kp.alpha(widget.value)
        if (data.mesh_kp.alpha() == 0.0):
            data.mesh_kp.toggle()

    def fold_kp_slider_func(self, widget, event):
        data.fold_kp.alpha(widget.value)
        if (data.fold_kp.alpha() == 0.0):
            data.fold_kp.toggle()

    def path_kp_slider_func(self, widget, event):
        data.path_kp.alpha(widget.value)
        if (data.path_kp.alpha() == 0.0):
            data.path_kp.toggle()


class Data():

    # Declare class variables.
    real_cell_mags = []
    real_cell_vectors = []
    real_cell_vertices = []
    real_cell_edges = []
    real_cell_faces = []
    recip_cell_mags = []
    recip_cell_vectors = []
    recip_cell_vertices = []
    recip_cell_edges = []
    recip_cell_faces = []
    bz_cell_faces = []
    bz_cell_edges = []
    bz_cell_vertices = []
    num_kpoints_abc = []
    mesh_kpoints = []
    num_mesh_kpoints = 0
    folded_kpoints = []
    num_folded_kpoints = 0
    kpoint_weights = []
    path_kpoints = []
    path_kp_mag = []
    num_total_BZ_path_KP = 0
    num_high_symmetry_BZ_paths = 0
    num_high_symmetry_BZ_points_per_path = ()
    high_symmetry_BZ_kpoints = []

    real_cell = v.Mesh()
    recip_cell = v.Mesh()
    bz_cell = v.Mesh()
    mesh_kp = v.Mesh()
    fold_kp = v.Mesh()
    path_kp = v.Mesh()

    def __init__(self, settings):
        # Read the class variable values from the settings file.
        exec(open("BZ." + f"{settings.bz_to_show}").read())
        Data.num_mesh_kpoints = m.prod(Data.num_kpoints_abc)
        Data.num_folded_kpoints = len(Data.kpoint_weights)


    def form_real_cell(self, settings):
        Data.real_cell = v.Mesh([Data.real_cell_vertices,
                                Data.real_cell_faces])
        Data.real_cell.alpha(settings.real_alpha)
        Data.real_cell.color('green')
        Data.real_cell.linecolor('black').linewidth(5)
        return Data.real_cell


    def form_recip_cell(self, settings):
        Data.recip_cell = v.Mesh([Data.recip_cell_vertices,
                            Data.recip_cell_faces])
        Data.recip_cell.alpha(settings.recip_alpha)
        Data.recip_cell.color('red')
        Data.recip_cell.linecolor('black').linewidth(5)
        return Data.recip_cell


    def form_bz_cell(self, settings):
        Data.bz_cell = v.Mesh([Data.bz_cell_vertices,
                          Data.bz_cell_faces])
        Data.bz_cell.alpha(settings.bz_alpha)
        Data.bz_cell.color('blue')
        Data.bz_cell.linecolor('black').linewidth(5)
        return Data.bz_cell


    def form_mesh_kp(self, settings):
        # Create a set of uniform spheres for the mesh kpoints.
        Data.mesh_kp = v.Spheres(v.Mesh([Data.mesh_kpoints]),
                r=min(np.array(Data.recip_cell_mags)
                      / (np.array(Data.num_kpoints_abc)+1)) / 4.0)
        Data.mesh_kp.color('red')
        return Data.mesh_kp


    def form_fold_kp(self, settings):
        # Create a set of glyphs instead of a set of uniform spheres. This
        #   way, each glyph can be individually scaled according to the
        #   weighting factor.
        s = v.Sphere()  # Use a sphere as the glyph.

        # Compute the scale factor such that the largest sphere equals the
        #   size of the spheres for the mesh kpoint and is representative of
        #   the most heavily weighted kpoint. The problem here is that for
        #   visual correctness, the mesh kpoints and the folded kpoint with
        #   the *smallest possible* weight should be the same size. Instead,
        #   we have the *largest* weight matching the size of the mesh
        #   kpoints. Unfortunately, if we make the radii of the spheres
        #   maintain correct proportions and use a linear scaling factor,
        #   then the most heavily weighted folded kpoints may be unplesantly
        #   and unhelpfully large. So, I picked this--more visually
        #   appealing--style instead.
        scale_factor = min(np.array(Data.recip_cell_mags) /
                (np.array(Data.num_kpoints_abc)+1)) / 4.0 / \
                max(Data.kpoint_weights)
        # The kpoints don't have an "orientation", but it seems that this
        #   variable cannot be an array of scalars for some reason. It must
        #   be a numpy array of three-vectors. Hence, I just take the
        #   scalar weighting factors and tack on two empty columns to make
        #   a set of three-vectors of magnitude "kpoint_weights".
        orientation = np.c_[np.array(Data.kpoint_weights) * scale_factor,
                            np.zeros(Data.num_folded_kpoints),
                            np.zeros(Data.num_folded_kpoints)]
        Data.fold_kp = v.Glyph(v.Mesh([Data.folded_kpoints]).points,
                s, orientation, scale_by_vector_size=True)
        Data.fold_kp.color('white')
        return Data.fold_kp


    def form_path_kp(self, settings):

        # Make space to hold the BZ path kpoints and their magnitudes.
        Data.path_kpoints = np.zeros([Data.num_total_BZ_path_KP, 3])
        Data.path_kp_mag = np.zeros([Data.num_total_BZ_path_KP])

        # Convert the high symmetry BZ kpoints from scaled fractional to
        #   direct reciprocal space coordinates.
        latt = np.array(Data.recip_cell_vectors)
        direct_high_symmetry_BZ_kpoints = []
        for kp in Data.high_symmetry_BZ_kpoints:
            kp = np.dot(np.array(kp), latt)
            direct_high_symmetry_BZ_kpoints.append(kp.tolist())

        # Store the high symmetry BZ kpoints in an array with indices based
        #   on the number of paths.
        Data.high_symmetry_BZ_kpoints = np.zeros([
                Data.num_high_symmetry_BZ_paths,
                max(Data.num_high_symmetry_BZ_points_per_path),3])
        kPointCounter = 0
        for i in range(Data.num_high_symmetry_BZ_paths):
            for j in range(Data.num_high_symmetry_BZ_points_per_path[i]):
                Data.high_symmetry_BZ_kpoints[i][j] = \
                        direct_high_symmetry_BZ_kpoints[kPointCounter]
                kPointCounter += 1

        # Compute the distances between consecutive kpoints on the paths.
        #   NOTE: when transitioning from one path to the next, the distance
        #   between those two kpoints should be zero.
        symKPDistVect = np.zeros((Data.num_high_symmetry_BZ_paths,
                max(Data.num_high_symmetry_BZ_points_per_path),3))
        symKPDistMag = np.zeros((Data.num_high_symmetry_BZ_paths,
                max(Data.num_high_symmetry_BZ_points_per_path)))
        for i in range(Data.num_high_symmetry_BZ_paths):
            for j in range(Data.num_high_symmetry_BZ_points_per_path[i]):
                if ((i == 0) and (j == 0)): # First point on first path.
                    symKPDistVect[i][0] = [0.0, 0.0, 0.0]
                    symKPDistMag[i][0] = 0.0
                elif (j == 0): # First point on any other path.
                    # Needs to be equal in position and magnitude to the last
                    # point on the previous path. This will signify later that
                    # we have started a new path.
                    symKPDistVect[i][0] = symKPDistVect[i-1][
                            Data.num_high_symmetry_BZ_points_per_path[i-1]-1]
                    symKPDistMag[i][0] = symKPDistMag[i-1][
                            Data.num_high_symmetry_BZ_points_per_path[i-1]-1]
                else:
                    symKPDistVect[i][j] = \
                            Data.high_symmetry_BZ_kpoints[i][j] \
                            - Data.high_symmetry_BZ_kpoints[i][j-1]
                    symKPDistMag[i][j] = symKPDistMag[i][j-1] \
                            + m.sqrt(sum(symKPDistVect[i][j]**2))


        # Consider the special case where the number of requested path kpoints
        #   is equal to the number of high symmetry kpoints provided.
        kPointCounter = 0
        if (Data.num_total_BZ_path_KP ==
            sum(Data.num_high_symmetry_BZ_points_per_path)):
            for i in range(Data.num_high_symmetry_BZ_paths):
                for j in range(Data.num_high_symmetry_BZ_points_per_path[i]):
                    Data.path_kpoints[kPointCounter] = \
                            Data.high_symmetry_BZ_kpoints[i][j]
                    Data.path_kp_mag[kPointCounter] = symKPDistMag[i][j]
                    kPointCounter += 1
        else:
            average_delta = symKPDistMag[
                    Data.num_high_symmetry_BZ_paths - 1,
                    Data.num_high_symmetry_BZ_points_per_path[
                            Data.num_high_symmetry_BZ_paths-1]-1] \
                    / (Data.num_total_BZ_path_KP - 1)

            for i in range(Data.num_high_symmetry_BZ_paths):
                for j in range(Data.num_high_symmetry_BZ_points_per_path[i]):
                    if (j != 0):

                        # Determine the number of kpoints to use for the
                        #   current segment between highly symmetric kpoints
                        #   (j) and (j-1) of this path.
                        numSegmentKPoints = int((symKPDistMag[i][j] -
                                symKPDistMag[i][j-1]) / average_delta)

                        # Check that this number of segment kpoints will not
                        #   take us over the limit. If it does, then reduce
                        #   the number of points in this segment.
                        if (kPointCounter + numSegmentKPoints >
                                Data.num_total_BZ_path_KP):
                            numSegmentKPoints = Data.num_total_BZ_path_KP \
                                    - kPointCounter

                        # If there are not going to be any segment kpoints for
                        #   this segment then cycle to the next symmetric
                        #   kpoint. This is most commonly the case when we
                        #   start a new path.
                        if (numSegmentKPoints == 0):
                            continue
                        else:
                            # Get the size of the x,y,z deltas for this
                            #   segment.
                            segmentDelta = symKPDistVect[i][j] \
                                    / numSegmentKPoints

                            for k in range(numSegmentKPoints):

                                # Store the position of the next kpoint on
                                #   this segment by adding the above
                                #   determined delta to the last known path
                                #   kpoint.
                                Data.path_kpoints[kPointCounter] = \
                                        Data.path_kpoints[kPointCounter - 1] \
                                        + segmentDelta

                                # Calculate the magnitude of the vector for
                                #   the above path kpoint.
                                Data.path_kp_mag[kPointCounter] = \
                                        Data.path_kp_mag[kPointCounter - 1] \
                                        + m.sqrt(sum(segmentDelta**2))

                                # Increment the kpoint counter.
                                kPointCounter += 1
                    else:

                        # Now, record the j-th high symmetry kpoint. For the
                        #   first kpoint (where i==1 and j==1) this will
                        #   happen first, before the code segment above.

                        # Store the position of the high symmetry kpoint in
                        #   the path array.
                        Data.path_kpoints[kPointCounter] = \
                                Data.high_symmetry_BZ_kpoints[i][j]

                        # Calculate the magnitude of the vector for this high
                        #   symmetry kpoint.
                        if ((i == 0) and (j == 0)):
                            Data.path_kp_mag[kPointCounter] = 0.0
                        elif (j == 0):
                            Data.path_kp_mag[kPointCounter] = \
                                    Data.path_kp_mag[kPointCounter - 1]
                        else:
                            Data.path_kp_mag[kPointCounter] = \
                                    Data.path_kp_mag[kPointCounter - 1] \
                                    + m.sqrt(sum(segmentDelta**2))

                        # Increment the kpoint counter.
                        kPointCounter += 1

        # Correct the number of path kpoints.
        Data.total_num_BZ_path_KP = kPointCounter

        Data.path_kp = v.Spheres(v.Mesh([Data.path_kpoints]),
                r=min(np.array(Data.recip_cell_mags)) \
                        * len(direct_high_symmetry_BZ_kpoints)
                        / (np.array(Data.num_total_BZ_path_KP)) / 4.0)
        Data.path_kp.color('green')
        return Data.path_kp


def get_data(settings):

    # Read in the data.
    global data
    data = Data(settings)

    # Compute implicit data.


def plot_data(settings):

    # Modify default vedo settings.
    v.settings.multi_samples = 1 # Very weak (no) antialiasing
    v.settings.use_parallel_projection = True # Orthographic projection

    # Create the plotter.
    global plt
    plt = v.Plotter(axes=2)  # Axis arrows at 0,0,0.

    # Set up the real space cell.
    plt.add(data.form_real_cell(settings))

    # Set up the reciprocal space cell.
    plt.add(data.form_recip_cell(settings))

    # Set up the Brillouin Zone cell.
    plt.add(data.form_bz_cell(settings))

    # Set up the mesh kpoints.
    plt.add(data.form_mesh_kp(settings))

    # Set up the fold kpoints.
    plt.add(data.form_fold_kp(settings))

    # Set up the path kpoints.
    plt.add(data.form_path_kp(settings))

def main():

    # Get script settings from a combination of the resource control file
    #   and parameters given by the user on the command line.
    settings = ScriptSettings()

    # Start executing the main activities of the program.

    # Get the data.
    get_data(settings)

    # Initialize the vedo plotter.
    plot_data(settings)

    # Create user interface elements for the render window.
    ui = UserInterface(settings, data, plt)

    # Define the camera.
    cam = dict(
        pos=(4.0000, 4.0000, 25.000),
        focal_point = (0, 0, 0),
        viewup=(0, 1.00000, 0),
        roll=0,
        parallel_scale=63.5785,
        clipping_range=(169.941, 341.540),
    )
    
    # Show the scene.
    plt.show(camera=cam).close()


    # Finalize the program activities and quit.


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. The purpose of this is to allow another
    #   python program to import *this* script and call its functions
    #   internally.
    main()
