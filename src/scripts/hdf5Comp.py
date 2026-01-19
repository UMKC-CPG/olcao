#!/usr/bin/env python3

import argparse as ap
import os
import sys
import math as m
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime


# Define the main class that holds script data structures and settings.
class ScriptSettings():
    """The instance variables of this object are the user settings that
       control the program. The variable values are pulled from a list
       that is created within a resource control file and that are then
       reconciled with command line parameters."""


    def __init__(self):
        """Define default values for the graph parameters by pulling them
        from the resource control file in the default location:
        $OLCAO_RC/hdf5Comprc.py or from the current working directory if a
        local copy of hdf5Comprc.py is present."""

        # Read default variables from the resource control file.
        sys.path.insert(1, os.getenv('OLCAO_RC'))
        from hdf5Comprc import parameters_and_defaults
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

        # The first group of default variables defines the hdf5 file(s) to
        #   work on and where to find them.

        # The first file name.
        self.file1 = default_rc["file1"]

        # The second file name.
        self.file2 = default_rc["file2"]

        # The directory of the first file.
        self.dir1 = default_rc["dir1"]

        # The directory of the second file.
        self.dir2 = default_rc["dir2"]

        # The second group of default variables defines the group in the
        #   hdf5 file to access and the dataset to use.

        # The first group name.
        self.group1 = default_rc["group1"]

        # The second group name.
        self.group2 = default_rc["group2"]

        # The name of the first dataset.
        self.dataset1 = default_rc["dataset1"]

        # The name of the second dataset.
        self.dataset2 = default_rc["dataset2"]

        # The next group of default variables defines the specific way that
        #   the datasets should be treated.

        # Flag indicating that the dataset should be unpacked.
        self.unpack = default_rc["unpack"]

        # Flag indicating that the imaginary part of the packed matrix should
        #   be taken as the dataset.
        self.imaginary = default_rc["imaginary"]

        # Flag indicating that the absolute magnitude of the packed matrix
        #   should be taken as the dataset. (I.e., sqrt(z*z) for z=x+iy.)
        self.magnitude = default_rc["magnitude"]


    def parse_command_line(self):
    
        # Create the parser tool.
        prog_name = "HDF5 Matrix Comparison Tool"

        description_text = """
Version 0.1: 

The purpose of the program is to allow direct comparison of the matrices in
one or more OLCAO HDF5 files.
"""

        epilog_text = """
Please contact Paul Rulis (rulisp@umkc.edu) regarding questions.
Defaults are given in ./hdf5Comprc.py or $OLCAO_RC/hdf5Comprc.py.
"""

        parser = ap.ArgumentParser(prog = prog_name,
                formatter_class = ap.RawDescriptionHelpFormatter,
                description = description_text, epilog = epilog_text)
    
        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()


    def add_parser_arguments(self, parser):

        # Initialize local variables.
        store_action = False
    
        # Define the name of file 1.
        parser.add_argument('-f1', '--file1', dest='file1',
                            type=str, default=self.file1,
                            help=f'Name of file one. Default: {self.file1}')

        # Define the name of file 2.
        parser.add_argument('-f2', '--file2', dest='file2',
                            type=str, default=self.file2,
                            help=f'Name of file two. Default: {self.file2}')

        # Define the name of the directory to find file 1 in.
        parser.add_argument('-d1', '--dir1', dest='dir1',
                            type=str, default=self.dir1,
                            help=f'Name of directory one. Default: '
                            + '{self.dir1}')

        # Define the name of the directory to find file 2 in.
        parser.add_argument('-d2', '--dir2', dest='dir2',
                            type=str, default=self.dir2,
                            help=f'Name of directory two. Default: '
                            + '{self.dir2}')
    
        # Define the name of group 1.
        parser.add_argument('-g1', '--group1', dest='group1',
                            type=str, default=self.group1,
                            help=f'Name of group one. Default: {self.group1}')

        # Define the name of group 2.
        parser.add_argument('-g2', '--group2', dest='group2',
                            type=str, default=self.group2,
                            help=f'Name of group two. Default: {self.group2}')

        # Define the name of dataset 1.
        parser.add_argument('-s1', '--dataset1', dest='dataset1',
                            type=str, default=self.dataset1,
                            help=f'Name of dataset one. Default: '
                            + f'{self.dataset1}')

        # Define the name of dataset 2.
        parser.add_argument('-s2', '--dataset2', dest='dataset2',
                            type=str, default=self.dataset2,
                            help=f'Name of dataset two. Default: '
                            + f'{self.dataset2}')

        # Flag stating that the dataset should be unpacked.
        if (self.unpack == False): # False default means flag stores true.
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-u', '--unpack', action=store_action,
                            dest='unpack', default=self.unpack,
                            help=f'Flag to request that the datasets be '
                            + f'unpacked. Default: {self.unpack}')

        # Flag stating that the imaginary part should be taken.
        if (self.imaginary == False): # False default means flag stores true.
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-i', '--imaginary', action=store_action,
                            dest='imaginary', default=self.imaginary,
                            help=f'Flag to request the imaginary part of the '
                            + f'packed dataset. Default: {self.imaginary}')

        # Flag stating that the absolute magnitude should be taken.
        if (self.magnitude == False): # False default means flag stores true.
            store_action = "store_true"
        else:
            store_action = "store_false"
        parser.add_argument('-m', '--magnitude', action=store_action,
                            dest='magnitude', default=self.magnitude,
                            help=f'Flag to request the absolute magnitude of '
                            + f'the packed dataset. Default: {self.magnitude}')


    def reconcile(self, args):
        self.file1 = args.file1
        self.file2 = args.file2
        self.dir1 = args.dir1
        self.dir2 = args.dir2
        self.group1 = args.group1
        self.group2 = args.group2
        self.dataset1 = args.dataset1
        self.dataset2 = args.dataset2
        self.unpack = args.unpack
        self.imaginary = args.imaginary
        self.magnitude = args.magnitude


    def recordCLP(self):
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write(f"Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


class h5Data():

    # Declare class variables.
    file1_fid = 0
    file2_fid = 0
    data1_did = 0
    data2_did = 0
    data1 = []
    data2 = []
    diff = 0

    def __init__(self, settings):
        self.file1_fid = h5.File(f'{settings.dir1}/{settings.file1}')
        self.file2_fid = h5.File(f'{settings.dir2}/{settings.file2}')

        # The matrices are stored in one of a few formats:
        # (1) Packed form of two columns. Each column is the packed upper
        #   triangle of a symmetric or Hermitian matrix. The first column
        #   contains the real values, and the second column contains the
        #   imaginary values.
        # (2) Not packed. A full matrix of either real or imaginary values.
        if (settings.unpack):
            temp_data1 = self.file1_fid[
                    f'{settings.group1}/{settings.dataset1}'][()]
            temp_data2 = self.file2_fid[
                    f'{settings.group2}/{settings.dataset2}'][()]
            self.dimension_a = round((-1 + m.sqrt(1+8*len(temp_data1))) / 2)
            self.dimension_b = self.dimension_a
            self.data1 = np.empty([self.dimension_a, self.dimension_b])
            self.data2 = np.empty([self.dimension_a, self.dimension_b])

            # Determine if we are unpacking the real or imaginary matrix,
            #   or if we are taking the absolute magnitude.
            if ((not settings.imaginary) and (not settings.magnitude)):
                print("Got here 1")
                # Initialize the array index counter.
                array_index = 0

                # Unpack the symmetric real matrix.
                for i in range(self.dimension_a):
                    for j in range(i+1):
                        self.data1[j, i] = temp_data1[array_index, 0]
                        self.data1[i, j] = temp_data1[array_index, 0]
                        self.data2[j, i] = temp_data2[array_index, 0]
                        self.data2[i, j] = temp_data2[array_index, 0]
                        array_index += 1
            elif ((settings.imaginary) and (not settings.magnitude)):
                print("Got here 2")

                # Initialize the array index counter.
                array_index = 0

                # Unpack the imaginary part of the hermitian matrix.
                for i in range(self.dimension_a):
                    for j in range(i+1):
                        self.data1[j, i] = temp_data1[array_index, 1]
                        self.data1[i, j] = -temp_data1[array_index, 1]
                        self.data2[j, i] = temp_data2[array_index, 1]
                        self.data2[i, j] = -temp_data2[array_index, 1]
                        array_index += 1

            else: # magnitude
                print("Got here 3")

                # Initialize the array index counter.
                array_index = 0

                # Unpack the hermitian matrix.
                for i in range(self.dimension_a):
                    for j in range(i+1):
                        self.data1[j, i] = m.sqrt(temp_data1[array_index,0]**2\
                                + temp_data1[array_index, 1]**2)
                        self.data1[i, j] = self.data1[j, i]
                        self.data2[j, i] = m.sqrt(temp_data2[array_index,0]**2\
                                + temp_data2[array_index, 1]**2)
                        self.data2[i, j] = self.data2[j, i]
                        array_index += 1

        elif (not settings.magnitude): 
            # A full matrix (real OR imaginary) is explicitly stored.
            self.data1_did = \
                    self.file1_fid[f'{settings.group1}/{settings.dataset1}']
            self.data2_did = \
                    self.file2_fid[f'{settings.group2}/{settings.dataset2}']
            self.data1 = self.data1_did[()]
            self.data2 = self.data2_did[()]
            self.dimension_a = round((-1 + m.sqrt(1+8*len(self.data1))) / 2)
            self.dimension_b = round((-1 + m.sqrt(1+8*len(self.data1[0]))) / 2)

        else:
            # A full matrix magnitude is explicitly stored.
            temp_real = self.file1_fid[
                    f'{settings.group1}/real{settings.dataset1}'][()]
            temp_imag = self.file1_fid[
                    f'{settings.group1}/imag{settings.dataset1}'][()]
            self.data1 = np.sqrt(temp_real*temp_real + temp_imag*temp_imag)
            temp_real = self.file2_fid[
                    f'{settings.group2}/real{settings.dataset2}'][()]
            temp_imag = self.file2_fid[
                    f'{settings.group2}/imag{settings.dataset2}'][()]
            self.data2 = np.sqrt(temp_real*temp_real + temp_imag*temp_imag)
            self.dimension_a = len(self.data1)
            self.dimension_b = len(self.data1[0])


    def computeDiff(self):
        self.diff = self.data2 - self.data1
        print(self.data1)
        print(self.data2)
        print(self.diff)


    def show(self, settings):
        fig, ax = plt.subplots()
        #cax = ax.imshow(self.diff, origin="lower", cmap="cividis",
        #        extent=[0, self.dimension_a, 0, self.dimension_b])
        cax = ax.imshow(self.diff, cmap="cividis")
        cbar = fig.colorbar(cax)
        plt.show()


def main():

    # Get script settings from a combination of the resource control file
    #   and parameters given by the user on the command line.
    settings = ScriptSettings()

    # Start executing the main activities of the program.

    # Open and populate the hdf5Data.
    data = h5Data(settings)

    # Compute the difference.
    data.computeDiff()

    # Create the visualization.
    data.show(settings)


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. The purpose of this is to allow another
    #   python program to import *this* script and call its functions
    #   internally.
    main()
