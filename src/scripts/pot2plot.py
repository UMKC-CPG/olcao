#!/usr/bin/env python3

import argparse as ap
import os
import sys
import math as m
import numpy as np

def print_help():
    output = """
#Script to convert a OLCAO scfV and the associated structure.dat file into a
#   set of files that each contain a plottable set of data columns. Each file
#   is associated with one specific element. Each row is associated with
#   one potential type of that element. Each column is associated with one
#   potential function coefficient.
#Last updated July, 24, 2023
######################################################################
# Note that default parameters are define by the file
#   "$OLCAO_RC/pot2plotrc".
#
#---------------------------------------------------------------------
#
#USAGE:  pot2plot [-i scfV_file] [-s struct_file] [-d]
#OR:     pot2plot -help
#
#If the "-i" option is given then the user can select which scf potential
#   file to use as input.
#If the "-s" option is given then the user can select which OLCAO structure
#   file to use as input.
#If the "-d" option is given then the second half of the potential file
#   is used which corresponds to the spin-down case.
#
#DEFAULTS: Check the pot2plotrc file in either the local directory or in
#   the $OLCAO_RC directory for total confirmation of default values.
#For -i: The default is probably gs_scfV-fb.dat.
#For -s: The default is probably structure.dat.
#For -d: The default is probably "False".
#
#---------------------------------------------------------------------
#
#REQUIREMENTS
#
#You need to have some OLCAO scf potential file and a structure.dat file.
#
"""
    print (output)


# Define the main class that holds script settings.
class ScriptSettings():
    """The instance variables of this object are the user settings that
       control the program. The variable values are pulled from a list
       that is created within a resource control file and that are then
       reconciled with command line parameters."""


    def __init__(self):
        """Define default values for the graph parameters by pulling them
        from the resource control file in the default location:
        $OLCAO_RC/pot2plotrc.py or from the current working directory if a
        local copy of pot2plotrc.py is present."""

        # Read default variables from the resource control file.
        sys.path.insert(1, os.getenv('OLCAO_RC'))
        from pot2plotrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults file.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the arguments with the rc file.
        self.reconcile(args)


    def assign_rc_defaults(self, default_rc):

        # Define the default potential function file.
        self.scfV_file = default_rc[0]

        # Define the default OLCAO structure file.
        self.struct_file = default_rc[1]

        # Establish default boolean value for the -d option.
        self.do_down = default_rc[2]


    def parse_command_line(self):
    
        # Create the parser tool.
        prog_name = "pot2plot"
        description_text = """
Script to convert a OLCAO scfV and the associated structure.dat file into a
   set of files that each contain a plottable set of data columns. Each file
   is associated with one specific element. Each row is associated with
   one potential type of that element. Each column is associated with one
   potential function coefficient.
Last updated July, 23, 2023
"""
        epilog_text = """
Note that default parameters are define by the file "$OLCAO_RC/pot2plotrc".
"""
        parser = ap.ArgumentParser(prog = prog_name,
                                   description = description_text,
                                   epilog = epilog_text)
    
        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()


    def add_parser_arguments(self, parser):
    
        # Define the scfV_file input file argument.
        parser.add_argument('-i', '--infile', dest='scfV_file',
                            type=str, default=self.scfV_file,
                            help='Input scfV file. Default: ' +
                            f'{self.scfV_file}')

        # Define the structure input file argument.
        parser.add_argument('-s', '--struct', dest='struct_file',
                            type=str, default=self.struct_file,
                            help='Input structure file. Default: ' +
                            f'{self.struct_file}')
    
        # Define the do_down spin argument.
        parser.add_argument('-d', '--down_spin', dest='do_down',
                            type=bool, default=self.do_down,
                            help='Use second part of scfV file. ' +
                            f'Default: {self.do_down}')
    
        # Execute the argument parsing.
        args = parser.parse_args()
    
        # Return the results.
        return(args)


    def reconcile(self, args):

        self.scfV_file = args.scfV_file
        self.struct_file = args.struct_file
        self.do_down = args.do_down


# Define the class that holds and operates on potential function data.
class PotData():

    def __init__(self):
        pass


    def read_data(self, settings):

        # Read the structure file and extract useful data.
        self.read_structure(settings)

        # Read the potential file and extract useful data.
        self.read_scfV(settings)


    def read_structure(self, settings):

        # Read the structure file to get the element of each type.
        with open (settings.struct_file) as f:
            lines = f.readlines()

        # Determine the number of atomic sites.
        num_atom_sites = int(lines[5])

        # Determine the number of potential sites.
        num_pot_sites = int(lines[8+num_atom_sites])

        # Create an array to hold the element name of each type.
        self.type_element = []

        # Create a set of the unique elements in the system.
        self.elements = set()

        # Make a counter for the number of types.
        num_types = 0

        # Visit each of the lines defining a potential site and assign the
        #   element name for the associated type.
        for line in(lines[10+num_atom_sites:10+num_atom_sites+num_pot_sites]):

            # Get each item in the line.
            values = line.split()

            # Determine if we have encountered a new potential type.
            if (int(values[1]) > num_types):
                num_types += 1

                # Store the element ID for this type.
                self.type_element.append(values[5])

                # If this element is not already present, then add the element
                #   ID to the set of elements.
                if (not values[5] in self.elements):
                    self.elements.add(values[5])


    def read_scfV(self, settings):

        # Read the scf potential file to get
        with open (settings.scfV_file) as f:
            lines = f.readlines()

        # Get the number of potential types.
        self.num_pot_types = int(lines[0].split()[1])

        # Initialize a counter to track which line number we are on.
        current_line = 1

        # Create an array to store the number of terms for each type.
        self.num_terms = []

        # Create an array to store the list of coefficients for all types.
        self.coeffs = []

        # Create an array to store the list of alphas for all types.
        self.alphas = []

        # For each potential type, read in the term coefficients.
        for pot_type in range(self.num_pot_types):

            # Move the line position counter.
            current_line += 1
            
            # Get the number of terms for this type.
            self.num_terms.append(int(lines[current_line].split()[0]))

            # Create the nested array to hold the coeffs for this type.
            self.coeffs.append([])

            # Create the nested array to hold the alphas for this type.
            self.alphas.append([])

            # Store each term.
            for pot_term in range(self.num_terms[pot_type]):

                # Move the line position counter.
                current_line += 1

                # Store this coefficient and alpha.
                self.coeffs[pot_type].append(
                        float(lines[current_line].split()[0]))
                self.alphas[pot_type].append(
                        float(lines[current_line].split()[1]))


    def rewrite_pot(self):

        # Consider each element.
        for element in (self.elements):

            # Open an element specific output file.
            with open(f"pot_coeffs_{element}.plot",'w') as f:

                # Initialize a counter of the number of types for this elem.
                num_elem_types = 0

                # Visit all types.
                for pot_type in range(self.num_pot_types):

                    # Only write terms from types of the same element.
                    if (self.type_element[pot_type] == element):

                        # Increment the number of types of this element.
                        num_elem_types += 1

                        # Write the potential coeffs of this type of this elem
                        f.write(f"{num_elem_types}")
                        for coeff in self.coeffs[pot_type]:
                            f.write(f" {coeff}")
                        f.write("\n")


    def evaluate_pot(self):

        # Consider each element.
        for element in (self.elements):

            # Open an element specific output file.
            with open(f"pot_eval_{element}.plot",'w') as f:

                # Evaluate the potential function of each set of types for
                #   this element at "x".
                for x in np.linspace(0, 5, 1000):
                    f.write(f"{x}")

                    # Visit all types.
                    for pot_type in range(self.num_pot_types):

                        # Only evaluate terms from types of the same element.
                        if (self.type_element[pot_type] == element):

                            # Consider all terms in this type and evaluate
                            #   them at x.
                            summation = 0
                            for term in range(len(self.coeffs[pot_type])):
                                summation += self.coeffs[pot_type][term] * \
                                        m.exp(-self.alphas[pot_type][term] * 
                                                x*x)

                            f.write(f" {summation}")

                    f.write(f"\n")


def main():

    # Set default values from the user resource control file.
    settings = ScriptSettings()

    # Create the object that will hold potential and structural data.
    pot_data = PotData()

    # Read in the potential function and structure data.
    pot_data.read_data(settings)

    # Print the potential function data in the new format.
    pot_data.rewrite_pot()

    # Evaluate all the potential functions.
    pot_data.evaluate_pot()


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. The purpose of this is to allow another
    #   python program to import *this* script and call its functions
    #   internally.
    main()
