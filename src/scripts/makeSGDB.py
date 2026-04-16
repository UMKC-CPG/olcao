#!/usr/bin/env python3
"""
PROGRAM:  makeSGDB.py
PURPOSE:  This program is used to create a database of space groups
          and their symmetry operations for use in reading the
          olcao.dat input files. The database is consumed by the
          applySpaceGroup Fortran program.
UPDATED:  Mar. 18, 2026

The program reads a specially formatted plain-text data file that
lists every crystallographic space group together with its symmetry
operations (rotations, reflections, inversions, and fractional
translations).  Each space group may have one or more subgroups
(alternate settings / origins).  For every subgroup the program:
  1. Parses the symmetry operations expressed in compact x,y,z form
     (where x,y,z actually refer to the possibly non-orthogonal
     a,b,c lattice axes — the convention used by "A Hypertext Book
     of Crystallographic Space Group Diagrams and Tables" from
     Birkbeck College, University of London).
  2. Decomposes each operation into a 3x3 rotation/reflection matrix
     and a 3-component fractional shift vector.
  3. Handles non-primitive lattices by applying shifted repetitions
     of the base symmetry operations (e.g., body-centred,
     face-centred, etc.).
  4. Writes the processed data to individual files named after the
     space group (one file per subgroup).
  5. Creates symbolic links so that each space group can be
     referenced by its numeric ID (with letter suffixes for
     multiple subgroups).

USAGE:  makeSGDB.py -i <spaceGroupInFile>

The -i option is used to specify which data file should be read and
converted into the form readable by the applySpaceGroup program.
This is a required argument and has no default value.

See additional documentation at the end of this script regarding
the origins of the data used here (uniform ASCII notation for
space groups, point groups, and crystal systems).
"""

import argparse as ap
import os
import re
import sys
from datetime import datetime


class ScriptSettings():
    """The instance variables of this object hold the user settings
       that control the program. The variable values are reconciled
       from command line parameters."""

    def __init__(self):
        """Parse the command line, validate inputs, and record the
        command that was used."""

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments.
        self.reconcile(args)

        # Record the command line parameters that were used.
        self.record_clp()

    def parse_command_line(self):
        """Set up the argument parser and parse sys.argv."""

        prog_name = "makeSGDB.py"

        description_text = """
Create a database of space groups and their symmetry operations.

Reads a specially formatted plain-text data file listing every
crystallographic space group together with its symmetry operations
and produces individual database files suitable for the
applySpaceGroup Fortran program.

See the end of this script for full documentation on the uniform
ASCII notation for space groups.
"""

        epilog_text = """
Please contact the OLCAO maintainers regarding questions.
"""

        parser = ap.ArgumentParser(
                prog=prog_name,
                formatter_class=ap.RawDescriptionHelpFormatter,
                description=description_text,
                epilog=epilog_text)

        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()

    def add_parser_arguments(self, parser):
        """Define all command line arguments."""

        # The -i option specifies the input data file containing
        #   space group definitions.  This is a required argument
        #   and has no default value.
        parser.add_argument(
                '-i', '--input', dest='sg_in_file',
                type=str, required=True,
                help='Input data file of space group '
                     'definitions (required).')

    def reconcile(self, args):
        """Transfer parsed arguments into instance variables."""
        self.sg_in_file = args.sg_in_file

    def record_clp(self):
        """Append the command line invocation to the 'command'
        file for reproducibility."""
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write("Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


class SpaceGroupDB():
    """Holds all data structures needed to build the space group
    database, and provides methods that correspond to the logical
    steps of the conversion process.

    Instance variables
    ------------------
    sg_in_file : str
        Path to the input data file.
    num_symmetry_ops : int
        Total number of symmetry operations for the current space
        group subgroup (including shifted repetitions).
    symmetry_ops : list
        3-D list indexed as [axis_prime][axis_contrib][op_number].
        Stores the rotation/reflection coefficients (+1, -1, or 0)
        that define how each axis of the new position depends on
        the original x, y, z (a, b, c) coordinates.
        Note: although the files contain x,y,z the meaning is
        understood to be possibly non-orthogonal a,b,c axes.
        This was the convention used by "A Hypertext Book of
        Crystallographic Space Group Diagrams and Tables" by
        Birkbeck College, University of London, which was used
        as the source for the definition of the space groups.
    symmetry_shifts : list
        2-D list indexed as [axis][op_number].
        Stores the fractional translation component for each
        symmetry operation along each axis.
    num_shifts : int
        Number of shifted repetitions of original symmetry ops
        (1 for primitive cells, 2+ for centred cells).
    output_file : str
        Name for each space group output file.
    sub_group_name : str
        Name for the current subgroup of the current space group.
    space_lattice : str
        Character ID for the type of lattice
        (P, A, B, C, I, F, R, or H).
    space_group_tag : str
        String identifying the space group and its subgroup,
        expressed using a uniform ASCII form.  See documentation
        at the end of this script.
    """

    def __init__(self, settings):
        """Initialize the space group database builder.

        Parameters
        ----------
        settings : ScriptSettings
            Parsed and validated program settings.
        """
        self.sg_in_file = settings.sg_in_file
        self.num_symmetry_ops = 0
        self.symmetry_ops = []
        self.symmetry_shifts = []
        self.num_shifts = 0
        self.output_file = ""
        self.sub_group_name = ""
        self.space_lattice = ""
        self.space_group_tag = ""

    def make_db(self):
        """Read the data file and create the full database.

        This is the top-level driver that reads the input file,
        iterates over every space group and its subgroups, and
        produces one output file per subgroup along with
        convenience symbolic links.
        """

        # Open the data file for reading.
        with open(self.sg_in_file, "r") as infile:
            for line in infile:
                # Prepare the line components.
                values = line.split()

                # Determine the space group number and how
                #   many variations there are for this space
                #   group.
                group_number = values[0]
                num_sub_groups = int(values[1])

                # Read the empty line that follows the space
                #   group number and number of subgroups line.
                next(infile)

                # Read each of the subgroups and convert them.
                for sub_group in range(
                        1, num_sub_groups + 1):

                    # Obtain all naming information.
                    self.get_sub_group_name(
                            sub_group, infile)

                    # Open the file for this space group
                    #   subgroup for writing.
                    with open(self.output_file, "w") \
                            as outfile:

                        # Write the title for this space
                        #   group.
                        outfile.write(
                                self.space_group_tag)

                        # Write the root space group number
                        #   and subgroup for this space group.
                        outfile.write(
                                f"{group_number}"
                                f" {sub_group}\n")

                        # Read all the symmetry operations.
                        self.read_symmetry_components(
                                infile)

                        # Record the symmetry operations.
                        self.write_symmetry_components(
                                outfile)

                    # Create a few links for easy reference
                    #   in the olcao.skl files.  Note that
                    #   this will modify sub_group_name to
                    #   include escaped chars.
                    self.make_soft_links(
                            group_number,
                            num_sub_groups,
                            sub_group)

                # Read the empty line that follows the space
                #   group symmetry operations.
                next(infile)

    def get_sub_group_name(self, sub_group, infile):
        """Obtain all naming information for a subgroup.

        Reads the subgroup title line from the input file, extracts
        the group name, optional extension name (for alternate
        settings), and the space lattice type character.  Prepends
        the lattice info to the space group tag so that the Fortran
        program can determine if the cell can be reduced to a
        primitive cell.

        Parameters
        ----------
        sub_group : int
            The current subgroup number (1-based).
        infile : file object
            Open input data file positioned at the subgroup title.
        """

        # Read the subgroup title.
        line = next(infile)
        self.space_group_tag = line
        values = line.split()

        # Extract the name.
        current_group_name = values[0]

        # Extract the subName and the reason for it if the
        #   subname exists.  Prepend a "_" for the filename
        #   for convenience.
        if len(values) > 1:
            current_ext_name = "_" + values[1]
            # values[2] would be the sub_name_comment
            #   (reason for the subname), if present.
        else:
            current_ext_name = ""

        # Extract the status of the cell lattice (primitive or
        #   non-primitive) and if it is non-primitive, then what
        #   type is it?
        # The first character of the group name indicates the
        #   lattice type: P, A, B, C, I, F, R, or H.
        self.space_lattice = current_group_name[0]

        # Prepend the lattice info to the space group tag so
        #   that the FORTRAN program that uses this later can
        #   know if it can reduce this cell to a primitive cell.
        self.space_group_tag = (
                self.space_lattice + " "
                + self.space_group_tag)

        # Record the name of this group for reference when
        #   creating a link to it.
        self.sub_group_name = (
                current_group_name + current_ext_name)

        # Create the output file name.
        self.output_file = (
                current_group_name + current_ext_name)

    def read_symmetry_components(self, infile):
        """Read all symmetry operations for the current subgroup.

        Reads symmetry components line by line until a blank line
        is found or end-of-group.  Lines containing '(' indicate
        shifted repetitions (non-primitive lattices) and are
        handled by apply_shifted_repetition().  All other lines
        are individual symmetry operations parsed by
        save_symmetry_operation().

        Initializes and populates symmetry_ops and
        symmetry_shifts.  Also computes num_shifts (the ratio
        of total operations to the base set).

        Parameters
        ----------
        infile : file object
            Open input data file positioned at the first
            symmetry operation line.
        """

        # Initialize the number of symmetry ops currently read
        #   in, and the number that is the total for this space
        #   group.  These two numbers are the same except when
        #   there are shifted repetitions.  Then the
        #   current_num_ops represents the initial group, and
        #   the num_symmetry_ops represents the total number
        #   (original + shifted as they are applied).
        current_num_ops = 0
        self.num_symmetry_ops = 0

        # Reset the symmetry data structures for each new
        #   subgroup.  The op-axis is 1-indexed with a
        #   [None] sentinel at slot 0 to match the Perl
        #   convention where $numSymmetryOps is pre-
        #   incremented before the first store.
        self.symmetry_ops = [
                [[None], [None], [None]]
                for _ in range(3)]
        self.symmetry_shifts = [
                [None] for _ in range(3)]

        # Read symmetry components and save each one until we
        #   find either the end (signified by a blank line) or
        #   we find and call for shifted repetition.  This
        #   takes the form of a line with something like
        #   +(1/2 0 1/2) on it.  This means to repeat all the
        #   symmetry operations with the shifts given on the
        #   line applied to each of them.
        for line in infile:
            if line.strip() == "":
                break
            elif "(" in line:
                self.apply_shifted_repetition(
                        line, current_num_ops)
            else:
                current_num_ops += 1
                self.num_symmetry_ops += 1
                self.save_symmetry_operation(
                        line, self.num_symmetry_ops)

        # Compute the number of shifted sets of symmetry
        #   operations that exist for this space group.
        self.num_shifts = (
                self.num_symmetry_ops // current_num_ops)

    def save_symmetry_operation(self, line, op_num):
        """Parse and store a single symmetry operation.

        Each symmetry operation is given as three space-separated
        expressions (one per axis).  Each expression is decomposed
        by get_axis_op_components() into rotation/reflection
        coefficients and a fractional shift.

        Parameters
        ----------
        line : str
            A line containing three axis operation expressions
            separated by whitespace (e.g. "x -y z+1/2").
        op_num : int
            The 1-based index at which to store this operation.
        """

        values = line.split()

        # Deal with each axis separately.
        for axis in range(3):
            # Get the list of components to the current axis
            #   operation.
            component_list = self.get_axis_op_components(
                    values[axis])

            self.symmetry_ops[axis][0].append(
                    component_list[0])
            self.symmetry_ops[axis][1].append(
                    component_list[1])
            self.symmetry_ops[axis][2].append(
                    component_list[2])
            self.symmetry_shifts[axis].append(
                    component_list[3])

    def get_axis_op_components(self, components):
        """Decompose one axis operation expression into its parts.

        Parses a string like "x", "-y+1/2", "x-z", etc. character
        by character to extract:
          - The coefficient of x (a): +1, -1, or 0
          - The coefficient of y (b): +1, -1, or 0
          - The coefficient of z (c): +1, -1, or 0
          - The fractional shift (e.g. 1/2 → 0.5)

        The parser tracks sign state across characters: a '+' or
        '-' sets the sign for the next component; 'x', 'y', 'z'
        contribute the current sign to the respective coefficient;
        a digit begins a fraction (numerator/denominator).

        Parameters
        ----------
        components : str
            A single axis operation expression (e.g. "-x+1/2").

        Returns
        -------
        list of float
            [coeff_x, coeff_y, coeff_z, shift] where each
            coefficient is +1.0, -1.0, or 0.0 and shift is
            a float (0.0 if no fractional translation).
        """

        # Get every character in the current axis operation.
        chars = list(components)

        # Assume a positive sign for the first term.
        sign_factor = 1.0

        # Initialize the component list to be returned.
        component_list = [0.0, 0.0, 0.0, 0.0]

        # Consider each character in this axis operation.
        char_idx = 0
        while char_idx < len(chars):
            c = chars[char_idx]
            if c == '+':
                # Positive sign for next component.
                sign_factor = 1.0
                char_idx += 1
            elif c == '-':
                # Negative sign for next component.
                sign_factor = -1.0
                char_idx += 1
            elif c == 'x':
                component_list[0] = sign_factor
                char_idx += 1
            elif c == 'y':
                component_list[1] = sign_factor
                char_idx += 1
            elif c == 'z':
                component_list[2] = sign_factor
                char_idx += 1
            elif c.isdigit():
                numerator = int(chars[char_idx])
                # Skip the '/' separator (char_idx+1) and
                #   read the denominator (char_idx+2).
                denominator = int(chars[char_idx + 2])
                component_list[3] = (
                        (sign_factor * numerator)
                        / denominator)
                char_idx += 3
            else:
                char_idx += 1

        return component_list

    def apply_shifted_repetition(self, line,
                                 current_num_ops):
        """Handle shifted repetitions for non-primitive lattices.

        When a line like "+(1/2 0 1/2)" is encountered, all
        symmetry operations recorded so far are duplicated with
        the specified shift vector added to their fractional
        translations.  This is how body-centred (I), face-centred
        (F), and other non-primitive lattice types generate their
        full set of symmetry operations from the primitive ones.

        Parameters
        ----------
        line : str
            The shift specification line, e.g. "+(1/2 0 1/2)".
        current_num_ops : int
            The number of symmetry operations recorded so far
            (before this shifted repetition).
        """

        # Prepare the shift vector.
        # Split on '+(' or ')' or whitespace to extract the
        #   three fraction components.
        values = re.split(r'\+\(|\)|\s', line.strip())
        # Filter out empty strings from the split.
        values = [v for v in values if v]

        shift_vector = [0.0, 0.0, 0.0]
        for axis in range(3):
            parts = values[axis].split('/')
            if parts[0].strip() != '0' and \
                    parts[0].strip().isdigit():
                shift_vector[axis] = (
                        float(parts[0]) / float(parts[1]))
            else:
                shift_vector[axis] = 0.0

        # Copy all the symmetry operations recorded so far
        #   and apply the shift_vector to the symmetry_shifts
        #   as they are copied.
        for sym_op in range(1, current_num_ops + 1):
            self.num_symmetry_ops += 1
            for axis in range(3):
                for axis2 in range(3):
                    # Copy the old operation to the new one.
                    self.symmetry_ops[axis][axis2].append(
                            self.symmetry_ops[axis]
                            [axis2][sym_op])

                # Copy the old shift to the new one and apply
                #   the shift_vector.
                self.symmetry_shifts[axis].append(
                        self.symmetry_shifts[axis][sym_op]
                        + shift_vector[axis])

    def write_symmetry_components(self, outfile):
        """Write the processed symmetry data to an output file.

        Each symmetry operation is written as a set of four
        triplets:
          - The first three lines define the 3x3 rotation/
            reflection matrix (contributions to a, b, c
            components respectively).
          - The fourth line defines the fractional translation
            shift for each axis (a, b, c).

        The file header contains the total number of symmetry
        operations and the number of shifted repetitions.

        Parameters
        ----------
        outfile : file object
            Open output file for writing.
        """

        # Write the total number of symmetry operations and
        #   the number of shifted repetitions that exist in
        #   this data.
        outfile.write(
                f"{self.num_symmetry_ops}"
                f" {self.num_shifts}\n")

        # Each symmetry operation is actually a set of four
        #   triplets.  The first defines the contributions to
        #   the a component, the second to the b, and the
        #   third to the c.  The last triplet defines the
        #   amount of fractional shift to apply to each
        #   axis (a,b,c).
        for op in range(1, self.num_symmetry_ops + 1):
            outfile.write("\n")
            for axis in range(3):
                a = self.symmetry_ops[axis][0][op]
                b = self.symmetry_ops[axis][1][op]
                c = self.symmetry_ops[axis][2][op]
                outfile.write(
                        f"{a:12.8f}"
                        f" {b:12.8f}"
                        f" {c:12.8f}\n")

            # Print the shift triplet.
            sa = self.symmetry_shifts[0][op]
            sb = self.symmetry_shifts[1][op]
            sc = self.symmetry_shifts[2][op]
            outfile.write(
                    f"{sa:12.8f}"
                    f" {sb:12.8f}"
                    f" {sc:12.8f}\n")

    def make_soft_links(self, group_number, num_sub_groups,
                        sub_group):
        """Create symbolic links for easy reference.

        When there is only one subgroup, a simple numeric link
        is created (e.g. "1" → "P1").  When there are multiple
        subgroups, each gets a letter-suffixed link (e.g.
        "48_a" → "Pnnn_1", "48_b" → "Pnnn_2").

        Special characters in the subgroup name (backslashes,
        apostrophes) are escaped before creating the link.

        Parameters
        ----------
        group_number : str
            The numeric space group identifier (1–230).
        num_sub_groups : int
            Total number of subgroups for this space group.
        sub_group : int
            The current subgroup number (1-based).
        """

        # Obtain the ordinal value of the character "a"
        #   (typically in the ASCII set).
        ord_a = ord("a")

        # Modify the default name to include characters that
        #   need to be escaped when making a soft link to the
        #   file.
        escaped_name = self.sub_group_name.replace(
                "\\", "\\\\")
        escaped_name = escaped_name.replace("'", "\\'")

        # In the case that there is only one subgroup, then
        #   life is easy.  We simply make a link to it with
        #   the associated group number.  If there are
        #   multiple subgroups, we make a link to this one
        #   where the link name is the associated group
        #   number including the "_a" or "_b" etc extension.
        if num_sub_groups == 1:
            # Create a soft link to the first subgroup
            #   using the space group number.
            link_name = group_number
            if os.path.exists(link_name):
                os.remove(link_name)
            os.symlink(escaped_name, link_name)
        else:
            # Create a link name that is the group number
            #   plus a char extension.
            link_name = (
                    group_number + "_"
                    + chr(sub_group - 1 + ord_a))

            # Remove the old link if it exists and then
            #   create the new link.
            if os.path.exists(link_name):
                os.remove(link_name)
            os.symlink(escaped_name, link_name)


def main():
    """Main entry point for the makeSGDB program."""

    # Get script settings from the command line.
    settings = ScriptSettings()

    # Create the space group database builder and run it.
    db = SpaceGroupDB(settings)
    db.make_db()


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or
    #   a request to import information from external modules.
    #   Only now do we actually start running the program.  The
    #   purpose of this is to allow another python program to
    #   import *this* script and call its functions internally.
    main()


# SPECIAL DOCUMENTATION:
# This information is also expressed in the structure_control.py
#   module.
#
# The space group can be given as either a number (which will use
#   the standard default origin and unique axis), or you can
#   specify the exact space group name using the uniform ASCII
#   form for specifying space groups which was derived from the
#   list given here:  "Uniform ASCII symbols for space groups,
#   point groups, and crystal systems.  By P. SUSSE, Poster
#   presented at the 17th General Meeting of the International
#   Mineralogical Association, Toronto, Canada, Aug. 9 - 14,
#   1998. s. Abstracts  page A62.  It is also available at the
#   following web site:
#   http://www.saint-hilaire.ca/main/space.htm. This site is no
#   longer directly accessible, but it was archived by the
#   Internet Archive - WayBackMachine.  It can be accessed at
#   this https URL:
#   web.archive.org/web/20080723040057/
#     www.saint-hilaire.ca/main/space.htm
#   However, just in case that goes away too, I have copied the
#   content below.
# The space groups and alternate origins actually available are
#   derived from the book: A Hypertext Book of Crystallographic
#   Space Group Diagrams and Tables Copyright 1997-1999.
#   Birkbeck College, University of London.
#   Available at the web site:
#   http://img.chem.ucl.ac.uk/sgp/mainmenu.htm.
#   All names are from the book. Notational scheme is from the
#   reference above.
#
# New notation symbols for space groups, point groups and
#   crystal systems
#
# The "Uniform ASCII symbols for space groups, point groups,
#   and crystal systems" was presented at the 17th General
#   Meeting of the International Mineralogical Association,
#   Toronto, Canada, Aug. 9 - 14, 1998 by Doctor Peter Susse.
#
# The use of subscript and superscript characters in space group
#   and point group symbols is quite disadvantageous for printing
#   as well as for electronic data processing. While preparing
#   the 6th edition of the mineral database MINABS, new symbols
#   were designed consisting of a string of 5 ASCII characters
#   for SCHOENFLIES style, and of up to 7 characters for the
#   symbols after HERMANN-MAUGUIN.  It proved to be very useful
#   to extend the new symbols to crystal systems as well.
#
# The symbol for the crystal system is made taking the
#   SCHOENFLIES symbol of the holoedric point group, bringing
#   the subscript characters in line, and filling the remainder
#   of the 5 spaces with asterisks (ASCII 42), while in the
#   point group symbol the remaining spaces are filled with
#   periods (ASCII 46).  The space group symbol is made adding
#   the superscript number of SCHOENFLIES.
#
# The ASCII format for the HERMANN-MAUGUIN style symbols is
#   achieved by denoting inversion axes with the tilde " ~ "
#   character (ASCII 126).  The neutral screw axes 21, 42, and
#   63 are indicated by an apostrophe " ' " (ASCII 39), while
#   the (right handed) screw axes 31, 41, and 61 are indicated
#   by a right bracket " ] " (ASCII 93), or by a right brace
#   " } " (ASCII 125) in case of 62.  The (left handed) screw
#   axes 32, 43, and 65 are denoted using the left bracket
#   " [ " (ASCII 91), or the left brace " { " (ASCII 123) in
#   case of 64.
#
# For retrieving purposes in databases, the SCHOENFLIES symbols
#   are useful because of their hierarchical nature and of being
#   independent of crystal settings.  The advantage of HERMANN-
#   MAUGUIN symbols is the orientational information they
#   contain.  The combination of both, using the new ASCII
#   symbols, is being successfully applied in the database
#   MINABS.
#
# If you require more information on this new notation, you may
#   contact Dr. Peter Susse or visit his website.
#
# 1. New symbols for the seven crystal systems.
#
# ci***   triclinic
# c2h**   monoclinic
# d2h**   orthorhombic
# d4h**   tetragonal
# d3d**   trigonal
# d6h**   hexagonal
# oh***   isometric
#
#
# 2. New symbols for the 32 point groups.
# Left column - new symbol Schoenflies type
# Right column - Hermann-Mauguin type
#
# c1...   1
# ci...   1~
# c2...   2
# cs...   m
# c2h..   2/m
# d2...   2 2 2
# c2v..   m m 2
# d2h..   m m m
# c4...   4
# s4...   4~
# c4h..   4/m
# d4...   4 2 2
# c4v..   4 m m
# d2d..   4~2 m
# d4h..   4/m m m
# c3...   3
#
# c3i..   3~
# c3v..   3 m
# d3...   3 2
# d3d..   3~m
# c6...   6
# c3h..   6~
# c6h..   6/m
# d6...   6 2 2
# c6v..   6 m m
# d3h..   6~2 m
# d6h..   6 m m
# t....   2 3
# th...   m 3
# o....   4 3 2
# td...   4~3 m
# oh...   m 3 m
#
#
# 3. New symbols for the 230 space groups.
# 1st column  - Space group
# 2nd column  - New symbol Schoenflies type
# 3rd column  - New symbol Hermann-Mauguin type
# 4th column  - Other Orientations
#
# Space  Sch. type   H-M type   Other
# 1 -    c1..1   P1
# 2 -    ci...1       P1~
# 3 -    c2...1       P2 B2
# 4 -    c2..2   P2'      B2'
# 5 -    c2..3   C2       A2 I2 F2
# 6 -    cs..1   Pm       Bm
# 7 -    cs..2   Pc       Pa Pn Bd
# 8 -    cs..3   Cm       Am Im Fm
# 9 -    cs..4   Cc       Aa Ia Fd
# 10 -   c2h.1   P2/m     B2/m
# 11 -   c2h.2   P2'/m    B2'/m
# 12 -   c2h.3   C2/m     A2/m I2/m F2/m
# 13 -   c2h.4   P2/c     P2/a P2/n B2/d
# 14 -   c2h.5   P2'/c    P2'/a P2'/n B2'/d
# 15 -   c2h.6   C2/c     A2/a I2/a F2/d
# 16 -   d2..1   P222
# 17 -   d2..2   P222'    P2'22 P22'2
# 18 -   d2..3   P2'2'2   P22'2' P2'22'
# 19 -   d2..4   P2'2'2'
# 20 -   d2..5   C222'    A2'22 B22'2
# 21 -   d2..6   C222     A222 B222
# 22 -   d2..7   F222
# 23 -   d2..8   I222
# 24 -   d2..9   I2'2'2'
# 25 -   c2v.1   Pmm2     P2mm Pm2m
# 26 -   c2v.2   Pmc2'    P2'ma Pb2'm Pm2'b Pcm2' P2'am
# 27 -   c2v.3   Pcc2     P2aa Pb2b
# 28 -   c2v.4   Pma2     P2mb Pc2m Pm2a Pbm2 P2cm
# 29 -   c2v.5   Pca2'    P2'ab Pc2'b Pb2'a Pbc2' P2'ca
# 30 -   c2v.6   Pnc2     P2na Pb2n Pn2b Pcn2 P2an
# 31 -   c2v.7   Pmn2'    P2'mn Pn2'm Pm2'n Pnm2' P2'nm
# 32 -   c2v.8   Pba2     P2cb Pc2a
# 33 -   c2v.9   Pna2     P2'nb Pc2'n Pn2'a Pbn2' P2'cn
# 34 -   c2v10   Pnn2     P2nn Pn2n
# 35 -   c2v11   Cmm2     A2mm Bm2m
# 36 -   c2v12   Cmc2'    A2'ma Bb2'm Bm2'b Ccm2' A2'am
# 37 -   c2v13   Ccc2     A2aa Bb2b
# 38 -   c2v14   Amm2     B2mm Cm2m Am2m Bmm2 C2mm
# 39 -   c2v15   Abm2     B2cm Cm2a Ac2m Bma2 C2mb
# 40 -   c2v16   Ama2     B2mb Cc2m Am2a Bbm2 C2cm
# 41 -   c2v17   Aba2     B2cb Cc2a Ac2a Bba2 C2cb
# 42 -   c2v18   Fmm2     F2mm Fm2m
# 43 -   c2v19   Fdd2     F2dd Fd2d
# 44 -   c2v20   Imm2     I2mm Im2m
# 45 -   c2v21   Iba2     I2cb Ic2a
# 46 -   c2v22   Ima2     I2mb Ic2m Im2a Ibm2 I2cm
# 47 -   d2h.1   Pmmm
# 48 -   d2h.2   Pnnn
# 49 -   d2h.3   Pccm     Pmaa Pbmb
# 50 -   d2h.4   Pban     Pncb Pcna
# 51 -   d2h.5   Pmma     Pbmm Pmcm Pmam Pmmb Pcmm
# 52 -   d2h.6   Pnna     Pbnn Pncn Pnan Pnnb Pcnn
# 53 -   d2h.7   Pmna     Pbmn Pncm Pman Pnmb Pcnm
# 54 -   d2h.8   Pcca     Pbaa Pbcb Pbab Pccb Pcaa
# 55 -   d2h.9   Pbam     Pmcb Pcma
# 56 -   d2h10   Pccn     Pnaa Pbnb
# 57 -   d2h11   Pbcm     Pmca Pbma Pcmb Pcam Pmab
# 58 -   d2h12   Pnnm     Pmnn Pnmn
# 59 -   d2h13   Pmmn     Pnmm Pmnm
# 60 -   d2h14   Pbcn     Pnca Pbna Pcnb Pcan Pnab
# 61 -   d2h15   Pbca     Pcab
# 62 -   d2h16   Pnma     Pbnm Pmcn Pnam Pmnb Pcmn
# 63 -   d2h17   Cmcm     Amma Bbmm Bmmb Ccmm Amam
# 64 -   d2h18   Cmca     Abma Bbcm Bmab Ccmb Acam
# 65 -   d2h19   Cmmm     Bmmm
# 66 -   d2h20   Cccm     Amaa Bbmb
# 67 -   d2h21   Cmma     Abmm Bmcm Bmam Cmmb Acmm
# 68 -   d2h22   Ccca     Abaa Bbcb Bbab Cccb Acaa
# 69 -   d2h23   Fmmm
# 70 -   d2h24   Fddd
# 71 -   d2h25   Immm
# 72 -   d2h26   Ibam     Imcb Icma
# 73 -   d2h27   Ibca     Icab
# 74 -   d2h28   Imma     Ibmm Imcm Imam Immb Icmm
# 75 -   c4..1   P4       C4
# 76 -   c4..2   P4]      C4]
# 77 -   c4..3   P4'      C4'
# 78 -   c4..4   P4[      C4[
# 79 -   c4..5   I4       F4
# 80 -   c4..6   I4]      F4]
# 81 -   s4..1   P4~      C4~
# 82 -   s4..2   I4~      F4~
# 83 -   c4h.1   P4/m     C4/m
# 84 -   c4h.2   P4'/m    C4'/m
# 85 -   c4h.3   P4/n     C4/a
# 86 -   c4h.4   P4'/n    C4'/a
# 87 -   c4h.5   I4/m     F4/m
# 88 -   c4h.6   I4]/a    F4]/d
# 89 -   d4..1   P422     C422
# 90 -   d4..2   P42'2    C422'
# 91 -   d4..3   P4]22    C4]22
# 92 -   d4..4   P4]2'2   C4]22'
# 93 -   d4..5   P4'22    C4'22
# 94 -   d4..6   P4'2'2   C4'22'
# 95 -   d4..7   P4[22    C4[22
# 96 -   d4..8   P4[2'2   C4[22'
# 97 -   d4..9   I422     F422
# 98 -   d4.10   I4]22    F4]22
# 99 -   c4v.1   P4mm     C4mm
# 100 -  c4v.2   P4bm     C4mb
# 101 -  c4v.3   P4'cm    C4'mc
# 102 -  c4v.4   P4'nm    C4'mn
# 103 -  c4v.5   P4cc     C4cc
# 104 -  c4v.6   P4nc     C4cn
# 105 -  c4v.7   P4'mc    C4'cm
# 106 -  c4v.8   P4'bc    C4'cb
# 107 -  c4v.9   I4mm     F4mm
# 108 -  c4v10   I4cm     F4mc
# 109 -  c4v11   I4]md    F4]dm
# 110 -  c4v12   I4]cd    F4]dc
# 111 -  d2d.1   P4~2m    C4~m2
# 112 -  d2d.2   P4~2c    C4~c2
# 113 -  d2d.3   P4~2'm   C4~m2'
# 114 -  d2d.4   P4~2'c   C4~c2'
# 115 -  d2d.4   P4~2'c   C4~c2'
# 116 -  d2d.5   P4~m2    C4~2m
# 117 -  d2d.6   P4~c2    C4~2c
# 118 -  d2d.7   P4~b2    C4~2b
# 119 -  d2d.8   P4~n2    C4~2n
# 120 -  d2d.9   I4~m2    F4~2m
# 121 -  d2d10   I4~c2    F4~2c
# 122 -  d2d11   I4~2m    F4~m2
# 123 -  d2d12   I4~2d    F4~d2
# 124 -  d4h.1   P4/mmm   C4/mmm
# 125 -  d4h.2   P4/mcc   C4/mcc
# 126 -  d4h.3   P4/nbm   C4/amb
# 127 -  d4h.4   P4/nnc   C4/acn
# 128 -  d4h.5   P4/mbm   C4/mmb
# 129 -  d4h.6   P4/mnc   C4/mcn
# 130 -  d4h.7   P4/nmm   C4/amm
# 131 -  d4h.8   P4/ncc   C4/acc
# 132 -  d4h.9   P4'/mmc  C4'/mcm
# 133 -  d4h10   P4'/mcm  C4'/mmc
# 134 -  d4h11   P4'/nbc  C4'/acb
# 135 -  d4h12   P4'/nnm  C4'/amn
# 136 -  d4h13   P4'/mbc  C4'/mcb
# 137 -  d4h14   P4'/mnm  C4'/mmn
# 138 -  d4h15   P4'/nmc  C4'/acm
# 139 -  d4h16   P4'/ncm  C4'/amc
# 140 -  d4h17   I4/mmm   F4/mmm
# 141 -  d4h18   I4/mcm   F4/mmc
# 142 -  d4h19   I4]/amd  F4]/ddm
# 143 -  d4h20   I4]/acd  F4]/ddc
# 144 -  c3..1   P3
# 145 -  c3..2   P3]
# 146 -  c3..3   P3[      R3
# 147 -  c3i.1   1 P3~
# 148 -  c3i.2   R3~
# 149 -  d3..1   P312
# 150 -  d3..2   P321
# 151 -  d3..3   P3]12
# 152 -  d3..4   P3]21
# 153 -  d3..5   P3[12
# 154 -  d3..6   P3[21
# 155 -  d3..7   R32
# 156 -  c3v.1   P3m1
# 157 -  c3v.2   P31m
# 158 -  c3v.3   P3c1
# 159 -  c3v.4   P31c
# 160 -  c3v.5   R3m
# 161 -  c3v.6   R3c
# 162 -  d3d.1   P3~1m
# 163 -  d3d.2   P3~1c
# 164 -  d3d.3   P3~m1
# 165 -  d3d.4   P3~c1
# 166 -  d3d.5   R3~m
# 167 -  d3d.6   R3~c
# 168 -  c6..1   P6
# 169 -  c6..2   P6]
# 170 -  c6..3   P6[
# 171 -  c6..4   P6}
# 172 -  c6..5   P6{
# 173 -  c6..6   P6'
# 174 -  c3h.1   P6~
# 175 -  c6h.1   P6/m
# 176 -  c6h.2   P6'/m
# 177 -  d6..1   P622
# 178 -  d6..2   P6]22
# 179 -  d6..3   P6[22
# 180 -  d6..4   P6}22
# 181 -  d6..5   P6{22
# 182 -  d6..6   P6'22
# 183 -  c6v.1   P6mm
# 184 -  c6v.2   P6cc
# 185 -  c6v.3   P6'cm
# 186 -  c6v.4   P6'mc
# 187 -  d3h.1   P6~m2
# 188 -  d3h.2   P6~c2
# 189 -  d3h.3   P6~2m
# 190 -  d3h.4   P6~2c
# 191 -  d6h.1   P6/mmm
# 192 -  d6h.2   P6/mcc
# 193 -  d6h.3   P6'/mcm
# 194 -  d6h.4   P6'/mmc
# 195 -  t...1   P23
# 196 -  t...2   F23
# 197 -  t...3   I23
# 198 -  t...4   P2'3
# 199 -  t...5   I2'3
# 200 -  th..1   Pm3
# 201 -  th..2   Pn3
# 202 -  th..3   Fm3
# 203 -  th..4   Fd3
# 204 -  th..5   Im3
# 205 -  th..6   Pa3
# 206 -  th..7   Ia3
# 207 -  o...1   P432
# 208 -  o...2   P4'32
# 209 -  o...3   F432
# 210 -  o...4   F4]32
# 211 -  o...5   I432
# 212 -  o...6   P4[32
# 213 -  o...7   P4]32
# 214 -  o...8   I4]32
# 215 -  td..1   P4~3m
# 216 -  td..2   F4~3m
# 217 -  td..3   I4~3m
# 218 -  td..4   P4~3n
# 219 -  td..5   F4~3c
# 220 -  td..6   I4~3d
# 221 -  oh..1   Pm3m
# 222 -  oh..2   Pn3n
# 223 -  oh..3   Pm3n
# 224 -  oh..4   Pn3m
# 225 -  oh..5   Fm3m
# 226 -  oh..6   Fm3c
# 227 -  oh..7   Fd3m
# 228 -  oh..8   Fd3c
# 229 -  oh..9   Im3m
# 230 -  oh.10   Ia3d
#
#
# If you have any comments or suggestions, please let me know.
