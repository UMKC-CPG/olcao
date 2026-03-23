#!/usr/bin/env python3

"""dump2skl.py -- Convert a LAMMPS dump file to an OLCAO skeleton file.

PROGRAM:  dump2skl.py
PURPOSE:  This program takes a LAMMPS xyz coordinate dump file as
   input and writes an OLCAO skeleton file (olcao.skl) of a
   specified timestep.
UPDATED:  Mar. 20, 2026
USAGE:    This program takes a LAMMPS dump file and an associated
   data file as input (specified with the -d or --dump command
   line switch and the -a or --data command line switch) and
   writes a file called olcao.skl.

   The user should specify either a timestep number (using -t
   or --timestep) or a frame number (using -f or --frame) to
   be used.  Timesteps correspond to those values listed in the
   dump file under "ITEM: TIMESTEP" and frame numbers are
   integers from 1 through the number of frames (timesteps
   printed in the file).  You can use a negative frame number
   to specify a frame counted backwards from the final frame,
   e.g. -f -1 to mean the last frame.

   COMMAND LINE OPTIONS:
     -d | --dump         Path to the LAMMPS dump file.
     -a | --data         Path to the LAMMPS data file.
     -f | --frame        Frame number (1 to N, or negative
                         to count from the end).
     -t | --timestep     Timestep number matching a value in
                         the dump file.
     -n | --name         System name for the skeleton file
                         header.  Default: "name".
     --np | --nonPeriodic  Add padding around the system to
                         make it act like a non-periodic
                         system.

WORKFLOW:
   1. Read the LAMMPS data file to determine which chemical
      element corresponds to each LAMMPS atom type (by
      comparing masses).
   2. Scan the dump file to count frames and record timestep
      values.
   3. Resolve the user's frame/timestep selection.
   4. Read box geometry and atom coordinates for the selected
      frame.
   5. Process coordinates: wrap into bounds, optionally pad,
      shift to origin, convert to fractional.
   6. Write the olcao.skl file.
"""

import argparse as ap
import os
import sys
from datetime import datetime

from element_data import ElementData


# ---------------------------------------------------------- #
#                     Script Settings                          #
# ---------------------------------------------------------- #

class ScriptSettings():
    """The instance variables of this object are the user
    settings that control the program.  The variable values are
    pulled from a list that is created within a resource control
    file and that are then reconciled with command line
    parameters."""

    def __init__(self):
        """Define default values for the parameters by pulling
        them from the resource control file in the default
        location: $OLCAO_RC/dump2sklrc.py or from the current
        working directory if a local copy of dump2sklrc.py is
        present."""

        # Read default variables from the resource control
        # file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit(
                "Error: $OLCAO_RC is not set. "
                "See instructions.")
        sys.path.insert(1, rc_dir)
        from dump2sklrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc
        # file.
        self.reconcile(args)

        # At this point, the command line parameters are set
        # and accepted.  When this initialization subroutine
        # returns the script will start running.  So, we use
        # this as a good spot to record the command line
        # parameters that were used.
        self.recordCLP()

    def assign_rc_defaults(self, default_rc):
        """Pull every setting from the rc dictionary into
        instance attributes."""

        # Path to the LAMMPS dump file.
        self.dump_file = default_rc["dump_file"]

        # Path to the LAMMPS data file.
        self.data_file = default_rc["data_file"]

        # Frame number (1 to N, or negative to count from
        # end).
        self.frame = default_rc["frame"]

        # Timestep number from the dump file.
        self.timestep = default_rc["timestep"]

        # System name for the skeleton file header.
        self.name = default_rc["name"]

        # Flag to add padding for non-periodic behavior.
        self.non_periodic = default_rc["non_periodic"]

    def parse_command_line(self):
        """Build and run the argument parser."""

        # Create the parser tool.
        prog_name = "dump2skl.py"

        description_text = """\
Convert a LAMMPS dump file to an OLCAO skeleton (.skl) file.

This program reads a LAMMPS dump file and an associated LAMMPS
data file, extracts the atomic coordinates and box geometry for
a user-specified timestep (or frame), determines element
identities by comparing atomic masses from the data file
against the OLCAO element database, and writes an olcao.skl
file suitable for use with the OLCAO code.

The program handles both orthogonal and non-orthogonal
(triclinic) simulation boxes.  Coordinates are converted to
fractional form for the skeleton file.

USAGE NOTES:
  - Either a frame number (-f) or a timestep number (-t) must
    be given (but not both).
  - Frame numbers range from 1 to the total number of frames
    in the dump file.
  - Negative frame numbers count backwards from the last frame
    (e.g. -f -1 means the last frame, -f -2 means second to
    last).
  - Timestep numbers must match an "ITEM: TIMESTEP" value that
    actually appears in the dump file."""

        epilog_text = """\
Defaults are given in ./dump2sklrc.py or $OLCAO_RC/dump2sklrc.py.
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

        # Path to the LAMMPS dump file.
        parser.add_argument(
            '-d', '--dump',
            dest='dump_file',
            type=str,
            default=self.dump_file,
            help='Path to the LAMMPS dump file. '
                 f'Default: "{self.dump_file}"')

        # Path to the LAMMPS data file.
        parser.add_argument(
            '-a', '--data',
            dest='data_file',
            type=str,
            default=self.data_file,
            help='Path to the LAMMPS data file. '
                 f'Default: "{self.data_file}"')

        # Frame number.  Positive values count from the
        # beginning (1 = first frame); negative values count
        # from the end (-1 = last frame).
        parser.add_argument(
            '-f', '--frame',
            dest='frame',
            type=int,
            default=self.frame,
            help='Frame number (1 to N, or negative to '
                 'count from end). '
                 f'Default: {self.frame}')

        # Timestep number.  Must match an "ITEM: TIMESTEP"
        # value in the dump file.
        parser.add_argument(
            '-t', '--timestep',
            dest='timestep',
            type=int,
            default=self.timestep,
            help='Timestep number matching an '
                 'ITEM: TIMESTEP value in the dump file. '
                 f'Default: {self.timestep}')

        # System name written into the skeleton file header.
        parser.add_argument(
            '-n', '--name',
            dest='name',
            type=str,
            default=self.name,
            help='Name for the system (written to the skl '
                 'file header). '
                 f'Default: "{self.name}"')

        # Non-periodic flag.  When set, 100 Angstroms of
        # padding are added on every side of the cell to
        # prevent interactions due to periodicity.
        # NOTE: The Perl version used -np as a short flag.
        # In argparse, -np would conflict with -n, so we
        # use --np (long form) instead.
        parser.add_argument(
            '--np', '--nonPeriodic',
            dest='non_periodic',
            action='store_true',
            default=self.non_periodic,
            help='Add 100 Ang padding on each side to '
                 'suppress periodic interactions. '
                 f'Default: {self.non_periodic}')

    def reconcile(self, args):
        """Reconcile command line arguments with rc file
        defaults."""

        self.dump_file = args.dump_file
        self.data_file = args.data_file
        self.frame = args.frame
        self.timestep = args.timestep
        self.name = args.name
        self.non_periodic = args.non_periodic

    def recordCLP(self):
        """Record the command line invocation to the 'command'
        file."""

        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime(
                "%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write("Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


# ---------------------------------------------------------- #
#                        Subroutines                           #
# ---------------------------------------------------------- #

def get_lammps_types(data_file):
    """Read the LAMMPS data file to determine the element for
    each LAMMPS atom type.

    The LAMMPS data file contains a "Masses" section that lists
    the atomic mass for each numbered atom type.  This function
    reads those masses and compares each one against the known
    atomic masses of all elements in the OLCAO element database.
    The element whose database mass is closest to the LAMMPS
    mass is assigned to that atom type.

    Parameters
    ----------
    data_file : str
        Path to the LAMMPS data file.

    Returns
    -------
    elements : dict
        Dictionary mapping LAMMPS atom type numbers
        (1-indexed int) to element name strings
        (e.g. {1: 'si', 2: 'o'}).
    num_atom_types : int
        Total number of atom types found in the data file.
    """

    # Initialize the element database so we can look up atomic
    # masses for comparison with the LAMMPS masses.
    ed = ElementData()
    ed.init_element_data()
    num_elements = ed.num_elements
    atomic_masses = ed.atomic_masses      # 1-indexed
    element_names = ed.element_names      # 1-indexed

    # Open the data file for reading.
    try:
        f = open(data_file, 'r')
    except FileNotFoundError:
        sys.exit(
            f"Could not open data file {data_file} "
            "for reading.")

    num_atom_types = 0
    elements = {}

    # Read through the file looking for the number of atom
    # types and the Masses section.
    for line in f:
        values = line.split()

        # Skip blank lines.
        if not values:
            continue

        # If we find the number of atomic types, record it.
        # The line format is:  "N atom types"
        if (len(values) > 2
                and values[1] == "atom"
                and values[2] == "types"):
            num_atom_types = int(values[0])

        # If we find the "Masses" tag, read the mass data.
        # Assume that there is a blank line between the
        # "Masses" tag and the first line of data.
        if values[0] == "Masses":
            f.readline()  # Skip the blank line.

            for atom_type in range(1, num_atom_types + 1):
                # Get the mass of this atom type.
                mass_values = f.readline().split()
                current_mass = float(mass_values[1])

                # Compare this mass to all masses in the
                # database and find the one with the smallest
                # difference.  This defines the element.
                min_diff = 1000.0
                min_diff_index = 0
                for z in range(1, num_elements + 1):
                    diff = abs(
                        current_mass - atomic_masses[z])
                    if diff < min_diff:
                        min_diff = diff
                        min_diff_index = z

                # Assign the element to this atom type.
                elements[atom_type] = (
                    element_names[min_diff_index])

            break  # Done reading masses.

    f.close()

    # Validate that we found what we needed.
    if num_atom_types == 0:
        sys.exit(
            "Error: Could not find 'atom types' in "
            f"the data file {data_file}.")
    if not elements:
        sys.exit(
            "Error: Could not find 'Masses' section "
            f"in the data file {data_file}.")

    return elements, num_atom_types


def in_bounds(atom_coords, box_bounds, tot_atoms):
    """Wrap atoms that have drifted outside the simulation box
    back inside.

    For each axis (x, y, z), if an atom's coordinate exceeds
    the upper box boundary, it is shifted back by the box
    length along that axis.  Similarly, if it is below the
    lower boundary, it is shifted forward.  This ensures all
    atoms are within the box for the subsequent coordinate
    processing.

    Parameters
    ----------
    atom_coords : list of list of float
        atom_coords[i] = [x, y, z] for atom i (0-indexed).
        Modified in place.
    box_bounds : list of list of float
        box_bounds[axis] = [lo, hi] for axis 0=x, 1=y, 2=z.
    tot_atoms : int
        Total number of atoms.
    """

    for axis in range(3):
        box_length = (box_bounds[axis][1]
                      - box_bounds[axis][0])

        for i in range(tot_atoms):
            # If the atom is beyond the box for this axis
            # in the positive direction, then we need to
            # shift the atom back by the length of the box
            # along this axis.
            if atom_coords[i][axis] > box_bounds[axis][1]:
                atom_coords[i][axis] -= box_length

            # Similarly, if the atom is beyond the box for
            # this axis in the negative direction, then we
            # need to shift the atom forward by the length
            # of the box along this axis.
            if atom_coords[i][axis] < box_bounds[axis][0]:
                atom_coords[i][axis] += box_length


def max_coord(atom_coords, tot_atoms, axis):
    """Find the maximum coordinate value along the specified
    axis.

    Takes an axis number as input:
      0 = x-axis
      1 = y-axis
      2 = z-axis

    Compares all the coordinate values in the column for that
    axis and finds the maximum.

    Parameters
    ----------
    atom_coords : list of list of float
        atom_coords[i] = [x, y, z] for atom i (0-indexed).
    tot_atoms : int
        Total number of atoms.
    axis : int
        0 for x-axis, 1 for y-axis, 2 for z-axis.

    Returns
    -------
    float
        The maximum coordinate value for the specified axis.
    """

    max_so_far = atom_coords[0][axis]
    for i in range(1, tot_atoms):
        if atom_coords[i][axis] > max_so_far:
            max_so_far = atom_coords[i][axis]
    return max_so_far


def min_coord(atom_coords, tot_atoms, axis):
    """Find the minimum coordinate value along the specified
    axis.

    Takes an axis number as input:
      0 = x-axis
      1 = y-axis
      2 = z-axis

    Compares all the coordinate values in the column for that
    axis and finds the minimum.

    Parameters
    ----------
    atom_coords : list of list of float
        atom_coords[i] = [x, y, z] for atom i (0-indexed).
    tot_atoms : int
        Total number of atoms.
    axis : int
        0 for x-axis, 1 for y-axis, 2 for z-axis.

    Returns
    -------
    float
        The minimum coordinate value for the specified axis.
    """

    min_so_far = atom_coords[0][axis]
    for i in range(1, tot_atoms):
        if atom_coords[i][axis] < min_so_far:
            min_so_far = atom_coords[i][axis]
    return min_so_far


def count_frames(dump_file):
    """Scan the dump file and catalog all frames and their
    timestep values.

    Reads through the entire LAMMPS dump file counting
    occurrences of "ITEM: TIMESTEP" to determine the total
    number of frames.  The timestep value that follows each
    such header is recorded in a list indexed by frame number
    (1-based).

    Parameters
    ----------
    dump_file : str
        Path to the LAMMPS dump file.

    Returns
    -------
    tot_frames : int
        Total number of frames found in the dump file.
    timesteps : list of int
        1-indexed list of timestep values.
        timesteps[i] is the timestep value for frame i
        (i = 1..tot_frames).  timesteps[0] is None (unused
        sentinel to preserve 1-based indexing).
    """

    tot_frames = 0
    # The timesteps list is 1-indexed: the frame number is
    # the index and the timestep number is the value.
    timesteps = [None]  # Index 0 unused.
    saw_timestep = False

    try:
        f = open(dump_file, 'r')
    except FileNotFoundError:
        sys.exit(
            f"Could not open dump file {dump_file} "
            "for reading.")

    for line in f:
        if "TIMESTEP" in line:
            tot_frames += 1
            saw_timestep = True
        elif saw_timestep:
            timesteps.append(int(line.strip()))
            saw_timestep = False

    f.close()

    if tot_frames == 0:
        sys.exit(
            "Error: No frames (TIMESTEP entries) found "
            f"in {dump_file}.")

    return tot_frames, timesteps


def resolve_frame_and_timestep(frame, timestep, timesteps,
                               tot_frames):
    """Validate and resolve the user's frame or timestep
    selection.

    The user specifies either a frame number or a timestep
    number.  This function ensures the choice is valid and
    returns both the resolved frame number and the
    corresponding timestep number.

    Frame numbers are 1-indexed (1 through tot_frames).
    Negative frame numbers count backwards from the end:
      -1 = last frame
      -2 = second to last
      etc.

    Parameters
    ----------
    frame : int or None
        User-requested frame number.
    timestep : int or None
        User-requested timestep number.
    timesteps : list of int
        1-indexed list of timestep values (from
        count_frames).
    tot_frames : int
        Total number of frames.

    Returns
    -------
    frame : int
        Resolved frame number (1-indexed, positive).
    timestep : int
        Resolved timestep number.
    """

    if frame is not None:
        # Check that the frame is within valid range.
        if (frame > tot_frames
                or frame < -tot_frames
                or frame == 0):
            sys.exit(
                "Frame out of range. Valid range: "
                f"1 to {tot_frames} "
                f"(or -1 to -{tot_frames}).")

        # If a negative frame number is given, convert it
        # to a positive frame number.  -1 means the last
        # frame, -2 means second to last, etc.
        if frame < 0:
            frame = tot_frames + (frame + 1)

        # Look up the corresponding timestep.
        timestep = timesteps[frame]

    elif timestep is not None:
        # Search for the frame that matches this timestep.
        found = False
        for step in range(1, tot_frames + 1):
            if timesteps[step] == timestep:
                frame = step
                found = True
                break
        if not found:
            sys.exit(
                f"Timestep {timestep} not found in "
                "the dump file.")

    else:
        sys.exit(
            "Must specify either a timestep (-t) or a "
            "frame number (-f).")

    return frame, timestep


def read_frame(dump_file, target_timestep):
    """Read box geometry and atom data for the target timestep.

    Seeks through the dump file to the frame matching the
    target timestep, then reads:

    1. The total number of atoms.
    2. The box bounds (orthogonal) or lattice vectors
       (triclinic) from the BOX BOUNDS section.
    3. The column layout from the ITEM: ATOMS header line.
    4. All atom coordinates and type numbers.

    LAMMPS dump files can have atoms listed in arbitrary order
    and with non-contiguous IDs.  This function reads all atom
    entries, handles gaps in the ID sequence by using a
    dictionary, and returns compacted (contiguous) lists of
    coordinates and types sorted by atom ID.

    Box type detection:
      - If the BOX BOUNDS header contains "abc origin", the
        box is triclinic (non-orthogonal) and the next three
        lines give lattice vectors directly.
      - Otherwise the box is orthogonal and the next three
        lines give lo/hi bounds for x, y, z.  Lattice
        vectors are derived as a diagonal matrix.

    Column detection:
      The ITEM: ATOMS header lists the column names, e.g.:
        ITEM: ATOMS id type x y z
      The first two tokens are "ITEM:" and "ATOMS", so the
      data column indices are offset by -2 from the header
      token positions.

    Parameters
    ----------
    dump_file : str
        Path to the LAMMPS dump file.
    target_timestep : int
        The timestep value to read.

    Returns
    -------
    tot_atoms : int
        Total number of atoms in this frame.
    box_bounds : list of list of float, or None
        box_bounds[axis] = [lo, hi] for orthogonal boxes.
        None for triclinic boxes.
    lattice : list of list of float
        3x3 lattice vectors.  For orthogonal boxes, this is
        a diagonal matrix derived from the box bounds.
    atom_coords : list of list of float
        Compacted [x, y, z] coordinates (0-indexed).
    atom_types : list of int
        Compacted LAMMPS atom type numbers (0-indexed).
    is_triclinic : bool
        True if the box is triclinic (non-orthogonal).
    """

    try:
        f = open(dump_file, 'r')
    except FileNotFoundError:
        sys.exit(
            f"Could not open dump file {dump_file} "
            "for reading.")

    # Move the read cursor to the correct timestep in the
    # dump file.  We search for "ITEM: TIMESTEP" and then
    # check whether the following line matches the target.
    found = False
    for line in f:
        if "TIMESTEP" in line:
            ts_line = f.readline().strip()
            if ts_line == str(target_timestep):
                found = True
                break

    if not found:
        sys.exit(
            f"Could not find timestep {target_timestep} "
            "in the dump file.")

    # --------------------------------------------------
    # Read the total number of atoms for this frame.
    # --------------------------------------------------
    # Next line is "ITEM: NUMBER OF ATOMS".
    f.readline()
    tot_atoms = int(f.readline().strip())

    # --------------------------------------------------
    # Read the box bounds / lattice vectors.
    # --------------------------------------------------
    # The next line is the BOX BOUNDS header which tells
    # us whether the box is orthogonal or triclinic.
    box_header = f.readline().strip()

    box_bounds = None
    lattice = [[0.0, 0.0, 0.0],
               [0.0, 0.0, 0.0],
               [0.0, 0.0, 0.0]]
    is_triclinic = False

    if "abc origin" in box_header:
        # ---- Non-orthogonal (triclinic) box ----
        # The next three lines give the lattice vectors
        # directly (a_x a_y a_z, b_x b_y b_z, c_x c_y c_z).
        is_triclinic = True
        for i in range(3):
            vals = f.readline().split()
            lattice[i][0] = float(vals[0])
            lattice[i][1] = float(vals[1])
            lattice[i][2] = float(vals[2])
    else:
        # ---- Orthogonal box ----
        # The next three lines give lo/hi bounds for the
        # x, y, and z axes respectively.
        box_bounds = [[0.0, 0.0],
                      [0.0, 0.0],
                      [0.0, 0.0]]

        for axis in range(3):
            vals = f.readline().split()
            box_bounds[axis][0] = float(vals[0])
            box_bounds[axis][1] = float(vals[1])

        # Convert orthogonal box bounds to lattice vectors.
        # For an orthogonal box, the lattice is a diagonal
        # matrix with the box dimensions on the diagonal.
        lattice[0][0] = (box_bounds[0][1]
                         - box_bounds[0][0])
        lattice[1][1] = (box_bounds[1][1]
                         - box_bounds[1][0])
        lattice[2][2] = (box_bounds[2][1]
                         - box_bounds[2][0])

    # --------------------------------------------------
    # Parse the ITEM: ATOMS header to find column layout.
    # --------------------------------------------------
    # The ITEM: ATOMS header line lists column names, e.g.:
    #   ITEM: ATOMS id type x y z
    # The data column indices are offset by -2 from the
    # header token positions because the header starts
    # with "ITEM:" and "ATOMS".
    while True:
        line = f.readline()
        if "ITEM: ATOMS" in line:
            header_tokens = line.split()
            break

    # Determine column indices for the data we need.
    xyz_idx = [None, None, None]
    type_idx = None

    for idx, token in enumerate(header_tokens):
        if token == "x":
            xyz_idx[0] = idx - 2
        elif token == "y":
            xyz_idx[1] = idx - 2
        elif token == "z":
            xyz_idx[2] = idx - 2
        elif token == "type":
            type_idx = idx - 2

    # --------------------------------------------------
    # Read atom coordinates and types for this frame.
    # --------------------------------------------------
    # LAMMPS dump files may list atoms in arbitrary order
    # and with non-contiguous IDs.  We store them by ID
    # in dictionaries and compact into contiguous lists
    # afterwards.
    raw_types = {}    # atom_id -> lammps_type
    raw_coords = {}   # atom_id -> [x, y, z]

    for line in f:
        line = line.strip()

        # End of this frame: either we hit the next
        # ITEM: TIMESTEP or the end of the file.
        if not line or "TIMESTEP" in line:
            break

        vals = line.split()
        atom_id = int(vals[0])
        raw_types[atom_id] = int(vals[type_idx])
        raw_coords[atom_id] = [
            float(vals[xyz_idx[0]]),
            float(vals[xyz_idx[1]]),
            float(vals[xyz_idx[2]])]

    f.close()

    # Compact the atom data into contiguous 0-indexed
    # lists, preserving the natural order of atom IDs.
    # This handles any gaps in the ID sequence.
    sorted_ids = sorted(raw_types.keys())
    atom_types = []
    atom_coords = []
    for aid in sorted_ids:
        atom_types.append(raw_types[aid])
        atom_coords.append(raw_coords[aid])

    return (tot_atoms, box_bounds, lattice,
            atom_coords, atom_types, is_triclinic)


def process_coordinates(atom_coords, box_bounds,
                        tot_atoms, non_periodic):
    """Process Cartesian coordinates into fractional form.

    This function performs the coordinate processing pipeline
    for orthogonal boxes:

    1. Wrap atoms into box bounds (in_bounds).
    2. Find the actual min/max extent of atom coordinates.
    3. Compare with LAMMPS-defined box bounds and take the
       larger extent (the box should encompass both).
    4. Optionally add 100 Angstrom padding on each side for
       non-periodic systems (to prevent interactions due to
       periodicity).
    5. Shift all coordinates so the box starts at the origin
       (all lo values become 0).
    6. Convert Cartesian coordinates to fractional by
       dividing by the box dimensions.

    Parameters
    ----------
    atom_coords : list of list of float
        atom_coords[i] = [x, y, z] in Cartesian coordinates.
        Modified in place to fractional coordinates.
    box_bounds : list of list of float
        box_bounds[axis] = [lo, hi] for axis 0=x, 1=y, 2=z.
    tot_atoms : int
        Total number of atoms.
    non_periodic : bool
        If True, add 100 Ang padding on each side.

    Returns
    -------
    cell_dims : list of float
        [a, b, c] cell dimensions in Angstroms (after any
        padding and shifting).
    """

    # Step 1: Make sure the box bounds completely encompass
    # all the atom coordinates.
    in_bounds(atom_coords, box_bounds, tot_atoms)

    # Step 2: Find the actual extent of the atom
    # coordinates along each axis.
    hi = [max_coord(atom_coords, tot_atoms, ax)
          for ax in range(3)]
    lo = [min_coord(atom_coords, tot_atoms, ax)
          for ax in range(3)]

    # Step 3: Compare with LAMMPS-defined box bounds and
    # redefine if needed.  We take the larger extent so the
    # cell fully encompasses the atoms and the declared box.
    for ax in range(3):
        if box_bounds[ax][1] > hi[ax]:
            hi[ax] = box_bounds[ax][1]
        if box_bounds[ax][0] < lo[ax]:
            lo[ax] = box_bounds[ax][0]

    # Step 4: If the non-periodic command line argument has
    # been selected, add padding space to the cell to
    # prevent interactions due to periodicity.
    if non_periodic:
        for ax in range(3):
            hi[ax] += 100.0
            lo[ax] -= 100.0

    # Step 5: Shift coordinates so that all box boundaries
    # go from 0 to some positive number.
    shift = [-lo[ax] for ax in range(3)]
    for i in range(tot_atoms):
        atom_coords[i][0] += shift[0]
        atom_coords[i][1] += shift[1]
        atom_coords[i][2] += shift[2]

    # Redefine max coordinates after shifting.  The min
    # coordinates are now all zero.
    for ax in range(3):
        hi[ax] += shift[ax]

    # Step 6: Rewrite atom coordinates in fractional form
    # by dividing each coordinate by the corresponding
    # box dimension.
    for i in range(tot_atoms):
        atom_coords[i][0] /= hi[0]
        atom_coords[i][1] /= hi[1]
        atom_coords[i][2] /= hi[2]

    # Return the cell dimensions (box lengths after shift).
    return hi


def write_skl(name, timestep, cell_dims, angles,
              tot_atoms, atom_coords, atom_types,
              elements):
    """Write the OLCAO skeleton (.skl) file.

    The skeleton file is the primary input geometry file for
    the OLCAO code.  It contains:
      - A title section with the system name, timestep, and
        date.
      - Cell parameters (lattice lengths and angles).
      - Fractional atomic coordinates with element labels.
      - Space group and supercell information.

    The file is written as 'olcao.skl' in the current
    directory.

    Parameters
    ----------
    name : str
        System name for the title section.
    timestep : int
        Timestep number for the title section.
    cell_dims : list of float
        [a, b, c] cell dimensions in Angstroms.
    angles : list of float
        [alpha, beta, gamma] cell angles in degrees.
    tot_atoms : int
        Total number of atoms.
    atom_coords : list of list of float
        Fractional coordinates [x, y, z] for each atom.
    atom_types : list of int
        LAMMPS atom type for each atom.
    elements : dict
        Mapping from LAMMPS atom type number to element
        name string.
    """

    with open('olcao.skl', 'w') as skl:
        # Get the current date.
        date = datetime.now().strftime(
            "%a %b %d %H:%M:%S %Y")

        # Write header information.
        skl.write("title\n")
        skl.write(f"{name}\n")
        skl.write(f"timestep {timestep}\n")
        skl.write(f"{date}\n")
        skl.write("end\n")

        # Write cell information.
        skl.write("cell\n")
        skl.write(
            f"{cell_dims[0]:.4f} "
            f"{cell_dims[1]:.4f} "
            f"{cell_dims[2]:.4f} "
            f"{angles[0]:.4f} "
            f"{angles[1]:.4f} "
            f"{angles[2]:.4f}\n")

        # Write atom information in fractional coordinates.
        skl.write(f"fractional {tot_atoms}\n")

        for i in range(tot_atoms):
            elem = elements[atom_types[i]]
            skl.write(
                f"{elem:<10s}"
                f"{atom_coords[i][0]:<15.10f}"
                f"{atom_coords[i][1]:<15.10f}"
                f"{atom_coords[i][2]:<15.10f}\n")

        # Print space group, supercell, and full information.
        skl.write(
            "space 1_a\nsupercell 1 1 1\nfull\n")


# ---------------------------------------------------------- #
# *********************************************************  #
# ---------------------------------------------------------- #
#                      Main Program                            #
# ---------------------------------------------------------- #
# *********************************************************  #
# ---------------------------------------------------------- #

def main():
    """Main program execution.

    Flow:
      1. Parse command line arguments and load settings.
      2. Read the LAMMPS data file to map atom types to
         chemical elements.
      3. Scan the dump file to count frames and record
         timestep values.
      4. Resolve the user's frame/timestep selection.
      5. Read the box geometry and atom coordinates for the
         selected frame.
      6. Process coordinates and write the olcao.skl file.
    """

    # Step 1: Get script settings from a combination of the
    # resource control file and parameters given by the user
    # on the command line.
    settings = ScriptSettings()

    # Step 2: Read the LAMMPS data file to determine which
    # chemical element corresponds to each LAMMPS atom type.
    # This is done by comparing the atomic masses listed in
    # the data file to the known masses in the OLCAO element
    # database.
    if not settings.data_file:
        sys.exit(
            "Error: No data file specified. "
            "Use -a or --data.")
    elements, num_atom_types = get_lammps_types(
        settings.data_file)

    # Step 3: Scan the dump file to find the total number of
    # frames and the timestep value associated with each
    # frame.
    if not settings.dump_file:
        sys.exit(
            "Error: No dump file specified. "
            "Use -d or --dump.")
    tot_frames, timesteps = count_frames(
        settings.dump_file)

    # Step 4: Check the requested timestep or frame number
    # and put either in terms of both the frame number and
    # the timestep (i.e. resolve both values).
    frame, timestep = resolve_frame_and_timestep(
        settings.frame, settings.timestep,
        timesteps, tot_frames)
    settings.frame = frame
    settings.timestep = timestep

    # Step 5: Read the box geometry and atom coordinates
    # for the selected frame from the dump file.
    (tot_atoms, box_bounds, lattice,
     atom_coords, atom_types, is_triclinic) = read_frame(
        settings.dump_file, timestep)

    # Step 6: Process coordinates and write the skeleton
    # file.  The processing differs for orthogonal vs.
    # triclinic boxes.
    if not is_triclinic:
        # Orthogonal box: full processing pipeline
        # (wrap, pad, shift, convert to fractional).
        cell_dims = process_coordinates(
            atom_coords, box_bounds,
            tot_atoms, settings.non_periodic)
    else:
        # Triclinic box: the lattice vectors define the
        # cell directly.  A full fractional coordinate
        # conversion for non-orthogonal cells would require
        # a matrix inversion of the lattice vectors.
        # NOTE: The original Perl version did not fully
        # support triclinic boxes in the coordinate
        # processing pipeline either.  For now, compute
        # the cell dimensions from the lattice vector
        # magnitudes and pass coordinates through as-is.
        cell_dims = [
            (lattice[0][0]**2
             + lattice[0][1]**2
             + lattice[0][2]**2) ** 0.5,
            (lattice[1][0]**2
             + lattice[1][1]**2
             + lattice[1][2]**2) ** 0.5,
            (lattice[2][0]**2
             + lattice[2][1]**2
             + lattice[2][2]**2) ** 0.5]
        print(
            "WARNING: Triclinic box detected. "
            "Coordinate processing is limited. "
            "Verify the output carefully.",
            file=sys.stderr)

    # Cell angles default to 90 degrees (orthogonal).
    angles = [90.0, 90.0, 90.0]

    # Write the skeleton file.
    write_skl(
        settings.name, timestep, cell_dims, angles,
        tot_atoms, atom_coords, atom_types, elements)


if __name__ == '__main__':
    # Everything before this point was a subroutine
    # definition or a request to import information from
    # external modules.  Only now do we actually start
    # running the program.  The purpose of this is to allow
    # another Python program to import this script and call
    # its functions internally.
    main()
