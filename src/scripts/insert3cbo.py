#!/usr/bin/env python3

"""insert3cbo -- Merge two- and three-center bond orders.

PROGRAM:  insert3cbo.py
PURPOSE:  To create new bond order raw and structure files
   such that they include new information describing
   three center bonds and such that information describing
   two center bonds is removed (except for real two center
   bonds which are identified by the fact that there are
   no three center bonds with them in them).

FIX: This program is basically broken until it can be
   remerged with the makeBOND program. The makeBOND
   program now uses olcao.skl to populate the atom and
   lattice data. This program worked by directly modifying
   the structure.dat file instead. So, perhaps this
   program can create another file that will be read in
   by makeBOND and used to modify the data obtained from
   the olcao.skl file in the future. Until then, this
   program is broken.

It is assumed that the three center bond order information
   is stored in the gs_bond-mb.3c.raw file. This default
   can be overridden with the -3cdata option.

It is assumed that the two center bond order information
   is stored in the gs_bond-mb.raw file. This default can
   be overridden with the -2cdata option.

It is assumed that the structure information will be found
   in structure.dat. This default can be overridden with
   the -sdata option.

It is assumed that the raw data output file will be
   gs_bond-mb.23cbo.raw. This is the original 2-center
   bond order file that has been manipulated to contain
   information for 2 and 3 center bond order. This
   default can be overridden with the -rawout option.

It is assumed that the structure output file will be
   structure.23cbo.dat. This is the original structure
   file modified to include the addition of fake atoms
   that represent the centroids of the three center
   bonds. This default can be overridden with the -sout
   option.

AUTHOR:  Paul Rulis (Perl original); Python port by
   Claude.
LAST MODIFIED:  Mar. 20, 2026

USAGE:
   insert3cbo.py [-3cdata DATA_3C] [-2cdata DATA_2C]
                  [-sdata STRUCT_DATA]
                  [-rawout RAW_23C_OUT]
                  [-sout STRUCT_23C_OUT]
                  [-h | --help]
"""

import argparse as ap
import math
import os
import sys
from datetime import datetime


# -----------------------------------------------------------
# Helper function
# -----------------------------------------------------------

def prep_line(file_handle):
    """Read one line, strip it, and split on whitespace.

    This is the Python equivalent of the Perl
    StructureControl::prepLine subroutine. It reads a
    single line from *file_handle*, strips leading/trailing
    whitespace, and splits on whitespace. Python's
    str.split() with no argument handles leading whitespace
    automatically, so no explicit shift of an empty leading
    token is needed.

    Parameters
    ----------
    file_handle : file object
        An open file from which to read one line.

    Returns
    -------
    list of str
        The whitespace-separated tokens from the line.
    """
    return file_handle.readline().strip().split()


# -----------------------------------------------------------
# ScriptSettings class
# -----------------------------------------------------------

class ScriptSettings():
    """The instance variables of this object are the user
    settings that control the program. The variable values
    are pulled from a list that is created within a
    resource control file and that are then reconciled
    with command line parameters."""

    def __init__(self):
        """Define default values for the parameters by
        pulling them from the resource control file in
        the default location:
        $OLCAO_RC/insert3cborc.py or from the current
        working directory if a local copy of
        insert3cborc.py is present."""

        # Read default variables from the resource
        # control file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit(
                "Error: $OLCAO_RC is not set. "
                "See instructions."
            )
        sys.path.insert(1, rc_dir)
        from insert3cborc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to settings from rc defaults.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the
        # rc file.
        self.reconcile(args)

        # Record the command line parameters used.
        self.record_clp()

    def assign_rc_defaults(self, default_rc):
        """Assign default values from the resource control
        file dictionary to instance variables."""

        # Input file names.
        self.data_3c = default_rc["data_3c"]
        self.data_2c = default_rc["data_2c"]
        self.struct_data = default_rc["struct_data"]

        # Output file names.
        self.raw_23c_out = default_rc["raw_23c_out"]
        self.struct_23c_out = default_rc["struct_23c_out"]

    def parse_command_line(self):
        """Set up and run the argument parser.

        Returns
        -------
        argparse.Namespace
            Parsed command line arguments.
        """

        # Create the parser tool.
        prog_name = "insert3cbo.py"

        description_text = """\
Merge two-center and three-center bond order data.

Creates new bond order raw and structure files that
include information describing three-center bonds.
Two-center bonds that are part of three-center bonds
are removed; genuine two-center bonds (those not
involved in any three-center bond) are retained.

FIX: This program is basically broken until it can be
remerged with the makeBOND program. The makeBOND
program now uses olcao.skl to populate the atom and
lattice data. This program worked by directly
modifying the structure.dat file instead."""

        epilog_text = """\
Please contact Paul Rulis regarding questions.
Defaults: ./insert3cborc.py or \
$OLCAO_RC/insert3cborc.py"""

        parser = ap.ArgumentParser(
            prog=prog_name,
            formatter_class=(
                ap.RawDescriptionHelpFormatter
            ),
            description=description_text,
            epilog=epilog_text,
        )

        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()

    def add_parser_arguments(self, parser):
        """Add all command line arguments to the parser.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            The parser to add arguments to.
        """

        # 3-center bond order data input file.
        parser.add_argument(
            '-3cdata', dest='data_3c',
            type=str, default=self.data_3c,
            help=(
                "Path to the 3-center bond order "
                "data file. "
                f"Default: {self.data_3c}"
            ),
        )

        # 2-center bond order data input file.
        parser.add_argument(
            '-2cdata', dest='data_2c',
            type=str, default=self.data_2c,
            help=(
                "Path to the 2-center bond order "
                "data file. "
                f"Default: {self.data_2c}"
            ),
        )

        # Structure data input file.
        parser.add_argument(
            '-sdata', dest='struct_data',
            type=str, default=self.struct_data,
            help=(
                "Path to the structure data file. "
                f"Default: {self.struct_data}"
            ),
        )

        # Combined raw output file.
        parser.add_argument(
            '-rawout', dest='raw_23c_out',
            type=str, default=self.raw_23c_out,
            help=(
                "Path to the combined 2+3 center "
                "bond order output file. "
                f"Default: {self.raw_23c_out}"
            ),
        )

        # Structure output file.
        parser.add_argument(
            '-sout', dest='struct_23c_out',
            type=str, default=self.struct_23c_out,
            help=(
                "Path to the output structure file "
                "with fake centroid atoms. "
                f"Default: {self.struct_23c_out}"
            ),
        )

    def reconcile(self, args):
        """Reconcile command line arguments with resource
        control file defaults.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed command line arguments.
        """
        self.data_3c = args.data_3c
        self.data_2c = args.data_2c
        self.struct_data = args.struct_data
        self.raw_23c_out = args.raw_23c_out
        self.struct_23c_out = args.struct_23c_out

    def record_clp(self):
        """Record the command line used to the 'command'
        file for provenance tracking."""
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime(
                "%b. %d, %Y: %H:%M:%S"
            )
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write("Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


# -----------------------------------------------------------
# Core processing functions
# -----------------------------------------------------------

def read_2cbo(settings):
    """Read the existing 2-center bond order file.

    Reads the two-center bond order data file into a
    dictionary of per-atom records. Each atom record
    contains identification data, orbital charge
    information, bonded atom lists, and bond angle data.

    The data structure uses 1-based atom numbering to
    match the conventions of the OLCAO bond order files.

    Parameters
    ----------
    settings : ScriptSettings
        The script settings containing file paths.

    Returns
    -------
    bo_data : dict
        Dictionary keyed by 1-based atom number. Each
        value is a dict with keys:
        - 'atom_num'             : int
        - 'system_num'           : int
        - 'element_name'         : str
        - 'element_id'           : int
        - 'species_id'           : int
        - 'type_id'              : int
        - 'atom_charge'          : float
        - 'num_orbital_charges'  : int
        - 'orbital_charges'      : list of [str, str]
        - 'num_bonded_atoms'     : int
        - 'bonded_atoms'         : list of
            [int id, float length, float bo]
        - 'num_bond_angles'      : int
        - 'bond_angles'          : list of str
    num_atoms : int
        Number of atoms in the system (before any
        3-center bond centroids are added).
    """

    bo_data = {}

    # Open the data file that contains the original
    # 2-center bond order data.
    with open(settings.data_2c, 'r') as f:

        # Get the number of atoms before three center
        # bonds are added.
        values = prep_line(f)
        num_atoms = int(values[-1])

        # For each atom we will read its information
        # into a dictionary structure.
        for atom in range(1, num_atoms + 1):
            rec = {}

            # ATOM_NUM (integer)
            values = prep_line(f)
            rec['atom_num'] = int(values[-1])

            # SYSTEM_NUM (integer)
            values = prep_line(f)
            rec['system_num'] = int(values[-1])

            # ELEMENT_NAME (character)
            values = prep_line(f)
            rec['element_name'] = values[-1]

            # ELEMENT_ID (integer)
            values = prep_line(f)
            rec['element_id'] = int(values[-1])

            # SPECIES_ID (integer)
            values = prep_line(f)
            rec['species_id'] = int(values[-1])

            # TYPE_ID (integer)
            values = prep_line(f)
            rec['type_id'] = int(values[-1])

            # ATOM_CHARGE (real)
            values = prep_line(f)
            rec['atom_charge'] = float(values[-1])

            # ATOM_ORBITAL_CHARGE count (integer)
            values = prep_line(f)
            num_orb = int(values[-1])
            rec['num_orbital_charges'] = num_orb

            # Read orbital charge entries. Each line
            # contains a label and a charge value.
            rec['orbital_charges'] = []
            for _ in range(num_orb):
                values = prep_line(f)
                rec['orbital_charges'].append(
                    [values[0], values[1]]
                )

            # NUM_BONDED_ATOMS (integer)
            values = prep_line(f)
            num_bonded = int(values[-1])
            rec['num_bonded_atoms'] = num_bonded

            # Read bonded atom entries. Each line:
            #   bonded_atom_ID  bond_length  bond_order
            rec['bonded_atoms'] = []
            for _ in range(num_bonded):
                values = prep_line(f)
                rec['bonded_atoms'].append([
                    int(values[0]),
                    float(values[1]),
                    float(values[2]),
                ])

            # NUM_BOND_ANGLES (integer)
            values = prep_line(f)
            num_angles = int(values[-1])
            rec['num_bond_angles'] = num_angles

            # Read bond angle entries (stored as raw
            # strings to preserve their exact format).
            rec['bond_angles'] = []
            for _ in range(num_angles):
                line = f.readline().rstrip('\n')
                rec['bond_angles'].append(line)

            bo_data[atom] = rec

    return bo_data, num_atoms


def _remove_bond(atom_rec, target_atom_id):
    """Remove a bond to *target_atom_id* from *atom_rec*.

    Searches the bonded_atoms list in *atom_rec* for a
    bond whose first element (atom ID) matches
    *target_atom_id*. If found, removes it and decrements
    the bonded atom count.

    This implements the same logic as the Perl original's
    loop-and-splice pattern for removing two-center bonds
    that are subsumed by a three-center bond.

    Parameters
    ----------
    atom_rec : dict
        A single atom's bond order data record.
    target_atom_id : int
        The atom ID to search for and remove.
    """
    for i, bond in enumerate(
        atom_rec['bonded_atoms']
    ):
        if bond[0] == target_atom_id:
            del atom_rec['bonded_atoms'][i]
            atom_rec['num_bonded_atoms'] -= 1
            break


def read_3cbo(settings, bo_data, num_atoms):
    """Read 3-center bond order data and modify 2CBO.

    Opens the three-center bond order file and reads its
    contents. As it is read we process the existing
    2-center bond order data structure by:

    1. Deleting 2-center bonds that are actually part
       of three-center bonds.
    2. Adding "2-center" bonds where the bonded "atom"
       is actually the centroid of the three-center
       bond.
    3. Extending the bond order data structure to
       include the centroids of three-center bonds as
       "atoms" with placeholder data.

    Additionally, the bond-weighted centroid positions
    (shifted inside the unit cell) are recorded for
    later insertion into the structure file.

    Parameters
    ----------
    settings : ScriptSettings
        The script settings containing file paths.
    bo_data : dict
        The 2-center bond order data dictionary, which
        will be modified in place.
    num_atoms : int
        Number of real atoms in the system.

    Returns
    -------
    bo_data : dict
        The modified bond order data dictionary, now
        including fake centroid "atoms".
    num_3c_bonds : int
        Number of three-center bonds found.
    bo3c_positions : dict
        Dictionary keyed by 1-based 3C bond index.
        Each value is [x, y, z] coordinates of the
        bond-weighted centroid (shifted inside the
        unit cell).
    """

    # Okay, this is stupid and at some point should be
    # improved.  We have to read the last line in the
    # file to get the number of three center bonds.
    with open(settings.data_3c, 'r') as f_tmp:
        for last_line in f_tmp:
            pass
    last_values = last_line.strip().split()
    num_3c_bonds = int(last_values[-1])

    bo3c_positions = {}

    # Open the three center bond order data file.
    with open(settings.data_3c, 'r') as f:

        # Read past the header.
        f.readline()

        # Read the lattice vectors. Each line has the
        # format:
        #   LATTICE_X  mag  angle  ax  ay  az
        # where indices 3, 4, 5 are the x, y, z
        # components of the lattice vector.
        lattice_v = [[0.0] * 4 for _ in range(4)]
        for vec in range(1, 4):
            values = prep_line(f)
            lattice_v[vec][1] = float(values[3])
            lattice_v[vec][2] = float(values[4])
            lattice_v[vec][3] = float(values[5])

        # Read information for each of the three center
        # bonds and modify the two center bond data
        # structure accordingly.
        for bo3c in range(1, num_3c_bonds + 1):

            # Read which atoms are involved in this
            # three center bond and call them atom1,
            # atom2, and atom3.
            values = prep_line(f)
            atom1 = int(values[1])
            atom2 = int(values[2])
            atom3 = int(values[3])

            # For the first atom, look in its list of
            # bonded atoms and remove from that list
            # any atom number that participates in this
            # three center bond. Note that we will have
            # to decrement the number of bonds that
            # this atom has for each bond that is
            # removed.
            _remove_bond(bo_data[atom1], atom2)
            _remove_bond(bo_data[atom1], atom3)

            # Now we have to do essentially the same
            # thing for the second atom but we only
            # have to search for the third atom because
            # the first atom will always have a lower
            # index number and thus will not be present
            # in the bond list.
            _remove_bond(bo_data[atom2], atom3)

            # Read the LATTICE_INDICES_1,2,3. The
            # LATTICE_INDICES are used to identify
            # which replicated cell a particular atom
            # should be in so as to achieve minimal
            # distances between the three atoms of a
            # possible three-center bond. The
            # LATTICE_INDICES of atom1 should always
            # be zero because we will force that atom
            # to always be inside the central cell.
            # The atom2 indices will position the atom
            # into an appropriate replicated cell so
            # as to be minimally distant from atom1.
            # The same is true for the atom3 indices.
            # (I.e., minimally distant from atom1.)
            lat_idx = [[0] * 4 for _ in range(4)]
            for a in range(1, 4):
                values = prep_line(f)
                lat_idx[a][1] = int(values[1])
                lat_idx[a][2] = int(values[2])
                lat_idx[a][3] = int(values[3])

            # Get the x,y,z coordinates for each atom
            # in the three center bond (au). Apply the
            # lattice index transformation to position
            # each atom in the correct replicated cell.
            coords = [
                [0.0] * 4 for _ in range(4)
            ]
            for a in range(1, 4):
                values = prep_line(f)
                for xyz in range(1, 4):
                    coords[xyz][a] = (
                        float(values[xyz])
                        + lat_idx[a][1]
                        * lattice_v[1][xyz]
                        + lat_idx[a][2]
                        * lattice_v[2][xyz]
                        + lat_idx[a][3]
                        * lattice_v[3][xyz]
                    )

            # Read past the geometric centroid.
            f.readline()

            # Get the bond weighted centroid coordinates
            # and compute the displacement between each
            # atom of the three center bond and the
            # weighted centroid.
            values = prep_line(f)
            cx = float(values[1])
            cy = float(values[2])
            cz = float(values[3])
            disp = [
                [0.0] * 4 for _ in range(4)
            ]
            for a in range(1, 4):
                disp[1][a] = coords[1][a] - cx
                disp[2][a] = coords[2][a] - cy
                disp[3][a] = coords[3][a] - cz

            # Compute the "bond length" between the
            # atom and the bond centroid.
            bond_len = [0.0] * 4
            for a in range(1, 4):
                bond_len[a] = math.sqrt(
                    disp[1][a] ** 2
                    + disp[2][a] ** 2
                    + disp[3][a] ** 2
                )

            # Read the position of the bond weighted
            # centroid that was shifted to be inside
            # the unit cell.
            values = prep_line(f)
            bo3c_positions[bo3c] = [
                float(values[1]),
                float(values[2]),
                float(values[3]),
            ]

            # Read past the 2C bond order pairs.
            f.readline()

            # Read in the 3C bond order.
            values = prep_line(f)
            current_3cbo = float(values[1])

            # Now, we are going to add a new bond to a
            # fake atom that represents the centroid of
            # the three center bond. Therefore, we will
            # first increment the number of bonds that
            # each atom has.
            next_3c_atom = num_atoms + bo3c

            # Record the data for the new bond from
            # each participating atom to the fake
            # centroid atom.
            for a_id, bl in [
                (atom1, bond_len[1]),
                (atom2, bond_len[2]),
                (atom3, bond_len[3]),
            ]:
                bo_data[a_id][
                    'bonded_atoms'
                ].append([
                    next_3c_atom, bl, current_3cbo
                ])
                bo_data[a_id][
                    'num_bonded_atoms'
                ] += 1

            # Finally, we have to add the "atoms" to
            # the list of bo_data. These are fake
            # entries representing the centroid of the
            # three center bond.
            bo_data[next_3c_atom] = {
                'atom_num': next_3c_atom,
                'system_num': 1,
                'element_name': "3c",
                'element_id': 0,
                'species_id': 0,
                'type_id': 0,
                'atom_charge': 0.0,
                'num_orbital_charges': 0,
                'orbital_charges': [],
                'num_bonded_atoms': 0,
                'bonded_atoms': [],
                'num_bond_angles': 0,
                'bond_angles': [],
            }

    return bo_data, num_3c_bonds, bo3c_positions


def write_2cbo(settings, bo_data, num_atoms,
               num_3c_bonds):
    """Write out the modified bond order data.

    Writes the combined two-center and three-center bond
    order data to a new file. This includes the original
    2-center bonds (with some removed), new bonds to
    three-center centroids, and fake "atom" entries for
    each centroid.

    For now we open a differently named file to avoid
    accidentally destroying the original data.

    Parameters
    ----------
    settings : ScriptSettings
        The script settings containing file paths.
    bo_data : dict
        The complete bond order data dictionary.
    num_atoms : int
        Number of real atoms.
    num_3c_bonds : int
        Number of three-center bonds (fake atoms).
    """

    # Compute the number of atoms (combined real and
    # fake).
    atom_number = num_atoms + num_3c_bonds

    # Open the new file.
    with open(settings.raw_23c_out, 'w') as f:

        # Print out the new boData.
        f.write(f"NUM_ATOMS {atom_number}\n")

        for atom in range(1, atom_number + 1):
            rec = bo_data[atom]

            f.write(
                "ATOM_NUM            "
                f"{rec['atom_num']}\n"
            )
            f.write(
                "SYSTEM_NUM          "
                f"{rec['system_num']}\n"
            )
            f.write(
                "ELEMENT_NAME        "
                f"{rec['element_name']}\n"
            )
            f.write(
                "ELEMENT_ID          "
                f"{rec['element_id']}\n"
            )
            f.write(
                "SPECIES_ID          "
                f"{rec['species_id']}\n"
            )
            f.write(
                "TYPE_ID             "
                f"{rec['type_id']}\n"
            )
            f.write(
                "ATOM_CHARGE         "
                f"{rec['atom_charge']}\n"
            )
            f.write(
                "ATOM_ORBITAL_CHARGE "
                f"{rec['num_orbital_charges']}\n"
            )

            for orb in rec['orbital_charges']:
                f.write(
                    f"{orb[0]} {orb[1]}\n"
                )

            f.write(
                "NUM_BONDED_ATOMS    "
                f"{rec['num_bonded_atoms']}\n"
            )

            for bond in rec['bonded_atoms']:
                f.write(
                    f"{bond[0]} "
                    f"{bond[1]} "
                    f"{bond[2]}\n"
                )

            f.write(
                "NUM_BOND_ANGLES     "
                f"{rec['num_bond_angles']}\n"
            )

            for angle in rec['bond_angles']:
                f.write(f"{angle}\n")


def add_to_struct(settings, num_atoms, num_3c_bonds,
                  bo3c_positions):
    """Add fake centroid atoms to the structure file.

    Now we have to add the fake atoms to the structure
    file so that when we try to plot a figure using
    makeBOND these "atoms" will be found and plotted.

    The fake atoms represent the centroids of the three-
    center bonds. They are added to both the atomic site
    list and the potential site list in the structure
    file.

    Parameters
    ----------
    settings : ScriptSettings
        The script settings containing file paths.
    num_atoms : int
        Number of real atoms.
    num_3c_bonds : int
        Number of three-center bonds (fake atoms to
        add).
    bo3c_positions : dict
        Centroid positions keyed by 1-based 3C bond
        index, each value being [x, y, z].
    """

    total_atoms = num_atoms + num_3c_bonds

    # Open the structure file for reading.
    with open(settings.struct_data, 'r') as fin:

        # Open a new file that will contain the
        # modified structure.
        with open(
            settings.struct_23c_out, 'w'
        ) as fout:

            # Read and reproduce the first few lines
            # (5 header lines: title, space group,
            # lattice vectors, etc.).
            for _ in range(5):
                fout.write(fin.readline())

            # Read the line containing the number of
            # atoms and replace with the new total
            # (including fake centroid atoms).
            fin.readline()
            fout.write(f"{total_atoms}\n")

            # Read past and reproduce the
            # num_type_x_y_z_elem tag line.
            fout.write(fin.readline())

            # Read and reproduce all of the atomic
            # sites.
            last_line = ""
            for _ in range(num_atoms):
                last_line = fin.readline()
                fout.write(last_line)

            # Extract the type number of the last
            # atom.
            last_vals = last_line.strip().split()
            last_type = int(last_vals[1])

            # Add the new "atoms" to the list. These
            # are fake entries for the centroids of
            # the three-center bonds.
            type_number = last_type + 1
            for a in range(1, num_3c_bonds + 1):
                atom_number = num_atoms + a
                pos = bo3c_positions[a]
                fout.write(
                    f" {atom_number}"
                    f"  {type_number}"
                    f" {pos[0]}"
                    f" {pos[1]}"
                    f" {pos[2]} 3c\n"
                )

            # Read and reproduce the potential site
            # tag.
            fout.write(fin.readline())

            # Print the new number of potential sites
            # (including the fake centroid "atoms").
            fin.readline()
            fout.write(f"{total_atoms}\n")

            # Read past and reproduce the
            # num_type_x_y_z_elem tag line.
            fout.write(fin.readline())

            # Read and reproduce all of the potential
            # sites.
            for _ in range(num_atoms):
                fout.write(fin.readline())

            # Add the new "atoms" to the potential
            # site list as well.
            for a in range(1, num_3c_bonds + 1):
                atom_number = num_atoms + a
                pos = bo3c_positions[a]
                fout.write(
                    f" {atom_number}"
                    f"  {type_number}"
                    f" {pos[0]}"
                    f" {pos[1]}"
                    f" {pos[2]} 3c\n"
                )


# -----------------------------------------------------------
# Main entry point
# -----------------------------------------------------------

def main():
    """Main program flow.

    1. Parse command line and set up file paths.
    2. Read the 2-center bond order data.
    3. Read the 3-center bond order data and modify the
       2-center data structure accordingly.
    4. Write the combined bond order data file.
    5. Write the modified structure file with fake
       centroid atoms.
    6. Print a recommended makeBOND command.
    """

    print("\n\nScript Executing.\n")

    # Define file names, initialize global variables,
    # and parse the command line parameters.
    now = datetime.now().strftime(
        "%b. %d, %Y: %H:%M:%S"
    )
    print(
        "\nInitializing the environment"
        f".................{now}."
    )
    now = datetime.now().strftime(
        "%b. %d, %Y: %H:%M:%S"
    )
    print(
        "\nParse the command line parameters"
        f"............{now}."
    )
    settings = ScriptSettings()

    # Read in the 2 center bond order data.
    now = datetime.now().strftime(
        "%b. %d, %Y: %H:%M:%S"
    )
    print(
        "\nReading in the 2-center bond order "
        f"data......{now}."
    )
    bo_data, num_atoms = read_2cbo(settings)

    # Read in the 3 center bond order data and modify
    # the data structures that were created to hold the
    # two center bond order information.
    now = datetime.now().strftime(
        "%b. %d, %Y: %H:%M:%S"
    )
    print(
        "\nReading in the 3-center bond order "
        f"data......{now}."
    )
    bo_data, num_3c_bonds, bo3c_positions = (
        read_3cbo(settings, bo_data, num_atoms)
    )

    # Print out the new bond order data file that
    # contains 2 and 3 center bonds. Note that many
    # 2-center bonds may have been removed and replaced
    # with 3-center bonds.
    now = datetime.now().strftime(
        "%b. %d, %Y: %H:%M:%S"
    )
    print(
        "\nWriting the new bond order data"
        f"..............{now}."
    )
    write_2cbo(
        settings, bo_data, num_atoms, num_3c_bonds
    )

    # Print out a new structure file with extra "atoms"
    # added that represent the centroids of the three
    # center bonds.
    now = datetime.now().strftime(
        "%b. %d, %Y: %H:%M:%S"
    )
    print(
        "\nWriting the new structure file"
        f"...............{now}."
    )
    add_to_struct(
        settings, num_atoms, num_3c_bonds,
        bo3c_positions,
    )

    # Print a recommended command line for makeBOND.
    print(
        "\n\nA possible command for the "
        "makeBOND script is:"
    )
    print(
        "makeBOND -pos structure.23cbo.dat "
        "-data gs_bond-mb.23cbo.raw"
    )
    print("-gc groupControl -minbo 0.0 -model")
    print(
        "Note that groupControl contains "
        "only \"SYSTEM_NUM\""
    )


if __name__ == '__main__':
    # Everything before this point was a subroutine
    # definition or a request to import information
    # from external modules. Only now do we actually
    # start running the program. The purpose of this
    # is to allow another python program to import
    # *this* script and call its functions internally.
    main()
