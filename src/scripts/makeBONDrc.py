#!/usr/bin/env python3

"""makeBONDrc.py -- Default parameters for the makeBOND.py script.

This resource control file defines the default values for all
parameters used by makeBOND.py.  Users can place a customised copy
in their $OLCAO_RC directory or in the current working directory
to override the built-in defaults without modifying the script
itself.

Each key in the dictionary returned by parameters_and_defaults()
corresponds to a command-line option or internal default in
makeBOND.py.  Comments beside each key explain the meaning and
acceptable values.
"""


def parameters_and_defaults():
    """Return a dictionary of all makeBOND default parameters.

    The keys are grouped into four categories:

    1. **Input file names** -- paths to the data and position
       files that makeBOND reads.
    2. **Control files** -- optional files that restrict which
       bonds are printed or how the output is grouped.
    3. **Output mode flags** -- which output format to produce
       (scatter, openDX model, VTK, mesh, bond profile).
    4. **Numerical thresholds and tuning parameters** -- bond
       order / bond length limits, radius factors, mesh
       resolution, etc.
    """
    param_dict = {

        # ---- Input file names ----

        # Path to the raw bond order data file produced by
        #   the OLCAO bond program.
        "data_file": "gs_bond-mb.raw",

        # Path to the structure file containing lattice
        #   vectors and atom positions (OLCAO format).
        "pos_file": "structure.dat",

        # ---- Control files ----

        # Path to a bond control file that filters bonds by
        #   element-pair-specific BO/BL ranges.  Empty string
        #   means no control file is used.
        "bond_control_file": "",

        # Path to a group control file that defines how atoms
        #   and bonds are tagged for output grouping.  Empty
        #   string means the default grouping (ELEMENT_NAME
        #   and SPECIES_ID) is used.
        "group_control_file": "",

        # Name of the bond profile output file.
        "bond_profile_file": "bondProfile.plot",

        # ---- Output mode flags ----
        #   Only one output mode is active at a time.  If none
        #   is explicitly requested the default is scatter.

        # Produce scatter-plot files grouped by atom/bond tag.
        "scatter_plot": False,

        # Compute and print dipole information.
        "dipole": False,

        # Produce an openDX structural model with coloured
        #   bonds and atoms.
        "system_model_dx": False,

        # Show three-centre bonds in the openDX model.
        "bo3c": False,

        # Produce a VTK Legacy format structural model.
        "system_model_vtk": False,

        # Produce a 3-D openDX mesh of bond order values.
        "mesh_dx": False,

        # Whether the openDX model is atom-based (True) or
        #   bond-based (False).
        "atom_based_dx": True,

        # Whether the openDX model is bond-based.
        "bond_based_dx": False,

        # Produce a bond-profile curve along a lattice axis.
        "bond_profile_plot": False,

        # The lattice axis for the bond profile ("a","b","c").
        "profile_axis": "",

        # ---- Bond order / bond length thresholds ----

        # Minimum bond order to include in output.
        "min_bo": -1000.0,

        # Maximum bond order to include in output.
        "max_bo": 1000.0,

        # Minimum bond length to include in output.
        "min_bl": 0.0,

        # Maximum bond length to include in output.
        "max_bl": 1000.0,

        # ---- Model tuning parameters ----

        # Multiplicative factor applied to covalent radii
        #   when drawing atom spheres in the openDX model.
        "radius_factor": 1.0,

        # Use grey-scale colours for atom spheres instead
        #   of the default element-based colour palette.
        "grey_scale": False,

        # ---- Mesh parameters ----

        # Atom index whose properties are subtracted from
        #   mesh values to create a difference mesh.  Zero
        #   means no subtraction.
        "comp_atom": 0,

        # Full width at half maximum of the Gaussian
        #   weighting function used when evaluating the
        #   mesh (in Angstroms).
        "fwhm": 2.0,

        # Number of mesh grid points along the a, b, c
        #   lattice directions.
        "num_mesh_points": [10, 10, 10],

        # Maximum distance (Angstroms) between a bond
        #   midpoint and a mesh point for the bond to
        #   contribute to that mesh point.
        "limit_dist": 4.0,
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
