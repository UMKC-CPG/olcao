#!/usr/bin/env python3

"""porosityrc.py -- Default parameters for the porosity.py script.

This resource control file defines the default values for all
parameters used by porosity.py.  Users can place a customised
copy in their $OLCAO_RC directory or in the current working
directory to override the built-in defaults without modifying
the script itself.

Each key in the dictionary returned by parameters_and_defaults()
corresponds to a command-line option or internal default in
porosity.py.  Comments beside each key explain the meaning and
acceptable values.
"""


def parameters_and_defaults():
    """Return a dictionary of all porosity default parameters.

    The keys are grouped into three categories:

    1. **Input file name** -- the skeleton structure file to
       read.
    2. **Mesh / resolution control** -- how the sampling grid
       is defined (either by explicit mesh counts or by a
       target resolution).
    3. **Physical parameters** -- the interaction radius and
       scale factor for covalent radii that define pore
       boundaries.
    """
    param_dict = {

        # ---- Input file ----

        # Path to the OLCAO skeleton (.skl) input file
        #   containing the structural description.
        "in_file": "olcao.skl",

        # ---- Mesh / resolution control ----

        # Number of sampling points along each a, b, c
        #   direction.  If given on the command line with
        #   -m, the resolution is computed from these.
        #   1-indexed: [None, na, nb, nc] to match the Perl
        #   @numMeshPoints[1..3] convention.  [None, 0, 0, 0]
        #   means "not yet specified".
        "num_mesh_points": [None, 0, 0, 0],

        # Target resolution (Angstroms) for the sampling
        #   mesh.  If given with -r, the mesh point counts
        #   are derived from this.  A value of 0.0 means
        #   "not yet specified".
        "resolution": 0.0,

        # ---- Physical parameters ----

        # Radius (Angstroms) of the sphere within which
        #   atoms can contribute to the pore map.  This
        #   defines how far from a mesh point to look for
        #   atoms when determining whether the point is
        #   inside or outside an atom.
        "limit_dist": 4.0,

        # Multiplicative scale factor applied to the
        #   covalent radii from the element database.
        #   The scaled atomic radius defines the pore
        #   boundary: mesh points within this radius of
        #   any atom are considered solid (not void).
        "scale_factor": 1.0,
    }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
