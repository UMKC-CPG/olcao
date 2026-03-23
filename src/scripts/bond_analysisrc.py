#!/usr/bin/env python3

"""bond_analysisrc.py -- Default parameters for bond_analysis.py.

This resource control file provides a single function that returns a
dictionary of all default parameters used by bond_analysis.  Users may
place a customized copy of this file in the current working directory
or in the directory pointed to by $OLCAO_RC to override any defaults.

The parameters control how geometric bond analysis is performed on an
atomic structure.  They set the input/output files, the analysis mode,
coordinate system, bonding criteria, histogram options, 3D printing
scales, Q^n distribution species, ring distribution bounds, and the
spatial bounding box.
"""


def parameters_and_defaults():

    param_dict = {

        # ---- Input / output files ----
        # Default input file (olcao.skl format).
        "in_file":             "olcao.skl",
        # Base name for output; a suffix is appended based on the
        # analysis operation (e.g. ".bl", ".ba", ".bs", etc.).
        "out_file":            "bondAnalysis",
        # Bond order raw data file used when -s flag is given.
        "bond_order_file":     "S1__S1.dat",

        # ---- Analysis operation ----
        # Which analysis mode to run.  Allowed values:
        #   1 = ball-and-stick     (-bs)
        #   2 = bond angles        (-ba)
        #   3 = bond lengths       (-bl)  [default]
        #   4 = bond orient. order (-boo)
        #   5 = OpenDX             (-dx)
        #   6 = statistics         (-st)
        #   7 = coordination       (-co)
        #   8 = 3D full model      (-3dfull)
        #   9 = 3D parts model     (-3dparts)
        #  11 = VTK                (-vtk)
        #  12 = Q^n distribution   (-qn)
        #  13 = ring distribution  (-rn)
        "operation":           3,

        # ---- Bonding criteria ----
        # Maximum distance (angstroms) beyond which atom pairs are
        # not considered.
        "limit_dist":          4.0,
        # Multiplicative factor applied to the sum of the covalent
        # radii to decide if two atoms are bonded.  A value of 1.1
        # means 10% more than the sum of the covalent radii.
        "bonding_factor":      1.1,
        # Whether to use OLCAO bond order data (True/False).
        "use_olcao_bo":        False,

        # ---- Coordinate type ----
        # 1 = direct x,y,z (cartesian)
        # 2 = direct a,b,c
        # 3 = fractional a,b,c  [default]
        "coord_type":          3,

        # ---- Bond orientational order (BOO) ----
        # The l quantum number for the Y_lm spherical harmonic used
        # in BOO analysis.  Allowed values: 0, 1, 6.
        "ylm_l":               6,

        # ---- 3D printing scale factors ----
        # Atom radius scale: 1 angstrom = a_scale mm.
        "a_scale":             10.0,
        # Overall crystal lattice size scale.
        "c_scale":             12.0,
        # Bond cylinder diameter scale (multiplied by 10% of the
        # largest covalent radius).
        "b_scale":             1.0,
        # Font size scale (multiplied by bond cylinder radius).
        "f_scale":             0.9,
        # Ratio of bond rod radius to covalent radius.
        "bond_radius_scale":   0.1,

        # ---- 3D printing volume (mm) ----
        # Printable volume dimensions for -3dparts. These are set
        # to 5% less than the max volume of a Prusa Mendel 3D
        # printer.
        "x_vol_3d":            150.0,
        "y_vol_3d":            150.0,
        "z_vol_3d":            150.0,

        # ---- Bounding box ----
        # Spatial box that constrains which atoms are considered.
        # Each axis has [min, max].  Use "max" as the max value to
        # automatically include the full extent of the cell.
        "box_min1":            0.0,
        "box_max1":            "max",
        "box_min2":            0.0,
        "box_max2":            "max",
        "box_min3":            0.0,
        "box_max3":            "max",
        # Zone: 1 = only atoms inside the box, 2 = only outside.
        "zone":                1,

        # ---- Histogram ----
        # Whether to produce a histogram alongside the output.
        "do_hist":             False,
        # Bin width for histograms.
        "hist_delta":          0.1,

        # ---- Q^n distribution ----
        # Bridging atoms (lower-case element names).
        "bridges":             ["o"],
        # Coordinated cations (lower-case element names).
        "ions":                ["si"],

        # ---- Ring distribution ----
        # Minimum and maximum ring lengths to consider.
        "min_ring_len":        3,
        "max_ring_len":        8,
    }

    return param_dict


if __name__ == "__main__":
    for k, v in parameters_and_defaults().items():
        print(f"  {k:24s} = {v!r}")
