#!/usr/bin/env python3

"""structure_control.py -- Structure data container and manipulation tools.

This module is the Python equivalent of the Perl StructureControl.pm module.
It is the central library for all OLCAO pre- and post-processing scripts.

The single class StructureControl holds the full description of an atomic
structure (lattice, atom positions, species, bonding, ...) and provides
methods for reading, writing, and analysing that structure.

TYPICAL USAGE::

    from structure_control import StructureControl

    sc = StructureControl()
    sc.read_olcao_skl()           # read olcao.skl in the current directory
    sc.assign_coval_radii()       # pull covalent radii from ElementData
    sc.create_bonding_list()      # compute bonds based on covalent radii
    sc.print_olcao()              # write the structure back out

ROLE IN THE OLCAO WORKFLOW:
    Scripts that use StructureControl fall into three categories:
    1. Input preparation   -- convert external formats to olcao.skl
    2. Execution control   -- set up and launch OLCAO Fortran calculations
    3. Results processing  -- read output files and compute derived quantities

INDEXING CONVENTION:
    All atom-indexed and element-indexed arrays are 1-indexed, matching the
    OLCAO Fortran source and the original Perl code.  Index 0 of every such
    list is a None placeholder and must not be used.  Loop ranges therefore
    read ``range(1, num_atoms + 1)``, which keeps the loop variable directly
    interpretable as the physical atom or element number.

NOTE ON GETTERS:
    The Perl module exposed ~80 getter subroutines that returned references to
    package-global arrays.  In Python, all of that data is stored as instance
    attributes on the StructureControl object and can be read directly
    (e.g. ``sc.num_atoms``, ``sc.fract_abc``).  Getter methods are therefore
    not reproduced here.  Setter methods are kept where they perform
    non-trivial initialisation beyond a simple assignment.
"""

import os
import math
import re
import random
import subprocess

from element_data import ElementData


# ---------------------------------------------------------------------------
# Module-level physical and numerical constants
# ---------------------------------------------------------------------------

PI = 3.1415926535897932384626433832795
EPSILON = 0.00001           # Numerical tolerance for coordinate comparisons
BIG_INT = 1_000_000_000     # Large integer sentinel value
BIG_REAL = 1_000_000_000.0  # Large float sentinel value
BOHR_RAD = 0.5291772180     # Bohr radius in Angstroms
HARTREE = 27.21138386       # Hartree energy in eV


# ---------------------------------------------------------------------------
# StructureControl class
# ---------------------------------------------------------------------------

class StructureControl:
    """Atomic structure container and manipulation library for OLCAO.

    All state describing the current structure is held as instance attributes,
    organised into the groups listed below.  Methods are arranged in matching
    sections further down the file.

    INSTANCE ATTRIBUTE GROUPS
    -------------------------

    Element data (from ElementData)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    element_data : ElementData
        The fully initialised periodic-table database object.

    File / path names
    ~~~~~~~~~~~~~~~~~
    olcao_skl : str
        Path to the OLCAO skeleton input file (default: 'olcao.skl').
    model_xyz : str
        Path to an XYZ structure file (default: 'model.xyz').
    model_pdb : str
        Path to a PDB structure file (default: 'model.pdb').
    model_hin : str
        Path to a HIN structure file (default: 'model.hin').
    space_group_db : str
        Path to the space group database ($OLCAO_DATA/spaceDB).
    atomic_bdb : str
        Path to the atomic basis set database ($OLCAO_DATA/atomicBDB).
    atomic_pdb : str
        Path to the atomic potential database ($OLCAO_DATA/atomicPDB).
    title : str
        Single-line title for the structure.
    system_title : list of str
        Multi-line title block read from the skeleton file.  1-indexed.

    Atom counts and identity
    ~~~~~~~~~~~~~~~~~~~~~~~~
    num_atoms : int
        Number of atoms in the central simulation cell.
    num_atoms_ext : int
        Number of atoms in the extended (replicated) cell used for PBC.
    num_elements : int
        Number of unique chemical elements in the system.
    num_atom_uniq_species : int
        Number of unique atomic species (element + local environment tag).
    atom_uniq_species : list
        Unique species index for each atom.  1-indexed.
    atom_element_name : list of str
        Lower-case element abbreviation for each atom.  1-indexed.
    atom_element_id : list of int
        Atomic number Z for each atom.  1-indexed.
    atomic_z : list of int
        Atomic number Z for each atom (synonym, kept for Perl parity).
        1-indexed.
    element_list : list of str
        Unique element abbreviations present in the system.  1-indexed.
    element_count : list of int
        Number of atoms of each unique element.  1-indexed.
    atom_species_id : list of int
        Species ID for each atom.  1-indexed.
    num_species : list of int
        Number of species for each element.  1-indexed.
    species_list : list
        Unique species tags for each element.  1-indexed.
    atom_tag : list of str
        Atom name/tag from HIN or PDB files.  1-indexed.
    connection_tag : list
        List of atoms connected in an undefined way to each atom.  1-indexed.

    Atom-to-atom mapping
    ~~~~~~~~~~~~~~~~~~~~
    dat_skl_map : list of int
        dat_skl_map[dat_atom] = skl_atom.  1-indexed.
    skl_dat_map : list of int
        skl_dat_map[skl_atom] = dat_atom.  1-indexed.
    central_to_ext_item_map : list
        Maps central cell atom index to its extended list entries.  1-indexed.
    ext_to_central_item_map : list of int
        Maps extended atom index to the central cell atom.  1-indexed.

    Atom coordinates
    ~~~~~~~~~~~~~~~~
    fract_abc : list of list
        fract_abc[atom] = [a, b, c] fractional coordinates.  1-indexed.
    direct_abc : list of list
        direct_abc[atom] = [a, b, c] direct-space coords along lattice.
        1-indexed.
    direct_xyz : list of list
        direct_xyz[atom] = [x, y, z] Cartesian coordinates (Angstroms).
        1-indexed.
    ext_direct_xyz : list
        Cartesian coordinates for atoms in each replicated cell, structured
        by cell index.  1-indexed.
    ext_direct_xyz_list : list of list
        Flat list of [x, y, z] for every atom in the extended system.
        1-indexed.
    ext_fract_abc_list : list of list
        Flat list of [a, b, c] fractional coords in the extended system.
        1-indexed.
    max_pos : list of float
        [x_max, y_max, z_max] bounding box maximum.
    min_pos : list of float
        [x_min, y_min, z_min] bounding box minimum.

    Per-atom elemental properties
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    coval_radii : list of float
        Covalent radius for each atom (drawn from ElementData).  1-indexed.
    color_vtk : list of list
        VTK [r, g, b, alpha] color for each atom.  1-indexed.

    Molecule / residue labels (PDB / HIN)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    molecule_name : list of str
        Molecule name for each atom.  1-indexed.
    molecule_seq_num : list of int
        Molecule sequence number for each atom.  1-indexed.
    residue_name : list of str
        Residue name for each atom.  1-indexed.
    residue_seq_num : list of int
        Residue sequence number for each atom.  1-indexed.

    Q^n network analysis
    ~~~~~~~~~~~~~~~~~~~~
    atom_qn : list of int
        Q^n value for each atom (Q^n = number of non-bridging oxygens).
        1-indexed.
    absolute_sys_qn : list of int
        Sum of each Q^n value over all atoms.  1-indexed.
    fractional_sys_qn : list of float
        absolute_sys_qn / num_atoms for each Q^n.  1-indexed.
    num_qn_atoms : int
        Number of atoms that are Q^n target atoms.
    net_conn_qn : float
        Network connectivity = sum_i(i * fractional_sys_qn[i]).

    Ring statistics
    ~~~~~~~~~~~~~~~
    atom_rings : list
        Ring indices that each atom belongs to.  1-indexed.
    rings : list
        Each entry is a ring (list of atom indices).  1-indexed.
    ring_counts : list of int
        Count of rings of each length.  1-indexed.

    Coordination
    ~~~~~~~~~~~~
    coordination : list
        Elemental coordination list for each atom.  1-indexed.
    coordination_summary : list
        Summary of elemental coordination for each element.  1-indexed.

    Lattice parameters
    ~~~~~~~~~~~~~~~~~~
    real_lattice : list of list
        real_lattice[abc_axis][xyz_axis]: both 1-indexed (1..3).
        real_lattice[1]=[None,ax,ay,az] (a vector), [2]=b, [3]=c.
        Index 0 of the outer list and index 0 of each inner list are None.
        Lattice vectors in Angstroms.
    real_lattice_inv : list of list
        Inverse of the real lattice matrix.  Same 1-indexed layout.
    recip_lattice : list of list
        Reciprocal lattice vectors (2*pi * real_lattice_inv^T).  1-indexed.
    mag : list of float
        1-indexed: mag[1]=|a|, mag[2]=|b|, mag[3]=|c| in Angstroms.
        mag[0] is None.
    angle : list of float
        1-indexed: angle[1]=alpha, [2]=beta, [3]=gamma in radians.
        alpha=angle(b,c); beta=angle(a,c); gamma=angle(a,b).  Index 0 None.
    angle_deg : list of float
        Same angles in degrees.  1-indexed; index 0 is None.
    sine_angle : list of float
        1-indexed: sine_angle[1]=sin(alpha), [2]=sin(beta), [3]=sin(gamma).
    full_cell_mag : list of float
        Conventional (full) cell lattice magnitudes.  1-indexed.
    full_cell_angle : list of float
        Conventional cell angles in radians.  1-indexed.
    full_cell_angle_deg : list of float
        Conventional cell angles in degrees.  1-indexed.
    mag_recip : list of float
        Magnitudes of the reciprocal lattice vectors.  1-indexed.
    angle_recip : list of float
        Reciprocal lattice angles in radians.  1-indexed.
    angle_deg_recip : list of float
        Reciprocal lattice angles in degrees.  1-indexed.
    real_cell_volume : float
        Volume of the real-space cell in Angstroms^3.
    recip_cell_volume : float
        Volume of the reciprocal-space cell.
    abc_order : list of int
        1-indexed: maps lattice vector index to sorted order.
    xyz_order : list of int
        1-indexed: maps Cartesian axis index to sorted order.
    sorted_lat_indices : list of int
        Indices used to sort lattice angles and magnitudes.

    Space group / symmetry
    ~~~~~~~~~~~~~~~~~~~~~~
    space_group : str
        Space group name or number as read from the skeleton file.
    space_group_num : int
        Numeric space group root number.
    space_group_sub_num : int
        Space group sub-number (a=1, b=2, ...).
    space_group_name : str
        Full Hermann-Mauguin space group name.
    lattice_type : str
        Bravais lattice type letter: F, I, P, A, B, C, R.
    do_full_cell : bool
        Use the conventional (full) cell rather than the primitive cell.
    do_not_cell_shift : bool
        Do not force atoms back into the central cell.
    supercell : list of int
        Supercell repeat counts along [a, b, c].
    sc_mirror : list of bool
        Flags requesting mirroring when building a supercell.

    Brillouin zone
    ~~~~~~~~~~~~~~
    b_zone_surfaces : list
        Surface normal vectors of the Brillouin zone faces.
    b_zone_vertices : list
        Vertices of the Brillouin zone.

    Replication / extension
    ~~~~~~~~~~~~~~~~~~~~~~~
    neg_bit : list of int
        Number of replicated cells needed in the -a, -b, -c directions.
    pos_bit : list of int
        Number of replicated cells needed in the +a, +b, +c directions.
    buffer : float
        Vacuum buffer (Angstroms) between the system and the simulation box
        for non-periodic systems (e.g. biomolecules).  Default: 10.

    Pair filtering (skip lists)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    skip_item_1 : list of bool
        Exclude item from first slot of a pair search.  1-indexed.
    skip_item_2 : list of bool
        Exclude item from second slot of a pair search.  1-indexed.
    skip_pair : list of list of bool
        Exclude a specific pair.  1-indexed on both axes.

    Distance and bonding
    ~~~~~~~~~~~~~~~~~~~~
    limit_dist : float
        Distance threshold for requested actions.
    limit_dist_sqrd : float
        limit_dist squared (cached for speed).
    min_dist : list of list of float
        Minimum distance between atom pairs accounting for PBC.  1-indexed.
    self_min_dist : list of float
        Minimum distance between an atom and its periodic image.  1-indexed.
    min_net_dist : list of list of float
        Minimum distance between atom pairs through a bonding network.
        1-indexed.
    num_bonds : list of int
        Number of bonds for each atom.  1-indexed.
    bonded : list of list
        Central-cell atom indices bonded to each atom.  1-indexed.
    bonded_ext : list of list
        Extended-cell atom indices bonded to each atom.  1-indexed.
    bond_length_ext : list of list
        Bond lengths to each bonded extended-cell partner.  1-indexed.
    bond_angles_ext : list of list
        Bond angle values and extended-cell atom indices.  1-indexed.
    num_bond_angles : list of int
        Number of bond angles centered on each atom.  1-indexed.
    num_bonds_total : int
        Total number of bonds in the system.
    bond_tag_id : list of int
        Unique bond-type ID for each bond.  1-indexed.
    num_unique_bond_tags : int
        Number of unique bond types (sorted alphanumerically).
    unique_bond_tags : list of str
        Tags for all unique bond types.  1-indexed.
    num_angles_total : int
        Total number of bond angles in the system.
    num_unique_angle_tags : int
        Number of unique bond-angle types.
    unique_angle_tags : list of str
        Tags for all unique angle types.  1-indexed.
    angle_tag_id : list of int
        Unique angle-type ID for each angle.  1-indexed.
    angle_bonded : list of list
        For each atom and angle, the pair of bonded atoms.  1-indexed.

    Bond orientational order
    ~~~~~~~~~~~~~~~~~~~~~~~~
    q_bar : list
        Complex multi-component orientational order parameter.
    q_order : list of float
        Real scalar bond orientational correlation value.
    ylm_coeff : list
        Coefficients to the order-6 spherical harmonics.
    ylm_l : int
        Order of q_bar.  Number of q_bar components = 2*ylm_l + 1.

    Selection box / border
    ~~~~~~~~~~~~~~~~~~~~~~
    border_zone : int
        0 = no constraint, 1 = atoms inside border, 2 = atoms outside.
    border_low : list of float
        Lower [x, y, z] (or [a, b, c]) border constraint.
    border_high : list of float
        Upper [x, y, z] (or [a, b, c]) border constraint.

    RPDF (radial pair distribution function)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rpdf_select : int
        0 = all pairs, 1 = restrict to rpdf_select_atoms pair only.
    rpdf_select_atoms : list of str
        Element names for the atom pair to use in RPDF calculation.
    rpdf : list of float
        RPDF result array.

    Rotation
    ~~~~~~~~
    rot_matrix : list of list
        3x3 rotation matrix.

    Spectral scan
    ~~~~~~~~~~~~~
    scan_element : str
        Element name for spectral scan filtering.
    num_scan_points : int
        Number of points in the spectral scan mesh.
    scan_points : list
        Scan mesh point coordinates.
    ext_scan_dist : list
        Extended-cell distances for each scan point.
    ext_pore_map : list
        Extended pore map for the scan.

    Implicit basis / potential parameters
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    These are derived from the ElementData database once the elemental
    composition of the system is known, and are set by
    compute_implicit_info().  In the Perl source they lived alongside
    explicit reference pointers into the ElementData package
    ($elementNames_ref, $covalRadii_ref, $valenceCharge_ref, etc.); in
    Python the ElementData object is accessed directly via
    self.element_data, so no reference pointers are needed.

    min_pot_alpha : float
        Smallest potential Gaussian exponent across all species in the
        system.  Initialised to 0.0; set to 1000.0 as a sentinel at the
        start of compute_implicit_info() before being reduced.
    max_num_pot_alphas : int
        Largest number of potential Gaussian terms across all species.
    max_num_atom_alphas : int
        Largest number of atom (wavefunction) Gaussian terms across all
        species.
    max_num_vale_states : int
        Largest number of valence states from any atom in the system.
    """

    def __init__(self):
        # Initialise the element database and all structure state.
        self.element_data = ElementData()
        self.element_data.init_element_data()
        self._init_state()

    def _init_state(self):
        """Set every instance attribute to its default (empty) value.

        Called by __init__() and reset().  Grouping all defaults here keeps
        __init__() short and makes reset() trivial.
        """
        data_dir = os.environ.get('OLCAO_DATA', '')

        # ------------------------------------------------------------------
        # File and path names
        # ------------------------------------------------------------------
        self.olcao_skl = 'olcao.skl'
        self.model_xyz = 'model.xyz'
        self.model_abc = 'model.abc'
        self.model_pdb = 'model.pdb'
        self.model_hin = 'model.hin'
        self.model_cif = 'model.cif'
        self.model_lmp = 'model.lmp'
        self.bond_analysis_bl = 'bondAnalysis.bl'
        self.bond_analysis_ba = 'bondAnalysis.ba'
        self.model_struct = 'model.struct'
        self.space_group_db = os.path.join(data_dir, 'spaceDB')
        self.atomic_bdb = os.path.join(data_dir, 'atomicBDB')
        self.atomic_pdb = os.path.join(data_dir, 'atomicPDB')
        self.title = ''
        self.system_title = [None]      # 1-indexed

        # ------------------------------------------------------------------
        # Atom counts and identity
        # ------------------------------------------------------------------
        self.num_atoms = 0
        self.num_atoms_ext = 0
        self.num_elements = 0
        self.num_atom_uniq_species = 0
        self.atom_uniq_species = [None]     # 1-indexed
        self.atom_element_name = [None]     # 1-indexed
        self.atom_element_id = [None]       # 1-indexed
        self.atomic_z = [None]              # 1-indexed (alias for atom_element_id)
        self.element_list = [None]          # 1-indexed
        self.element_count = [None]         # 1-indexed
        self.atom_species_id = [None]       # 1-indexed
        self.num_species = [None]           # 1-indexed
        self.species_list = [None]          # 1-indexed
        self.atom_tag = [None]              # 1-indexed
        self.connection_tag = [None]        # 1-indexed

        # ------------------------------------------------------------------
        # Atom-to-atom mapping
        # ------------------------------------------------------------------
        self.dat_skl_map = [None]           # 1-indexed
        self.skl_dat_map = [None]           # 1-indexed
        self.central_to_ext_item_map = [None]   # 1-indexed
        self.ext_to_central_item_map = [None]   # 1-indexed

        # ------------------------------------------------------------------
        # Atom coordinates
        # ------------------------------------------------------------------
        self.fract_abc = [None]             # 1-indexed; entries: [a, b, c]
        self.direct_abc = [None]            # 1-indexed; entries: [a, b, c]
        self.direct_xyz = [None]            # 1-indexed; entries: [x, y, z]
        self.ext_direct_xyz = [None]        # 1-indexed by cell
        self.ext_direct_xyz_list = [None]   # 1-indexed flat; entries: [x,y,z]
        self.ext_fract_abc_list = [None]    # 1-indexed flat; entries: [a,b,c]
        self.max_pos = [None, 0.0, 0.0, 0.0]  # 1-indexed; [None,x_max,y_max,z_max]
        self.min_pos = [None, 0.0, 0.0, 0.0]  # 1-indexed; [None,x_min,y_min,z_min]

        # ------------------------------------------------------------------
        # Per-atom elemental properties
        # ------------------------------------------------------------------
        self.coval_radii = [None]           # 1-indexed
        self.color_vtk = [None]             # 1-indexed; entries: [r,g,b,alpha]

        # ------------------------------------------------------------------
        # Molecule / residue labels (PDB / HIN)
        # ------------------------------------------------------------------
        self.molecule_name = [None]         # 1-indexed
        self.molecule_seq_num = [None]      # 1-indexed
        self.residue_name = [None]          # 1-indexed
        self.residue_seq_num = [None]       # 1-indexed

        # ------------------------------------------------------------------
        # Q^n network analysis
        # ------------------------------------------------------------------
        self.atom_qn = [None]              # 1-indexed
        # The Q^n histograms are 0-indexed (index 0 holds the Q^0 count)
        # to match Perl's @absoluteSysQn[0..8] convention.  n is
        # literally the index, so there is no sentinel slot.
        self.absolute_sys_qn = []
        self.fractional_sys_qn = []
        self.num_qn_atoms = 0
        self.net_conn_qn = 0.0

        # ------------------------------------------------------------------
        # Ring statistics
        # ------------------------------------------------------------------
        self.atom_rings = [None]           # 1-indexed
        self.rings = [None]                # 1-indexed
        self.ring_counts = [None]          # 1-indexed

        # ------------------------------------------------------------------
        # Coordination
        # ------------------------------------------------------------------
        self.coordination = [None]         # 1-indexed
        self.coordination_summary = [None] # 1-indexed

        # ------------------------------------------------------------------
        # Lattice parameters
        # ------------------------------------------------------------------
        # All lattice arrays are 1-indexed to match the Perl source and the
        # OLCAO Fortran code.  Index 0 of every list is a None placeholder.
        #
        # real_lattice[abc_axis][xyz_axis]: abc_axis in {1,2,3} for a,b,c;
        #   xyz_axis in {1,2,3} for x,y,z.
        #   real_lattice[1] = [None, ax, ay, az]  (a lattice vector)
        #   real_lattice[2] = [None, bx, by, bz]  (b lattice vector)
        #   real_lattice[3] = [None, cx, cy, cz]  (c lattice vector)
        self.real_lattice     = [[None, 0.0, 0.0, 0.0] for _ in range(4)]
        self.real_lattice_inv = [[None, 0.0, 0.0, 0.0] for _ in range(4)]
        self.recip_lattice    = [[None, 0.0, 0.0, 0.0] for _ in range(4)]
        self.real_lattice[0]     = None  # index 0 unused
        self.real_lattice_inv[0] = None
        self.recip_lattice[0]    = None

        # 1-indexed lattice magnitudes and angles.
        # mag[1]=|a|, mag[2]=|b|, mag[3]=|c|;  index 0 is None.
        self.mag          = [None, 0.0, 0.0, 0.0]  # Angstroms
        self.angle        = [None, 0.0, 0.0, 0.0]  # alpha, beta, gamma (radians)
        self.angle_deg    = [None, 0.0, 0.0, 0.0]  # alpha, beta, gamma (degrees)
        self.sine_angle   = [None, 0.0, 0.0, 0.0]  # sin(alpha,beta,gamma)
        self.full_cell_mag       = [None, 0.0, 0.0, 0.0]
        self.full_cell_angle     = [None, 0.0, 0.0, 0.0]
        self.full_cell_angle_deg = [None, 0.0, 0.0, 0.0]
        self.mag_recip       = [None, 0.0, 0.0, 0.0]
        self.angle_recip     = [None, 0.0, 0.0, 0.0]
        self.angle_deg_recip = [None, 0.0, 0.0, 0.0]
        self.real_cell_volume  = 0.0
        self.recip_cell_volume = 0.0
        self.abc_order = [None, 1, 2, 3]     # 1-indexed
        self.xyz_order = [None, 1, 2, 3]     # 1-indexed
        self.sorted_lat_indices = []

        # ------------------------------------------------------------------
        # Space group / symmetry
        # ------------------------------------------------------------------
        self.space_group = ''
        self.space_group_num = 0
        self.space_group_sub_num = 0
        self.space_group_name = ''
        self.lattice_type = 'P'
        self.do_full_cell = False
        self.do_not_cell_shift = False
        self.supercell = [None, 1, 1, 1]        # 1-indexed; repeat along a,b,c
        self.sc_mirror = [None, False, False, False]  # 1-indexed

        # ------------------------------------------------------------------
        # Brillouin zone
        # ------------------------------------------------------------------
        self.b_zone_surfaces = []
        self.b_zone_vertices = []

        # ------------------------------------------------------------------
        # Replication / extension
        # ------------------------------------------------------------------
        self.neg_bit = [None, 0, 0, 0]      # 1-indexed; cells needed in -a,-b,-c
        self.pos_bit = [None, 0, 0, 0]      # 1-indexed; cells needed in +a,+b,+c
        self.buffer = 0.0                  # Angstroms; for non-periodic systems

        # ------------------------------------------------------------------
        # Pair filtering (skip lists)
        # ------------------------------------------------------------------
        self.skip_item_1 = [None]           # 1-indexed
        self.skip_item_2 = [None]           # 1-indexed
        self.skip_pair = [None]             # 1-indexed on both axes

        # ------------------------------------------------------------------
        # Distance and bonding
        # ------------------------------------------------------------------
        self.limit_dist = 0.0
        self.limit_dist_sqrd = 0.0
        self.min_dist = [None]              # 1-indexed; entries: [dist, ...]
        self.self_min_dist = [None]         # 1-indexed
        self.min_net_dist = [None]          # 1-indexed
        self.num_bonds = [None]             # 1-indexed
        self.bonded = [None]                # 1-indexed; entries: list of ints
        self.bonded_ext = [None]            # 1-indexed; entries: list of ints
        self.bond_length_ext = [None]       # 1-indexed
        self.bond_angles_ext = [None]       # 1-indexed
        self.num_bond_angles = [None]       # 1-indexed
        self.num_bonds_total = 0
        self.bond_tag_id = [None]           # 1-indexed
        self.num_unique_bond_tags = 0
        self.unique_bond_tags = [None]      # 1-indexed
        self.num_angles_total = 0
        self.num_unique_angle_tags = 0
        self.unique_angle_tags = [None]     # 1-indexed
        self.angle_tag_id = [None]          # 1-indexed
        self.angle_bonded = [None]          # 1-indexed

        # ------------------------------------------------------------------
        # Bond orientational order
        # ------------------------------------------------------------------
        self.q_bar = []
        self.q_order = []
        self.ylm_coeff = []
        self.ylm_l = 0

        # ------------------------------------------------------------------
        # Selection box / border
        # ------------------------------------------------------------------
        # border_zone controls which atoms are included in calculations:
        #   0 = no constraint applied (all atoms)
        #   1 = only atoms within the border are included
        #   2 = only atoms outside the border are included
        # In the Perl original, $borderCoords_ref and $borderCoordsExt_ref
        # were reference pointers selecting between the central-cell and
        # extended-cell coordinate arrays at runtime.  Python accesses
        # those arrays directly, so no equivalent reference pointers are
        # needed.
        self.border_zone = 0                        # 0=none, 1=inside, 2=outside
        self.border_low  = [None, 0.0, 0.0, 0.0]   # 1-indexed lower bounds
        self.border_high = [None, 0.0, 0.0, 0.0]   # 1-indexed upper bounds

        # ------------------------------------------------------------------
        # RPDF (radial pair distribution function)
        # ------------------------------------------------------------------
        # rpdf_select controls which atom pairs are included:
        #   0 = all pairs included
        #   1 = only the pair named in rpdf_select_atoms is included
        self.rpdf_select = 0            # 0=all pairs, 1=restrict to rpdf_select_atoms
        self.rpdf_select_atoms = []     # element name pair, e.g. ['si', 'o']
        self.rpdf = []

        # ------------------------------------------------------------------
        # Rotation
        # ------------------------------------------------------------------
        # rot_matrix[row][col]: both 1-indexed (1..3).  Index 0 is None.
        self.rot_matrix = [[None, 0.0, 0.0, 0.0] for _ in range(4)]
        self.rot_matrix[0] = None  # index 0 unused

        # ------------------------------------------------------------------
        # Spectral scan
        # ------------------------------------------------------------------
        self.scan_element = ''
        self.num_scan_points = 0
        self.scan_points = []
        self.ext_scan_dist = []
        self.ext_pore_map = []

        # ------------------------------------------------------------------
        # Implicit basis / potential parameters
        # (set by compute_implicit_info() after species are known)
        # ------------------------------------------------------------------
        # In the Perl source these lived alongside explicit reference
        # pointers into the ElementData package ($elementNames_ref,
        # $covalRadii_ref, $valenceCharge_ref, etc.).  Python accesses
        # ElementData directly via self.element_data, so no reference
        # pointers are needed here.
        self.min_pot_alpha       = 0.0  # smallest potential Gaussian exponent
        self.max_num_pot_alphas  = 0    # max # of potential Gaussians (any species)
        self.max_num_atom_alphas = 0    # max # of atom Gaussians (any species)
        self.max_num_vale_states = 0    # max # of valence states from any atom

    # ======================================================================
    # Section: Reset
    # ======================================================================

    def reset(self):
        """Reset all data structures so that another file can be read in and
        manipulated.

        Every attribute set by _init_state() is returned to its default:
        structure-specific arrays (atoms, lattice, bonding, rings, …),
        default file names (olcao.skl, model.xyz, …), database paths, and
        numeric parameters (buffer, rpdf_select, …).

        The loaded ElementData object is retained (Python-specific behaviour;
        the Perl original cleared its element-database references on reset).
        Call this method between successive file reads on a single instance.
        """
        self._init_state()

    # ======================================================================
    # Section: Setters
    # These perform non-trivial initialisation beyond a simple assignment.
    # ======================================================================

    def set_num_atoms(self, n):
        """Set the number of atoms in the central (simulation) cell.

        This records only the atom count for the central cell.  It is distinct
        from ``num_atoms_ext``, which counts atoms in the periodically-extended
        neighbourhood used during bond/interaction searches.  Per-atom data
        arrays (fract_abc, direct_xyz, atom_element_id, …) are populated
        lazily by the file-reader that calls this setter, so only the count is
        stored here.

        Parameters
        ----------
        n : int
            Number of atoms in the central cell (positive integer).
        """
        self.num_atoms = n

    def set_lattice_from_mag_angle(self, mags, angles_deg):
        """Populate real_lattice from lattice magnitudes and angles.

        Parameters
        ----------
        mags : sequence of float
            1-indexed: [None, |a|, |b|, |c|] lattice vector
            magnitudes in Angstroms (mags[1] = |a|, etc.).
        angles_deg : sequence of float
            1-indexed: [None, alpha, beta, gamma] in degrees.
        """
        # Both input sequences use the 1-indexed convention
        # that is standard throughout the codebase. Angles
        # are stored in both degrees (for output) and radians
        # (for math).
        for i in range(1, 4):
            self.mag[i] = mags[i]
            self.angle_deg[i] = angles_deg[i]
            self.angle[i] = angles_deg[i] * PI / 180.0

        # Build the Cartesian lattice vectors from magnitudes and angles, then
        # derive the direct-space inverse and the reciprocal lattice matrix.
        self.get_abc_vectors()
        self.make_inv_or_recip_lattice(make_recip=False)  # real-space inverse
        self.make_inv_or_recip_lattice(make_recip=True)   # reciprocal lattice

    def set_limit_dist(self, dist):
        """Set the distance threshold and cache its square.

        Parameters
        ----------
        dist : float
            Distance threshold in Angstroms.
        """
        self.limit_dist = dist
        # Cache the squared threshold so that hot-path distance comparisons can
        # compare squared distances directly and avoid a sqrt call.
        self.limit_dist_sqrd = dist * dist

    def set_abc_xyz_assignment_order(self, abc_order, xyz_order):
        """Copy the abc/xyz axis assignment order arrays into this instance.

        Parameters
        ----------
        abc_order : sequence
            1-indexed sequence of 3 abc axis labels.
        xyz_order : sequence
            1-indexed sequence of 3 xyz axis labels.
        """
        # Both sequences are already 1-indexed; copy elements 1–3 directly.
        for axis in range(1, 4):
            self.abc_order[axis] = abc_order[axis]
            self.xyz_order[axis] = xyz_order[axis]

    def set_buffer(self, buf):
        """Set the vacuum buffer size for non-periodic (e.g. bio) systems.

        Parameters
        ----------
        buf : float
            Buffer in Angstroms between the system and the simulation box.
        """
        # Used by compute_crystal_parameters when building the enclosing
        # orthorhombic box for non-crystalline (molecule / bio) structures.
        self.buffer = buf

    def set_border(self, zone, low, high, coord_system='xyz'):
        """Define an atom-selection border region.

        Perl's ``setBorder`` received the six boundary values as
        positional arguments and stored them into ``@borderLow[1..3]``
        and ``@borderHigh[1..3]``.  The Python port groups each triple
        into a list but keeps the 1-indexed layout used everywhere else
        in this module, so axis 1 reads ``low[1]`` just like the stored
        ``self.border_low[1]`` does — no offset bookkeeping at the
        interface.

        Parameters
        ----------
        zone : int
            0 = no constraint, 1 = atoms inside border only,
            2 = atoms outside border only.
        low : sequence of float
            Lower boundary as ``[None, v1, v2, v3]`` (1-indexed,
            ``v1..v3`` match the a/b/c or x/y/z axes depending on
            ``coord_system``).
        high : sequence of float
            Upper boundary in the same 1-indexed layout as ``low``.
        coord_system : str
            'xyz' for direct Cartesian, 'abc' for direct abc,
            'frac' for fractional abc.
        """
        self.border_zone = zone
        # Both input sequences are already 1-indexed [None, v1, v2, v3];
        # copy them into the 1-indexed border storage slot-by-slot.
        self.border_low  = [None, low[1],  low[2],  low[3]]
        self.border_high = [None, high[1], high[2], high[3]]
        # Record which coordinate system the boundary values are expressed in
        # so that cut_block / apply_filter can convert atom positions correctly.
        self.border_coord_type = coord_system

    def set_rpdf_select_atoms(self, elem1, elem2):
        """Restrict RPDF calculation to a specific element pair.

        Parameters
        ----------
        elem1 : str
            Element abbreviation for the first atom in the pair.
        elem2 : str
            Element abbreviation for the second atom in the pair.
        """
        # Flag value 1 signals compute_rpdf to filter by element pair rather
        # than summing over all pairs; 0 means use all pairs.
        self.rpdf_select = 1
        # 1-indexed pair list: index 1 = first element, index 2 = second.
        self.rpdf_select_atoms = [None, elem1, elem2]

    def set_ylm_l(self, l):
        """Set the order l for bond orientational order calculations.

        Parameters
        ----------
        l : int
            Angular momentum order (typically 4 or 6).
        """
        # l is used by init_ylm_coeff and accumulate_q_bar; it must be set
        # before any bond orientational order (BOO / Steinhardt Q_l) routine.
        self.ylm_l = l

    def set_xyz_mesh_points(self, nx, ny, nz):
        """Define a Cartesian mesh for scan calculations.

        Parameters
        ----------
        nx, ny, nz : int
            Number of mesh points along each Cartesian axis (1-indexed in use).
        """
        # Store counts as a 1-indexed list [None, nx, ny, nz] so axis indices
        # 1/2/3 map directly to x/y/z without an offset adjustment.
        self.num_mesh_points = [None, nx, ny, nz]

        # Pre-allocate the pore map as a nested dict keyed by (x, y, z) in the
        # range [1..n*].  Values are initialised to 0 (unoccupied / open pore).
        self.ext_pore_map = {}
        for x in range(1, nx + 1):
            self.ext_pore_map[x] = {}
            for y in range(1, ny + 1):
                self.ext_pore_map[x][y] = {z: 0 for z in range(1, nz + 1)}

    def set_scan_points(self, do_mesh, nx=None, ny=None, nz=None,
                        num_points=None, start_xyz=None, end_xyz=None):
        """Compute and store scan point coordinates for a line or mesh scan.

        Determine if the points list is defined by a line or a mesh (do_mesh).

        Line mode (do_mesh == 0): the scan is defined by two Cartesian XYZ
        end points (start_xyz, end_xyz) and a total point count (num_points).
        Points are placed at equal intervals from start to end, inclusive.

        Mesh mode (do_mesh != 0): the scan covers the entire unit cell as a
        3-D lattice mesh of nx × ny × nz points.  The grid spans from the
        origin to the far corner of the cell along each lattice axis, with
        one step per (n-1) intervals.  Points are enumerated in a-major,
        b-middle, c-minor order (i.e. c varies fastest).

        scan_points is 1-indexed; each entry is [None, x, y, z]
        (inner coordinates are also 1-indexed).

        Parameters
        ----------
        do_mesh : int
            0 = line scan between two Cartesian end points.
            Non-zero = 3-D lattice mesh.
        nx, ny, nz : int, optional
            Number of mesh points along the a, b, c lattice axes (mesh only).
        num_points : int, optional
            Number of evenly spaced points along the scan line (line only).
        start_xyz : list of float, optional
            Cartesian start coordinates ``[None, x, y, z]``, 1-indexed
            (line only).
        end_xyz : list of float, optional
            Cartesian end coordinates ``[None, x, y, z]``, 1-indexed
            (line only).
        """
        if do_mesh == 0:
            # Obtain the variables describing the end points and the number of
            # points to use between the end points.
            self.num_scan_points = num_points

            # Compute the delta necessary in each direction to move from start
            # to end in exactly num_points steps.  Delta, start, and end are
            # all 1-indexed [None, x, y, z] so the axis loop reads the same
            # slot from each of them without any offset.
            delta = [None, 0.0, 0.0, 0.0]
            for axis in range(1, 4):
                delta[axis] = (end_xyz[axis] - start_xyz[axis]) / (num_points - 1)

            # Compute the position of each point along the path.
            # Build 1-indexed list; inner coords are also 1-indexed [None,x,y,z].
            self.scan_points = [None]
            for point in range(1, num_points + 1):
                coords = [None, 0.0, 0.0, 0.0]
                for axis in range(1, 4):
                    coords[axis] = (start_xyz[axis] +
                                    (point - 1) * delta[axis])
                self.scan_points.append(coords)

        else:
            # Obtain the variables describing the number of mesh points to use
            # along each axis.
            num_mesh = [None, nx, ny, nz]  # 1-indexed for axis arithmetic

            # Compute the total number of points in this scan.
            self.num_scan_points = nx * ny * nz

            # Compute the x,y,z delta necessary along each a,b,c axis.
            # delta[abc][xyz]: contribution of one step along abc to xyz.
            delta = [None,
                     [None, 0.0, 0.0, 0.0],
                     [None, 0.0, 0.0, 0.0],
                     [None, 0.0, 0.0, 0.0]]
            for abc in range(1, 4):
                for xyz in range(1, 4):
                    # Compute abc-axis contribution to x,y,z delta.
                    delta[abc][xyz] = (self.real_lattice[abc][xyz]
                                       / (num_mesh[abc] - 1))

            # Compute the position of each point along the path.
            # Enumerate a-major, b-middle, c-minor (c varies fastest).
            self.scan_points = [None]
            for a_pt in range(0, nx):
                for b_pt in range(0, ny):
                    for c_pt in range(0, nz):
                        coords = [None]
                        for xyz in range(1, 4):
                            val = (delta[1][xyz] * a_pt +
                                   delta[2][xyz] * b_pt +
                                   delta[3][xyz] * c_pt)
                            coords.append(val)
                        self.scan_points.append(coords)

    def set_scan_element(self, element):
        """Restrict spectral scan to atoms of the given element.

        Parameters
        ----------
        element : str
            Element abbreviation (e.g. 'si').  Empty string = all elements.
        """
        # Normalise to lowercase to match atom_element_name storage convention.
        self.scan_element = element.lower()

    # ======================================================================
    # Section: File I/O -- Read
    # ======================================================================

    def read_olcao_skl(self, filename=None, use_file_species=True):
        """Read the OLCAO skeleton (.skl) input file.

        The skeleton file is the primary OLCAO structure format.  It
        specifies the title, space group, lattice parameters, and all atom
        fractional coordinates with their species assignments.

        File format (all keywords are matched case-insensitively):

        ``title`` ... ``end``
            The title block.  The word ``title`` on a line by itself opens
            the block; the word ``end`` on a line by itself closes it.
            Between the two the user may include any free-form text (e.g.
            citations, previous calculations, pertinent notes).

        ``cell`` [``fixed``]  *or*  ``cellxyz`` [``fixed``]
            Specifies the lattice.  If ``cell`` is used, the *next* line must
            contain six real numbers: a b c alpha beta gamma (Angstroms and
            degrees).  The a-axis is then placed colinear with the Cartesian
            x-axis, the b-axis in the x-y plane, and c out of the x-y plane
            (following Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299, 2010).
            If ``cellxyz`` is used, the *next three* lines give the x,y,z
            components of the a, b, and c vectors respectively (one per line,
            Angstroms); the axes are used exactly as given.  The optional
            ``fixed`` modifier suppresses automatic reordering of the lattice
            parameters.  Without ``fixed`` the lattice is automatically
            arranged so that a ≤ b ≤ c and alpha ≤ beta ≤ gamma where allowed
            by the crystal system.

        ``frac N`` [do_not_shift]  *or*  ``cart N`` [do_not_shift]
            Introduces the atom list.  N is the number of symmetry-equivalent
            atoms that follow.  The optional third integer, if nonzero,
            inhibits shifting atoms into the central unit cell.  Each of the
            N subsequent lines has the form::

                <tag>  <coord1>  <coord2>  <coord3>  [<connection_tag>]

            The tag is an element symbol (case-insensitive) followed by an
            integer species ID, e.g. ``si1`` or ``O2``.  If no digit is
            present a ``1`` is appended automatically.  For ``frac`` the
            coordinates are fractional a,b,c; for ``cart`` they are
            Cartesian x,y,z (Angstroms).  At least eight decimal places are
            recommended to ensure correct application of symmetry operations.
            The optional ``connection_tag`` is a colon-separated element list
            (e.g. ``H:H``); if found on an O-atom line it creates a water
            molecule at that site.
            Note (Thoreau's fourth theory): if N is less than the number of
            atom lines actually present in the file, only the first N are
            read and the rest are silently ignored.

        ``space <group>``
            Space group identifier on a single line.  May be given as a
            plain integer (uses the standard default origin and unique axis)
            or as the exact ASCII name derived from the Susse (1998) uniform
            notation scheme.  The value must match a file (or symlink) in the
            space-group database directory ($OLCAO_DATA/spaceDB).  The
            database file format is: line 1 = ``<lattice_type> <name>
            [<alt_suffix>]``; line 2 = ``<sg_num> <sg_sub_num>``.  When an
            alternative origin is available the file name is extended by
            ``_<letter>`` (e.g. 227 = ``Fd3~m_a``).

        ``supercell n1 n2 n3`` [m1 m2 m3]
            Number of times to duplicate the cell along each of the a, b, c
            axes (use ``1 1 1`` for a single unit cell).  The optional three
            binary mirror flags (0 or 1) activate mirroring along the
            corresponding axis; they default to 0 if omitted.

        ``full``  *or*  ``prim``
            Request a full (non-primitive) or primitive cell.  The space
            group operations can convert full → prim, but not the reverse.

        Parameters
        ----------
        filename : str, optional
            Path to the skeleton file.  Defaults to self.olcao_skl
            ('olcao.skl' in the current directory).
        use_file_species : bool, optional
            When True (default), numeric suffixes in atom tags (e.g. 'si2')
            distinguish different species of the same element.  When False,
            all atoms of a given element are collapsed to species 1.
        """
        if filename is None:
            filename = self.olcao_skl

        with open(filename) as fh:
            skeleton = fh.readlines()

        def tok(line):
            return line.strip().split()

        # ----------------------------------------------------------------
        # Extract the title (between 'title' and 'end' keywords)
        # ----------------------------------------------------------------
        # We read a line and check 'end' *before* checking 'title' so that
        # neither keyword is captured as part of the title text itself.
        # Lines that consist solely of a semicolon are also excluded.
        getting_title = False
        found_title = False
        for raw in skeleton:
            vals = tok(raw)
            if vals and vals[0] == 'end' and len(vals) == 1:
                getting_title = False
                found_title = True
                break
            if getting_title and raw.strip() != ';':
                self.system_title.append(raw)
            if vals and vals[0] == 'title' and len(vals) == 1:
                getting_title = True

        if not found_title:
            raise ValueError(
                "Finished reading the skeleton file but did not successfully "
                "read the title.  The title must have the word 'title' on a "
                "line by itself followed by any number of lines describing the "
                "system (e.g. citations, previous calculations, other pertinent "
                "information).  The title is closed by giving the word 'end' on "
                "a line by itself.  Please adjust your input accordingly."
            )

        # ----------------------------------------------------------------
        # Extract the space group
        # ----------------------------------------------------------------
        # The value after 'space' must exactly match the name of a file (or
        # symlink) in the space-group database (usually $OLCAO_DATA/spaceDB).
        # It may be a plain integer (standard default origin) or the full
        # ASCII name from the Susse (1998) uniform notation scheme.
        #
        # Space-group database file format:
        #   Line 1: <lattice_type> <sg_name> [<alt_suffix>]
        #     lattice_type: primitive or non-primitive designation
        #     sg_name: full space group name in ASCII notation
        #     alt_suffix: lowercase letter present when an alternative origin
        #       is defined; appended as '_<letter>' to form the complete name
        #       (e.g. 227 = Fd3~m_a).
        #   Line 2: <sg_num> <sg_sub_num>
        #     Root space group number and sub-number index.
        found_space = False
        for raw in skeleton:
            vals = tok(raw)
            if vals and vals[0] == 'space' and len(vals) == 2:
                self.space_group = vals[1]
                sg_path = os.path.join(self.space_group_db, self.space_group)
                with open(sg_path) as sg_fh:
                    line1 = sg_fh.readline().split()
                    self.lattice_type = line1[0]
                    self.space_group_name = line1[1]
                    # If an alternative-origin suffix (lowercase letter) is
                    # present, append it to form the complete unique name.
                    if len(line1) > 2 and re.search(r'[a-z]', line1[2]):
                        self.space_group_name += '_' + line1[2]
                    line2 = sg_fh.readline().split()
                    self.space_group_num = int(line2[0])
                    self.space_group_sub_num = int(line2[1])
                found_space = True
                break

        if not found_space:
            raise ValueError(
                "Finished reading the skeleton file but did not successfully "
                "read the space group information.  Use 'space <group>' on a "
                "single line where <group> is one of the available options "
                "from the spaceDB directory.  Please adjust your input "
                "accordingly."
            )

        # ----------------------------------------------------------------
        # Extract the supercell
        # ----------------------------------------------------------------
        # Three integers define the number of duplicate cells along each of
        # the a, b, c axes (use 1 1 1 for a single unit cell).  An optional
        # second trio of binary values (0 or 1) activates mirroring along
        # the corresponding axis; they default to 0 (no mirroring) if
        # omitted.
        found_supercell = False
        for raw in skeleton:
            vals = tok(raw)
            if vals and vals[0] == 'supercell' and len(vals) in (4, 7):
                self.supercell[1] = int(vals[1])
                self.supercell[2] = int(vals[2])
                self.supercell[3] = int(vals[3])
                if len(vals) == 7:
                    self.sc_mirror[1] = int(vals[4]) != 0
                    self.sc_mirror[2] = int(vals[5]) != 0
                    self.sc_mirror[3] = int(vals[6]) != 0
                else:
                    self.sc_mirror[1] = False
                    self.sc_mirror[2] = False
                    self.sc_mirror[3] = False
                found_supercell = True
                break

        if not found_supercell:
            raise ValueError(
                "Finished reading the skeleton file but did not successfully "
                "read the supercell information.  Use 'supercell n1 n2 n3' on "
                "a single line where n1 n2 n3 are the number of duplicate cells "
                "along each of the a, b, c axes.  Optionally append three "
                "binary mirror flags (0 or 1).  Please adjust your input "
                "accordingly."
            )

        # ----------------------------------------------------------------
        # Extract the full/prim flag
        # ----------------------------------------------------------------
        # The lattice parameters are expected to describe the conventional
        # (full) cell of the associated space group.  'full' keeps the
        # conventional cell; 'prim' requests reduction to the primitive cell.
        # Note that this can only convert full → prim; the reverse is not
        # possible because the additional atoms needed for a full cell are
        # generated by the space group operations.
        found_cell_type = False
        for raw in skeleton:
            vals = tok(raw)
            if not vals:
                continue
            if vals[0] == 'full' and len(vals) == 1:
                self.do_full_cell = True
                found_cell_type = True
                break
            elif vals[0] == 'prim' and len(vals) == 1:
                self.do_full_cell = False
                found_cell_type = True
                break

        if not found_cell_type:
            raise ValueError(
                "Finished reading the skeleton file but did not successfully "
                "read the full/prim flag.  Use either 'full' or 'prim' on a "
                "line by itself to indicate whether to use the full "
                "(non-primitive) or primitive cell.  Please adjust your input "
                "accordingly."
            )

        # ----------------------------------------------------------------
        # Extract lattice parameters ('cell' or 'cellxyz')
        # ----------------------------------------------------------------
        # 'cell':    next line = a b c alpha beta gamma (Å and degrees).
        #   The a-axis is placed colinear with Cartesian x, the b-axis in
        #   the x-y plane, and c out of the x-y plane.  Reference:
        #   W. Setyawan & S. Curtarolo, Comp. Mat. Sci. 49, 299-312 (2010)
        #   extended via the Birkbeck crystallographic space group tables.
        # 'cellxyz': next 3 lines = x,y,z vector for a, then b, then c
        #   (Å); axes are used exactly as given.
        # In both cases the optional 'fixed' modifier (second token on the
        # 'cell'/'cellxyz' line) suppresses automatic reordering.
        #
        # Note on full_cell storage: mag/angle/angle_deg serve as both the
        # full-cell and primitive-cell parameters (the same variables are
        # overwritten when the primitive cell is created).  To preserve the
        # original values we store them separately in full_cell_mag,
        # full_cell_angle, and full_cell_angle_deg.
        mag_angle_form = False
        found_cell = False
        for i, raw in enumerate(skeleton):
            vals = tok(raw)
            if vals and vals[0] in ('cell', 'cellxyz') and len(vals) <= 2:
                if vals[0] == 'cell':
                    mag_angle_form = True
                    nv = tok(skeleton[i + 1])
                    self.mag[1] = float(nv[0])
                    self.mag[2] = float(nv[1])
                    self.mag[3] = float(nv[2])
                    self.angle_deg[1] = float(nv[3])
                    self.angle_deg[2] = float(nv[4])
                    self.angle_deg[3] = float(nv[5])
                    self.angle[1] = PI / 180.0 * self.angle_deg[1]
                    self.angle[2] = PI / 180.0 * self.angle_deg[2]
                    self.angle[3] = PI / 180.0 * self.angle_deg[3]
                    # Preserve the original full-cell parameters separately
                    # because mag/angle are later overwritten for 'prim' cells.
                    self.full_cell_mag[1]       = self.mag[1]
                    self.full_cell_mag[2]       = self.mag[2]
                    self.full_cell_mag[3]       = self.mag[3]
                    self.full_cell_angle[1]     = self.angle[1]
                    self.full_cell_angle[2]     = self.angle[2]
                    self.full_cell_angle[3]     = self.angle[3]
                    self.full_cell_angle_deg[1] = self.angle_deg[1]
                    self.full_cell_angle_deg[2] = self.angle_deg[2]
                    self.full_cell_angle_deg[3] = self.angle_deg[3]
                else:  # cellxyz: next 3 lines are the a, b, c vectors
                    for abc_axis in range(1, 4):
                        xv = tok(skeleton[i + abc_axis])
                        self.real_lattice[abc_axis][1] = float(xv[0])
                        self.real_lattice[abc_axis][2] = float(xv[1])
                        self.real_lattice[abc_axis][3] = float(xv[2])
                found_cell = True
                break

        if not found_cell:
            raise ValueError(
                "Finished reading the skeleton file but did not successfully "
                "read the cell parameter information.  Use 'cell' followed by "
                "a line with 6 real numbers (a b c alpha beta gamma) or "
                "'cellxyz' followed by three lines each giving the x,y,z "
                "components of the a, b, c vectors.  Units are Angstroms "
                "(and degrees for angles).  Please adjust your input "
                "accordingly."
            )

        # If magnitudes and angles were given, compute the Cartesian lattice
        # vectors.  Otherwise, we received the vectors directly and need to
        # derive the magnitudes and angles from them.
        if mag_angle_form:
            self.get_abc_vectors()
        else:
            self.abc_alpha_beta_gamma()

        # Compute the inverse of the real lattice now so that we can convert
        # Cartesian atom coordinates (cart) to fractional (abc) during atom
        # reading.  It will be recomputed after applying space group,
        # supercell, and any lattice reordering.
        self.make_inv_or_recip_lattice(make_recip=False)
        self.get_angle_sine()

        # ----------------------------------------------------------------
        # Extract atomic positions
        # ----------------------------------------------------------------
        # The 'frac' or 'cart' keyword introduces the atom list.  N (the
        # second token) is the number of symmetry-equivalent atoms that
        # follow.  An optional third integer, if present and nonzero,
        # suppresses shifting atoms into the central unit cell.
        #
        # Thoreau's fourth theory of adaptation: if N is less than the
        # number of atom lines actually written in the file, only the first
        # N are read and the remainder are silently ignored.  The print
        # statement below gives the user a chance to notice any discrepancy.
        #
        # Each atom line: <tag> <c1> <c2> <c3> [<connection_tag>]
        #   tag: element symbol + integer species ID (e.g. 'si1', 'O2');
        #     if no digit is present a '1' is appended automatically.
        #   c1..c3: fractional a,b,c ('frac') or Cartesian x,y,z Å ('cart').
        #     Eight or more decimal places are recommended for correct
        #     application of symmetry operations.
        #   connection_tag (optional): colon-separated element list
        #     (e.g. 'H:H'); used to add covalently bonded atoms.  For
        #     example, placing 'H:H' on an O-atom line creates water.
        self.num_atoms = 0
        found_atoms = False
        for i, raw in enumerate(skeleton):
            vals = tok(raw)
            if (vals and
                    (vals[0].startswith('frac') or vals[0].startswith('cart'))
                    and len(vals) >= 2):
                coord_type = vals[0]
                self.num_atoms = int(vals[1])
                # Optional third value: if nonzero, atoms stay outside the
                # central cell rather than being folded back in.
                self.do_not_cell_shift = (int(vals[2]) != 0) if len(vals) > 2 else False
                print(f"Reading in {self.num_atoms} atoms.")

                for atom in range(1, self.num_atoms + 1):
                    av = tok(skeleton[i + atom])
                    tag = av[0].lower()
                    # Append species ID '1' if the tag has no trailing digit.
                    if not re.search(r'[0-9]', tag):
                        tag += '1'
                    self.atom_tag.append(tag)
                    # Connection tag is the fifth field if present (index 4).
                    self.connection_tag.append(av[4].lower() if len(av) > 4 else '')

                    if coord_type.startswith('frac'):
                        fa = float(av[1])
                        fb = float(av[2])
                        fc = float(av[3])
                        self.fract_abc.append([None, fa, fb, fc])
                        # Compute direct-space a,b,c (fractional × magnitude).
                        self.direct_abc.append([None,
                                                fa * self.mag[1],
                                                fb * self.mag[2],
                                                fc * self.mag[3]])
                        self.direct_xyz.append([None, 0.0, 0.0, 0.0])
                        self.get_direct_xyz(atom)
                    else:  # cart: given in Cartesian x,y,z Å
                        dx = float(av[1])
                        dy = float(av[2])
                        dz = float(av[3])
                        self.direct_xyz.append([None, dx, dy, dz])
                        self.fract_abc.append([None, 0.0, 0.0, 0.0])
                        self.direct_abc.append([None, 0.0, 0.0, 0.0])
                        # Convert Cartesian → direct abc → fractional abc.
                        self.get_direct_abc(atom)
                        self.get_fract_abc(atom)

                found_atoms = True
                break

        if not found_atoms:
            raise ValueError(
                "Finished reading the skeleton file but did not successfully "
                "read the atomic coordinate information.  Use 'frac N' or "
                "'cart N' (optionally followed by a do-not-shift flag) on a "
                "single line, followed by N atom lines each giving a species "
                "tag and three coordinates.  Coordinates must be fractional "
                "a,b,c for 'frac' or Cartesian x,y,z (Angstroms) for 'cart'. "
                "At least eight decimal places are recommended.  Please adjust "
                "your input accordingly."
            )

        # ----------------------------------------------------------------
        # Post-read processing
        # ----------------------------------------------------------------
        # Add any atoms implied by connection tags (e.g. H:H on an O → water).
        self.add_connection_atoms()
        # Build element/species bookkeeping structures.
        self.create_element_list()
        self.create_species_data(use_file_species)
        # Expand the asymmetric unit using space group symmetry operations,
        # then tile the cell according to the supercell request.
        self.apply_space_group()
        self.apply_supercell()
        # Note: automatic reordering of lattice parameters (a≤b≤c,
        # alpha≤beta≤gamma) is currently disabled; the 'fixed' keyword
        # documents the intent but reordering is not yet implemented.

    def read_olcao_species(self, filename=None):
        """Read species assignments from inputs/olcao.fract-mi.

        This file is generated by the ``makeinput`` script and records the
        OLCAO atom positions in fractional coordinates together with the
        species tag assigned to each atom.

        File layout:
          - Lines 1-7: lattice parameters and other header records (skipped).
          - Lines 8+:  one atom per line in dat order.  Each line begins with
            a species tag such as ``si2``, where the letter part is the element
            symbol and the trailing integer is the species ID.
          - After num_atoms atom lines the index number resets, marking the
            start of the potential-position section; we stop reading there.

        The species ID is extracted from the tag (e.g. ``si2`` -> ``2``) and
        stored in ``atom_species_id`` indexed by the corresponding skeleton
        atom number, obtained via ``dat_skl_map``.

        Parameters
        ----------
        filename : str, optional
            Path to the fract-mi file. Defaults to 'inputs/olcao.fract-mi'.
        """
        if filename is None:
            filename = 'inputs/olcao.fract-mi'

        # This subroutine requires the dat-to-skl atom index map.
        self.read_dat_skl_map()

        with open(filename, 'r') as f:
            lines = f.readlines()

        # Read past the lattice parameters and other unnecessary lines.
        data_lines = [l for l in lines[7:] if l.strip()]

        import re
        for dat_atom in range(1, self.num_atoms + 1):
            tag = data_lines[dat_atom - 1].split()[0]
            # Extract the species number from the tag, e.g. 'si2' -> ['', '2'].
            parts = re.split(r'[a-zA-Z]+', tag)
            species_id = int(parts[-1])
            skl_atom = self.dat_skl_map[dat_atom]
            self.atom_species_id[skl_atom] = species_id

    def read_dat_skl_map(self, filename=None):
        """Read the dat-to-skl atom index mapping file.

        Populates two complementary, 1-indexed lookup arrays:

        - ``dat_skl_map[dat_atom]`` — given a *dat* atom number (the index
          used in the ``.dat`` / ``.fract-mi`` output files), returns the
          corresponding *skeleton* atom number from the original ``.skl``
          file.
        - ``skl_dat_map[skl_atom]`` — the reverse: given a skeleton atom
          number, returns the dat atom number.

        File format (``inputs/datSkl.map`` by default)::

            <header line — skipped>
            <dat_atom_num>  <skl_atom_num>
            ...  (one line per atom, 1-indexed)

        Parameters
        ----------
        filename : str, optional
            Path to the mapping file. Defaults to 'inputs/datSkl.map'.
        """
        if filename is None:
            filename = 'inputs/datSkl.map'

        self.dat_skl_map = [None] * (self.num_atoms + 1)
        self.skl_dat_map = [None] * (self.num_atoms + 1)

        with open(filename, 'r') as f:
            lines = f.readlines()

        # Skip header line; remaining lines are atom records.
        data_lines = [l for l in lines[1:] if l.strip()]

        for dat_atom in range(1, self.num_atoms + 1):
            values = data_lines[dat_atom - 1].split()
            # values[1] is the skeleton atom number corresponding to this dat
            # atom.  Store both directions so callers can translate freely
            # between the two numbering schemes.
            skl_atom = int(values[1])
            self.dat_skl_map[dat_atom] = skl_atom
            self.skl_dat_map[skl_atom] = dat_atom

    def read_input_file(self, filename=None, use_file_species=True):
        """Read a structure file in any supported format and populate this instance.

        Dispatches to the appropriate format-specific reader based on substring
        matching on the filename (case-insensitive), then calls
        map_element_number() and compute_implicit_info() to finish setting up
        derived data.

        The dispatch order matters: 'olcao' and 'skl' both map to the OLCAO
        skeleton reader.  'dat' maps to the legacy binary-style struct reader —
        note that this is a poor naming convention (the Perl source flags it as
        "Not ideal naming convention. FIX") because 'dat' is too generic; prefer
        more specific extensions when possible.

        Parameters
        ----------
        filename : str
            Path to the input file.  The extension (or any substring of the
            name) is used to determine the format:
              'olcao' or 'skl' → read_olcao_skl
              'pdb'            → read_pdb
              'xyz'            → read_xyz
              'dat'            → read_struct  (legacy; non-ideal extension)
              'abc'            → read_abc
              'hin'            → read_hin
              'cif'            → read_cif
              'lmp'            → read_lmp
        use_file_species : bool, optional
            If True (default), use the species tags present in the input file
            to distinguish atoms of the same element into different species.
            If False, collapse all atoms of each element into a single
            species (species 1).  The distinction matters for calculations
            where the user wants to explicitly control potential types via
            the species assignments in the skeleton file, versus cases where
            all atoms of an element should share one potential type that may
            later be refined by grouping methods (reduce, target, block).
        """
        if filename is None:
            raise ValueError('read_input_file: filename must be provided')

        # Initialise element-property arrays before any format-specific reader
        # populates atom lists (no-op in the Python port; kept for API parity).
        self.setup_data_base_ref()

        # Dispatch to the appropriate reader by substring matching on the
        # lowercased filename.  This mirrors the Perl readInputFile heuristic.
        # The use_file_species flag is forwarded to each reader so that
        # create_species_data is called with the correct setting.
        fname = filename.lower()
        if 'olcao' in fname or 'skl' in fname:
            self.read_olcao_skl(filename, use_file_species)
        elif 'pdb' in fname:
            self.read_pdb(filename, use_file_species)
        elif 'xyz' in fname:
            self.read_xyz(filename, use_file_species)
        elif 'dat' in fname:
            # Not ideal naming convention — 'dat' is too generic. (Perl: "FIX")
            self.read_struct(filename, use_file_species)
        elif 'abc' in fname:
            self.read_abc(filename, use_file_species)
        elif 'hin' in fname:
            self.read_hin(filename, use_file_species)
        elif 'cif' in fname:
            self.read_cif(filename, use_file_species)
        elif 'lmp' in fname:
            self.read_lmp(filename, use_file_species)
        else:
            raise ValueError(f'read_input_file: unknown file type for "{filename}"')

        # Post-read: assign atomic-number-based element IDs, then derive the
        # basis/potential sizing parameters from species counts.
        self.map_element_number()
        self.compute_implicit_info()

    def setup_data_base_ref(self):
        """Pull element property arrays from ElementData into this instance.

        In the Python port the ElementData object is already initialised in
        __init__, so this is a no-op kept for API compatibility with the Perl
        version.
        """
        pass

    def read_pdb(self, filename=None, use_file_species=True):
        """Read a Protein Data Bank (.pdb) structure file.

        Parameters
        ----------
        filename : str, optional
            Path to the PDB file.  Defaults to self.model_pdb.
        use_file_species : bool, optional
            If True (default), species are assigned from the uniqueness of
            PDB atom-name fields within each element.  If False, all atoms
            of a given element are collapsed to species 1.

        PDB format notes
        ----------------
        Records are identified by a 6-character name in columns 1-6 (0-indexed
        cols 0-5).  Parsing stops at the first TER record.

        CRYST1 record — unit-cell parameters (fixed-width, Angstroms/degrees):
          cols  6-14  a (9 chars)
          cols 15-23  b (9 chars)
          cols 24-32  c (9 chars)
          cols 33-39  alpha (7 chars, degrees → converted to radians here)
          cols 40-46  beta  (7 chars, degrees → converted to radians here)
          cols 47-53  gamma (7 chars, degrees → converted to radians here)
        After reading, get_abc_vectors() converts a,b,c + angles into the
        3×3 real-space lattice vector representation used throughout OLCAO.

        ATOM / HETATM records — one atom per line:
          cols 11-15  PDB atom name (5 chars; see caveat below)
          cols 17-19  residue name (3 chars)
          cols 22-25  residue sequence number (4 chars)
          cols 30-37  x coordinate (8 chars, orthogonal Angstroms)
          cols 38-45  y coordinate (8 chars, orthogonal Angstroms)
          cols 46-53  z coordinate (8 chars, orthogonal Angstroms)
          cols 76-77  element symbol (2 chars, preferred source of element name)

        Note: the PDB atom name spans cols 12-15 in the official standard, but
        some files misuse col 12 (which belongs to the preceding serial-number
        field) and bleed a digit there.  We therefore read cols 11-15 (5 chars)
        to capture these cases; well-formed files will have a leading space in
        col 11 that gets stripped harmlessly.

        Element name resolution:
          1. If cols 76-77 are present and non-empty, strip digits/spaces and
             use that as the element symbol (preferred).
          2. Otherwise, guess from the PDB atom-name field by taking the leading
             alphabetic run after stripping digits.  This works reliably only
             for single-letter element symbols.  Two-letter elements such as Ca
             cannot be distinguished from "Carbon type A" using this heuristic.

        Coordinate convention:
          PDB stores atom positions in *orthogonal* Cartesian x,y,z (Angstroms),
          not fractional a,b,c.  We therefore convert direct_xyz → direct_abc →
          fract_abc after loading all atoms.

        Species assignment strategy:
          The pdb_atom_tag (raw PDB atom name) is fundamentally different from
          the element+IDnumber atomTag used elsewhere in OLCAO — it is never a
          valid element+species pair on its own.  Rather than modifying
          create_species_data() to accept a pdb_atom_tag array as an edge case,
          we replicate the uniqueness-detection logic here to build valid
          atom_tag values (element + numeric species ID), then call
          create_species_data() normally.
        """
        if filename is None:
            filename = self.model_pdb

        explicit_abc = False
        self.num_atoms = 0

        # Local parallel array: PDB-specific atom-name tag (1-indexed).
        # This is an analog to atom_tag but holds the raw PDB atom name,
        # which is NOT the element+species-number combination used by OLCAO.
        pdb_atom_tag = [None]

        with open(filename, 'r') as fh:
            for raw in fh:
                line = raw.rstrip('\n')
                if len(line) < 6:
                    continue
                record = line[:6]

                if 'TER' in record:
                    break

                elif 'CRYST1' in record:
                    # Mark that crystal parameters are given explicitly so we
                    # skip the fallback compute_crystal_parameters() call below.
                    explicit_abc = True
                    # Extract a,b,c magnitudes and alpha,beta,gamma angles.
                    # Angles arrive in degrees; store in radians.
                    self.mag[1] = float(line[6:15])
                    self.mag[2] = float(line[15:24])
                    self.mag[3] = float(line[24:33])
                    self.angle[1] = PI / 180.0 * float(line[33:40])
                    self.angle[2] = PI / 180.0 * float(line[40:47])
                    self.angle[3] = PI / 180.0 * float(line[47:54])
                    # Convert a,b,c + angles into x,y,z vector lattice format.
                    self.get_abc_vectors()
                    self.get_angle_sine()

                elif 'ATOM' in record or 'HETATM' in record:
                    self.num_atoms += 1
                    atom = self.num_atoms

                    # Grow per-atom lists.
                    self.atom_tag.append(None)
                    self.atom_element_name.append(None)
                    self.atom_element_id.append(None)
                    self.atom_species_id.append(None)
                    self.fract_abc.append([None, 0.0, 0.0, 0.0])
                    self.direct_abc.append([None, 0.0, 0.0, 0.0])
                    self.direct_xyz.append([None, 0.0, 0.0, 0.0])
                    self.residue_name.append(None)
                    self.residue_seq_num.append(None)
                    self.molecule_seq_num.append(None)
                    pdb_atom_tag.append(None)

                    # Read the PDB atom name from cols 11-15 (5 chars).  We
                    # grab one more character than the strict standard (cols
                    # 12-15) because some files rudely use col 12 when they
                    # are not supposed to.  Strip leading/trailing whitespace
                    # and take the first token as the canonical PDB atom tag.
                    # NOTE: pdb_atom_tag is NOT an element+species-ID string.
                    temp_tag = line[11:16].lower()
                    parts = temp_tag.split()
                    pdb_atom_tag[atom] = parts[0] if parts else ''

                    # Preferred: element symbol from cols 76-77.  Strip digits
                    # and whitespace; result is the lowercase element symbol.
                    # Fallback: guess from the PDB atom-name field.  Split on
                    # digits and take the first non-empty alphabetic run.
                    # WARNING: this fallback is only reliable for single-letter
                    # elements.  For two-letter elements (e.g. Ca) the tag
                    # might be "CA" and we cannot distinguish Calcium from
                    # "Carbon type A".  We presently assume no such ambiguous
                    # atoms appear in the model.
                    if len(line) >= 77:
                        elem = re.sub(r'[0-9]', '', line[76:78]).strip().lower()
                    else:
                        elem = ''

                    if elem:
                        self.atom_element_name[atom] = elem
                    else:
                        # Split pdb_atom_tag on digits; first non-empty run
                        # of letters is the element name guess.
                        parts2 = re.split(r'[0-9]', pdb_atom_tag[atom])
                        letters = next(
                            (p.strip() for p in parts2 if p.strip()), '')
                        self.atom_element_name[atom] = letters

                    # Residue name, cols 17-19.
                    self.residue_name[atom] = line[17:20].strip()

                    # Residue sequence number, cols 22-25.
                    raw_seq = line[22:26].strip()
                    self.residue_seq_num[atom] = int(raw_seq) if raw_seq else 0

                    # Assume a single molecule is present.
                    self.molecule_seq_num[atom] = 1

                    # PDB coordinates are orthogonal Cartesian x,y,z in
                    # Angstroms (not fractional a,b,c).  Cols 30-37/38-45/46-53.
                    self.direct_xyz[atom][1] = float(line[30:38])
                    self.direct_xyz[atom][2] = float(line[38:46])
                    self.direct_xyz[atom][3] = float(line[46:54])

                    # Assume species 1 for now; refined below if use_file_species.
                    self.atom_species_id[atom] = 1
                    self.atom_tag[atom] = self.atom_element_name[atom] + '1'

        # If no CRYST1 record was found, derive a bounding-box lattice from
        # the atom positions (builds an orthorhombic P1 enclosing cell).
        if not explicit_abc:
            self.compute_crystal_parameters()

        # Convert direct_xyz → direct_abc → fract_abc for every atom.
        for atom in range(1, self.num_atoms + 1):
            self.get_direct_abc(atom)
            self.get_fract_abc(atom)

        # Build element list and element IDs (uses atom_tag).
        self.create_element_list()

        # Assign species from PDB atom-name uniqueness (per element).
        # We replicate the createSpeciesData logic here rather than adding a
        # pdb_atom_tag edge-case to that subroutine.  The result is a valid
        # atom_tag array of element+IDnumber strings so that create_species_data
        # can then be called normally.
        if use_file_species:
            pdb_species_list = {}   # {element_id: {species_id: pdb_tag}}
            num_species_pdb = {}    # {element_id: count}

            for atom in range(1, self.num_atoms + 1):
                curr_eid = self.atom_element_id[atom]
                if curr_eid not in pdb_species_list:
                    pdb_species_list[curr_eid] = {}
                    num_species_pdb[curr_eid] = 0

                found = 0
                for sid, tag in pdb_species_list[curr_eid].items():
                    if pdb_atom_tag[atom] == tag:
                        found = sid
                        break

                if found == 0:
                    num_species_pdb[curr_eid] += 1
                    new_sid = num_species_pdb[curr_eid]
                    pdb_species_list[curr_eid][new_sid] = pdb_atom_tag[atom]
                    self.atom_species_id[atom] = new_sid
                    self.atom_tag[atom] = (
                        self.atom_element_name[atom] + str(new_sid))
                else:
                    self.atom_species_id[atom] = found
                    self.atom_tag[atom] = (
                        self.atom_element_name[atom] + str(found))

        self.create_species_data(use_file_species)

    def read_xyz(self, filename=None, use_file_species=True):
        """Read an XYZ structure file.

        XYZ file layout::

            <num_atoms>
            <title/comment line>
            <Element>  <x>  <y>  <z>
            ...

        Line 1: integer atom count (first whitespace-delimited token).
        Line 2: free-form comment; the first token is stored as self.title.
        Lines 3..N+2: one atom per line — a 1- or 2-letter element symbol
        (case-insensitive) followed by three Cartesian coordinates in Angstroms.

        The XYZ format carries *no* lattice-parameter information and no species
        distinctions — explicitABC is always False.  A bounding-box orthorhombic
        P1 cell is constructed automatically via compute_crystal_parameters().
        All atoms default to species 1; create_species_data() then refines
        species assignments from the atom_tag values.

        Residue and molecule fields (residue_name, residue_seq_num,
        molecule_seq_num) are initialised to the same neutral defaults used by
        read_pdb so that downstream code that inspects those arrays works
        uniformly regardless of the source format.

        Parameters
        ----------
        filename : str, optional
            Path to the XYZ file.  Defaults to self.model_xyz.
        use_file_species : bool, optional
            If True (default), species are assigned from uniqueness of atom
            tags within each element.  If False, all atoms of a given element
            collapse to species 1.
        """
        if filename is None:
            filename = self.model_xyz

        with open(filename, 'r') as fh:
            lines = fh.readlines()

        # Line 1: total atom count.
        self.num_atoms = int(lines[0].split()[0])

        # Line 2: comment / title — take only the first whitespace token.
        self.title = lines[1].split()[0] if lines[1].split() else ''

        for atom in range(1, self.num_atoms + 1):
            # Atom lines begin at index 2 (0-based): lines[1 + atom] maps
            # atom 1 → index 2, atom 2 → index 3, etc.
            parts = lines[1 + atom].split()

            # Element symbol: assumed to be a 1- or 2-letter abbreviation;
            # lowercased for consistency with the rest of the module.
            elem = parts[0].lower()

            # XYZ carries no species information; all atoms default to species 1.
            # atom_tag is set to the bare element name here; create_species_data
            # will append the numeric species suffix later.
            self.atom_tag.append(elem)
            self.atom_element_name.append(elem)
            self.atom_element_id.append(None)
            self.atom_species_id.append(1)

            # XYZ coordinates are Cartesian (Angstroms), not fractional.
            self.direct_xyz.append([None, float(parts[1]),
                                          float(parts[2]),
                                          float(parts[3])])
            self.direct_abc.append([None, 0.0, 0.0, 0.0])
            self.fract_abc.append([None, 0.0, 0.0, 0.0])

            # PDB-compatible placeholder fields: unassigned residue, single
            # molecule.  Set to the same defaults as read_pdb so that
            # format-agnostic downstream code does not need special cases.
            self.residue_name.append('-')
            self.residue_seq_num.append(0)
            self.molecule_seq_num.append(1)

        # XYZ provides no lattice — derive an enclosing orthorhombic P1 cell
        # from the atom positions (compute_crystal_parameters handles this).
        self.compute_crystal_parameters()

        # Convert Cartesian positions to direct-abc and then fractional coords.
        for atom in range(1, self.num_atoms + 1):
            self.get_direct_abc(atom)
            self.get_fract_abc(atom)

        self.create_element_list()
        self.create_species_data(use_file_species)

    def read_struct(self, filename=None, use_file_species=True):
        """Read an OLCAO internal struct file (CELL_VECTORS / NUM_ATOM_SITES format).

        The struct format stores lattice vectors in Bohr and Cartesian atom
        coordinates in Bohr.  This reader converts both to Angstroms.

        File layout::

            CELL_VECTORS
            ax  ay  az       # a-vector components (Bohr)
            bx  by  bz
            cx  cy  cz
            NUM_ATOM_SITES
            N
            NUM_TYPE_X_Y_Z_ELEM
            num  type  x  y  z  elem   # one line per atom; x/y/z in Bohr

        Columns in the per-atom line:
          - num  : sequential atom index (1-based; ignored on read)
          - type : species ID within the element (used when use_file_species=True)
          - x y z: Cartesian coordinates in Bohr
          - elem : element symbol (case-insensitive)

        The struct format carries no title; title is left empty.

        Parameters
        ----------
        filename : str, optional
            Path to the struct file.  Defaults to self.model_struct.
        use_file_species : bool, optional
            If True (default), species IDs are taken from the type column of the
            atom lines.  If False, all atoms of a given element are collapsed
            to species 1.
        """
        if filename is None:
            filename = self.model_struct

        # The struct format carries no title; assume an empty title.
        self.title = ''

        with open(filename, 'r') as fh:
            lines = fh.readlines()

        idx = 0

        # Read past the CELL_VECTORS header.
        idx += 1

        # Read the abc lattice parameters in vector xyz format (Bohr → Angstrom).
        for axis in range(1, 4):
            vals = lines[idx].split()
            idx += 1
            self.real_lattice[axis][1] = float(vals[0]) * BOHR_RAD
            self.real_lattice[axis][2] = float(vals[1]) * BOHR_RAD
            self.real_lattice[axis][3] = float(vals[2]) * BOHR_RAD

        # Read past the NUM_ATOM_SITES header.
        idx += 1

        # Read and extract the number of atoms.
        self.num_atoms = int(lines[idx].split()[0])
        idx += 1

        # Read past the NUM_TYPE_X_Y_Z_ELEM header.
        idx += 1

        # Read all the type identification of each atom, the atomic coordinates,
        # and the elemental designation of each atom.
        for atom in range(1, self.num_atoms + 1):
            vals = lines[idx].split()
            idx += 1

            elem = vals[5].lower()

            # Store the type (species) number for this atom.
            species_id = int(vals[1]) if use_file_species else 1

            self.atom_element_name.append(elem)
            self.atom_element_id.append(None)
            self.atom_species_id.append(species_id)
            # Construct the atom tag from element name + species ID.
            self.atom_tag.append(elem + str(species_id))
            self.direct_xyz.append([None,
                                     float(vals[2]) * BOHR_RAD,
                                     float(vals[3]) * BOHR_RAD,
                                     float(vals[4]) * BOHR_RAD])
            self.direct_abc.append([None, 0.0, 0.0, 0.0])
            self.fract_abc.append([None, 0.0, 0.0, 0.0])
            # Assume an unassigned residue name.
            self.residue_name.append('-')
            # Assume an unassigned residue sequence number.
            self.residue_seq_num.append(0)
            # Assume only 1 molecule is present.
            self.molecule_seq_num.append(1)

        # Obtain the lattice vectors in a,b,c,alpha,beta,gamma form.
        self.abc_alpha_beta_gamma()

        # Obtain the inverse of the real lattice.  This must be done now so
        # that we can get the abc coordinates of atoms when given xyz.
        self.make_inv_or_recip_lattice(make_recip=False)

        # Convert direct_xyz → direct_abc → fract_abc for every atom.
        for atom in range(1, self.num_atoms + 1):
            self.get_direct_abc(atom)
            self.get_fract_abc(atom)

        # Create a set of data structures that contain information about the
        # relationship between atom numbers, elements, and species.  Also make
        # lists of the elements in the system and species for each element.
        self.create_element_list()
        self.create_species_data(use_file_species)

    def read_abc(self, filename=None, use_file_species=True):
        """Read a structure file in direct abc (fractional) coordinate format.

        The abc format is the simplest of the OLCAO readers: no title line,
        no metadata, just the bare lattice matrix followed by one atom per
        line.  It contains no explicit space-group or supercell directives
        (P1 is assumed) and no per-atom species information.

        File format
        -----------
        Line 1 (lattice matrix):
            9 whitespace-separated floats in row-major order:
                ax ay az  bx by bz  cx cy cz
            i.e. all three components of the a-vector first, then b, then c.
        Lines 2…N (one per atom):
            <tag>  <fa>  <fb>  <fc>
            where <tag> is the atom label (element symbol + optional numeric
            suffix, e.g. "si1", "o3") and fa/fb/fc are the fractional
            (crystallographic) coordinates along the a, b, c axes.

        Atom tag convention
        -------------------
        Tags follow the same element-symbol-plus-suffix convention used in
        olcao.skl.  The bare element name is recovered by stripping any
        trailing digits (e.g. "si1" → "si").  In the Perl original this
        extraction happened inline during the read loop via
            split(/[0-9]/, $atomTag[$numAtoms])
        Here we defer it to create_element_list(), which applies the same
        rule uniformly across all formats.

        All atoms default to species 1; create_species_data() may refine
        species assignments if use_file_species allows it.

        Parameters
        ----------
        filename : str, optional
            Path to the abc file.  Defaults to self.model_abc.
        use_file_species : bool, optional
            Passed through to create_species_data.
        """
        if filename is None:
            filename = self.model_abc

        with open(filename, 'r') as fh:
            lines = fh.readlines()

        # First line: 9 whitespace-separated values forming the 3×3 lattice
        # matrix, stored row-major (a_x a_y a_z  b_x b_y b_z  c_x c_y c_z).
        # The index expression maps (abc_axis, xyz_axis) → flat position:
        #   abc_axis=1,xyz_axis=1 → vals[0]  (a_x)
        #   abc_axis=1,xyz_axis=3 → vals[2]  (a_z)
        #   abc_axis=3,xyz_axis=3 → vals[8]  (c_z)
        vals = lines[0].split()
        for abc_axis in range(1, 4):
            for xyz_axis in range(1, 4):
                self.real_lattice[abc_axis][xyz_axis] = float(
                    vals[(abc_axis - 1) * 3 + (xyz_axis - 1)])

        # Derive magnitudes, angles, and sine of angles from the vectors.
        self.abc_alpha_beta_gamma()

        # Remaining lines: one atom per non-blank line (tag fa fb fc).
        # The abc format carries no species information; all atoms default to
        # species 1.  Element name is inferred later by create_element_list.
        self.num_atoms = 0
        for raw in lines[1:]:
            vals = raw.split()
            if not vals:
                continue

            self.num_atoms += 1
            atom = self.num_atoms

            tag = vals[0].lower()
            self.atom_tag.append(tag)
            self.atom_element_id.append(None)
            self.atom_species_id.append(1)

            fa = float(vals[1])
            fb = float(vals[2])
            fc = float(vals[3])
            self.fract_abc.append([None, fa, fb, fc])
            # direct_abc[axis] = fract_abc[axis] * mag[axis]  (Bohr along each
            # lattice direction, not Cartesian — see get_direct_xyz for xyz).
            self.direct_abc.append([None,
                                    fa * self.mag[1],
                                    fb * self.mag[2],
                                    fc * self.mag[3]])
            self.direct_xyz.append([None, 0.0, 0.0, 0.0])
            self.atom_element_name.append(None)
            self.get_direct_xyz(atom)

        self.create_element_list()
        self.create_species_data(use_file_species)

    def read_cif(self, filename=None, use_file_species=True):
        """Read a Crystallographic Information File (.cif).

        .. note::
            At present this will only read CIF files that are symmetry P1 and
            in the format created by the OLCAO ``makeinput`` program.  Perhaps
            in the future this can be made more general.

        The expected file layout is::

            data_<title>
            <5 header lines>
            _cell_length_a   <value>    # lattice magnitude a, in Angstroms
            _cell_length_b   <value>    # lattice magnitude b, in Angstroms
            _cell_length_c   <value>    # lattice magnitude c, in Angstroms
            _cell_angle_alpha  <value>  # angle alpha, in degrees
            _cell_angle_beta   <value>  # angle beta,  in degrees
            _cell_angle_gamma  <value>  # angle gamma, in degrees
            <6 loop-header lines>
            <one atom per line: label  type_symbol  fract_x  fract_y  fract_z>
            ...

        The ``type_symbol`` column (column 1) is present in the CIF atom loop
        but is not used here; atom identity is derived entirely from the
        ``label`` column (column 0) instead.

        Atom label convention
        ---------------------
        CIF atom labels produced by ``makeinput`` follow the pattern
        ``<element><digit>_<rest>``, e.g. ``Si1_3``.  The part before the
        first ``_`` is extracted and lowercased to form the OLCAO atom tag
        (e.g. ``si1``).  ``create_species_data`` later separates the element
        symbol from any trailing species digit in that tag.

        Parameters
        ----------
        filename : str, optional
            Path to the CIF file.  Defaults to self.model_cif.
        use_file_species : bool, optional
            If True (default), species are taken from atom label prefixes
            (the part before the first ``_``).  If False, all atoms of a
            given element are collapsed to species 1.
        """
        if filename is None:
            filename = self.model_cif

        # CIF files from makeinput always use P1 symmetry.
        self.space_group = '1_a'
        self.num_atoms = 0

        with open(filename, 'r') as fh:
            # Title line: "data_<title>"
            line = fh.readline()
            parts = line.strip().split('data_', 1)
            self.title = parts[1].strip() if len(parts) > 1 else ''

            # Skip 5 header lines.
            for _ in range(5):
                fh.readline()

            # Cell parameters: "_cell_length_a VALUE" etc. (6 lines).
            self.mag[1] = float(fh.readline().split()[1])
            self.mag[2] = float(fh.readline().split()[1])
            self.mag[3] = float(fh.readline().split()[1])
            self.angle[1] = PI / 180.0 * float(fh.readline().split()[1])
            self.angle[2] = PI / 180.0 * float(fh.readline().split()[1])
            self.angle[3] = PI / 180.0 * float(fh.readline().split()[1])

            # Build lattice vectors and populate sine_angle.
            self.get_abc_vectors()
            self.get_angle_sine()

            # Skip 6 atom-loop header lines.
            for _ in range(6):
                fh.readline()

            # Atom data: one atom per line until EOF.
            # Columns: label  type_symbol  fract_x  fract_y  fract_z
            for raw in fh:
                values = raw.strip().split()
                if not values:
                    continue

                self.num_atoms += 1
                atom = self.num_atoms

                # Label like "Si1_3": take part before first "_" as the
                # OLCAO atom tag (element symbol + optional species digit).
                # values[1] is CIF type_symbol (element symbol); we skip it
                # and rely on the label alone for species information.
                tag = values[0].split('_')[0].lower()
                fa = float(values[2])
                fb = float(values[3])
                fc = float(values[4])

                self.atom_tag.append(tag)
                self.atom_element_name.append(None)
                self.atom_element_id.append(None)
                self.atom_species_id.append(None)
                self.fract_abc.append([None, fa, fb, fc])
                # direct_abc = fract * mag: position along each lattice axis
                # in Angstroms (not Cartesian XYZ).  get_direct_xyz converts
                # to Cartesian coordinates below.
                self.direct_abc.append([None,
                                        fa * self.mag[1],
                                        fb * self.mag[2],
                                        fc * self.mag[3]])
                self.direct_xyz.append([None, 0.0, 0.0, 0.0])
                self.residue_name.append('-')
                self.residue_seq_num.append(0)
                self.molecule_seq_num.append(1)

                self.get_direct_xyz(atom)

        self.create_element_list()
        self.create_species_data(use_file_species)

    def read_hin(self, filename=None, use_file_species=True):
        """Read a HyperChem HIN structure file.

        HIN is the native format of the HyperChem molecular modeling package.
        The file is keyword-driven and organises atoms into a hierarchy of
        molecules and residues.  There is no concept of a periodic simulation
        cell; one is derived after reading using compute_crystal_parameters.

        File structure
        --------------
        Records are identified by the first whitespace-delimited token
        (case-insensitive).  The tokens of interest are:

        ``mol <n> <name>``
            Opens a new molecule block.  Increments the molecule counter and
            resets the residue counter to zero.
        ``endmol <n>``
            Closes the current molecule block (contains "mol" *and* "end", so
            it is distinguished by the "end" guard).
        ``res <n> <resname> <...>``
            Opens a new residue within the current molecule.  Increments the
            residue counter and records the residue name from field [2].
        ``endres <n>``
            Closes the current residue block (distinguished the same way as
            endmol).
        ``atom <serial> <name> <element> <type> <formal_chg> <partial_chg> <x> <y> <z> <num_bonds> [<bonded_atom> <bond_order> ...]``
            Describes one atom.  The fixed layout is:
              [0]  "atom"
              [1]  atom serial number within the molecule
              [2]  atom name (e.g. "C1", "N3") — stored as atom_tag (lowercased)
              [3]  element symbol (e.g. "C", "N", "Lp") — stored as atom_element_name
              [4]  Sybyl/HyperChem atom type (e.g. "C.ar", "N.am")
              [5]  formal charge
              [6]  partial charge
              [7]  x coordinate (Angstroms, Cartesian)
              [8]  y coordinate (Angstroms, Cartesian)
              [9]  z coordinate (Angstroms, Cartesian)
              [10] number of bonds (followed by bond-pair fields)
            Fields [7..9] are always at these fixed offsets regardless of the
            bond list that follows.

        Lone pairs
        ----------
        HyperChem uses the token "Lp" in field [3] to denote lone-pair
        pseudo-atoms.  These are not real atoms but occupy a numbered slot in
        the molecule's atom list.  They are silently skipped (num_atoms is
        not incremented) so that the resulting index sequence matches only
        genuine atoms.

        Cell derivation
        ---------------
        Because HIN files carry no periodic boundary information,
        compute_crystal_parameters is called after parsing to construct a
        minimal orthogonal P1 box that encloses all atom positions with a
        small buffer.  Coordinates are then converted from direct_xyz to
        direct_abc and fract_abc.

        Parameters
        ----------
        filename : str, optional
            Path to the HIN file.  Defaults to self.model_hin.
        use_file_species : bool, optional
            If True (default), species are assigned from uniqueness of atom
            tags within each element.  If False, all atoms of a given element
            collapse to species 1.
        """
        if filename is None:
            filename = self.model_hin

        self.num_atoms = 0
        num_molecules = 0
        num_residues = 0
        current_res_name = ''

        with open(filename, 'r') as fh:
            for raw in fh:
                values = raw.split()
                if not values:
                    continue

                line_id = values[0].lower()

                # "mol" opens a new molecule; "endmol" contains both substrings
                # so the "end" guard distinguishes them.
                if 'mol' in line_id and 'end' not in line_id:
                    num_molecules += 1
                    num_residues = 0

                # "res" opens a new residue; same end-guard logic as mol.
                elif 'res' in line_id and 'end' not in line_id:
                    num_residues += 1
                    current_res_name = values[2].lower() if len(values) > 2 else ''

                elif line_id == 'atom':
                    # Some atoms are actually lone pairs ("Lp" in field [3]).
                    # They are not real atoms but do occupy a numbered slot in
                    # the molecule.  Skip them without incrementing num_atoms.
                    if len(values) > 3 and values[3].lower() == 'lp':
                        continue

                    self.num_atoms += 1

                    # Grow per-atom lists (1-indexed; index 0 stays None).
                    self.atom_tag.append(values[2].lower())
                    self.atom_element_name.append(values[3].lower())
                    self.atom_element_id.append(None)
                    self.atom_species_id.append(1)
                    self.fract_abc.append([None, 0.0, 0.0, 0.0])
                    self.direct_abc.append([None, 0.0, 0.0, 0.0])
                    # Cartesian coordinates are always at fields [7], [8], [9]
                    # regardless of how many bond pairs follow in the record.
                    self.direct_xyz.append([None,
                                            float(values[7]),
                                            float(values[8]),
                                            float(values[9])])
                    self.residue_name.append(current_res_name)
                    self.residue_seq_num.append(num_residues)
                    self.molecule_seq_num.append(num_molecules)

        # HIN files carry no periodic cell; derive one from atom positions.
        self.compute_crystal_parameters()

        # Convert direct_xyz → direct_abc → fract_abc for every atom.
        for atom in range(1, self.num_atoms + 1):
            self.get_direct_abc(atom)
            self.get_fract_abc(atom)

        self.create_element_list()
        self.create_species_data(use_file_species)

    def read_lmp(self, filename=None, use_file_species=True):
        """Read a LAMMPS data file (as written by print_lmp / printLMP).

        This reads the LAMMPS *data* file format written by OLCAO's
        print_lmp / printLMP.  Note: the original Perl readLMP comment said
        "read a lammps dump file while also pulling information from an
        associated lammps data file", but that sub was never completed (it
        only read the title line).  The actual OLCAO–LAMMPS interchange
        format is the data file written by printLMP, and that is what this
        method reads.

        **File structure (as produced by print_lmp)**::

            <filename>                          ← title (first non-blank line)
                                                ← blank line
            N atoms
            M atom types
            0.000 lx xlo xhi                    ← box bounds; lo is always 0
            0.000 ly ylo yhi
            0.000 lz zlo zhi
            [xy xz yz xy xz yz]                 ← only for non-orthogonal cells
                                                ← blank before Masses keyword
            Masses

            <type_id> <mass> # <elem> <sid>     ← one line per unique (elem, species)

            Atoms

            <id> <type_id> <x> <y> <z> # <elem> <sid>   ← Cartesian angstroms

        **Box-bounds convention**: printLMP writes the lattice *magnitudes*
        (mag[1], mag[2], mag[3]) for lx/ly/lz rather than the true LAMMPS
        reduced-cell dimensions (a_x, b_y, c_z).  read_lmp reads them back
        the same way, so the OLCAO round-trip is self-consistent, but the
        resulting file uses a non-standard convention for triclinic cells and
        may not load correctly into LAMMPS itself.

        **Tilt factors** (written only when any angle ≠ 90°; see the LAMMPS
        triclinic documentation)::

            xy = |b| * cos(gamma)
            xz = |c| * cos(beta)
            yz = (|b|·|c|·cos(alpha) − xy·xz) / sqrt(|b|² − xy²)

        **Lattice recovery**: real_lattice is rebuilt from the six box
        parameters using the standard LAMMPS triclinic box-vector convention::

            a = (lx,  0,  0)
            b = (xy, ly,  0)
            c = (xz, yz, lz)

        **Masses section**: one entry per unique OLCAO (element, species)
        pair, in atom-traversal order (first occurrence wins, same ordering as
        printLMP).  The ``# elem sid`` trailing comment carries the element
        name and integer species ID used to populate the internal type_map
        when use_file_species is True.  If the comment is absent (e.g., a
        file not produced by OLCAO), elem is set to 'x' and species to 1.

        **Atoms section**: each atom's Cartesian (x, y, z) in angstroms.
        After all atoms are read, get_direct_abc and get_fract_abc are called
        for every atom to fill in the fractional and direct-ABC coordinates.

        **Perl bug note**: printLMP's Atoms section used ``$currentSpecies``
        (left over from the Masses loop) for every atom instead of the
        per-atom species ID, so the trailing comment was wrong for all atoms
        after the last Masses entry.  print_lmp and read_lmp both use the
        correct per-atom species ID.

        Parameters
        ----------
        filename : str, optional
            Path to the LAMMPS data file.  Defaults to self.model_lmp.
        use_file_species : bool, optional
            If True (default), species tags are read from the ``# elem
            species`` comment in the Masses section.  If False, all atoms of
            a given element collapse to species 1.
        """
        if filename is None:
            filename = self.model_lmp

        self.space_group = '1_a'
        self.num_atoms = 0

        # type_id (1-indexed int) → (element_name str, species_id int)
        type_map = {}

        # Box geometry from header
        lx = ly = lz = 0.0
        xy = xz = yz = 0.0

        with open(filename, 'r') as fh:
            lines = fh.readlines()

        # ------------------------------------------------------------------
        # Pass 1: title (first non-blank line)
        # ------------------------------------------------------------------
        body_start = 0
        for i, raw in enumerate(lines):
            if raw.strip():
                self.title = raw.strip()
                body_start = i + 1
                break

        # ------------------------------------------------------------------
        # Pass 2: header keywords (atoms, atom types, box bounds, tilt)
        # ------------------------------------------------------------------
        for raw in lines[body_start:]:
            values = raw.split()
            if not values:
                continue
            joined = raw.lower()
            if 'atoms' in joined and 'atom types' not in joined and len(values) >= 2:
                # "N atoms"
                try:
                    self.num_atoms = int(values[0])
                except ValueError:
                    pass
            elif 'xlo' in joined and 'xhi' in joined:
                lx = float(values[1]) - float(values[0])
            elif 'ylo' in joined and 'yhi' in joined:
                ly = float(values[1]) - float(values[0])
            elif 'zlo' in joined and 'zhi' in joined:
                lz = float(values[1]) - float(values[0])
            elif ('xy' in joined and 'xz' in joined and 'yz' in joined
                  and len(values) >= 3):
                # "xy xz yz xy xz yz"  (values after tilt keyword)
                xy = float(values[0])
                xz = float(values[1])
                yz = float(values[2])
            elif values[0].lower() == 'masses':
                break

        # ------------------------------------------------------------------
        # Build real_lattice from LAMMPS box vectors:
        #   a = (lx,  0,  0)
        #   b = (xy, ly,  0)
        #   c = (xz, yz, lz)
        # This is identical to the triclinic P1 convention used in OLCAO.
        # ------------------------------------------------------------------
        self.real_lattice[1] = [None, lx,  0.0, 0.0]
        self.real_lattice[2] = [None, xy,  ly,  0.0]
        self.real_lattice[3] = [None, xz,  yz,  lz ]
        self.abc_alpha_beta_gamma()
        self.make_inv_or_recip_lattice(make_recip=False)
        self.make_inv_or_recip_lattice(make_recip=True)

        # ------------------------------------------------------------------
        # Pass 3: Masses and Atoms sections
        # ------------------------------------------------------------------
        section = None
        atom_counter = 0
        for raw in lines:
            values = raw.split()
            if not values:
                continue
            keyword = values[0].lower()

            if keyword == 'masses':
                section = 'masses'
                continue
            if keyword == 'atoms':
                section = 'atoms'
                continue

            if section == 'masses':
                # <type_id> <mass> # <elem_name> <species_id>
                try:
                    tid = int(values[0])
                except ValueError:
                    continue
                if '#' in raw:
                    comment = raw.split('#', 1)[1].split()
                    elem_name  = comment[0].lower() if comment else 'x'
                    species_id = int(comment[1]) if len(comment) > 1 else 1
                else:
                    elem_name  = 'x'
                    species_id = 1
                type_map[tid] = (elem_name, species_id)

            elif section == 'atoms':
                # <atom_id> <type_id> <x> <y> <z> [# ...]
                if len(values) < 5:
                    continue
                try:
                    tid = int(values[1])
                    x   = float(values[2])
                    y   = float(values[3])
                    z   = float(values[4])
                except ValueError:
                    continue

                elem_name, species_id = type_map.get(tid, ('x', 1))
                tag = f'{elem_name}{species_id}'

                atom_counter += 1
                atom = atom_counter

                self.atom_tag.append(tag)
                self.atom_element_name.append(None)
                self.atom_element_id.append(None)
                self.atom_species_id.append(None)
                self.direct_xyz.append([None, x, y, z])
                self.fract_abc.append([None, 0.0, 0.0, 0.0])
                self.direct_abc.append([None, 0.0, 0.0, 0.0])
                self.residue_name.append('-')
                self.residue_seq_num.append(0)
                self.molecule_seq_num.append(1)

                self.get_direct_abc(atom)
                self.get_fract_abc(atom)

        # Honour the count from the header (atoms may be written out of order
        # in extended formats; atom_counter is the authoritative tally here).
        self.num_atoms = atom_counter

        self.create_element_list()
        self.create_species_data(use_file_species)

    def read_bond_analysis_bl(self, filename=None):
        """Read a bond-length analysis file (bondAnalysis.bl).

        This method does not assume that any other structural information has
        been previously loaded (e.g. no associated olcao.skl file is required).
        It reads BOND LENGTHS only; bond angles are in a separate .ba file.

        File format
        -----------
        The file contains one block per atom.  Each block begins with a
        single header line followed by one or more bond-data rows:

        Header line (one per atom):
            <elemTag>_<atomNum>  <something>  <numBonds>  ...
            e.g.  "fe_42  fe  6  ..."
            - elemTag   : element symbol (e.g. 'b', 'c', 'fe') — same tag
                          style used elsewhere in the bond-analysis output
            - atomNum   : 1-based global atom index
            - numBonds  : total number of bonds for this atom (field index 2)

        Bond-data rows (ceil(numBonds / 4) rows follow the header):
            Up to 4 bond pairs per row, each pair being:
                <elemTag>_<bondAtomNum>  <bondLength>
            e.g.  "o_3  1.625  o_7  1.618  o_12  1.631  o_19  1.640"
            - bondAtomNum : 1-based global index of the bonded atom
            - bondLength  : interatomic distance in Angstroms

        Double-listing: every bond A→B also appears as B→A, so the file
        contains each unique bond twice.  num_bonds_total is halved at the
        end to give the true count of unique bonds.

        Populates
        ---------
        num_bonds[atom]              : number of bonds for atom (1-indexed)
        bonded[atom][bondIndex]      : global index of bonded partner
        bond_length_ext[atom][bondIndex] : bond length in Angstroms
        bond_tag_id[atom][bondIndex] : index into unique_bond_tags (1-indexed)
        unique_bond_tags[tagID]      : canonical "elemA elemB" pair string
        num_unique_bond_tags         : count of distinct bond-type pairs
        num_bonds_total              : unique bond count (double-listing halved)

        Parameters
        ----------
        filename : str, optional
            Path to the bond-length file.  Defaults to self.bond_analysis_bl.
        """
        if filename is None:
            filename = self.bond_analysis_bl

        # Assume that all atoms have zero bonds initially.
        self.num_bonds           = [None] + [0     for _ in range(self.num_atoms)]
        self.bonded              = [None] + [[None] for _ in range(self.num_atoms)]
        self.bond_length_ext     = [None] + [[None] for _ in range(self.num_atoms)]
        self.bond_tag_id         = [None] + [[None] for _ in range(self.num_atoms)]
        # Initialise a count of the total number of bonds in the system.
        self.num_bonds_total     = 0
        # Initialise a count of the number of unique bond-type tags.
        self.num_unique_bond_tags = 0
        self.unique_bond_tags    = [None]  # 1-indexed; [0] is None placeholder

        with open(filename) as f:
            while True:
                line = f.readline()
                if not line:
                    break
                values = line.split()
                if not values:
                    continue

                # Get the atom element tag, atom number, and bond count from
                # the header line.  Field 0 is "elemTag_atomNumber";
                # field 2 is the number of bonds for this atom.
                parts       = values[0].split('_')
                atom_tag    = parts[0]   # e.g. 'b', 'c', 'fe'
                atom_number = int(parts[1])
                num_bonds_atom = int(values[2])
                self.num_bonds[atom_number] = num_bonds_atom

                # Bonds are written 4 per row; compute how many rows follow.
                num_bond_rows = math.ceil(num_bonds_atom / 4.0)
                for _ in range(num_bond_rows):
                    # Compute the number of bonds in this row (may be 4 or
                    # fewer for the last row).
                    row_values       = f.readline().split()
                    num_bonds_in_row = len(row_values) // 2
                    for bond in range(num_bonds_in_row):
                        self.num_bonds_total += 1

                        # Read the bonded atom's element tag and atom number.
                        bond_parts    = row_values[2 * bond].split('_')
                        bond_atom_tag = bond_parts[0]
                        bond_atom_num = int(bond_parts[1])

                        self.bonded[atom_number].append(bond_atom_num)
                        self.bond_length_ext[atom_number].append(
                            float(row_values[2 * bond + 1]))

                        # Build a canonical bond-type tag by placing the
                        # lexicographically smaller element tag first, so each
                        # bond pair (A-B and B-A) maps to the same unique key.
                        if atom_tag > bond_atom_tag:
                            bond_tag = f'{bond_atom_tag} {atom_tag}'
                        else:
                            bond_tag = f'{atom_tag} {bond_atom_tag}'

                        found = 0
                        for t in range(1, self.num_unique_bond_tags + 1):
                            if self.unique_bond_tags[t] == bond_tag:
                                found = t
                                break
                        # If this tag wasn't already known, register it.
                        if found == 0:
                            self.num_unique_bond_tags += 1
                            found = self.num_unique_bond_tags
                            self.unique_bond_tags.append(bond_tag)

                        # Store the unique bond-type ID for this bond.
                        self.bond_tag_id[atom_number].append(found)

        # Cut the bond count in half because all bonds are double-listed in
        # the bondAnalysis.bl file (each A→B bond also appears as B→A).
        self.num_bonds_total //= 2

    def read_bond_analysis_ba(self, filename=None):
        """Read a bond-angle analysis file (bondAnalysis.ba).

        Populates num_bond_angles, bond_angles_ext, angle_tag_id,
        angle_bonded, unique_angle_tags, num_unique_angle_tags, and
        num_angles_total.

        No prior system knowledge is assumed: this method can be called
        independently of any olcao.skl file.  self.num_atoms must already
        be set (e.g. from a prior read_bond_analysis_bl call).

        File format (written by the bondAnalysis script with -ba):
          Header line per atom:
            "<atom> Num bond angles:  <count>"
          Detail line per bond angle:
            "<tag1> <tag2> <tag3> <bonded1> <vertex> <bonded2> <angle_deg>"
          where tag2/vertex is the central (vertex) atom and tag1/tag3 are the
          two outer atoms.  Angles are not double-listed (unlike bonds in .bl).

        Parameters
        ----------
        filename : str, optional
            Path to the bond-angle file.  Defaults to self.bond_analysis_ba.
        """
        if filename is None:
            filename = self.bond_analysis_ba

        # Assume that all atoms have zero bond angles.
        self.num_bond_angles       = [None] + [0     for _ in range(self.num_atoms)]
        self.bond_angles_ext       = [None] + [[None] for _ in range(self.num_atoms)]
        self.angle_tag_id          = [None] + [[None] for _ in range(self.num_atoms)]
        self.angle_bonded          = [None] + [[None] for _ in range(self.num_atoms)]

        # Initialize a count of the total number of bond angles in the system.
        self.num_angles_total      = 0

        # Initialize a count of the total number of unique bond angle
        # configurations defined by the combination of element-species tags
        # and bond angles.
        self.num_unique_angle_tags = 0
        self.unique_angle_tags     = [None]

        with open(filename) as f:
            # Read each set of bond angle information in the given bondAnalysis file.
            while True:
                line = f.readline()
                if not line:
                    break
                values = line.split()
                if not values:
                    continue

                # Get the number of bond angles where this atom serves as a
                # vertex.  Also increment the total number of bond angles.
                # Header line: "<atomNumber> Num bond angles:  <count>"
                atom_number     = int(values[0])
                num_angles_atom = int(values[4])
                self.num_bond_angles[atom_number] = num_angles_atom
                self.num_angles_total += num_angles_atom

                # Read information about each bond angle.
                for _ in range(num_angles_atom):
                    # Detail line: tag1 tag2 tag3 bonded1 vertex bonded2 angle
                    vals    = f.readline().split()
                    tag1    = vals[0]
                    tag2    = vals[1]
                    tag3    = vals[2]
                    bonded1 = int(vals[3])
                    atom    = int(vals[4])   # vertex atom
                    bonded2 = int(vals[5])
                    raw     = float(vals[-1])

                    # Round the bond angle to the nearest integer or
                    # half-integer.  Breakpoints (fractional part):
                    #   [0.00, 0.25) → floor
                    #   [0.25, 0.75) → floor + 0.5
                    #   [0.75, 1.00) → floor + 1
                    decimal = raw - int(raw)
                    if decimal > 0.75:
                        current_angle = float(int(raw) + 1)
                    elif decimal > 0.25:
                        current_angle = int(raw) + 0.5
                    else:
                        current_angle = float(int(raw))

                    self.bond_angles_ext[atom].append(current_angle)

                    # Classify the bond angle type as a unique element-species
                    # combination.  The 1st and 3rd tags can appear in any order
                    # (physically the angle is symmetric), so we sort them
                    # lexicographically to produce a single canonical key.
                    # The 2nd tag must remain in the middle: it identifies the
                    # vertex atom of the angle.
                    if tag1 <= tag3:
                        angle_tag = f'{tag1} {tag2} {tag3} {current_angle}'
                    else:
                        angle_tag = f'{tag3} {tag2} {tag1} {current_angle}'

                    # Search the list of known tags for the one we just made.
                    # When unique_angle_tags is empty, num_unique_angle_tags is
                    # 0 and the loop body is never entered.
                    found = 0
                    for t in range(1, self.num_unique_angle_tags + 1):
                        if self.unique_angle_tags[t] == angle_tag:
                            found = t
                            break

                    # If the tag was not found, increase the number of known
                    # tags and store the new tag.
                    if found == 0:
                        self.num_unique_angle_tags += 1
                        found = self.num_unique_angle_tags
                        self.unique_angle_tags.append(angle_tag)

                    # For the current angle, store the associated unique tag number.
                    self.angle_tag_id[atom].append(found)

                    # Store the two atoms associated with this bond angle vertex.
                    # [None, bonded1, bonded2] gives 1-indexed access matching
                    # Perl's $angleBonded[$a][$i][1/2].
                    self.angle_bonded[atom].append([None, bonded1, bonded2])

    # ======================================================================
    # Section: File I/O -- Write / Print
    # ======================================================================

    def print_olcao(self, filename=None, title='', style='frac'):
        """Write the current structure in OLCAO skeleton (.skl) format.

        File format
        -----------
        The .skl file has four sections:

        1. Title block::

               title
               <title arg, if non-empty>
               <system_title lines, one per line>
               end

           ``system_title`` holds the raw title lines accumulated during input
           parsing (e.g. from a previous .skl read).  An extra caller-supplied
           ``title`` line is prepended when given.

        2. Cell block::

               cell
               <a> <b> <c> <alpha> <beta> <gamma>

           Six values on one line: the three lattice magnitudes (Angstroms) then
           the three angles (degrees).  ``self.angle_deg`` holds the pre-converted
           values; the Perl original converted from radians inline.

        3. Atom count / coordinate style line::

               <style> <num_atoms> 1

           The trailing ``1`` is a flag meaning "force atoms into the unit cell"
           (wrap fractional coordinates back into [0, 1)).

           ``style`` is the literal keyword written to the file.  The Perl
           original used ``"fractional"`` as the canonical token; this Python
           port defaults to ``'frac'`` for brevity but accepts any string
           (including ``'fractional'``) since the reader matches on the first
           few characters.

        4. Atom list::

               <tag>  <coord1>  <coord2>  <coord3>

           ``tag`` = element_name + species_id (e.g. ``si2``).
           Coordinates are fractional a,b,c when ``style`` starts with
           ``'frac'``, otherwise Cartesian x,y,z (Angstroms).

        5. Footer (always written, regardless of input source)::

               space 1_a
               supercell 1 1 1
               full

           ``space 1_a`` declares P1 symmetry (space group 1, setting a) — no
           symmetry operations will be applied on re-read.
           ``supercell 1 1 1`` declares no supercell replication.
           ``full`` signals that the structure is fully expanded and does not
           need further generation from a primitive cell.

        Parameters
        ----------
        filename : str, optional
            Output path.  Defaults to self.olcao_skl.
        title : str, optional
            Extra title line inserted after 'title' and before the stored
            system_title block.  Omitted when empty (default).
        style : str, optional
            Coordinate style keyword written to the file.  ``'frac'`` (default)
            writes fractional a,b,c; anything else (e.g. ``'cart'``) writes
            Cartesian x,y,z.
        """
        if filename is None:
            filename = self.olcao_skl

        # Build per-atom tag: element_name + species_id  (e.g. 'si2')
        atoms = [None]
        for atom in range(1, self.num_atoms + 1):
            atoms.append(self.atom_element_name[atom]
                         + str(self.atom_species_id[atom]))

        with open(filename, 'w') as fh:
            fh.write('title\n')
            if title:
                fh.write(title + '\n')
            # system_title entries are raw lines that already carry their '\n'
            for raw in self.system_title[1:]:
                fh.write(raw)
            fh.write('end\n')

            fh.write('cell\n')
            for axis in range(1, 4):
                fh.write(f'{self.mag[axis]:13.8f}')
            for axis in range(1, 4):
                # angle_deg holds pre-converted degree values (Perl converted
                # from self.angle radians inline: angle * 180/pi)
                fh.write(f'{self.angle_deg[axis]:13.8f}')
            fh.write('\n')

            # Trailing '1' = force atoms into cell on re-read
            fh.write(f'{style} {self.num_atoms} 1\n')

            if style.startswith('frac'):
                for atom in range(1, self.num_atoms + 1):
                    fa = self.fract_abc[atom]
                    fh.write(f'{atoms[atom]:>6s} {fa[1]:12.8f} {fa[2]:12.8f} {fa[3]:12.8f}\n')
            else:
                for atom in range(1, self.num_atoms + 1):
                    xyz = self.direct_xyz[atom]
                    fh.write(f'{atoms[atom]:>6s} {xyz[1]:12.8f} {xyz[2]:12.8f} {xyz[3]:12.8f}\n')

            # Always write P1 / no-supercell / full footer regardless of how
            # the structure was originally read (space group, primitive cell, etc.)
            fh.write('space 1_a\n')
            fh.write('supercell 1 1 1\n')
            fh.write('full\n')

    def print_olcao_map(self, filename=None):
        """Write the dat-to-skl atom index map file.

        This file records the identity of every atom in a fully-expanded OLCAO
        .dat structure so that other scripts can map back to the original
        skeleton (.skl) atom numbers.  It is the companion to
        `read_dat_skl_map`, which reads an equivalent file (typically saved as
        ``inputs/datSkl.map``) and builds two cross-reference arrays:

        - ``datSklMap[dat_atom]``  → corresponding skeleton atom number
        - ``sklDatMap[skl_atom]``  → corresponding dat atom number

        File format
        -----------
        Line 1 — fixed header::

            OATOM#   MOL#   RES#   RESNAME   ELEMENT   TYPE

        Lines 2..N+1 — one line per atom, formatted as
        ``%6d %7d %7d %10s %10s %7s``::

            OATOM#  : (int, w=6)  OLCAO atom index, 1..num_atoms
            MOL#    : (int, w=7)  molecule sequence number (from PDB/HIN input,
                                  or 1 for crystal structures)
            RES#    : (int, w=7)  residue sequence number (from PDB/HIN input,
                                  or 1 for crystal structures)
            RESNAME : (str, w=10) residue name (from PDB/HIN input, or element
                                  symbol for crystal structures)
            ELEMENT : (str, w=10) lowercase element symbol (e.g. 'si')
            TYPE    : (str, w=7)  OLCAO atom tag, i.e. the species identifier
                                  used in the .skl file (e.g. 'si_tb2')

        Parameters
        ----------
        filename : str, optional
            Output path.  Defaults to ``'olcao.map'``.
        """
        if filename is None:
            filename = 'olcao.map'

        with open(filename, 'w') as fh:
            fh.write('OATOM#   MOL#   RES#   RESNAME   ELEMENT   TYPE\n')
            for atom in range(1, self.num_atoms + 1):
                fh.write(
                    f'{atom:6d}'
                    f'{self.molecule_seq_num[atom]:7d}'
                    f'{self.residue_seq_num[atom]:7d}'
                    f'{self.residue_name[atom]:>10s}'
                    f'{self.atom_element_name[atom]:>10s}'
                    f'{self.atom_tag[atom]:>7s}\n'
                )

    def print_vasp(self, filename='POSCAR'):
        """Write the current structure in VASP POSCAR format.

        The Perl ``printVASP`` subroutine wrote four files in one call:
        POSCAR, POTCAR, KPOINTS, and INCAR.  This Python port only writes the
        POSCAR.  The POTCAR logic was not ported because it relies on
        shell-level decompression (``zcat``) and an external pseudopotential
        database referenced through environment variables
        (``$VASPPOT_DIR``, ``$VASPPOT_USPP_LDA``, ``$VASPPOT_PAW_PBE``,
        etc.).  If POTCAR/KPOINTS/INCAR generation is needed it should be
        handled as a separate step.

        POSCAR format
        -------------
        The POSCAR is written entirely in OLCAO's internal Bohr units with a
        universal scaling factor of ``1.000`` (i.e. VASP is told the vectors
        are already in the final units, and the factor is a no-op).

        Line 1 — title: ``" System <system_title[0]>"``

        Line 2 — scaling factor: ``1.000``

        Lines 3–5 — lattice vectors, one row per ``abc`` axis, three columns
        per ``xyz`` component, formatted as ``%16.8e``.

        Line 6 — space-separated atom count for each element, in the order
        produced by ``create_element_list`` / ``count_element_atoms``.

        Line 7 — coordinate mode keyword.  VASP uses the word ``direct`` to
        mean fractional coordinates, which is the opposite of OLCAO's
        convention.  The mapping is:

        - VASP ``direct``    == OLCAO ``fract``  (fractions of lattice vectors)
        - VASP ``Cartesian`` == OLCAO ``direct``  (Cartesian/XYZ in Bohr)

        Lines 8..N+7 — one fractional-coordinate line per atom, three columns
        formatted as ``%16.8e``.

        Perl POTCAR notes (not ported)
        ------------------------------
        The Perl code assembled a POTCAR by concatenating one pseudopotential
        per element from a directory tree rooted at ``$VASPPOT_DIR``.
        Potential types were selected by a ``potType`` argument:
        ``potLDA``, ``potGGA``, ``potpawLDA``, ``potpawGGA``, ``potpawPBE``,
        ``potpawLDA5x``, ``potpawPBE5x``.  For PAW-type potentials the element
        directory names use ``ucfirst`` capitalisation (e.g. ``Si``).  An
        optional ``subPotType`` suffix allowed selecting a sub-directory
        variant (e.g. ``_sv``, ``_pv``).

        Perl KPOINTS notes (not ported)
        --------------------------------
        If ``doGamma == 1`` a single Gamma-only KPOINTS file was written.
        Otherwise a Monkhorst-Pack mesh was generated with sizes chosen by
        cell dimension: if the lattice parameter > 10 Å → 1 k-point;
        > 5 Å → 2; otherwise 4.  The shift was always ``0 0 0``.

        Perl INCAR notes (not ported)
        ------------------------------
        The smearing method (``ISMEAR``) was set to ``-5`` (linear analytic
        tetrahedron) when the total number of k-points exceeded 4, otherwise
        ``0`` (Gaussian).  The relaxation flag ``ISIF`` was controlled by a
        ``jobType`` argument: ``relaxfull`` → 3, ``relaxion1`` → 2,
        ``relaxion2`` → 4, ``relaxvol`` → 7.  A canned INCAR template with
        ``ENCUT = 600 eV``, ``EDIFF = 1.0E-5``, ``EDIFFG = -1.0E-3``, etc.
        was then written.

        Parameters
        ----------
        filename : str, optional
            Output path.  Defaults to ``'POSCAR'``.
        """
        self.create_element_list()
        self.count_element_atoms()

        with open(filename, 'w') as f:
            # Title line — uses the first line of the system title block.
            title = self.system_title[0] if self.system_title else ''
            f.write(f' System {title}\n')
            # Universal scaling factor.  Raw Bohr units are written, so the
            # factor is 1.000 (a no-op in VASP).
            f.write('1.000\n')
            # Lattice vectors: one row per abc axis, three xyz columns each.
            for abc in range(1, 4):
                f.write(''.join(f'{self.real_lattice[abc][xyz]:16.8e}'
                                for xyz in range(1, 4)) + '\n')
            # Atom counts per element, in element_list order.
            counts = [str(self.element_count[e])
                      for e in range(1, self.num_elements + 1)]
            f.write(' '.join(counts) + '\n')
            # "direct" is VASP's word for fractional coordinates (== OLCAO
            # "fract").  Do not confuse with OLCAO "direct" (== VASP
            # "Cartesian").
            f.write('direct\n')
            for atom in range(1, self.num_atoms + 1):
                f.write(''.join(f'{self.fract_abc[atom][abc]:16.8e}'
                                for abc in range(1, 4)) + '\n')

    def print_pdb(self, filename=None):
        """Write the current structure in PDB format.

        Parameters
        ----------
        filename : str, optional
            Output path.  Defaults to self.model_pdb.

        PDB output format notes
        -----------------------
        The file opens with a HEADER record and three REMARK 1 lines, then a
        CRYST1 cell record, ATOM records for every atom, and a TER terminator.

        Perl counterpart (printPDB) accepted a second argument $pdbFormat and
        only wrote the CRYST1 record when $pdbFormat eq "crystal".  This
        Python version always writes CRYST1, which is appropriate for the
        OLCAO use-case where cell parameters are always defined.

        CRYST1 record layout (values in Bohr — OLCAO internal units):
          cols  1- 6   "CRYST1"
          cols  7-15   a (9.3f)
          cols 16-24   b (9.3f)
          cols 25-33   c (9.3f)
          cols 34-40   alpha (7.2f, degrees)
          cols 41-47   beta  (7.2f, degrees)
          cols 48-54   gamma (7.2f, degrees)
          cols 55-64   space-group name (right-justified in 10 chars)

        ATOM record layout (see:
        www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM):
          cols  1- 6   "ATOM  " (record type, 6s)
          cols  7-11   atom serial number (5d)
          cols 12-14   element symbol uppercase (3s, right-justified)
          cols 15-16   species ID (2d, left-justified)
          cols 17-22   residue name + chain: " MOL A" (6s)
          cols 23-26   residue sequence number: 1 (4d)
          cols 27-28   insertion code: "  " (2s, blank)
          cols 29-36   x (8.3f, Bohr)
          cols 37-44   y (8.3f, Bohr)
          cols 45-52   z (8.3f, Bohr)
          cols 53-58   occupancy: 1.00 (6.2f)
          cols 59-64   temperature factor: 0.00 (6.2f)
          cols 65-74   segment identifier: blank (10s)
          cols 75-76   element symbol uppercase (2s, right-justified)

        Coordinate units: OLCAO stores atom positions in Bohr internally.
        Unlike the PDB standard (Angstroms), output coordinates here are in
        Bohr.  Applications reading this file should apply the Bohr→Angstrom
        conversion as needed.

        Connectivity (CONECT) records after TER are not written; the Perl
        original also left this section empty as a future placeholder.
        """
        if filename is None:
            filename = self.model_pdb

        with open(filename, 'w') as f:
            # Print the HEADER and header-like REMARKS.
            f.write('HEADER\n')
            f.write('REMARK   1\n')
            f.write('REMARK   1 Created by print_pdb.\n')
            f.write('REMARK   1\n')

            # Print crystal cell parameters (CRYST1 record).  The Perl version
            # made this conditional on pdbFormat == "crystal"; we always write it.
            f.write(
                f'{"CRYST1":6s}'
                f'{self.mag[1]:9.3f}{self.mag[2]:9.3f}{self.mag[3]:9.3f}'
                f'{self.angle_deg[1]:7.2f}{self.angle_deg[2]:7.2f}'
                f'{self.angle_deg[3]:7.2f}'
                f'{self.space_group_name:>10s}\n'
            )

            # Print each atom in PDB ATOM record format.
            for atom in range(1, self.num_atoms + 1):
                elem_upper = self.atom_element_name[atom].upper()
                f.write(
                    f'{"ATOM  ":6s}'
                    f'{atom:5d}'
                    f'{elem_upper:>3s}'
                    f'{self.atom_species_id[atom]:<2d}'
                    f'{"":6s}'
                    f'{"MOL A":4s}'
                    f'{1:4d}'
                    f'{"  ":2s}'
                    f'{self.direct_xyz[atom][1]:8.3f}'
                    f'{self.direct_xyz[atom][2]:8.3f}'
                    f'{self.direct_xyz[atom][3]:8.3f}'
                    f'{1.0:6.2f}'
                    f'{0.0:6.2f}'
                    f'{"":10s}'
                    f'{elem_upper:>2s}\n'
                )

            # Terminate the ATOM section.
            # Connectivity (CONECT) records are not yet implemented.
            f.write('TER\n')

    def print_cif(self, filename=None):
        """Write the current structure in CIF (Crystallographic Information File) format.

        CIF is a standard text format for crystallographic data.  The file is
        structured as a series of key/value pairs and ``loop_`` blocks.

        Data-block name
        ---------------
        The CIF ``data_`` block identifier is derived from the filename stem
        (the filename with the ``.cif`` suffix stripped).  In the Perl original
        this name was passed explicitly as ``$currentName`` by the caller; here
        we derive it automatically.

        Symmetry / space group
        ----------------------
        OLCAO always writes fully-expanded P1 structures — no symmetry
        reduction is applied before output.  Consequently the cell setting is
        hardcoded to ``triclinic`` and the space group to ``'P 1'`` with the
        single symmetry operation ``X,Y,Z`` (the identity).  The Perl version
        accepted a ``$cellName`` argument so the caller could supply a
        descriptive name (e.g. ``cubic``), but the space group was still P 1.

        Cell parameters
        ---------------
        OLCAO stores lattice magnitudes internally in Bohr radii.  CIF requires
        Angstroms, so each magnitude is multiplied by BOHR_RAD before output.
        The inter-axial angles are stored in ``self.angle_deg`` (degrees) and
        are written directly.  The ordering is alpha (bc), beta (ac),
        gamma (ab) matching the standard CIF tags ``_cell_angle_alpha/beta/gamma``.

        Atom loop
        ---------
        Five fields are written per atom:

        ``_atom_site_label``
            A unique site identifier built as ``{element}{species_id}_{atom_index}``,
            e.g. ``si1_3`` for the third atom of silicon species 1.  In the
            Perl original the trailing integer came from a ``sortedAtomTypeID``
            array passed in by the caller; here we use the sequential atom
            index since atoms are written in their natural order.

        ``_atom_site_type_symbol``
            The chemical-element symbol with the first letter capitalised
            (Perl: ``ucfirst($element)``; Python: ``elem.capitalize()``),
            e.g. ``Si``, ``B``, ``Na``.

        ``_atom_site_fract_x/y/z``
            Fractional coordinates, formatted to 12.6f.

        Parameters
        ----------
        filename : str, optional
            Output path.  Defaults to self.model_cif.
        """
        if filename is None:
            filename = self.model_cif

        # The CIF data-block identifier is the filename stem (no .cif suffix).
        stem = filename.replace('.cif', '')

        # CIF requires Angstroms; OLCAO stores magnitudes in Bohr internally.
        a_ang = self.mag[1] * BOHR_RAD
        b_ang = self.mag[2] * BOHR_RAD
        c_ang = self.mag[3] * BOHR_RAD

        with open(filename, 'w') as f:
            # Header: data-block name and symmetry information.
            # Always P 1 (triclinic) — OLCAO writes fully-expanded structures.
            f.write(f'data_{stem}\n')
            f.write('_symmetry_cell_setting triclinic\n')
            f.write("_symmetry_space_group_name_H-M 'P 1'\n")
            f.write('loop_\n')
            f.write('   _symmetry_equiv_pos_as_xyz\n')
            f.write('              X,Y,Z\n')  # Identity — the only symmetry op.

            # Cell parameters: lengths in Angstroms, angles already in degrees.
            f.write(f'_cell_length_a {a_ang}\n')
            f.write(f'_cell_length_b {b_ang}\n')
            f.write(f'_cell_length_c {c_ang}\n')
            f.write(f'_cell_angle_alpha {self.angle_deg[1]}\n')
            f.write(f'_cell_angle_beta  {self.angle_deg[2]}\n')
            f.write(f'_cell_angle_gamma {self.angle_deg[3]}\n')

            # Atom-site loop header.
            f.write('loop_\n')
            f.write('   _atom_site_label\n')
            f.write('   _atom_site_type_symbol\n')
            f.write('   _atom_site_fract_x\n')
            f.write('   _atom_site_fract_y\n')
            f.write('   _atom_site_fract_z\n')

            # One line per atom: label, type symbol, fractional coordinates.
            # Label format: {elem}{species_id}_{atom_index}  e.g. "si1_3"
            # type_symbol:  elem.capitalize()  (ucfirst in Perl)  e.g. "Si"
            for atom in range(1, self.num_atoms + 1):
                eid  = self.atom_element_id[atom]
                elem = self.element_list[eid]
                sid  = self.atom_species_id[atom]
                f.write(
                    f'   {elem:2s}{sid:<1d}{"_"}{atom:<3d} '
                    f'{elem.capitalize():2s} '
                    f'{self.fract_abc[atom][1]:12.6f} '
                    f'{self.fract_abc[atom][2]:12.6f} '
                    f'{self.fract_abc[atom][3]:12.6f}\n'
                )

    def print_lmp(self, filename=None):
        """Write the current structure as a LAMMPS data file.

        Produces a LAMMPS "data" file with the following sections in order:
          1. Header line (filename) followed by a blank line.
          2. Counts: "N atoms" and "M atom types".
          3. Box bounds: orthogonal extents as "0.000 <mag> xlo xhi" for each
             of the three axes (a→x, b→y, c→z).
          4. Tilt factors (xy xz yz) — written only when the cell is non-
             orthogonal (any lattice angle != 90°).  The formulas follow the
             LAMMPS triclinic convention:
               xy = b * cos(gamma)
               xz = c * cos(beta)
               yz = (b*c*cos(alpha) - xy*xz) / sqrt(b^2 - xy^2)
             See: https://lammps.sandia.gov/doc/Howto_triclinic.html
          5. "Masses" section: one line per unique atom type (sequential integer
             starting at 1), with the atomic mass and a comment giving the
             element symbol and species ID.
          6. "Atoms" section: one line per atom with
               atom_id  type_id  x  y  z  # element_name species_id
             Coordinates are Cartesian (direct_xyz), in Angstroms.

        LAMMPS "atom type" corresponds to a unique (element, species) pair.
        The type ID is assigned by walking the atom list in order and
        incrementing the counter each time the element or species changes.
        Hence the atoms must already be sorted by element/species for the
        type IDs to be contiguous — which is guaranteed by the OLCAO atom
        ordering conventions.

        Parameters
        ----------
        filename : str, optional
            Output path.  Defaults to self.model_lmp.
        """
        if filename is None:
            filename = self.model_lmp

        # Walk atoms once to assign each atom its unique-species type ID.
        # A new type is recorded whenever the element or species changes
        # relative to the previous atom (only types, not individual atoms,
        # are printed in the Masses section).
        num_unique = 0
        current_element = 0
        current_species  = 0
        atom_uniq = [None] * (self.num_atoms + 1)
        for atom in range(1, self.num_atoms + 1):
            eid = self.atom_element_id[atom]
            sid = self.atom_species_id[atom]
            if eid != current_element:
                current_element = eid
                current_species  = sid
                num_unique += 1
            elif sid != current_species:
                current_species = sid
                num_unique += 1
            atom_uniq[atom] = num_unique
        num_types = num_unique

        with open(filename, 'w') as f:
            # Header: filename followed by a required blank line.
            f.write(f'{filename}\n\n')

            # Counts section.
            f.write(f'{self.num_atoms} atoms\n')
            f.write(f'{num_types} atom types\n')

            # Orthogonal box bounds: lower bound is always 0.
            f.write(f'0.000 {self.mag[1]} xlo xhi\n')
            f.write(f'0.000 {self.mag[2]} ylo yhi\n')
            f.write(f'0.000 {self.mag[3]} zlo zhi\n')

            # Tilt factors are only needed for non-orthogonal (triclinic) cells.
            # angle[1]=alpha (b^c), angle[2]=beta (a^c), angle[3]=gamma (a^b).
            tol = 1e-10
            if (abs(self.angle_deg[1] - 90.0) > tol or
                    abs(self.angle_deg[2] - 90.0) > tol or
                    abs(self.angle_deg[3] - 90.0) > tol):
                xy = self.mag[2] * math.cos(self.angle[3])
                xz = self.mag[3] * math.cos(self.angle[2])
                yz = (self.mag[2]*self.mag[3]*math.cos(self.angle[1]) - xy*xz) \
                     / math.sqrt(self.mag[2]**2 - xy**2)
                f.write(f'{xy} {xz} {yz} xy xz yz\n\n')
            else:
                f.write('\n')

            # Masses section: one line per unique type (element+species pair).
            # The same change-detection logic as above is used so that the
            # type IDs here exactly match those assigned in atom_uniq[].
            f.write('Masses\n\n')
            num_unique = 0
            current_element = 0
            current_species  = 0
            for atom in range(1, self.num_atoms + 1):
                eid = self.atom_element_id[atom]
                sid = self.atom_species_id[atom]
                if eid != current_element:
                    current_element = eid
                    current_species  = sid
                    num_unique += 1
                    z = self.atomic_z[atom]
                    mass = self.element_data.atomic_masses[z]
                    f.write(f'{num_unique} {mass} '
                            f'# {self.atom_element_name[atom]} {sid}\n')
                elif sid != current_species:
                    current_species = sid
                    num_unique += 1
                    z = self.atomic_z[atom]
                    mass = self.element_data.atomic_masses[z]
                    f.write(f'{num_unique} {mass} '
                            f'# {self.atom_element_name[atom]} {sid}\n')

            # Atoms section: atom_id, type_id, Cartesian x y z, comment.
            f.write('\nAtoms\n\n')
            for atom in range(1, self.num_atoms + 1):
                f.write(
                    f'{atom} {atom_uniq[atom]} '
                    f'{self.direct_xyz[atom][1]} '
                    f'{self.direct_xyz[atom][2]} '
                    f'{self.direct_xyz[atom][3]} '
                    f'# {self.atom_element_name[atom]} '
                    f'{self.atom_species_id[atom]}\n'
                )

    def print_isaacs(self, filename=None):
        """Write the current structure in ISAACS XML + Chem3D format.

        ISAACS (Interactive Structure Analysis of Amorphous and Crystal
        Systems) is a tool for analysing structural properties of both
        crystalline and disordered materials.  It uses a two-file scheme:

          - <filename>.ipf  : ISAACS XML project file (I.S.A.A.C.S. v1.1)
          - <filename>.c3d  : Chem3D coordinate file referenced by the IPF

        The .ipf file is an XML document with the following major sections:

        <data>
            Names the companion C3D file using an absolute path
            (cwd + '/' + c3d_file), so the IPF and C3D should be used
            from the same working directory in which they were generated.

        <chemistry>
            <atoms>     — total atom count
            <species>   — element list; each <label id="N"> uses 0-based
                          element IDs (element index − 1)
            <element>   — per-element block with:
                <name>    full element name (e.g. "Silicon")
                <z>       atomic number Z
                <mass>    atomic mass (from element database)
                <rad>     always 0 (not used by ISAACS in this workflow)
                <radius>  atomic radius (from element database)
                <nscatt>  neutron scattering length (from element database)
                <xscatt>  X-ray scattering factor; written as Z (integer
                          approximation sufficient for ISAACS)

        <box>
            Lattice edge magnitudes (a, b, c in Angstroms), angles
            (alpha, beta, gamma in degrees), and the full 3×3 Cartesian
            lattice vector matrix (a.x, a.y, a.z, b.x, …, c.z).

        <pbc>
            Periodic boundary conditions are always enabled (apply=TRUE).
            Coordinates in the C3D file are Cartesian, not fractional
            (fractional=FALSE, fractype=0).

        <cutoffs>
            <total>    — minimum over all element pairs of
                         (r_cov_i + r_cov_j) * 2.0
            <partials> — one <Ei-Ej> tag per ordered pair with the
                         pair-specific value (r_cov_i + r_cov_j) * 2.0,
                         formatted to 6 decimal places.
            The factor of 2.0 gives a generous outer shell around the
            sum-of-covalent-radii bond criterion.

        <project>
            Always TRUE; tells ISAACS to apply the projection when the
            project is loaded.

        The .c3d (Chem3D) file format is minimal:
            Line 1:  <num_atoms>
            Lines 2…: <ElementSymbol> <atom_index> <x> <y> <z>
        where coordinates are Cartesian (Angstroms) and atom_index is the
        1-based atom number.  Element symbols are title-cased (e.g. 'Si').

        Parameters
        ----------
        filename : str, optional
            Base name (without extension).  Defaults to 'model'.
        """
        if filename is None:
            filename = 'model'

        ipf_file = filename + '.ipf'
        c3d_file = filename + '.c3d'

        ed = self.element_data
        import os
        pwd = os.getcwd()

        # Z-number per element (1-indexed by element position in element_list).
        z_of_elem = [None] * (self.num_elements + 1)
        for e in range(1, self.num_elements + 1):
            z_of_elem[e] = ed.get_element_z(self.element_list[e])

        with open(ipf_file, 'w') as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write('<!-- I.S.A.A.C.S. v1.1 XML file -->\n')
            f.write('<isaacs-xml>\n')
            f.write(' <!-- Data format and file containing the configurations(s) -->\n')
            f.write(' <data>\n')
            f.write('  <type>Chem3D file</type>\n')
            f.write(f'  <file>{pwd}/{c3d_file}</file>\n')
            f.write(' </data>\n')
            f.write(' <!-- Chemistry information -->\n')
            f.write(' <chemistry>\n')
            f.write(f'  <atoms>{self.num_atoms}</atoms>\n')
            f.write(f'  <species number="{self.num_elements}">\n')
            for e in range(1, self.num_elements + 1):
                elem_cap = self.element_list[e].capitalize()
                f.write(f'   <label id="{e-1}">{elem_cap:<2s}</label>\n')
            f.write('  </species>\n')
            for e in range(1, self.num_elements + 1):
                elem_cap = self.element_list[e].capitalize()
                z = z_of_elem[e]
                full_name = ed.element_full_names[z].capitalize()
                f.write(f'  <element symbol="{elem_cap:<2s}">\n')
                f.write(f'   <name>{full_name:<15s}</name>\n')
                f.write(f'   <z>{z}</z>\n')
                f.write(f'   <mass>{ed.atomic_masses[z]}</mass>\n')
                f.write('   <rad>0</rad>\n')
                f.write(f'   <radius>{ed.atomic_radii[z]}</radius>\n')
                f.write(f'   <nscatt>{ed.neut_scatt[z]}</nscatt>\n')
                f.write(f'   <xscatt>{z}</xscatt>\n')
                f.write('  </element>\n')
            f.write(' </chemistry>\n')
            f.write(' <!-- Box information -->\n')
            f.write(' <box>\n')
            f.write('  <edges>\n')
            f.write(f'   <a>{self.mag[1]}</a>\n')
            f.write(f'   <b>{self.mag[2]}</b>\n')
            f.write(f'   <c>{self.mag[3]}</c>\n')
            f.write('  </edges>\n')
            f.write('  <angles>\n')
            f.write(f'   <alpha>{self.angle_deg[1]}</alpha>\n')
            f.write(f'   <beta>{self.angle_deg[2]}</beta>\n')
            f.write(f'   <gamma>{self.angle_deg[3]}</gamma>\n')
            f.write('  </angles>\n')
            f.write('  <vectors>\n')
            f.write(f'   <a.x>{self.real_lattice[1][1]}</a.x>\n')
            f.write(f'   <a.y>{self.real_lattice[1][2]}</a.y>\n')
            f.write(f'   <a.z>{self.real_lattice[1][3]}</a.z>\n')
            f.write(f'   <b.x>{self.real_lattice[2][1]}</b.x>\n')
            f.write(f'   <b.y>{self.real_lattice[2][2]}</b.y>\n')
            f.write(f'   <b.z>{self.real_lattice[2][3]}</b.z>\n')
            f.write(f'   <c.x>{self.real_lattice[3][1]}</c.x>\n')
            f.write(f'   <c.y>{self.real_lattice[3][2]}</c.y>\n')
            f.write(f'   <c.z>{self.real_lattice[3][3]}</c.z>\n')
            f.write('  </vectors>\n')
            f.write(' </box>\n')
            f.write(' <!-- PBC information -->\n')
            f.write(' <pbc>\n')
            f.write('  <apply>TRUE</apply>\n')
            f.write('  <fractional>FALSE</fractional>\n')
            f.write('  <fractype>0</fractype>\n')
            f.write(' </pbc>\n')
            f.write(' <!-- Cutoffs information -->\n')
            f.write(' <cutoffs>\n')
            min_cutoff = 1e8
            for e1 in range(1, self.num_elements + 1):
                for e2 in range(1, self.num_elements + 1):
                    curr = (ed.coval_radii[z_of_elem[e1]] +
                            ed.coval_radii[z_of_elem[e2]]) * 2.0
                    if curr < min_cutoff:
                        min_cutoff = curr
            f.write(f'  <total>{min_cutoff}</total>\n')
            f.write('  <partials>\n')
            for e1 in range(1, self.num_elements + 1):
                n1 = self.element_list[e1].capitalize()
                for e2 in range(1, self.num_elements + 1):
                    n2 = self.element_list[e2].capitalize()
                    val = (ed.coval_radii[z_of_elem[e1]] +
                           ed.coval_radii[z_of_elem[e2]]) * 2.0
                    f.write(f'   <{n1}-{n2}>{val:.6f}</{n1}-{n2}>\n')
            f.write('  </partials>\n')
            f.write(' </cutoffs>\n')
            f.write(' <!-- Apply project -->\n')
            f.write(' <project>TRUE</project>\n')
            f.write('</isaacs-xml>\n')

        with open(c3d_file, 'w') as f:
            f.write(f'{self.num_atoms}\n')
            for atom in range(1, self.num_atoms + 1):
                elem_cap = self.atom_element_name[atom].capitalize()
                x = self.direct_xyz[atom][1]
                y = self.direct_xyz[atom][2]
                z = self.direct_xyz[atom][3]
                f.write(f'{elem_cap} {atom}  {x} {y} {z}\n')

    # ======================================================================
    # Section: Lattice and coordinate geometry
    # ======================================================================

    def abc_alpha_beta_gamma(self):
        """Compute magnitudes and angles between vectors in a lattice.

        Originally written for real lattices only, but later adapted to serve
        reciprocal lattices as well.  For reciprocal lattices the same logic
        applies; the caller simply passes k_a/k_b/k_c and the results land in
        the reciprocal-space mag/angle/angle_deg arrays.  The inline notes
        below are written in terms of the real lattice, but everything is
        identical for the reciprocal case.

        Angles are defined as:
            alpha  -- angle between b and c  (angle[1])
            beta   -- angle between a and c  (angle[2])
            gamma  -- angle between a and b  (angle[3])

        Formula: angle = arccos(v1 · v2 / (|v1| |v2|)), result in radians.

        Populates: mag, angle, angle_deg, sine_angle.
        Note: sine_angle is only meaningful for real-space lattices.
        """
        # Compute the magnitude of each lattice vector (a, b, c).
        for axis in range(1, 4):
            self.mag[axis] = math.sqrt(
                sum(self.real_lattice[axis][k]**2 for k in range(1, 4)))

        # Dot products of each pair of lattice vectors used for angle calc.
        dot_bc = sum(self.real_lattice[2][k]*self.real_lattice[3][k] for k in range(1, 4))
        dot_ac = sum(self.real_lattice[1][k]*self.real_lattice[3][k] for k in range(1, 4))
        dot_ab = sum(self.real_lattice[1][k]*self.real_lattice[2][k] for k in range(1, 4))

        # Clamp the cosine argument to [-1, 1] to guard against floating-point
        # rounding that would make acos raise a domain error on cubic cells.
        self.angle[1] = math.acos(max(-1.0, min(1.0, dot_bc / (self.mag[2]*self.mag[3]))))
        self.angle[2] = math.acos(max(-1.0, min(1.0, dot_ac / (self.mag[1]*self.mag[3]))))
        self.angle[3] = math.acos(max(-1.0, min(1.0, dot_ab / (self.mag[1]*self.mag[2]))))

        for axis in range(1, 4):
            self.angle_deg[axis] = self.angle[axis] * 180.0 / PI

        # Compute the sine of each angle (real lattices only).
        self.get_angle_sine()

    def get_angle_sine(self):
        """Compute and store sine_angle[1..3] from the current angle[1..3] values.

        Angle/axis mapping (all values in radians):
          angle[1] = alpha  — angle between lattice vectors b and c
          angle[2] = beta   — angle between lattice vectors a and c
          angle[3] = gamma  — angle between lattice vectors a and b

        Note: presently this is only meaningful for real-space lattices.
        Reciprocal-lattice angles are handled separately and do not go through
        this routine.
        """
        for axis in range(1, 4):
            self.sine_angle[axis] = math.sin(self.angle[axis])

    def get_abc_vectors(self):
        """Compute Cartesian lattice vectors from magnitudes and angles.

        Given the space group, define the real-space lattice vectors (stored in
        real_lattice[abc][xyz]) according to the conventional-cell conventions
        of W. Setyawan and S. Curtarolo, Comp. Mat. Sci. v. 49, pp. 299-312
        (2010).

        If no space group has been set yet, space group 1 (triclinic P1) is
        assumed so that the general-triclinic path is always available as a
        fallback.

        The seven crystal systems are handled in order of their space-group
        number ranges:

        Triclinic   (sg 1-2):   alpha < beta < gamma, then a < b < c if possible.
                                a is co-axial with x; b is in the x,y plane;
                                c has components in all three directions.

        Monoclinic  (sg 3-15):  a <= b <= c; one angle < 90, the other two = 90.
                                Three unique-axis settings are encoded via the
                                sub-group number: y(b)-axis unique (most common),
                                z(c)-axis unique, or x(a)-axis unique.

        Orthorhombic (sg 16-74):  a < b < c; alpha = beta = gamma = 90.

        Tetragonal  (sg 75-142): a = b ≠ c; alpha = beta = gamma = 90.

        Trigonal    (sg 143-167): Two sub-cases depending on lattice type and
                                  sub-group number:
                                  - Hexagonal cell (P, H, or R with hexagonal
                                    axes, sg_sub == 2): a = b ≠ c;
                                    alpha = beta = 90, gamma = 120.  H-centered
                                    cells are made from three hexagonal cells;
                                    R with hexagonal axes from three primitive
                                    rhombohedral cells.
                                  - Primitive rhombohedral cell (R, sg_sub == 1):
                                    a = b = c; alpha = beta = gamma ≠ 90.

        Hexagonal   (sg 168-194): a = b ≠ c; alpha = beta = 90, gamma = 120.

        Cubic       (sg 195-230): a = b = c; alpha = beta = gamma = 90.

        After all components are set, values smaller than EPSILON are zeroed
        to remove floating-point noise from trigonometric calculations.
        """
        # If no space group has been defined yet, default to P1 (triclinic).
        if not self.space_group:
            self.space_group = '1'
        sg     = self.space_group_num
        sg_sub = self.space_group_sub_num
        lt     = self.lattice_type
        m      = self.mag
        a      = self.angle
        rl     = self.real_lattice

        if sg <= 2:  # Triclinic: alpha<beta<gamma then a<b<c if possible
            # a is co-axial with x.
            rl[1][1] = m[1]
            rl[1][2] = 0.0
            rl[1][3] = 0.0
            # b is in the x,y plane.
            rl[2][1] = m[2] * math.cos(a[3])   # b*cos(gamma)
            rl[2][2] = m[2] * math.sin(a[3])   # b*sin(gamma)
            rl[2][3] = 0.0
            # c has components in all three directions.
            rl[3][1] = m[3] * math.cos(a[2])   # c*cos(beta)
            rl[3][2] = m[3] * (math.cos(a[1]) - math.cos(a[3])*math.cos(a[2])) / math.sin(a[3])
            rl[3][3] = m[3] * math.sqrt(1.0 - math.cos(a[2])**2 - (rl[3][2]/m[3])**2)

        elif sg <= 15:  # Monoclinic: a<=b<=c; one angle<90; other two=90
            # Determine which axis is the unique (non-90-degree) axis from the
            # sub-group number.  y(b)-axis unique is the most common setting.
            if ((sg == 3  and sg_sub <=  2) or (sg == 4  and sg_sub <=  2) or
                (sg == 5  and sg_sub <=  4) or (sg == 6  and sg_sub <=  2) or
                (sg == 7  and sg_sub <=  5) or (sg == 8  and sg_sub <=  4) or
                (sg == 9  and sg_sub <=  4) or (sg == 10 and sg_sub <=  2) or
                (sg == 11 and sg_sub <=  2) or (sg == 12 and sg_sub <=  4) or
                (sg == 13 and sg_sub <=  5) or (sg == 14 and sg_sub <=  5) or
                (sg == 15 and sg_sub <=  7)):
                unique_axis = 2  # y(b)-axis unique
            elif ((sg == 3  and sg_sub <=  4) or (sg == 4  and sg_sub <=  4) or
                  (sg == 5  and sg_sub <=  8) or (sg == 6  and sg_sub <=  4) or
                  (sg == 7  and sg_sub <= 10) or (sg == 8  and sg_sub <=  8) or
                  (sg == 9  and sg_sub <=  8) or (sg == 10 and sg_sub <=  4) or
                  (sg == 11 and sg_sub <=  4) or (sg == 12 and sg_sub <=  8) or
                  (sg == 13 and sg_sub <= 10) or (sg == 14 and sg_sub <= 10) or
                  (sg == 15 and sg_sub <= 14)):
                unique_axis = 3  # z(c)-axis unique
            else:
                unique_axis = 1  # x(a)-axis unique

            if unique_axis == 1:  # x(a)-axis unique: beta is the oblique angle
                # a = (a*sin(beta), 0, a*cos(beta))
                rl[1][1] = m[1] * math.sin(a[2]);  rl[1][2] = 0.0;  rl[1][3] = m[1] * math.cos(a[2])
                # b = (0, b, 0)
                rl[2][1] = 0.0;                     rl[2][2] = m[2]; rl[2][3] = 0.0
                # c = (0, 0, c)
                rl[3][1] = 0.0;                     rl[3][2] = 0.0;  rl[3][3] = m[3]
            elif unique_axis == 2:  # y(b)-axis unique: gamma is the oblique angle
                # a = (a, 0, 0)
                rl[1][1] = m[1];                    rl[1][2] = 0.0;                  rl[1][3] = 0.0
                # b = (b*cos(gamma), b*sin(gamma), 0)
                rl[2][1] = m[2] * math.cos(a[3]);  rl[2][2] = m[2]*math.sin(a[3]); rl[2][3] = 0.0
                # c = (0, 0, c)
                rl[3][1] = 0.0;                     rl[3][2] = 0.0;                  rl[3][3] = m[3]
            else:  # unique_axis == 3: z(c)-axis unique: alpha is the oblique angle
                # a = (a, 0, 0)
                rl[1][1] = m[1];  rl[1][2] = 0.0;  rl[1][3] = 0.0
                # b = (0, b, 0)
                rl[2][1] = 0.0;   rl[2][2] = m[2]; rl[2][3] = 0.0
                # c = (0, c*cos(alpha), c*sin(alpha))
                rl[3][1] = 0.0;   rl[3][2] = m[3] * math.cos(a[1]); rl[3][3] = m[3] * math.sin(a[1])

        elif sg <= 74:  # Orthorhombic: a<b<c; alpha=beta=gamma=90
            rl[1][1] = m[1]; rl[1][2] = 0.0;  rl[1][3] = 0.0
            rl[2][1] = 0.0;  rl[2][2] = m[2]; rl[2][3] = 0.0
            rl[3][1] = 0.0;  rl[3][2] = 0.0;  rl[3][3] = m[3]

        elif sg <= 142:  # Tetragonal: a=b≠c; alpha=beta=gamma=90
            rl[1][1] = m[1]; rl[1][2] = 0.0;  rl[1][3] = 0.0
            rl[2][1] = 0.0;  rl[2][2] = m[2]; rl[2][3] = 0.0
            rl[3][1] = 0.0;  rl[3][2] = 0.0;  rl[3][3] = m[3]

        elif sg <= 167:  # Trigonal: two sub-cases (hexagonal vs rhombohedral cell)
            if lt == 'P' or lt == 'H' or (lt == 'R' and sg_sub == 2):
                # Hexagonal cell: a=b≠c; alpha=beta=90, gamma=120.
                # a = (a/2, -(a*sqrt(3))/2, 0)
                rl[1][1] =  m[1] / 2.0;  rl[1][2] = -(m[1]*math.sqrt(3.0)/2.0); rl[1][3] = 0.0
                # b = (b/2,  (b*sqrt(3))/2, 0)
                rl[2][1] =  m[2] / 2.0;  rl[2][2] =  (m[2]*math.sqrt(3.0)/2.0); rl[2][3] = 0.0
                # c = (0, 0, c)
                rl[3][1] =  0.0;          rl[3][2] =  0.0;                        rl[3][3] = m[3]
            else:  # Primitive rhombohedral cell: a=b=c; alpha=beta=gamma≠90
                # a = (a*cos(alpha/2), -a*sin(alpha/2), 0)
                rl[1][1] =  m[1] * math.cos(a[1]/2.0)
                rl[1][2] = -m[1] * math.sin(a[1]/2.0)
                rl[1][3] =  0.0
                # b = (a*cos(alpha/2), a*sin(alpha/2), 0)
                rl[2][1] =  m[1] * math.cos(a[1]/2.0)
                rl[2][2] =  m[1] * math.sin(a[1]/2.0)
                rl[2][3] =  0.0
                # c = (a*cos(alpha)/cos(alpha/2), 0,
                #      a*sqrt(1 - cos(alpha)^2 / cos(alpha/2)^2))
                rl[3][1] =  m[1] * math.cos(a[1]) / math.cos(a[1]/2.0)
                rl[3][2] =  0.0
                rl[3][3] =  m[1] * math.sqrt(1.0 - math.cos(a[1])**2 / math.cos(a[1]/2.0)**2)

        elif sg <= 194:  # Hexagonal: a=b≠c; alpha=beta=90, gamma=120
            # a = (a/2, -(a*sqrt(3))/2, 0)
            rl[1][1] =  m[1] / 2.0;  rl[1][2] = -(m[1]*math.sqrt(3.0)/2.0); rl[1][3] = 0.0
            # b = (b/2,  (b*sqrt(3))/2, 0)
            rl[2][1] =  m[2] / 2.0;  rl[2][2] =  (m[2]*math.sqrt(3.0)/2.0); rl[2][3] = 0.0
            # c = (0, 0, c)
            rl[3][1] =  0.0;          rl[3][2] =  0.0;                        rl[3][3] = m[3]

        elif sg <= 230:  # Cubic: a=b=c; alpha=beta=gamma=90
            rl[1][1] = m[1]; rl[1][2] = 0.0;  rl[1][3] = 0.0
            rl[2][1] = 0.0;  rl[2][2] = m[2]; rl[2][3] = 0.0
            rl[3][1] = 0.0;  rl[3][2] = 0.0;  rl[3][3] = m[3]

        # Zero out numerical noise from trigonometric calculations.
        for abc_axis in range(1, 4):
            for xyz_axis in range(1, 4):
                if abs(self.real_lattice[abc_axis][xyz_axis]) < EPSILON:
                    self.real_lattice[abc_axis][xyz_axis] = 0.0

    def make_inv_or_recip_lattice(self, lattice=None, lattice_inv=None,
                                   make_recip=True):
        """Compute the real-space lattice inverse and, optionally, the
        reciprocal lattice, and store the corresponding cell volumes.

        The method uses the cofactor / adjugate approach to invert a 3×3
        matrix without any external library.  The same cycling-index trick
        used in the Perl original is preserved so the two implementations
        can be compared directly:

            ic1 = (i % 3) + 1          # the "next" abc axis after i
            ic2 = ((i+1) % 3) + 1     # the "one after that" abc axis

        These select the two rows of the input lattice whose cross-product
        gives the cofactor column for output row i.  Accumulating
        sum_j(cofactor[i][j] * inLattice[j][i]) over the inner loop yields
        the determinant (cell volume) by cofactor expansion along column i.
        All output coefficients are in Angstroms (or Å⁻¹ for the inverse,
        2π/Å for the reciprocal).

        The Perl original also supported a third piFactor=-1 mode (scale by
        0.5/π) for inverting the reciprocal lattice, but that case is never
        called in practice and is not exposed here.

        An older, commented-out version of this Perl sub (also named
        makeInvOrRecipLattice) had a copy-paste bug where both inLattice_ref
        and outLattice_ref were set to $_[0], and it stored results into a
        transposed temporary array.  The current Perl version (and this port)
        stores directly into the output array as outLattice[i][j] with i as
        the output-row (abc_axis) index.

        Parameters
        ----------
        lattice : list or None
            1-indexed 4×4 lattice to invert.  Defaults to self.real_lattice.
        lattice_inv : list or None
            1-indexed 4×4 array to receive the plain inverse (piFactor=1).
            Defaults to self.real_lattice_inv.
        make_recip : bool
            If True (default), also compute the reciprocal lattice by
            repeating the cofactor expansion with piFactor = 2π, stored in
            self.recip_lattice.  The reciprocal cell volume is then
            (2π)³ / real_cell_volume.
            If False, only the plain inverse is computed and
            self.real_cell_volume is updated.
        """
        if lattice is None:
            lattice = self.real_lattice
        if lattice_inv is None:
            lattice_inv = self.real_lattice_inv
        latt = lattice
        inv  = lattice_inv
        last_vol = 0.0
        for i in range(1, 4):
            cell_vol = 0.0
            ic1 = (i % 3) + 1          # next abc axis (cycles 1→2→3→1)
            ic2 = ((i + 1) % 3) + 1   # one further (cycles 1→3→2→1)
            for j in range(1, 4):
                jc1 = (j % 3) + 1
                jc2 = ((j + 1) % 3) + 1
                # 2×2 minor: cross-product of the two "other" lattice rows
                inv[i][j] = (latt[jc1][ic1] * latt[jc2][ic2] -
                             latt[jc2][ic1] * latt[jc1][ic2])
                # Accumulate determinant by cofactor expansion along column i
                cell_vol += inv[i][j] * latt[j][i]
            for j in range(1, 4):
                inv[i][j] /= cell_vol  # units: Å⁻¹ (plain inverse)
            last_vol = cell_vol
        self.real_cell_volume = abs(last_vol)  # volume is always positive
        if make_recip:
            recip   = self.recip_lattice
            two_pi  = 2.0 * PI
            for i in range(1, 4):
                cell_vol = 0.0
                ic1 = (i % 3) + 1
                ic2 = ((i + 1) % 3) + 1
                for j in range(1, 4):
                    jc1 = (j % 3) + 1
                    jc2 = ((j + 1) % 3) + 1
                    recip[i][j] = (latt[jc1][ic1] * latt[jc2][ic2] -
                                   latt[jc2][ic1] * latt[jc1][ic2])
                    cell_vol += recip[i][j] * latt[j][i]
                for j in range(1, 4):
                    recip[i][j] = two_pi * recip[i][j] / cell_vol  # 2π/Å
            self.recip_cell_volume = two_pi**3 / self.real_cell_volume

    def get_min_max_xyz(self, atom1=None, atom2=None):
        """Scan atom Cartesian coordinates to populate max_pos and min_pos.

        Iterates over direct_xyz[atom][axis] for atoms atom1..atom2 and
        records the per-axis extrema in self.max_pos and self.min_pos
        (both 1-indexed: 1=x, 2=y, 3=z).

        BIG_REAL is used as a ±infinity sentinel so that every real coordinate
        will immediately displace the initial bounds on the first comparison.

        Parameters
        ----------
        atom1, atom2 : int, optional
            Inclusive 1-indexed atom range.  Defaults to all atoms.
        """
        if atom1 is None:
            atom1 = 1
        if atom2 is None:
            atom2 = self.num_atoms
        self.max_pos = [None, -BIG_REAL, -BIG_REAL, -BIG_REAL]  # 1=x 2=y 3=z
        self.min_pos = [None,  BIG_REAL,  BIG_REAL,  BIG_REAL]  # 1=x 2=y 3=z
        for atom in range(atom1, atom2 + 1):
            for axis in range(1, 4):
                if self.direct_xyz[atom][axis] > self.max_pos[axis]:
                    self.max_pos[axis] = self.direct_xyz[atom][axis]
                if self.direct_xyz[atom][axis] < self.min_pos[axis]:
                    self.min_pos[axis] = self.direct_xyz[atom][axis]

    def get_direct_xyz(self, atom):
        """Convert fractional abc to Cartesian xyz for a single atom.

        Reads self.fract_abc[atom][1..3] and self.real_lattice; writes the
        result into self.direct_xyz[atom][1..3].  Delegates the actual
        matrix multiply to get_direct_xyz_point().

        Parameters
        ----------
        atom : int
            Atom index (1-indexed).
        """
        xyz = self.get_direct_xyz_point(self.fract_abc[atom])
        for axis in range(1, 4):
            self.direct_xyz[atom][axis] = xyz[axis]

    def get_direct_xyz_point(self, fract):
        """Convert a single fractional coordinate to Cartesian via lattice matrix multiply.

        The transformation is:
            Px = Pa*ax + Pb*bx + Pc*cx
            Py = Pa*ay + Pb*by + Pc*cy
            Pz = Pa*az + Pb*bz + Pc*cz

        where Pa, Pb, Pc are the fractional (abc) coordinates and ax/ay/az,
        bx/by/bz, cx/cy/cz are the xyz-component coefficients of the a, b, c
        lattice vectors stored in self.real_lattice[abc_axis][xyz_axis].

        In other words: point[xyz] = sum over abc of ( fractABC[abc] * latt[abc][xyz] )

        Parameters
        ----------
        fract : sequence of float
            1-indexed [None, a, b, c] fractional coordinate.  Slot 0 is
            an unused sentinel matching the convention used everywhere in
            this module and mirrors Perl's ``@fractABC[1..3]``.

        Returns
        -------
        list of float
            1-indexed [None, x, y, z] Cartesian coordinate.  Slot 0 is a
            None sentinel so that callers index directly with axis values
            1, 2, 3 without an offset.
        """
        xyz = [None, 0.0, 0.0, 0.0]
        for xyz_axis in range(1, 4):
            for abc_axis in range(1, 4):
                xyz[xyz_axis] += (fract[abc_axis] *
                                  self.real_lattice[abc_axis][xyz_axis])
        return xyz

    def get_direct_abc(self, atom):
        """Compute direct-space abc coordinates from xyz for a single atom.

        Reads self.direct_xyz[atom][1..3] and self.real_lattice_inv; writes
        the result into self.direct_abc[atom][1..3].

        The computation proceeds in two steps:

          Step 1 — project onto lattice axes (fractional ABC):
              Pabc_frac = Pxyz · L⁻¹
          where Pxyz is the Cartesian position row-vector and L⁻¹ is the
          real-lattice inverse matrix (real_lattice_inv[xyz_axis][abc_axis]).
          After this step, direct_abc[atom][abc_axis] holds the fractional
          coordinate along that axis (range [0, 1) for atoms inside the cell).

          Step 2 — scale by lattice magnitude:
              Pabc = Pabc_frac · |a_axis|
          Multiplying by self.mag[abc_axis] converts fractional to direct-space
          (Ångström) coordinates along each lattice vector.

        The intermediate fractional ABC values are not stored separately; they
        exist only transiently inside the loop before being scaled in place.

        Parameters
        ----------
        atom : int
            Atom index (1-indexed).
        """
        for abc_axis in range(1, 4):
            self.direct_abc[atom][abc_axis] = 0.0
            for xyz_axis in range(1, 4):
                self.direct_abc[atom][abc_axis] += (self.direct_xyz[atom][xyz_axis] *
                                                    self.real_lattice_inv[xyz_axis][abc_axis])
            # At this point direct_abc holds the fractional coordinate; scale
            # by the lattice vector magnitude to get direct-space ABC (Å).
            self.direct_abc[atom][abc_axis] *= self.mag[abc_axis]

    def direct_xyz2fract_abc(self, xyz):
        """Convert Cartesian [x,y,z] to fractional [a,b,c].

        The transformation is a matrix-vector product:

            Pxyz = atom position in xyz orthogonal (Cartesian) coordinates.
            Pabc = atom position in abc fractional coordinates.
            L-1  = Real lattice matrix inverse.
            Pabc = Pxyz · L-1

        The inverse lattice matrix is stored as real_lattice_inv[xyz_axis][abc_axis]
        (note: this is the *transpose* of the real_lattice storage convention, which
        is real_lattice[abc_axis][xyz_axis]).  The dot product therefore sums over
        the xyz axis in the inner loop, accumulating into each abc component.

        Parameters
        ----------
        xyz : sequence of float
            1-indexed Cartesian coordinates in Angstroms:
            ``xyz[1]`` = x, ``xyz[2]`` = y, ``xyz[3]`` = z; slot 0 is
            an unused sentinel matching Perl's ``@Pxyz[1..3]``.

        Returns
        -------
        list of float
            1-indexed fractional (dimensionless) coordinates
            ``[None, a, b, c]`` so callers read ``result[1..3]``
            directly.
        """
        abc = [None, 0.0, 0.0, 0.0]
        for abc_axis in range(1, 4):
            for xyz_axis in range(1, 4):
                abc[abc_axis] += (xyz[xyz_axis] *
                                  self.real_lattice_inv[xyz_axis][abc_axis])
        return abc

    def direct_xyz2direct_abc(self, xyz):
        """Convert Cartesian [x,y,z] to direct-space abc coordinates.

        Two-step operation:
          1. Compute fractional ABC via the inverse-lattice matrix product
             (Pabc = Pxyz · L-1), which gives dimensionless coordinates in [0, 1).
          2. Scale each fractional component by the corresponding lattice vector
             magnitude (self.mag[1..3]) to obtain direct-space ABC in Angstroms.

        "Direct-space ABC" means the scalar projection of a position along each
        lattice vector direction, measured in Angstroms.  In contrast, fractional
        ABC coordinates are dimensionless (0.0 = origin, 1.0 = one full lattice
        repeat).  For orthogonal lattices the two are simply related by the vector
        lengths; for non-orthogonal lattices the L-1 matrix handles the angular
        factors.

        See direct_xyz2fract_abc for the full derivation of the L-1 matrix product.
        Perl equivalent: directXYZ2directABC.

        Parameters
        ----------
        xyz : sequence of float
            1-indexed Cartesian coordinates in Angstroms
            ``[None, x, y, z]``, matching the convention used by every
            other coordinate routine in this module.

        Returns
        -------
        list of float
            1-indexed direct-space abc coordinates in Angstroms
            ``[None, a, b, c]`` (projected along lattice vectors, NOT
            fractional).
        """
        abc = self.direct_xyz2fract_abc(xyz)
        # Scale fractional coords by lattice magnitudes to get direct-space ABC (Å).
        for abc_axis in range(1, 4):
            abc[abc_axis] *= self.mag[abc_axis]
        return abc

    def fract_abc2direct_xyz(self, fract):
        """Convert fractional ABC coordinates to orthogonal Cartesian XYZ.

        The transformation expresses a position given in fractional lattice
        coordinates (Pa, Pb, Pc) as an orthogonal Cartesian position (Px, Py,
        Pz) using the real-space lattice vectors:

            Px = Pa*ax + Pb*bx + Pc*cx
            Py = Pa*ay + Pb*by + Pc*cy
            Pz = Pa*az + Pb*bz + Pc*cz

        where ax/ay/az, bx/by/bz, cx/cy/cz are the xyz-component coefficients
        of the a, b, c lattice vectors stored in
        self.real_lattice[abc_axis][xyz_axis].

        This is the inverse of directXYZ2fractABC.  The actual matrix multiply
        is delegated to get_direct_xyz_point, which accepts and returns
        1-indexed sequences with slot 0 used as a sentinel.

        Parameters
        ----------
        fract : sequence of float
            1-indexed fractional coordinates ``[None, a, b, c]``.

        Returns
        -------
        list of float
            1-indexed Cartesian coordinates in Angstroms
            ``[None, x, y, z]``.
        """
        return self.get_direct_xyz_point(fract)

    def get_fract_abc(self, atom):
        """Compute fractional abc coordinates from direct abc for a single atom.

        Uses the existing direct-space ABC coordinates and the lattice vector
        magnitudes to produce fractional ABC coordinates:

            fract_abc[axis] = direct_abc[axis] / mag[axis]

        fract_abc[atom][0] is left as None (the Python equivalent of Perl's ""
        initialisation) as a placeholder that simplifies printing loops elsewhere.

        Parameters
        ----------
        atom : int
            Atom index (1-indexed).
        """
        # Use the existing direct-space ABC coordinates and the vector magnitudes
        # to get the fractional ABC coordinates.
        for axis in range(1, 4):
            self.fract_abc[atom][axis] = self.direct_abc[atom][axis] / self.mag[axis]

    def reorder_lattice_parameters(self):
        """Sort lattice vectors into a canonical ordering.

        The goal is a standardised form where the lattice vectors are ordered
        primarily by their inter-axial angles (smallest to largest) and
        secondarily, within any group of equal angles, by their vector
        magnitudes (smallest to largest).

        The conventional angle-to-axis relationship is:
          alpha  = angle between the b and c axes  →  axis a follows alpha's index
          beta   = angle between the a and c axes  →  axis b follows beta's index
          gamma  = angle between the a and b axes  →  axis c follows gamma's index
        So when the angles are permuted, the corresponding cell vector
        magnitudes must move with them to maintain that correspondence.

        After the initial angle-based sort a second pass handles tie-breaking:
        if all three angles are equal (e.g. cubic or rhombohedral), the three
        axes are re-sorted purely by magnitude.  If only a pair of angles is
        equal, up to three conditional swaps are applied to put the smaller
        magnitude first within each equal-angle pair.

        The Perl implementation contained a subtle indexing bug: after sorting
        @sortedLatIndices via a slice assignment ("@angle = @angle[@sortedLatIndices]")
        the array was collapsed from 1-indexed (0..3) to 0-indexed (0..2), making
        $angle[3] and $mag[3] undefined (→ 0) in every subsequent comparison.
        This made the "all angles equal" branch unreachable and corrupted the
        pair-swap comparisons in the else branch.  The Python port avoids the
        issue by working with a plain 0-indexed list `indices` (length 3) while
        keeping self.angle/mag 1-indexed throughout; `sorted_lat_indices` is
        stored as [None] + indices for callers that expect 1-indexed access.

        Postcondition: self.angle, self.angle_deg, self.mag, and
        self.real_lattice are all permuted consistently; self.sorted_lat_indices
        records the permutation for reorder_atom_coordinates().
        """
        # Sort axes 1..3 by angle ascending (Python sort is stable).
        indices = sorted([1, 2, 3], key=lambda i: self.angle[i])

        # Within equal-angle groups, also sort by magnitude.
        def swap_if_needed(idx, i, j):
            if (abs(self.angle[idx[i]] - self.angle[idx[j]]) < EPSILON and
                    self.mag[idx[i]] > self.mag[idx[j]]):
                idx[i], idx[j] = idx[j], idx[i]

        if (abs(self.angle[indices[0]] - self.angle[indices[1]]) < EPSILON and
                abs(self.angle[indices[1]] - self.angle[indices[2]]) < EPSILON):
            # All three angles equal: sort purely by magnitude.
            indices.sort(key=lambda i: self.mag[i])
        else:
            # Check all three possible equal-angle pairs and swap if the
            # smaller-magnitude axis does not come first.
            swap_if_needed(indices, 0, 2)
            swap_if_needed(indices, 0, 1)
            swap_if_needed(indices, 1, 2)

        # Store the permutation (1-indexed) for use by reorder_atom_coordinates.
        self.sorted_lat_indices = [None] + indices

        # Apply permutation to angle, angle_deg, mag.
        # The cell vector magnitudes must follow the same permutation as the
        # angles because alpha↔a, beta↔b, gamma↔c by convention.
        old_angle     = self.angle[:]
        old_angle_deg = self.angle_deg[:]
        old_mag       = self.mag[:]
        for new_axis, old_axis in enumerate(indices, start=1):
            self.angle[new_axis]     = old_angle[old_axis]
            self.angle_deg[new_axis] = old_angle_deg[old_axis]
            self.mag[new_axis]       = old_mag[old_axis]

        # Apply the same permutation to the (abc, xyz) lattice vector rows.
        old_lattice = [None] + [self.real_lattice[i][:] for i in range(1, 4)]
        for new_axis, old_axis in enumerate(indices, start=1):
            for xyz in range(1, 4):
                self.real_lattice[new_axis][xyz] = old_lattice[old_axis][xyz]

    def reorder_atom_coordinates(self):
        """Permute atom coordinates to match the reordered lattice vectors.

        Must be called after reorder_lattice_parameters() once the permutation
        stored in self.sorted_lat_indices is known.  Applies that same
        permutation to each atom's fractional (fract_abc) coordinates, then
        propagates the change to the direct-space abc (direct_abc) and
        Cartesian xyz (direct_xyz) representations.

        The Perl version (reorderAtomCoordinates) was called once per atom from
        the read_olcao_skl call chain.  The Python port loops over all atoms
        internally for simplicity.

        Coordinate update chain per atom:
          fract_abc  (permuted in place)
          → direct_abc  =  fract_abc * mag   (element-wise along abc axes)
          → direct_xyz  =  get_direct_xyz()  (matrix multiply with real_lattice)
        """
        for atom in range(1, self.num_atoms + 1):
            # Copy current fractional coordinates before overwriting.
            old_fract = self.fract_abc[atom][:]
            # Apply the lattice permutation to the fractional coordinates.
            for new_axis, old_axis in enumerate(self.sorted_lat_indices[1:], start=1):
                self.fract_abc[atom][new_axis] = old_fract[old_axis]
            # Propagate to direct-space abc coordinates.
            for axis in range(1, 4):
                self.direct_abc[atom][axis] = self.fract_abc[atom][axis] * self.mag[axis]
            # Propagate to Cartesian xyz coordinates.
            self.get_direct_xyz(atom)

    # ======================================================================
    # Section: Symmetry operations
    # ======================================================================

    def apply_space_group(self):
        """Generate all symmetry-equivalent atoms using the space group.

        Uses self.space_group (a name or number string) to locate the
        corresponding file in the space group database directory
        (self.space_group_db).  The DB file format is:
          line 1:   <lattice_type>  <sg_name>  [alt_letter]
          line 2:   <sg_num>  <sg_sub_num>
          lines 3+: symmetry operation records (passed verbatim to the binary)

        lattice_type indicates whether the space group uses a primitive ('P')
        or a non-primitive (e.g. 'F', 'I', 'A', ...) conventional cell.

        Some space groups have alternative definitions (different origin
        choices, axis settings, etc.).  When a file carries an alternative,
        the third token on line 1 is a lowercase letter ('a', 'b', ...) and
        the canonical name stored in self.space_group_name is extended with
        that letter, e.g. space group 227 becomes 'Fd3~m_a'.

        From these inputs this method assembles an 'sginput' control file,
        invokes the external applySpaceGroup Fortran binary, then reads the
        expanded structure back from 'sgoutput'.

        sginput layout (in order):
          lattice_type / sg_num sg_sub_num / symmetry ops / do_full_cell flag /
          3 lattice vectors / num_atoms / atom lines (eid sid fa fb fc) /
          do_not_cell_shift flag

        do_full_cell (1/0): request the full, possibly non-primitive,
          conventional cell rather than the primitive cell.
        do_not_cell_shift (1/0): suppress the automatic folding of all
          atoms back inside the unit cell after symmetry expansion.

        After reading sgoutput the per-atom arrays (atom_tag,
        atom_element_name, atom_element_id, atom_species_id, fract_abc,
        direct_xyz, direct_abc) are completely reset and rebuilt.
        Note: element_list and species_list are intentionally left untouched —
        space group operations never change which elements or species are
        present, only how many atoms of each there are.
        """
        # Use the space group ID (name or number) to open the space group
        # operations file from the space group database.
        sg_path = os.path.join(self.space_group_db, self.space_group)
        with open(sg_path) as fh:
            lines = fh.readlines()

        # Read the first line to get the cell lattice type (primitive or non-
        # primitive) and the space group name.
        values = lines[0].split()
        self.lattice_type    = values[0]
        self.space_group_name = values[1]

        # If there is an alternative definition for this space group it is
        # labeled in the space group file with an index starting at 'a'.  In
        # such a case the file name is extended by '_a' for the complete name
        # of the space group (i.e. not the number) (e.g. 227 = Fd3~m_a).
        if len(values) > 2 and values[2].islower():
            self.space_group_name = self.space_group_name + '_' + values[2]

        # Read the second line to get the root space group number.
        values = lines[1].split()
        self.space_group_num     = int(values[0])
        self.space_group_sub_num = int(values[1])

        # Read remaining lines and prepare input for the applySpaceGroup
        # Fortran program.  The symmetry op records from the DB file are
        # passed verbatim; the header (lattice_type / sg_num sg_sub_num) is
        # prepended, and the structure data are appended afterwards.
        sginput = []
        sginput.append(f'{self.lattice_type}\n')
        sginput.append(f'{self.space_group_num}  {self.space_group_sub_num}\n')
        sginput.extend(lines[2:])

        # Add the flag indicating whether to produce a full (possibly
        # non-primitive) conventional cell (1) or a primitive cell (0).
        sginput.append('1\n' if self.do_full_cell else '0\n')

        # Add lattice vectors to the sginput array (one row per a/b/c axis).
        for axis in range(1, 4):
            rl = self.real_lattice[axis]
            sginput.append(f'  {rl[1]}  {rl[2]}  {rl[3]}\n')

        # Add atomic data: total atom count then one line per atom with
        # element ID, species ID, and fractional abc coordinates.
        sginput.append(f'{self.num_atoms}\n')
        for atom in range(1, self.num_atoms + 1):
            fa  = self.fract_abc[atom]
            eid = self.atom_element_id[atom]
            sid = self.atom_species_id[atom]
            sginput.append(f'{eid} {sid}  {fa[1]}  {fa[2]}  {fa[3]}\n')

        # Add a flag indicating whether to prevent atoms from being shifted
        # back inside the unit cell (1 = do not shift, 0 = shift normally).
        sginput.append('1\n' if self.do_not_cell_shift else '0\n')

        # Write the assembled control file for the applySpaceGroup program.
        with open('sginput', 'w') as fh:
            fh.writelines(sginput)

        # Apply the space group to this structure.
        subprocess.run(['applySpaceGroup'], check=True)

        # Open the output file and process the expanded structure.
        with open('sgoutput') as fh:
            # Reassign the lattice parameters and propagate to other
            # representations (magnitudes, angles).
            for axis in range(1, 4):
                vals = fh.readline().split()
                self.real_lattice[axis][1] = float(vals[0])
                self.real_lattice[axis][2] = float(vals[1])
                self.real_lattice[axis][3] = float(vals[2])
            self.abc_alpha_beta_gamma()

            # Read the number of atoms in the newly expanded system.
            self.num_atoms = int(fh.readline().strip())

            # Reset all per-atom arrays — they are about to be rebuilt from
            # scratch.  element_list and species_list are intentionally not
            # reset here because space group expansion does not change which
            # elements or species are present in the model.
            self.atom_tag          = [None]
            self.atom_element_name = [None]
            self.atom_element_id   = [None]
            self.atom_species_id   = [None]
            self.fract_abc         = [None]
            self.direct_xyz        = [None]
            self.direct_abc        = [None]

            # Read each atom's element ID and species ID numbers and rebuild
            # the per-atom arrays.  atom_element_name and atom_tag are
            # recovered from the still-valid element_list and species_list.
            for atom in range(1, self.num_atoms + 1):
                vals = fh.readline().split()
                eid = int(vals[0])
                sid = int(vals[1])
                fa1, fa2, fa3 = float(vals[2]), float(vals[3]), float(vals[4])
                self.atom_element_id.append(eid)
                self.atom_species_id.append(sid)
                self.fract_abc.append([None, fa1, fa2, fa3])
                self.atom_element_name.append(self.element_list[eid])
                self.atom_tag.append(self.species_list[eid][sid])
                self.direct_xyz.append([None, 0.0, 0.0, 0.0])
                self.direct_abc.append([None, 0.0, 0.0, 0.0])

        # Propagate fractional coordinates to direct-space representations.
        for atom in range(1, self.num_atoms + 1):
            self.get_direct_xyz(atom)
            self.get_direct_abc(atom)

    def apply_supercell(self):
        """Replicate the unit cell to form a supercell.

        Reads self.supercell[1..3] (the number of cells along each a/b/c axis)
        and self.sc_mirror[1..3] (per-axis mirror flags) from the .skl
        "supercell" line, e.g.:
            supercell 2 2 1           -> tile 2×2×1, no mirroring
            supercell 2 2 1  1 0 0    -> tile 2×2×1, mirror alternate a-cells

        The optional three trailing integers on the supercell line are the
        sc_mirror flags.  If omitted they default to 0 (no mirroring).

        Algorithm (one axis at a time):
        1. Scale: multiply mag[abc] and real_lattice[abc][xyz] by n for each
           axis, so the new cell is n times larger in that direction.
        2. Replicate: for each axis with n > 1, loop over all current atoms
           and generate n copies.  For copy number `cell` (1-indexed):
             - Along the tiled axis the fractional coordinate is mapped into
               the sub-cell:
                 normal:  new_frac = frac/n + (cell-1)/n
                 mirror:  new_frac = (1-frac)/n + (cell-1)/n  (even cells only)
               The mirror formula reflects the atom about the sub-cell center,
               producing a structure that is a mirror image of its neighbour.
               Mirroring is applied only to even-numbered cells (cell % 2 == 0);
               odd-numbered cells always use the normal formula.
             - Along the other two axes the fractional coordinate is unchanged.
           After each axis the old per-atom arrays (fract_abc, atom_element_name,
           atom_element_id, atom_species_id, atom_tag) are discarded and replaced
           with the freshly built expanded versions (analogous to Perl's undef
           followed by array reassignment), so that the next axis loop always
           sees the already-tiled, up-to-date atom list.
        3. Recompute: after all axes are processed, re-derive the inverse and
           reciprocal lattice (make_inv_or_recip_lattice, formerly makeLatticeInv
           in older versions of the Perl code), the a/b/c magnitudes and angles
           for both the real and reciprocal lattice (abc_alpha_beta_gamma), and
           the Cartesian (direct_xyz) and direct_abc coordinates of every atom.
           Because num_atoms has grown, direct_xyz and direct_abc must be
           explicitly re-allocated before the per-atom coordinate calls (Python
           fixed-length lists, unlike Perl's auto-extending arrays).
        """
        # Compute the new lattice parameters in both a,b,c,alpha,beta,gamma and
        # [a,b,c][x,y,z] form.  Increase the dimensions of the model along each
        # axis by a factor equal to the number of cells requested in that
        # direction.
        for abc_axis in range(1, 4):
            n = self.supercell[abc_axis]
            self.mag[abc_axis] *= n
            # Increase each x,y,z component of each a,b,c vector in the same way.
            for xyz_axis in range(1, 4):
                self.real_lattice[abc_axis][xyz_axis] *= n

        # Apply the supercell requests along each axis to each atom.
        for abc_axis in range(1, 4):
            n = self.supercell[abc_axis]
            # There is no need to replicate atoms if the supercell request in
            # this direction is equal to one.
            if n == 1:
                continue
            mirror = self.sc_mirror[abc_axis]
            new_num_atoms        = 0
            new_fract_abc        = [None]
            new_atom_element_name = [None]
            new_atom_element_id  = [None]
            new_atom_species_id  = [None]
            new_atom_tag         = [None]

            # Loop over each atom and replicate them to fill the new cell.
            for atom in range(1, self.num_atoms + 1):
                for cell in range(1, n + 1):
                    new_num_atoms += 1
                    new_fa = [None, 0.0, 0.0, 0.0]
                    for atom_axis in range(1, 4):
                        if atom_axis == abc_axis:
                            # Normal: map frac into the cell-th sub-cell.
                            # Mirror: reflect about sub-cell center for even cells.
                            if not mirror or (cell % 2 == 1):
                                new_fa[atom_axis] = (
                                    self.fract_abc[atom][atom_axis] / n
                                    + (cell - 1.0) / n
                                )
                            else:
                                new_fa[atom_axis] = (
                                    (1.0 - self.fract_abc[atom][atom_axis]) / n
                                    + (cell - 1.0) / n
                                )
                        else:
                            # Non-tiled axes: fractional coordinate is unchanged.
                            new_fa[atom_axis] = self.fract_abc[atom][atom_axis]
                    new_fract_abc.append(new_fa)
                    new_atom_element_name.append(self.atom_element_name[atom])
                    new_atom_element_id.append(self.atom_element_id[atom])
                    new_atom_species_id.append(self.atom_species_id[atom])
                    new_atom_tag.append(self.atom_tag[atom])

            # Update the system information necessary for the next abc_axis.
            # Discard the old per-atom arrays and install the expanded ones so
            # the next axis loop operates on the already-tiled structure.
            # (In Perl: undef @array followed by @array = @newArray.)
            self.num_atoms         = new_num_atoms
            self.fract_abc         = new_fract_abc
            self.atom_element_name = new_atom_element_name
            self.atom_element_id   = new_atom_element_id
            self.atom_species_id   = new_atom_species_id
            self.atom_tag          = new_atom_tag

        # Reobtain the inverse lattice and the reciprocal lattice vectors.
        # (Perl formerly called makeLatticeInv twice; that routine was renamed
        # makeInvOrRecipLattice.  The flag 0 = inverse, 1 = reciprocal.)
        self.make_inv_or_recip_lattice(make_recip=True)
        # Recompute a/b/c magnitudes and angles for both the real and reciprocal
        # lattice so that mag, angle, mag_recip, angle_recip, etc. are
        # consistent with the newly scaled real_lattice vectors.
        self.abc_alpha_beta_gamma()

        # Obtain the direct_abc and direct_xyz coordinates of the atoms in the
        # model such that they match the fractional_abc.
        # Pre-allocate to the new (larger) atom count before filling; Python
        # lists do not auto-extend on indexed assignment the way Perl arrays do.
        self.direct_xyz = [None] + [[None, 0.0, 0.0, 0.0]
                                    for _ in range(self.num_atoms)]
        self.direct_abc = [None] + [[None, 0.0, 0.0, 0.0]
                                    for _ in range(self.num_atoms)]
        for atom in range(1, self.num_atoms + 1):
            self.get_direct_xyz(atom)
            self.get_direct_abc(atom)

    def prep_surface(self, hkl):
        """Orient and cut the cell to expose the (hkl) surface plane.

        We are given a reciprocal-space Miller index vector hkl = (h,k,l).
        This defines a family of planes in real space that intercepts the
        primitive cell at a/h, b/k, c/l.  An hkl value of zero means the
        real-space plane is parallel to that axis and does not intersect it
        (e.g. hkl=010 → plane parallel to the xz plane, no intersection
        with x or z).

        The algorithm proceeds in five stages:

        1.  Build the real-space normal vector (uvwNormal):
            Compute the least common multiple (LCM) of the non-zero hkl
            values.  Use it to construct two helper lattices, each
            [abc 1..3][xyz 1..3]:
              uvwNormalLattice — for non-zero hkl[i]:
                uvwNormalLattice[i] = (LCM/hkl[i]) * real_lattice[i]
              For zero hkl[i]: uvwNormalLattice[i] = 0 (plane parallel).
              uvwRealLattice — identical to uvwNormalLattice, except that
                hkl=0 rows use the original real_lattice vector instead of 0.
                This defines the actual replicated lattice for the search.
            uvwNormal = sum of all three uvwNormalLattice rows.  This is the
            real-space vector R = ua_1 + va_2 + wa_3 that is normal to the
            plane K = hb_1 + kb_2 + lb_3 (b_i are reciprocal lattice
            vectors).  It will become the new a-axis of the surface cell.

        2.  Find the minimal in-plane lattice vectors (new b and c axes):
            Execute a spiral search over integer multiples of uvwRealLattice
            (same algorithm as makeLattice in lattice.f90; note that here
            uvwRealLattice has [abc][xyz] ordering while the Fortran uses
            [xyz][abc]).  Search limits are maxRep=5 for hkl=0 axes (lattice
            unchanged in that direction) and maxRep=100 for non-zero axes.
            Sort resulting lattice points by distance from origin.  Scan in
            order of increasing distance for:
              - First point perpendicular to uvwNormal → new b-axis.
              - Next point perpendicular to uvwNormal AND non-colinear with
                b → new c-axis.
            (When one hkl=0 this is the easy case: the original lattice
            vector along that axis becomes the new a-axis unchanged, possibly
            after rotation.  See e.g. hkl=110 where old z → new |c| axis.)

        3.  Install the new lattice early:
            Store new a/b/c axes in real_lattice and recompute derived
            quantities (inverse, reciprocal, angles) before the atom-filling
            loop, because direct_xyz2fract_abc (used for the on-face test)
            requires real_lattice_inv.  Note: the code still references
            new_real_lattice directly in many places; real_lattice is kept
            in sync.
            For the inside-cell test build:
              - Center of the new parallelpiped cell.
              - Center of each of its six faces (three pairs of opposite faces
                sharing a common pair of lattice vectors).
              - Outward-pointing unit normals for each face (three primary
                normals via cross-product from the origin, three opposites by
                negation; then check each against the face-to-center vector:
                a positive dot product means the normal still points inward,
                so flip it).

        4.  Fill the new cell with atoms from a super-lattice:
            Iterate over all replicated cells (cells with lattice-point
            distance > diagonal_mag are skipped as a conservative cull).
            For each atom in each cell compute the trial Cartesian position
            and apply three filters:
              a. Outside test — if the dot product of the (face-center → atom)
                 vector with any outward face normal is positive the atom is
                 outside the cell; skip it.
              b. On-face test — if any fractional coordinate has |frac| ≈ 1
                 the atom sits on a cell face and is a periodic image of the
                 atom at frac=0 on the opposite face; skip to prevent
                 duplicates.
              c. Duplicate test — exact Cartesian match with an already-kept
                 atom; skip.

        5.  Rotate to canonical orientation:
            Rotation 1 — align uvwNormal (new a-axis) with Cartesian x-axis.
              Rotation axis = cross-product of x-axis and uvwNormal (i.e. the
              normal to the plane they span).  If uvwNormal has a positive
              y-component, invert the rotation axis so the rotation sweeps
              *back* toward x rather than away from it.  If uvwNormal is
              already aligned with x (normalizer=0), no rotation is needed.
            Rotation 2 — bring the b-axis into the Cartesian xy-plane.
              Rotation axis = x-axis (the now-fixed a-axis direction).
              Compute the rotation angle from the yz-projection of the b-axis,
              NOT b itself: using b directly would give a bogus non-zero angle
              if b were already in the xy-plane but not aligned with y.
              If the yz-projection has a positive z-component, invert the
              rotation axis so the sweep goes back toward y.
            After both rotations recompute all derived quantities and
            directABC/fractABC coordinates for every atom.

        Parameters
        ----------
        hkl : list
            1-indexed Miller index list [None, h, k, l].
        """
        # Pre-compute the diagonal of the current cell (used later as
        # a culling distance for replicated cells).  The diagonal is a
        # 1-indexed vector [None, x, y, z] to match the surrounding
        # coordinate convention.
        diagonal = [None, 0.0, 0.0, 0.0]
        for abc_axis in range(1, 4):
            for xyz_axis in range(1, 4):
                diagonal[xyz_axis] += self.real_lattice[abc_axis][xyz_axis]
        diagonal_mag = math.sqrt(diagonal[1]**2 +
                                 diagonal[2]**2 +
                                 diagonal[3]**2)

        # Compute the least common multiple of the non-zero hkl values.
        # An hkl value of zero means the real-space plane will not intersect
        # that axis, so we exclude it from the LCM calculation.  The upper
        # bound for the search is the product of all non-zero hkl values
        # (guaranteed to be divisible by each of them).
        product = 1
        for hkl_idx in range(1, 4):
            if hkl[hkl_idx] != 0:
                product *= hkl[hkl_idx]
        least_common_mult = 0
        for i in range(1, product + 1):
            found = True
            for hkl_idx in range(1, 4):
                if hkl[hkl_idx] != 0 and i % hkl[hkl_idx] != 0:
                    found = False
                    break
            if found:
                least_common_mult = i
                break

        # Build uvwNormalLattice and uvwRealLattice [abc 1..3][xyz 1..3].
        # uvwNormalLattice rows sum to the real-space vector normal to the hkl
        # plane: R = (LCM/h)*a1 + (LCM/k)*a2 + (LCM/l)*a3.  For a zero hkl
        # index the row is zeroed (that axis is parallel to the plane).
        # uvwRealLattice is identical except for hkl=0 axes: those rows use
        # the original real_lattice vector so that the resulting lattice has
        # valid (non-zero) periodicity in every direction and can be used for
        # the spiral search that finds the in-plane b and c axes.
        uvw_normal_lattice = [[None] * 4 for _ in range(4)]
        uvw_real_lattice   = [[None] * 4 for _ in range(4)]
        for abc_axis in range(1, 4):
            for xyz_axis in range(1, 4):
                if hkl[abc_axis] == 0:
                    uvw_normal_lattice[abc_axis][xyz_axis] = 0.0
                    uvw_real_lattice[abc_axis][xyz_axis]   = \
                            self.real_lattice[abc_axis][xyz_axis]
                else:
                    val = (least_common_mult *
                           self.real_lattice[abc_axis][xyz_axis] /
                           hkl[abc_axis])
                    uvw_normal_lattice[abc_axis][xyz_axis] = val
                    uvw_real_lattice[abc_axis][xyz_axis]   = val

        # uvwNormal = accumulated sum of the three uvwNormalLattice rows in xyz.
        # This is the real-space lattice vector that is normal to the plane
        # defined by hkl and that points to an actual lattice site (its
        # components are integral multiples of the real-lattice vectors).
        # It will become the new a-axis of the surface cell.  Stored
        # 1-indexed [None, x, y, z] to feed the 1-indexed vector helpers.
        uvw_normal = [None, 0.0, 0.0, 0.0]
        for abc_axis in range(1, 4):
            for xyz_axis in range(1, 4):
                uvw_normal[xyz_axis] += uvw_normal_lattice[abc_axis][xyz_axis]
        uvw_mag = math.sqrt(uvw_normal[1]**2 +
                            uvw_normal[2]**2 +
                            uvw_normal[3]**2)

        # Set spiral search limits per axis.  For a zero hkl index the lattice
        # vector in that direction is unchanged, so only a small neighbourhood
        # (maxRep=5) needs to be searched.  For non-zero indices we may need to
        # travel far to find in-plane lattice points, so maxRep=100 is used
        # (assumed sufficient for all practical hkl values).
        max_rep = [None, 0, 0, 0]
        for i in range(1, 4):
            max_rep[i] = 5 if hkl[i] == 0 else 100

        # Enumerate all lattice points of the uvwRealLattice within the spiral
        # search box.  Each point is the (distance-from-origin, 1-indexed
        # [None, x, y, z]) pair so it can be fed directly to dot_product
        # and cross_product_mag without any offset translation.
        # The loop centers each axis on zero by subtracting half the max range.
        cell_data = []
        for i in range(max_rep[1] + 2):
            rep1 = i - (1 + max_rep[1]) // 2
            for j in range(max_rep[2] + 2):
                rep2 = j - (1 + max_rep[2]) // 2
                for k in range(max_rep[3] + 2):
                    rep3 = k - (1 + max_rep[3]) // 2
                    reps = [None, rep1, rep2, rep3]
                    lp = [None, 0.0, 0.0, 0.0]
                    for l in range(1, 4):
                        for m in range(1, 4):
                            lp[m] += uvw_real_lattice[l][m] * reps[l]
                    dist = math.sqrt(lp[1]**2 + lp[2]**2 + lp[3]**2)
                    cell_data.append((dist, lp))

        # Sort lattice points by distance from the origin so we can scan in
        # order of increasing distance and find the *shortest* in-plane vectors.
        cell_data.sort(key=lambda x: x[0])
        num_cells   = len(cell_data)
        cell_dist   = [None] + [cd[0] for cd in cell_data]  # 1-indexed
        cell_points = [None] + [cd[1] for cd in cell_data]  # 1-indexed

        # Find the new b and c lattice vectors by scanning the distance-sorted
        # lattice points for vectors that lie in the plane perpendicular to
        # uvwNormal (dot product ≈ 0) and that are not co-linear with each
        # other (cross-product magnitude > epsilon).  Cell index 1 is always
        # (0,0,0) and is skipped.  The first qualifying point → b-axis
        # (new_real_lattice[2]); the next qualifying, non-colinear point →
        # c-axis (new_real_lattice[3]).
        new_mag          = [None, 0.0, 0.0, 0.0]
        new_real_lattice = [[None] * 4 for _ in range(4)]
        for cell in range(2, num_cells + 1):  # skip cell 1 = (0,0,0)
            tv = cell_points[cell]  # already [None, x, y, z]
            if abs(self.dot_product(tv, uvw_normal)) < EPSILON:
                if new_mag[2] == 0.0:
                    new_mag[2] = cell_dist[cell]
                    for axis in range(1, 4):
                        new_real_lattice[2][axis] = tv[axis]
                else:
                    rl2 = new_real_lattice[2]  # already 1-indexed row
                    if abs(self.cross_product_mag(tv, rl2)) > EPSILON:
                        new_mag[3] = cell_dist[cell]
                        for axis in range(1, 4):
                            new_real_lattice[3][axis] = tv[axis]
                        break

        # The a-axis is the surface normal (uvwNormal).
        new_mag[1] = uvw_mag
        for axis in range(1, 4):
            new_real_lattice[1][axis] = uvw_normal[axis]

        # Install the new lattice into self and recompute derived quantities now,
        # even though much of the code below still references new_real_lattice
        # directly.  We must do this here because direct_xyz2fract_abc (used
        # for the on-face duplicate check) requires real_lattice_inv, which is
        # built from real_lattice.  The angles will need to be recomputed again
        # after the rotation steps below.
        for abc_axis in range(1, 4):
            for xyz_axis in range(1, 4):
                self.real_lattice[abc_axis][xyz_axis] = \
                        new_real_lattice[abc_axis][xyz_axis]
            self.mag[abc_axis] = new_mag[abc_axis]
        self.make_inv_or_recip_lattice(make_recip=True)
        self.abc_alpha_beta_gamma()

        # --- Geometry for the inside-cell test ---
        # To determine whether an atom lies inside the new parallelpiped we:
        #   1. Compute the center of the cell (half-sum of the three lattice
        #      vectors).
        #   2. Compute the center of each of the six faces.  Faces come in
        #      three opposite pairs: (1,2) span axes 1&2; (3,4) span axes 1&3;
        #      (5,6) span axes 2&3.  Face 2 of each pair = face 1 + third vec.
        #   3. Build outward-pointing unit normals (see below).
        #   4. For each candidate atom, compute the vector from each face center
        #      to the atom.  If the dot product with any outward normal is
        #      positive (> epsilon) the atom is outside.
        # Every geometric vector in this block is 1-indexed
        # [None, x, y, z] to match the helper interface.
        new_cell_center = [None] + [
            (new_real_lattice[1][axis] +
             new_real_lattice[2][axis] +
             new_real_lattice[3][axis]) / 2.0
            for axis in range(1, 4)
        ]

        # Face centers (1-indexed list of 1-indexed [None, x, y, z] vectors).
        face_center    = [None] * 7
        face_center[1] = [None] + [
                (new_real_lattice[1][axis] + new_real_lattice[2][axis]) / 2.0
                for axis in range(1, 4)]
        face_center[2] = [None] + [
                face_center[1][axis] + new_real_lattice[3][axis]
                for axis in range(1, 4)]
        face_center[3] = [None] + [
                (new_real_lattice[1][axis] + new_real_lattice[3][axis]) / 2.0
                for axis in range(1, 4)]
        face_center[4] = [None] + [
                face_center[3][axis] + new_real_lattice[2][axis]
                for axis in range(1, 4)]
        face_center[5] = [None] + [
                (new_real_lattice[2][axis] + new_real_lattice[3][axis]) / 2.0
                for axis in range(1, 4)]
        face_center[6] = [None] + [
                face_center[5][axis] + new_real_lattice[1][axis]
                for axis in range(1, 4)]

        # Face normals: compute via cross-product from the origin for three
        # unique faces, negate for the three opposite faces.  Each normal is
        # then flipped if necessary so it points *outward* (away from the cell
        # center): build the vector from that face's center to the cell center;
        # if its dot product with the face normal is positive, both vectors
        # point in the same inward direction and the normal must be negated.
        # Finally normalise to unit length.
        origin_v = [None, 0.0, 0.0, 0.0]
        # Each new_real_lattice row is itself already 1-indexed, so they
        # can be passed directly as [None, x, y, z] vectors.
        rl1 = new_real_lattice[1]
        rl2 = new_real_lattice[2]
        rl3 = new_real_lattice[3]
        face_normal    = [None] * 7
        face_normal[1] = self.get_plane_normal(origin_v, rl1, rl2)
        face_normal[3] = self.get_plane_normal(origin_v, rl1, rl3)
        face_normal[5] = self.get_plane_normal(origin_v, rl2, rl3)
        face_normal[2] = [None] + [-face_normal[1][axis] for axis in range(1, 4)]
        face_normal[4] = [None] + [-face_normal[3][axis] for axis in range(1, 4)]
        face_normal[6] = [None] + [-face_normal[5][axis] for axis in range(1, 4)]

        # Ensure face normals point outward, then normalize to unit length.
        for face in range(1, 7):
            f2c = [None] + [new_cell_center[axis] - face_center[face][axis]
                            for axis in range(1, 4)]
            if self.dot_product(face_normal[face], f2c) > 0.0:
                face_normal[face] = [None] + [-face_normal[face][axis]
                                              for axis in range(1, 4)]
            norm = math.sqrt(face_normal[face][1]**2 +
                             face_normal[face][2]**2 +
                             face_normal[face][3]**2)
            face_normal[face] = [None] + [face_normal[face][axis] / norm
                                          for axis in range(1, 4)]

        # --- Fill the new cell with atoms from replicated original cells ---
        # Iterate over the distance-sorted list of super-lattice cells.  Cells
        # whose lattice-point distance from the origin exceeds diagonal_mag
        # (the body-diagonal of the original cell) cannot possibly contribute
        # atoms inside the new cell and are skipped.
        new_num_atoms         = 0
        new_direct_xyz        = [None]
        new_atom_element_name = [None]
        new_atom_element_id   = [None]
        new_atom_species_id   = [None]
        new_atom_tag          = [None]
        for cell in range(1, num_cells + 1):
            if cell_dist[cell] > diagonal_mag:
                continue
            cp = cell_points[cell]  # already [None, x, y, z]
            for atom in range(1, self.num_atoms + 1):
                new_xyz = [None,
                           self.direct_xyz[atom][1] + cp[1],
                           self.direct_xyz[atom][2] + cp[2],
                           self.direct_xyz[atom][3] + cp[3]]

                # Outside test: dot product of (face_center → atom) with each
                # outward unit normal.  A positive value (> epsilon) means the
                # atom is on the outer side of that face → outside the cell.
                outside = False
                for face in range(1, 7):
                    f2a = [None] + [new_xyz[axis] - face_center[face][axis]
                                    for axis in range(1, 4)]
                    if self.dot_product(face_normal[face], f2a) > EPSILON:
                        outside = True
                        break
                if outside:
                    continue

                # On-face test: an atom at fractional coordinate ±1 is the
                # periodic image of the atom at frac=0 on the opposite face;
                # keeping it would create a duplicate.  Eliminate it here.
                # direct_xyz2fract_abc now takes/returns [None, a, b, c].
                temp_fract = self.direct_xyz2fract_abc(new_xyz)
                if any(abs(abs(temp_fract[axis]) - 1.0) < EPSILON
                       for axis in range(1, 4)):
                    continue

                # Discard exact duplicates of already-kept atoms.
                same = False
                for prev in new_direct_xyz[1:]:
                    if (abs(prev[1] - new_xyz[1]) < EPSILON and
                            abs(prev[2] - new_xyz[2]) < EPSILON and
                            abs(prev[3] - new_xyz[3]) < EPSILON):
                        same = True
                        break
                if same:
                    continue

                # Keep this atom.
                new_num_atoms += 1
                new_direct_xyz.append(new_xyz)
                new_atom_element_name.append(self.atom_element_name[atom])
                new_atom_element_id.append(self.atom_element_id[atom])
                new_atom_species_id.append(self.atom_species_id[atom])
                new_atom_tag.append(self.atom_tag[atom])

        # Replace per-atom data with the new set.
        self.num_atoms         = new_num_atoms
        self.direct_xyz        = new_direct_xyz
        self.atom_element_name = new_atom_element_name
        self.atom_element_id   = new_atom_element_id
        self.atom_species_id   = new_atom_species_id
        self.atom_tag          = new_atom_tag
        self.direct_abc = [None] + [[None, 0.0, 0.0, 0.0]
                                    for _ in range(self.num_atoms)]
        self.fract_abc  = [None] + [[None, 0.0, 0.0, 0.0]
                                    for _ in range(self.num_atoms)]

        # --- Rotation 1: align uvwNormal (new a-axis) with Cartesian x-axis ---
        # The rotation axis is the normal to the plane spanned by x-axis and
        # uvwNormal (i.e. their cross-product, computed via get_plane_normal).
        # If the two vectors are already co-linear, get_plane_normal returns
        # (0,0,0): normaliser=0, so no rotation is needed and the angle is
        # set to zero (which produces the identity rotation matrix).  All
        # three vectors are stored as [None, x, y, z] to feed the
        # 1-indexed vector helpers directly.
        x_axis   = [None, 1.0, 0.0, 0.0]
        rot_axis = self.get_plane_normal(origin_v, x_axis, uvw_normal)
        norm     = math.sqrt(rot_axis[1]**2 + rot_axis[2]**2 + rot_axis[3]**2)
        if norm != 0.0:
            rot_axis = [None] + [rot_axis[axis] / norm
                                 for axis in range(1, 4)]
            rot_angle_deg = math.degrees(
                    self.get_vector_angle(uvw_normal, x_axis))
        else:
            rot_angle_deg = 0.0

        # If uvwNormal points into the positive-y half-space, the rotation as
        # computed would swing it *away* from x, so invert the rotation axis
        # to sweep it back toward the x-axis instead.
        if uvw_normal[2] > 0.0:
            rot_axis = [None] + [-rot_axis[axis] for axis in range(1, 4)]

        self.define_rot_matrix(rot_axis, rot_angle_deg)

        # Rotate the three lattice vectors in place (each row is already
        # 1-indexed, so the rotate_one_point output can be copied back
        # slot-by-slot without any offset bookkeeping).
        for abc_axis in range(1, 4):
            rpt = self.rotate_one_point(self.real_lattice[abc_axis])
            for xyz_axis in range(1, 4):
                self.real_lattice[abc_axis][xyz_axis] = rpt[xyz_axis]
        for abc_axis in range(1, 4):
            for xyz_axis in range(1, 4):
                if abs(self.real_lattice[abc_axis][xyz_axis]) < EPSILON:
                    self.real_lattice[abc_axis][xyz_axis] = 0.0

        # Rotate all atoms.
        for atom in range(1, self.num_atoms + 1):
            rpt = self.rotate_one_point(self.direct_xyz[atom])
            for xyz_axis in range(1, 4):
                self.direct_xyz[atom][xyz_axis] = rpt[xyz_axis]

        # --- Rotation 2: align b-axis into the Cartesian xy plane ---
        # Rotation axis is the current a-axis, which after rotation 1 is
        # coincident with the Cartesian x-axis.
        # We want to rotate the b-axis so it has zero z-component (lies in the
        # xy plane).  The key subtlety: we must NOT compute the angle between
        # b itself and the y-axis.  If b were already in the xy plane but not
        # aligned with y, that computation would yield a bogus non-zero angle
        # and we would incorrectly rotate b.  Instead, project b onto the yz
        # plane (yz_proj = [None, 0, b_y, b_z]) and compute the angle between
        # that projection and y.  This gives exactly the out-of-plane tilt
        # that needs to be removed.  If the yz-projection points in the
        # positive-z half-space, invert the rotation axis to sweep back toward
        # y.
        rot_axis = [None, 1.0, 0.0, 0.0]
        # Project the b-axis onto the yz plane; rotate this projection to y.
        # The 1-indexed slot layout is [None, x=0, y=b_y, z=b_z] so the
        # vector helper sees a conventional [None, x, y, z] vector with
        # x forced to zero.
        yz_proj = [None,
                   0.0,
                   self.real_lattice[2][2],
                   self.real_lattice[2][3]]
        y_axis        = [None, 0.0, 1.0, 0.0]
        rot_angle_deg = math.degrees(self.get_vector_angle(yz_proj, y_axis))

        # If the yz-projection points in the positive-y half-space, the naive
        # rotation would sweep away from y; invert the rotation axis to
        # correct.  (Perl's commentary referenced "positive z quadrant" but
        # the code checks the y-component; we preserve that behavior.)
        if yz_proj[2] > 0.0:
            rot_axis = [None, -1.0, 0.0, 0.0]

        self.define_rot_matrix(rot_axis, rot_angle_deg)

        # Rotate the three lattice vectors in place.
        for abc_axis in range(1, 4):
            rpt = self.rotate_one_point(self.real_lattice[abc_axis])
            for xyz_axis in range(1, 4):
                self.real_lattice[abc_axis][xyz_axis] = rpt[xyz_axis]
        for abc_axis in range(1, 4):
            for xyz_axis in range(1, 4):
                if abs(self.real_lattice[abc_axis][xyz_axis]) < EPSILON:
                    self.real_lattice[abc_axis][xyz_axis] = 0.0

        # Rotate all atoms.
        for atom in range(1, self.num_atoms + 1):
            rpt = self.rotate_one_point(self.direct_xyz[atom])
            for xyz_axis in range(1, 4):
                self.direct_xyz[atom][xyz_axis] = rpt[xyz_axis]

        # Recompute derived lattice quantities after the two rotations.
        self.make_inv_or_recip_lattice(make_recip=True)
        self.abc_alpha_beta_gamma()

        # Recompute direct_abc and fract_abc for all atoms.
        for atom in range(1, self.num_atoms + 1):
            self.get_direct_abc(atom)
            self.get_fract_abc(atom)

    # ======================================================================
    # Section: Structure manipulation
    # ======================================================================

    def add_connection_atoms(self):
        """Add covalently bonded connection atoms to tagged atoms.

        The connection_tag for an atom is a ':'-separated list of element
        symbols (e.g. 'H:H:O') naming atoms to be bonded to it.  Up to four
        connection atoms are supported; the geometry is tetrahedral (bond
        angles of 109.5°, beta = 1.91113553 rad).

        The element name of the central (host) atom is extracted from its
        atom_tag by splitting on digits, e.g. 'Si3' → 'si' (Perl:
        prepLine on the tag with regex '[0-9]+').

        Algorithm per atom (one axis at a time, all atoms frozen at their
        current num_atoms to avoid processing newly added connection atoms):

          1. Parse the ':'-separated connection_tag to get element names.
          2. Bond length = sum of covalent radii of central + connection atom
             (in angstroms, from coval_radii table).
          3. Atom 1: always added if we reach this point.  Choose a
             uniformly random direction: theta ∈ [0, π], phi ∈ [0, 2π].
             Place it at bond_length from the central atom in that direction:
               x = sin(θ)cos(φ)·L + x₀
               y = sin(θ)sin(φ)·L + y₀
               z = cos(θ)·L       + z₀
          4. Atom 2: position is constrained by atom 1.  Two-step process:
             (a) Compute a provisional position at the cell origin using
                 (theta+beta, phi) — i.e., the same azimuth as atom 1 but
                 tilted by the tetrahedral angle — with the appropriate
                 bond_length:
                   tmp = [sin(θ+β)cos(φ)·L, sin(θ+β)sin(φ)·L, cos(θ+β)·L]
             (b) Build the unit vector from the central atom toward atom 1:
                   conn_one_unit = (direct_xyz[atom1] − direct_xyz[central])
                                   / bond_length_1
                 (Perl computes this as if atom 1 were at the cell origin.)
             (c) Choose a random phi_prime ∈ [0, 2π) and rotate tmp about
                 conn_one_unit by phi_prime (rotate_arb_axis).
             (d) Translate the result by the central atom's position.
          5. Atom 3: identical two-step process (same tmp starting point as
             atom 2), but the rotation angle is phi_prime += 2.0944 rad (120°).
             Note: the Perl inline comment for atom 3 says "phiPrime + beta
             radians (109.5°)" — this is a documentation error in the Perl;
             the code uses 2.0944 rad = 120°, not beta = 109.5°.
             Also note: Perl accidentally looks up element index [2] instead
             of [3] for the covalent radius — likely a copy-paste bug.
          6. Atom 4: same tmp, phi_prime += 2.0944 (another 120°, total 240°
             from atom 2's phi_prime).  Perl comment says "phiPrime − beta"
             but the code again does += 2.0944; the comment is wrong.
             Same copy-paste bug: Perl uses element index [2] for radius.

        Each new atom is appended with atom_tag = lowercase_element + '1' and
        an empty connection_tag so it does not trigger further additions.
        Coordinates are stored in direct_xyz and then converted to direct_abc
        and fract_abc.

        Note: atoms 2–4 are not yet implemented (require rotate_arb_axis).
        """
        beta = 1.91113553  # radians; tetrahedral angle 109.5°

        # Snapshot the current atom count so newly added connection atoms are
        # never themselves processed as connection hosts in this same pass.
        initial_num_atoms = self.num_atoms
        for atom in range(1, initial_num_atoms + 1):
            # If the connection_tag list is empty, skip to the next atom.
            if not self.connection_tag[atom]:
                continue

            # Get the list of atoms connected to this atom.
            # connection_tag is ':'-separated element symbols, e.g. 'H:H:H:H'
            conn_elements = self.connection_tag[atom].split(':')
            num_conn = len(conn_elements)

            # If the number of connected atoms is greater than 4 we abort
            # because we don't know how to add such a pattern of atoms.
            if num_conn > 4:
                raise ValueError(
                    f"Cannot connect more than four atoms for atom number {atom}."
                )

            # Get the element name of the initial atom (the one being connected
            # to) by stripping any trailing digits from the atom_tag.
            init_elem   = re.split(r'[0-9]+', self.atom_tag[atom])[0].lower()
            # Get the covalent radius of the initial atom in angstroms.
            init_z      = self.element_data.get_element_z(init_elem)
            init_radius = self.element_data.coval_radii[init_z]

            # ------------------------------------------------------------------
            # Treat the first connection atom.  If we reached this point then
            # there is always at least one connection atom.  The first connection
            # atom is always added a distance away equal to the sum of the
            # covalent radii, and in a direction that is randomly determined in
            # spherical coordinates.
            # ------------------------------------------------------------------
            conn1_elem  = conn_elements[0].lower()
            conn1_z     = self.element_data.get_element_z(conn1_elem)
            # Find the bond distance as the sum of covalent radii of the
            # connection atom and the initial atom (the one being connected to).
            bond_length = init_radius + self.element_data.coval_radii[conn1_z]

            # Choose a random theta (0..pi) and phi (0..2pi).
            theta = random.uniform(0.0, PI)
            phi   = random.uniform(0.0, 2.0 * PI)

            # Increment the number of atoms in the system.
            self.num_atoms += 1
            # Record the atom tag for this atom (lowercase element + '1').
            self.atom_tag.append(conn1_elem + '1')
            # New atoms are never themselves connection hosts.
            self.connection_tag.append('')

            # Compute the x,y,z coordinates of the connected atom (including
            # the translation from the cell origin to the host atom position).
            dx = math.sin(theta) * math.cos(phi) * bond_length + self.direct_xyz[atom][1]
            dy = math.sin(theta) * math.sin(phi) * bond_length + self.direct_xyz[atom][2]
            dz = math.cos(theta) * bond_length                 + self.direct_xyz[atom][3]
            self.direct_xyz.append([None, dx, dy, dz])
            self.fract_abc.append([None, 0.0, 0.0, 0.0])
            self.direct_abc.append([None, 0.0, 0.0, 0.0])
            # Convert to direct a,b,c and then to fractional a,b,c.
            self.get_direct_abc(self.num_atoms)
            self.get_fract_abc(self.num_atoms)

            # If we don't have a second connection atom then skip to next atom.
            if num_conn < 2:
                continue

            # ------------------------------------------------------------------
            # Treat the second (and subsequent) connection atoms.  Each atom's
            # position is constrained by the position of the first.  Strategy:
            #   (a) Compute a provisional point at the cell origin using angle
            #       (theta+beta, phi) — same azimuth as atom 1, tilted by the
            #       tetrahedral angle — at the new bond_length.
            #   (b) Compute a unit vector from the central atom toward atom 1
            #       (as if atom 1 were at the cell origin) and rotate the
            #       provisional point about that axis:
            #         atom 2: by a random phi_prime ∈ [0, 2π)
            #         atom 3: by phi_prime += 2.0944 rad (120°)
            #         atom 4: by phi_prime += 2.0944 rad (another 120°, = 240°)
            #       (Note: Perl inline comments for atoms 3 and 4 incorrectly
            #       state "±beta"; the actual rotation increments are 120°.)
            #   (c) Translate the result by the central atom's position.
            # Atoms 2–4 require rotate_arb_axis, which is not yet implemented.
            # ------------------------------------------------------------------
            raise NotImplementedError(
                "Connection atoms beyond the first (2nd–4th) require "
                "rotate_arb_axis, which is not yet implemented."
            )

    def shift_xyz_center(self, atom1, atom2):
        """Translate all atoms in atom1..atom2 linearly along each orthogonal
        axis so that the atom system as a whole is centered in the simulation
        cell.

        This routine is only applicable to molecular (non-periodic) systems,
        and requires that max_pos / min_pos have been computed for the range.
        The Python version calls get_min_max_xyz internally for convenience;
        the Perl original required the caller to do so beforehand.

        Parameters
        ----------
        atom1, atom2 : int
            Inclusive 1-indexed atom range to shift.

        Algorithm (Perl-authoritative)
        --------------------------------
        1. AtomSysCenter[ax] = (max_pos[ax] + min_pos[ax]) / 2
           equivalently: Min + (Max - Min)/2 = (Max + Min)/2
        2. CellCenter[ax]    = mag[ax] / 2
        3. Diff[ax]          = AtomSysCenter[ax] - CellCenter[ax]
        4. For every atom in the range: direct_xyz[atom][ax] -= Diff[ax]
        5. Recompute direct_abc and fract_abc for each shifted atom.
        """
        # Refresh max_pos / min_pos over the requested range.
        self.get_min_max_xyz(atom1, atom2)

        # Atom system center: midpoint of the bounding box spanned by the atoms.
        # AtomSysCenter = Min + (Max - Min)/2 = (Max + Min)/2
        center_pos  = [None] + [(self.max_pos[ax] + self.min_pos[ax]) / 2.0
                                for ax in range(1, 4)]

        # Cell center: midpoint of the simulation cell along each axis.
        center_cell = [None, self.mag[1] / 2.0,
                             self.mag[2] / 2.0,
                             self.mag[3] / 2.0]

        # Diff = AtomSysCenter - CellCenter (how far the cloud is off-center).
        center_diff = [None] + [center_pos[ax] - center_cell[ax]
                                for ax in range(1, 4)]

        for atom in range(atom1, atom2 + 1):
            for axis in range(1, 4):
                # Shift each atom by the difference between the atom system
                # center and the cell center, moving the cloud to the middle.
                self.direct_xyz[atom][axis] -= center_diff[axis]

            # Keep derived coordinate arrays consistent.
            self.get_direct_abc(atom)
            self.get_fract_abc(atom)

    def translate_atoms(self, atom1, atom2, displacement):
        """Translate atoms atom1..atom2 by a fixed displacement vector.

        The displacement is given in Angstroms in the direct (Cartesian XYZ)
        coordinate frame, and is applied by adding it to direct_xyz for each
        atom in the range.  After the shift, direct_abc and fract_abc are
        recomputed from direct_xyz so that all three coordinate arrays remain
        consistent.

        Special case — displacement == [0, 0, 0]:
            Instead of a null translation, this acts as "center the atom
            cloud inside the simulation cell".  The midpoint of the
            axis-aligned bounding box of the selected atoms is computed and
            then shifted to coincide with the midpoint of the cell.  This
            behaviour is only meaningful for molecular (non-periodic) systems.
            In the Perl original the caller was required to call
            getMinMaxXYZ before shiftXYZCenter; here shift_xyz_center
            handles that internally.

        Note on displacement indexing:
            The displacement is 1-indexed: [None, dx, dy, dz], consistent
            with the Perl original and the rest of the codebase's coordinate
            arrays (direct_xyz, direct_abc, fract_abc).

        After either branch, check_bounding_box(atom1, atom2) is called to
        wrap any atom that drifted outside [0, 1) in fractional coordinates
        back into the simulation cell via periodic boundary conditions.  Only
        atoms in the translated range are wrapped, matching the Perl original.

        Parameters
        ----------
        atom1, atom2 : int
            Inclusive 1-indexed atom range to translate.
        displacement : list of float
            1-indexed translation in Angstroms:
            [None, dx, dy, dz].  Pass [None, 0.0, 0.0, 0.0]
            to center the atom cloud in the cell.
        """
        if (displacement[1] == 0.0
                and displacement[2] == 0.0
                and displacement[3] == 0.0):
            # Special 0 0 0 designation: center the atom
            # cloud in the cell.  shift_xyz_center computes
            # get_min_max_xyz internally, then shifts
            # direct_xyz and re-derives direct_abc /
            # fract_abc for each atom.
            self.shift_xyz_center(atom1, atom2)
        else:
            for atom in range(atom1, atom2 + 1):
                # Add the requested XYZ displacement to
                # each atom's direct Cartesian coordinates.
                for axis in range(1, 4):
                    self.direct_xyz[atom][axis] += (
                        displacement[axis]
                    )

                # Re-derive direct_abc and fract_abc from the updated direct_xyz
                # so all three coordinate representations stay in sync.
                self.get_direct_abc(atom)
                self.get_fract_abc(atom)

        # Wrap any atom that ended up outside the simulation
        # cell back inside using periodic boundary conditions
        # (fract_abc wrapping).  Only wrap the atoms that were
        # just translated — wrapping all atoms would disturb
        # molecules that were previously placed correctly.
        self.check_bounding_box(atom1, atom2)

    def insert_vacuum(self, vac_axis, vac_amt):
        """Insert a vacuum region by expanding a lattice vector magnitude.

        The vacuum is introduced purely by increasing the length of the chosen
        lattice vector (a, b, or c).  Atom positions in real space are
        completely unaffected — both the Cartesian (direct_xyz) and the
        direct-lattice Angstrom (direct_abc) coordinates remain the same.
        Only the fractional coordinates (fract_abc) must be recomputed,
        because they are expressed relative to the (now larger) lattice.

        Parameters
        ----------
        vac_axis : int
            Lattice vector to expand: 1 = a, 2 = b, 3 = c.
        vac_amt : float
            Length to add to the chosen lattice vector, in Angstroms.
        """
        # Increase the magnitude of the requested lattice direction.
        self.mag[vac_axis] += vac_amt

        # Rebuild ABC vectors (in XYZ form) from the new lattice magnitudes.
        self.get_abc_vectors()

        # Reobtain the inverse lattice and the reciprocal lattice vectors.
        self.make_inv_or_recip_lattice(make_recip=False)
        self.make_inv_or_recip_lattice(make_recip=True)

        self.abc_alpha_beta_gamma()

        # Recompute fractional ABC positions from the direct-space ABC
        # coordinates, which are still the same (atoms didn't move).
        # The direct-space XYZ coordinates are also still the same.
        for atom in range(1, self.num_atoms + 1):
            self.get_fract_abc(atom)

    def cut_block(self, zone, abcxyz_flag, block_border):
        """Remove a designated set of atoms from the currently stored model.

        The basic algorithm is to build a duplicate list containing only those
        atoms that should NOT be cut, then replace all relevant data structures
        with that filtered list.

        Parameters
        ----------
        zone : int
            Which region to cut:
              0 = cut atoms found *inside* the block (keep outside atoms);
              1 = cut atoms found *outside* the block (keep inside atoms).
        abcxyz_flag : int
            Coordinate system used for the block boundaries:
              0 = abc direct-lattice coordinates;
              1 = Cartesian xyz coordinates.
        block_border : list
            1-indexed list; block_border[axis] = [None, low, high] for
            axis in {1, 2, 3}.  The 'high' value may be the sentinel string
            'a', 'b', or 'c' (for axes 1, 2, 3 respectively), which is
            interpreted as a request for the full lattice magnitude and is
            replaced with self.mag[axis] before comparison.
        """
        # The block borders may use a letter sentinel in the "to" (high) slot
        #   to indicate a request for the maximum lattice extent.  Replace
        #   with the actual magnitude before any comparisons.
        sentinels = {1: 'a', 2: 'b', 3: 'c'}
        for axis in range(1, 4):
            if block_border[axis][2] == sentinels[axis]:
                block_border[axis][2] = self.mag[axis]

        # Build new parallel arrays that will hold only the atoms we keep.
        new_atom_element_name = [None]
        new_atom_element_id   = [None]
        new_atom_species_id   = [None]
        new_fract_abc         = [None]
        new_direct_abc        = [None]
        new_direct_xyz        = [None]
        new_atom_tag          = [None]

        # Loop over each atom and decide whether to keep it.
        for atom in range(1, self.num_atoms + 1):
            # Choose the coordinate set to compare against the block borders.
            if abcxyz_flag == 0:
                coords = self.direct_abc[atom]   # [None, a, b, c]
            else:
                coords = self.direct_xyz[atom]   # [None, x, y, z]

            # An atom is "inside" only if it falls within all three axis ranges.
            found_inside = True
            for axis in range(1, 4):
                if (coords[axis] < block_border[axis][1] or
                        coords[axis] > block_border[axis][2]):
                    found_inside = False
                    break

            # zone==0: cut inside atoms; zone==1: cut outside atoms.
            if zone == 0:
                cut = found_inside       # cutting inside
            else:
                cut = not found_inside   # cutting outside

            # Keep this atom by appending it to the new arrays.
            if not cut:
                new_atom_element_name.append(self.atom_element_name[atom])
                new_atom_element_id.append(self.atom_element_id[atom])
                new_atom_species_id.append(self.atom_species_id[atom])
                new_fract_abc.append(self.fract_abc[atom][:])
                new_direct_abc.append(self.direct_abc[atom][:])
                new_direct_xyz.append(self.direct_xyz[atom][:])
                new_atom_tag.append(self.atom_tag[atom])

        # Replace the system arrays with the filtered versions.
        self.num_atoms         = len(new_fract_abc) - 1
        self.atom_element_name = new_atom_element_name
        self.atom_element_id   = new_atom_element_id
        self.atom_species_id   = new_atom_species_id
        self.fract_abc         = new_fract_abc
        self.direct_abc        = new_direct_abc
        self.direct_xyz        = new_direct_xyz
        self.atom_tag          = new_atom_tag

    def cut_sphere(self, zone, abcxyz_flag, sphere_rad, target_atom, sphere_loc):
        """Remove atoms inside or outside a sphere centered at a given location.

        The algorithm builds a duplicate set of atom arrays containing only the
        atoms that should NOT be cut, then replaces the model's atom arrays with
        those survivors.  All coordinate representations (fract_abc, direct_abc,
        direct_xyz) and the identity arrays (element name/id, species id, tag)
        are rebuilt simultaneously.

        Parameters
        ----------
        zone : int
            Which region to cut:
              0 = cut atoms inside the sphere (keep the shell outside);
              1 = cut atoms outside the sphere (keep the core inside).
        abcxyz_flag : int
            Coordinate system for sphere_loc:
              0 = direct-lattice abc coordinates (Angstroms along a, b, c);
              1 = Cartesian xyz coordinates (Angstroms).
            Ignored and forced to 1 when target_atom is non-zero.
        sphere_rad : float
            Radius of the sphere in the same units as the coordinate system
            selected by abcxyz_flag (Angstroms).
        target_atom : int
            If non-zero, the sphere is centered on this atom's Cartesian xyz
            position.  sphere_loc is ignored and abcxyz_flag is forced to 1.
            Pass 0 to use sphere_loc explicitly.
        sphere_loc : list
            1-indexed center coordinates [None, coord1, coord2, coord3].
            Interpreted as abc or xyz depending on abcxyz_flag (unless
            overridden by target_atom).
        """
        # In the case that a particular atom site was given, define the sphere
        # location in XYZ coordinates according to that atom's position and
        # force the coordinate flag to xyz mode.
        if target_atom != 0:
            sphere_loc = [None,
                          self.direct_xyz[target_atom][1],
                          self.direct_xyz[target_atom][2],
                          self.direct_xyz[target_atom][3]]
            abcxyz_flag = 1

        # Build new (survivor-only) versions of every per-atom array.
        new_atom_element_name = [None]
        new_atom_element_id   = [None]
        new_atom_species_id   = [None]
        new_fract_abc         = [None]
        new_direct_abc        = [None]
        new_direct_xyz        = [None]
        new_atom_tag          = [None]

        # Loop over each atom and decide whether it should be kept or cut.
        for atom in range(1, self.num_atoms + 1):

            # Select the coordinate representation to use for the distance
            # test: direct abc (along lattice vectors) or Cartesian xyz.
            if abcxyz_flag == 0:
                coords = self.direct_abc[atom]   # [None, a, b, c]
            else:
                coords = self.direct_xyz[atom]   # [None, x, y, z]

            # Determine whether the atom is inside the sphere by computing its
            # distance from the sphere center and comparing to the radius.
            distance = ((coords[1] - sphere_loc[1])**2 +
                        (coords[2] - sphere_loc[2])**2 +
                        (coords[3] - sphere_loc[3])**2) ** 0.5
            found_inside = distance <= sphere_rad

            # zone==0: cut atoms found inside; zone==1: cut atoms found outside.
            # An atom is "cut" when it falls in the designated zone.
            if zone == 0:
                cut = found_inside       # cutting inside, found inside → cut
            else:
                cut = not found_inside   # cutting outside, found outside → cut

            # Accumulate all coordinate and identity data for atoms we keep.
            if not cut:
                new_atom_element_name.append(self.atom_element_name[atom])
                new_atom_element_id.append(self.atom_element_id[atom])
                new_atom_species_id.append(self.atom_species_id[atom])
                new_fract_abc.append(self.fract_abc[atom][:])
                new_direct_abc.append(self.direct_abc[atom][:])
                new_direct_xyz.append(self.direct_xyz[atom][:])
                new_atom_tag.append(self.atom_tag[atom])

        # Replace the model's atom arrays with the survivor-only versions.
        self.num_atoms         = len(new_fract_abc) - 1
        self.atom_element_name = new_atom_element_name
        self.atom_element_id   = new_atom_element_id
        self.atom_species_id   = new_atom_species_id
        self.fract_abc         = new_fract_abc
        self.direct_abc        = new_direct_abc
        self.direct_xyz        = new_direct_xyz
        self.atom_tag          = new_atom_tag

    def make_ortho(self):
        """Derive an orthorhombic cell from the current (non-orthorhombic) lattice.

        The goal is to take the crystal in its given lattice and deduce a set of
        orthorhombic lattice parameters that retains periodic boundary conditions.
        Essentially we make a=x, b=y, and c=z while retaining the magnitudes of
        the ax, by, and cz diagonal components for a, b, and c respectively.
        All off-diagonal lattice components are zeroed and all angles are set to
        90 degrees (pi/2 radians).

        This can only be done correctly for certain lattice types — hexagonal
        cells in particular — that have *already* been doubled (via the -sc
        supercell option) in the direction of the lattice vector that will be
        modified.  In other words, you have to plan ahead to use this; it will
        not do everything automatically.

        Atom Cartesian positions (direct_xyz) are kept fixed.  direct_abc and
        fract_abc are recomputed from those xyz positions using the newly
        defined orthogonal lattice, and then any atoms that have drifted outside
        [0, 1) are wrapped back into the simulation box.
        """
        # Clearly, all angles must be 90 degrees (= pi/2 radians) for an
        # orthorhombic cell.
        for axis in range(1, 4):
            self.angle[axis]     = math.pi / 2.0
            self.angle_deg[axis] = 90.0

        # Further, the magnitudes of the ABC vectors will be the values
        # currently stored in the diagonal ax, by, and cz lattice elements.
        for abc_axis in range(1, 4):
            self.mag[abc_axis] = self.real_lattice[abc_axis][abc_axis]

        # Eliminate components that are off-axis (i.e. not ax, by, or cz).
        for abc_axis in range(1, 4):
            for xyz_axis in range(1, 4):
                # The ABC axis and XYZ axis *should* be collinear for the
                # diagonal element, so no adjustment is needed there.
                if abc_axis == xyz_axis:
                    continue
                # It is necessary to make this ABC/XYZ axis pair orthogonal.
                self.real_lattice[abc_axis][xyz_axis] = 0.0

        # Obtain the inverse of the real lattice.  This must be done now so
        # that we can convert atom positions from xyz to abc using the newly
        # defined orthogonal lattice.  (It will be needed again later if a
        # space group or supercell is subsequently applied.)
        self.make_inv_or_recip_lattice(make_recip=False)

        # Finally, it is necessary to make sure that all atoms in the model
        # are within the simulation box.  This is done by first recomputing
        # direct_abc and fract_abc for every atom from their current direct_xyz
        # Cartesian positions under the new lattice definition, and then
        # calling check_bounding_box to wrap any out-of-range atoms back in.
        for atom in range(1, self.num_atoms + 1):
            self.get_direct_abc(atom)
            self.get_fract_abc(atom)

        self.check_bounding_box()

    def apply_perturbation(self, magnitude):
        """Apply a random displacement to each atom in the model.

        For each atom, a random displacement vector is computed using spherical
        coordinates (R, Theta, Phi) and then added to the atom's Cartesian
        (direct_xyz) position.  The magnitude parameter is the maximum allowed
        displacement radius R, specified in Angstroms and typically supplied
        from the command line.

        After each atom is displaced, the updated Cartesian coordinates are
        propagated to the other representations stored in the module
        (direct_abc and fract_abc), so all coordinate arrays remain consistent.

        Finally, periodicity is enforced: any atom that has drifted outside the
        simulation box is wrapped back inside via check_bounding_box.

        Note: the Perl original had a latent bug — it declared `my $numAtoms`
        as a local (uninitialized) variable, which shadowed the package-level
        $numAtoms.  The range `1..$numAtoms` therefore collapsed to `1..0`
        (empty), so the loop never executed.  This Python port correctly
        iterates over `range(1, self.num_atoms + 1)`.

        Parameters
        ----------
        magnitude : float
            Maximum perturbation radius R in Angstroms.  Each atom receives a
            random displacement with R in [0, magnitude), Theta in [0, pi),
            and Phi in [0, 2*pi).
        """
        import random
        import math

        # For each atom, apply the necessary translation.  The displacement is
        # chosen in spherical coordinates (R, Theta, Phi) and converted to XYZ.
        for atom in range(1, self.num_atoms + 1):
            # Compute the XYZ perturbation for this atom based on random R,
            # Theta, and Phi values.  The maximum R value is maxPertMag
            # (magnitude), given in Angstroms.
            r     = random.uniform(0, magnitude)
            theta = random.uniform(0, math.pi)
            phi   = random.uniform(0, 2.0 * math.pi)

            # Convert spherical displacement to Cartesian components.
            dx = r * math.sin(theta) * math.cos(phi)
            dy = r * math.sin(theta) * math.sin(phi)
            dz = r * math.cos(theta)

            self.direct_xyz[atom][1] += dx
            self.direct_xyz[atom][2] += dy
            self.direct_xyz[atom][3] += dz

            # Propagate the updated XYZ position to the other representations
            # (directABC and fractABC) so all coordinate arrays stay consistent.
            self.get_direct_abc(atom)
            self.get_fract_abc(atom)

        # Make sure that all atoms are inside the simulation box (adjusted for
        # periodicity).
        self.check_bounding_box()

    def apply_filter(self, min_dist_filter):
        """Remove atoms that are too close to any already-accepted atom.

        Iterates over atoms in order; each atom that has not been rejected is
        kept.  Any atom in the extended cell that is within min_dist_filter of
        a kept atom is added to the rejected set so it will be skipped.  After
        the pass, the atom list is compacted to contain only the kept atoms.

        set_limit_dist is called first so that create_min_dist_matrix builds
        the distance matrix using min_dist_filter as the periodic-cell
        interaction cutoff.

        Note: the maps and some extended-position arrays are NOT preserved by
        this routine and may need to be recomputed afterwards if required.

        Parameters
        ----------
        min_dist_filter : float
            Minimum allowed interatomic distance (same units as direct_xyz,
            i.e. Angstroms).  Pairs closer than or equal to this value cause
            the second atom to be rejected.
        """
        # Set the limit_dist value for periodic cell interaction using
        # min_dist_filter as the interaction factor for this operation.
        self.set_limit_dist(min_dist_filter)

        # Create the minimal distance matrix.
        self.create_min_dist_matrix()

        # Local copies of per-atom data for the surviving atoms.
        fract_abc_local         = [None]
        atom_element_name_local = [None]
        atom_element_id_local   = [None]
        atom_species_id_local   = [None]

        new_num_atoms = 0
        rejected = set()

        # Consider each atom in turn and determine if it has any neighbours
        # that are too close.  Any neighbour that is too close is recorded and
        # will be removed from the system.  Print a running progress count so
        # the user does not get too bored waiting; this process can take a
        # while for large systems.
        for atom in range(1, self.num_atoms + 1):
            if atom % 10 == 0:
                print('|', end='', flush=True)
            else:
                print('.', end='', flush=True)
            if atom % 50 == 0:
                print(f' {atom}')

            # If this atom has already been marked for rejection, skip it.
            if atom in rejected:
                continue

            # Keep this atom: record its data in the local survivor arrays.
            new_num_atoms += 1
            fract_abc_local.append([None,
                                     self.fract_abc[atom][1],
                                     self.fract_abc[atom][2],
                                     self.fract_abc[atom][3]])
            atom_element_name_local.append(self.atom_element_name[atom])
            atom_element_id_local.append(self.atom_element_id[atom])
            atom_species_id_local.append(self.atom_species_id[atom])

            # Mark any extended-cell neighbour that is too close as rejected.
            for neighbor_atom in range(atom + 1, self.num_atoms_ext + 1):
                if self.min_dist[atom][neighbor_atom] <= min_dist_filter:
                    central = self.ext_to_central_item_map[neighbor_atom]
                    print(f'Rejecting {central}')
                    rejected.add(central)

        if self.num_atoms % 50 != 0:
            print()

        self.set_num_atoms(new_num_atoms)

        # Save the new atom positions from the local survivor arrays and
        # propagate them to the other coordinate representations (direct_xyz,
        # direct_abc).  Also record the element and species naming data.
        for atom in range(1, self.num_atoms + 1):
            for axis in range(1, 4):
                self.fract_abc[atom][axis] = fract_abc_local[atom][axis]
            self.get_direct_xyz(atom)
            self.get_direct_abc(atom)
            self.atom_element_name[atom] = atom_element_name_local[atom]
            self.atom_element_id[atom]   = atom_element_id_local[atom]
            self.atom_species_id[atom]   = atom_species_id_local[atom]

    def check_bounding_box(self, atom1=None, atom2=None):
        """Wrap atoms outside [0, 1) in fractional coordinates
        back into the simulation cell via periodic boundary
        conditions.

        Parameters
        ----------
        atom1, atom2 : int, optional
            Inclusive 1-indexed atom range to check.  When
            omitted, all atoms (1..num_atoms) are processed.

        This is the translation-style wrapping (shiftStyle=0
        in the Perl original).  The Perl subroutine also
        supported a rotation style (shiftStyle=1) that is not
        needed here because rotate_all_atoms handles that case
        differently.

        Algorithm (translation style):
          For each atom, check all three fractional axes:
            - fract_abc > 1.0 -> subtract 1 (off +face)
            - fract_abc < 0.0 -> add 1      (off -face)
          If any axis was shifted, re-propagate the full
          coordinate triplet so that direct_xyz and direct_abc
          stay consistent with the corrected fract_abc.
        """
        if atom1 is None:
            atom1 = 1
        if atom2 is None:
            atom2 = self.num_atoms
        for atom in range(atom1, atom2 + 1):
            moved = False
            for axis in range(1, 4):
                if self.fract_abc[atom][axis] > 1.0:
                    self.fract_abc[atom][axis] -= 1.0
                    moved = True
                elif self.fract_abc[atom][axis] < 0.0:
                    self.fract_abc[atom][axis] += 1.0
                    moved = True
            if moved:
                self.get_direct_xyz(atom)
                self.get_direct_abc(atom)
                self.get_fract_abc(atom)

    def sort_atoms(self):
        """Sort atoms alphabetically by element name.

        NOTICE: This must be called before any calculations are performed,
        such as the creation of the minimal distance matrix or the extended
        atom positions that account for periodic boundary conditions.  It
        should only be used immediately after reading a data set.  Otherwise
        a data structure could be built from the unsorted data, the data
        would then be sorted, and the mapping between the data structure and
        the atom arrays would be lost.

        The sort key is atom_element_name (lowercase).  Python's built-in
        ``sorted`` is stable (Timsort) so the (key, original_index) pair
        technique used by Perl's ``stableSort`` is implemented here with
        a plain call — then the resulting permutation is passed to
        ``apply_sort``.
        """
        # Build 0-indexed list of element names for sorting.
        element_names = [self.atom_element_name[atom].lower()
                         for atom in range(1, self.num_atoms + 1)]

        # Sort 0-indexed positions by element name, then convert to a
        # 1-indexed order array where order[new_pos] = old_pos.
        sorted_0 = sorted(range(len(element_names)),
                          key=lambda i: element_names[i])
        order = [None] + [idx + 1 for idx in sorted_0]

        self.apply_sort(order)

    def apply_sort(self, order):
        """Reorder all per-atom arrays in-place according to a 1-indexed
        permutation array.

        Each per-atom array uses 1-based indexing: index 0 is a None
        placeholder and is never touched.  Only indices 1..num_atoms are
        moved.  The algorithm is: copy the list to a temporary array first,
        then write back using the permutation so that reads and writes don't
        interfere.

        Parameters
        ----------
        order : list
            1-indexed permutation where order[new_pos] = old_pos (1-based).
            order[0] is unused (None or empty).
        """
        def _reorder(lst):
            # Index 0 is the None placeholder, so only act when the list has
            # at least one real element (len > 1).
            if len(lst) <= 1:
                return
            # Copy the list to a temporary array so that reads below always
            # come from the original ordering.
            temp = lst[:]
            # Apply the sorted indices: position i in the result takes its
            # value from position order[i] in the original.
            for i in range(1, self.num_atoms + 1):
                lst[i] = temp[order[i]]

        _reorder(self.atom_element_name)
        _reorder(self.atom_element_id)
        _reorder(self.atom_species_id)
        _reorder(self.atomic_z)
        _reorder(self.atom_tag)
        _reorder(self.dat_skl_map)
        _reorder(self.skl_dat_map)
        _reorder(self.molecule_name)
        _reorder(self.molecule_seq_num)
        _reorder(self.residue_name)
        _reorder(self.residue_seq_num)
        _reorder(self.fract_abc)
        _reorder(self.direct_abc)
        _reorder(self.direct_xyz)

    # ======================================================================
    # Section: Rotation
    # ======================================================================

    def rotate_arb_axis(self, axis, angle_deg):
        """Rotate the entire structure about an arbitrary axis through the origin.

        In the Perl module, ``rotateArbAxis`` was a low-level single-point
        helper (it rotated one passed-in point and returned nothing).  That
        role is filled here by :meth:`rotate_one_point`.  The Python method
        with this name is the public, full-structure entry point: it
        normalises the axis, builds the rotation matrix, and then rotates
        every atom via the two-pass PBC-safe :meth:`rotate_all_atoms`.

        Perl's ``rotateArbAxis`` accepted the rotation angle in **radians**;
        this method accepts it in **degrees** for consistency with the rest
        of the Python rotation interface.

        Parameters
        ----------
        axis : sequence of float
            1-indexed rotation axis vector ``[None, x, y, z]`` (need
            not be unit length; normalisation is done internally).
        angle_deg : float
            Rotation angle in degrees (positive = right-hand rule about
            the axis direction).
        """
        # Magnitude sums over slots 1..3 (slot 0 is the sentinel).
        mag = math.sqrt(axis[1]**2 + axis[2]**2 + axis[3]**2)
        unit = [None] + [axis[axis_idx] / mag for axis_idx in range(1, 4)]
        self.define_rot_matrix(unit, angle_deg)
        self.rotate_all_atoms()

    def define_rot_matrix(self, axis, angle_deg):
        """Build the 3×3 rotation matrix for the given axis and angle.

        Uses Rodrigues' rotation formula, stored in rot_matrix[row][col]
        (both 1-indexed), and applied via the row-vector form:
          result[k] = p[1]*M[1][k] + p[2]*M[2][k] + p[3]*M[3][k]

        The Perl ``defineRotMatrix`` references these sources:
          - http://mathworld.wolfram.com/RotationFormula.html
          - http://www.gamedev.net/reference/articles/article1199.asp
              (Note: there is a known error in the final matrix of that
              article — row 1 col 3 should be ``txz - sy``, NOT ``txz - xy``.
              The derivation steps in the article are correct; only the
              summary matrix has the typo.)
          - http://www.fastgraph.com/makegames/3drotation/
          - http://en.wikipedia.org/wiki/Rotation_matrix

        Variables: c = cos(angle), s = sin(angle), t = 1 - cos(angle).

        The Perl module also contained a commented-out implementation using
        the ``Math::MatrixReal`` library; the direct array version shown
        here (and in Perl) is the one actually used.

        Parameters
        ----------
        axis : sequence of float
            1-indexed unit rotation axis ``[None, x, y, z]``.  Must
            already be normalised; call :meth:`rotate_arb_axis` for
            automatic normalisation.
        angle_deg : float
            Rotation angle in degrees.
        """
        angle_rad = angle_deg * PI / 180.0
        c = math.cos(angle_rad)
        s = math.sin(angle_rad)
        t = 1.0 - c
        # Read the axis components from the 1-indexed storage slots so the
        # interface matches every other coordinate-accepting method.
        x, y, z = axis[1], axis[2], axis[3]

        # rot_matrix is already allocated in __init__ as a 4×4 list
        self.rot_matrix[1][1] = t*x*x + c
        self.rot_matrix[1][2] = t*x*y + s*z
        self.rot_matrix[1][3] = t*x*z - s*y    # gamedev typo: -sy not -xy

        self.rot_matrix[2][1] = t*x*y - s*z
        self.rot_matrix[2][2] = t*y*y + c
        self.rot_matrix[2][3] = t*y*z + s*x

        self.rot_matrix[3][1] = t*x*z + s*y
        self.rot_matrix[3][2] = t*y*z - s*x
        self.rot_matrix[3][3] = t*z*z + c

    def rotate_one_point(self, point):
        """Apply rot_matrix to a single Cartesian point about the origin.

        Implements the row-vector form from Perl ``rotateOnePoint``:
          pointVector = point - origin               # translate to origin
          result[k]   = pV[1]*M[1][k] + pV[2]*M[2][k] + pV[3]*M[3][k]
                                                     # (pV_x,pV_y,pV_z) x rM(3x3)
          output      = result + origin              # translate back

        The Perl version accepted an explicit ``orig_ref`` parameter so it
        could rotate about any point in space.  Here the origin is always
        (0, 0, 0) (as in Perl's ``rotateArbAxis`` helper), so the
        translate-to-origin and translate-back steps are no-ops and are
        omitted from the body for clarity.

        Perl's ``rotateOnePoint`` mutated the point in-place via a
        reference.  This Python version returns a new list instead so
        callers do not need to worry about aliasing surprises — a
        Python-side non-mutation improvement.

        The Perl module also had a commented-out ``Math::MatrixReal``
        version (``$rotPointVector = $pointVector->multiply($rotMatrix)``);
        only the direct index form is implemented here.

        Parameters
        ----------
        point : sequence of float
            1-indexed Cartesian point ``[None, x, y, z]``.

        Returns
        -------
        list of float
            1-indexed rotated point ``[None, x', y', z']``.
        """
        # Origin is (0,0,0), so pointVector = point (translate step is a no-op).
        # Read directly from the 1-indexed slots to mirror the surrounding code.
        px, py, pz = point[1], point[2], point[3]

        # Apply rotation matrix: (pV_x, pV_y, pV_z) x rM(3x3).
        # Both the axis loop and rot_matrix access use the 1-indexed
        # convention directly, eliminating the previous +1 offset bridge.
        result = [None, 0.0, 0.0, 0.0]
        for k in range(1, 4):
            result[k] = (px * self.rot_matrix[1][k] +
                         py * self.rot_matrix[2][k] +
                         pz * self.rot_matrix[3][k])

        # Origin is (0,0,0), so translate-back step is also a no-op.
        return result

    def rotate_all_atoms(self):
        """Apply rot_matrix to the Cartesian coordinates of every atom.

        Implements the two-pass PBC-safe algorithm from Perl
        ``rotateAllAtoms``.  After any rotation it is almost certain that
        some atoms will lie outside the simulation box.  These cannot simply
        be wrapped to the other side (naive PBC) because that would change
        their physical position relative to the rotated structure.  The
        Perl comment calls the correct algorithm "ugly":

          Save a copy of the original direct_xyz positions (done before
          pass 1).  Apply the rotation to all representations — direct_abc,
          fract_abc, direct_xyz, but not the saved copy.  Then restore
          *only* direct_xyz to its pre-rotation value.  Check which atoms
          are outside the box ***according to the now-rotated fract_abc***.
          For any offending atom, use the restored direct_xyz to recover the
          original fract_abc and shift it by ±1 along the offending axis so
          that when the rotation is re-applied the atom stays in the box.
          Finally re-apply the rotation to all atoms.

        Each pass propagates in the order:
          direct_xyz  →  get_direct_abc  →  get_fract_abc
        (matching Perl's "Propogate [sic] to other representations" calls).

        In Perl, step (1)+(2) — restore and shift — are handled by a call to
        ``checkBoundingBox(1, numAtoms, 1)`` where the third argument is a
        flag meaning "do the shift".  Python inlines that logic directly
        rather than delegating to a separate method.

        Note: Perl's ``defineRotMatrix`` takes ``($rotAngle, $rot_ref)``
        — angle first, then axis ref — the opposite of Python's
        ``define_rot_matrix(axis, angle_deg)``.  Both ultimately build the
        same rotation matrix; only the call-site argument order differs.

        Requires :meth:`define_rot_matrix` to have been called first.

        Post-condition (from Perl): "At this point there should be no atoms
        outside the simulation box at all."
        """
        # --- Pass 1: rotate, collect direct_xyz_copy and do_move ----------
        direct_xyz_copy = [[None, 0.0, 0.0, 0.0]
                           for _ in range(self.num_atoms + 1)]
        direct_xyz_copy[0] = None

        for atom in range(1, self.num_atoms + 1):
            for ax in range(1, 4):
                direct_xyz_copy[atom][ax] = self.direct_xyz[atom][ax]

        for atom in range(1, self.num_atoms + 1):
            rpt = self.rotate_one_point(self.direct_xyz[atom])
            for ax in range(1, 4):
                self.direct_xyz[atom][ax] = rpt[ax]
            self.get_direct_abc(atom)
            self.get_fract_abc(atom)

        # Record which axes went outside [0, 1) after rotation
        do_move = [[0, 0, 0, 0] for _ in range(self.num_atoms + 1)]
        do_move[0] = None
        for atom in range(1, self.num_atoms + 1):
            for ax in range(1, 4):
                if self.fract_abc[atom][ax] > 1.0:
                    do_move[atom][ax] = 1
                elif self.fract_abc[atom][ax] < 0.0:
                    do_move[atom][ax] = -1

        # --- Restore original direct_xyz ----------------------------------
        for atom in range(1, self.num_atoms + 1):
            for ax in range(1, 4):
                self.direct_xyz[atom][ax] = direct_xyz_copy[atom][ax]

        # --- Pre-shift offending atoms ------------------------------------
        for atom in range(1, self.num_atoms + 1):
            if any(do_move[atom][ax] != 0 for ax in range(1, 4)):
                # Recompute fract_abc from original (restored) direct_xyz
                self.get_direct_abc(atom)
                self.get_fract_abc(atom)
                # Shift fract_abc so the upcoming rotation stays in-box
                for ax in range(1, 4):
                    if do_move[atom][ax] == 1:
                        self.fract_abc[atom][ax] -= 1.0
                    elif do_move[atom][ax] == -1:
                        self.fract_abc[atom][ax] += 1.0
                # Propagate shifted fractional coords to Cartesian
                self.get_direct_xyz(atom)
                self.get_direct_abc(atom)
                self.get_fract_abc(atom)

        # --- Pass 2: re-rotate all atoms ----------------------------------
        for atom in range(1, self.num_atoms + 1):
            rpt = self.rotate_one_point(self.direct_xyz[atom])
            for ax in range(1, 4):
                self.direct_xyz[atom][ax] = rpt[ax]
            self.get_direct_abc(atom)
            self.get_fract_abc(atom)

    # ======================================================================
    # Section: Element and atom utilities
    # ======================================================================

    def create_element_list(self):
        """Build the list of unique elements present and assign element IDs.

        For each atom, the element name is extracted from atom_tag by
        stripping all digit characters (Perl: prepLine with pattern '[0-9]+';
        Python: re.split(r'[0-9]+', tag)[0]).  For example, 'si2' → 'si'.

        element_list is grown in order of first appearance; elements are
        never reordered.  element_list[0] is None (index-0 placeholder,
        matching Perl's $elementList[0] = "").  Valid entries start at
        index 1.

        Sets:
            num_elements         -- count of distinct elements found
            element_list[1..]    -- unique element names in appearance order
            atom_element_name[]  -- element name for each atom (1-indexed)
            atom_element_id[]    -- index into element_list for each atom
                                    (1-indexed; value is 1..num_elements)

        Called by read_olcao_skl() immediately after the atom list is read,
        before apply_space_group() expands the symmetry.
        """
        # Initialize the element count and element list.
        self.num_elements     = 0
        self.element_list     = [None]   # index-0 placeholder; entries at 1..
        self.atom_element_name = [None] * (self.num_atoms + 1)
        self.atom_element_id  = [None] * (self.num_atoms + 1)

        for atom in range(1, self.num_atoms + 1):
            # Define the element name for this atom by stripping any trailing
            # digits from the atom tag (e.g. 'si2' → 'si', 'o1' → 'o').
            element_name = re.split(r'[0-9]+', self.atom_tag[atom])[0]
            self.atom_element_name[atom] = element_name

            # Find the element name of this atom in the list of unique elements.
            found = 0
            for element in range(1, self.num_elements + 1):
                if element_name == self.element_list[element]:
                    found = element
                    break

            # Add a new element to the list if it was not found.
            if found == 0:
                self.num_elements += 1
                self.element_list.append(element_name)
                self.atom_element_id[atom] = self.num_elements
            else:
                self.atom_element_id[atom] = found

    def create_species_data(self, use_file_species=True):
        """Assign species IDs to every atom.

        A "species" is a named sub-type of an element distinguished by the
        tag string stored in atom_tag (e.g. 'si1', 'si2').  Two atoms of the
        same element are the same species if and only if their atom_tag
        strings are identical.

        species_list[element_id][species_id] holds the raw atom_tag string
        (e.g. 'si2') for that species.  Both indices are 1-based; index 0 is
        a None placeholder.  Species IDs are assigned in first-encounter
        order as atoms are traversed from 1..num_atoms.

        If use_file_species is True (default), the tag strings read from the
        input file (stored in atom_tag) are used directly to distinguish
        species.  This is the normal path: atoms tagged 'si1' and 'si2' in
        the SKL file will end up as two separate species of silicon.

        If use_file_species is False, species distinctions are discarded.
        Every element is given exactly one species (species 1), and its
        species_list entry is set to the bare element name concatenated with
        '1' (e.g. 'si' → 'si1').  This collapses all tag-based distinctions
        so that every silicon atom, regardless of its original tag, maps to
        species 1.

        Populates: atom_species_id, num_species, species_list.

        Parameters
        ----------
        use_file_species : bool
            True  -- respect species tags from the input file (default).
            False -- treat all atoms of the same element as one species.

        Called by read_olcao_skl() immediately after create_element_list().
        """
        # species_list[element_id][species_id] = tag string; both 1-indexed.
        self.species_list    = [None] + [[None] for _ in range(self.num_elements)]
        self.num_species     = [None] + [0] * self.num_elements
        self.atom_species_id = [None] * (self.num_atoms + 1)

        # If we are asked NOT to use the species designated in the file, set
        # the number of species of each element to 1, assign all atoms of
        # each element to species 1, and store the species tag as the bare
        # element name concatenated with '1' (e.g. 'si' → 'si1'), then return.
        if not use_file_species:
            for element in range(1, self.num_elements + 1):
                self.num_species[element] = 1
            for atom in range(1, self.num_atoms + 1):
                eid = self.atom_element_id[atom]
                self.species_list[eid] = [None, self.atom_element_name[atom] + '1']
                self.atom_species_id[atom] = 1
            return

        # Proceed only if we will use the species designated in the input
        # file, which are stored in the atom_tag array.

        # Initialize the count of the number of unique species for each
        # element to zero before scanning the atom list.
        # (num_species was already zeroed during array construction above.)

        for atom in range(1, self.num_atoms + 1):
            eid = self.atom_element_id[atom]

            # Determine if the tag associated with this atom's element
            # already exists in the species list for that element.
            found = 0
            for species in range(1, self.num_species[eid] + 1):
                if self.atom_tag[atom] == self.species_list[eid][species]:
                    found = species
                    break

            if found == 0:
                # New species encountered: extend the list and assign the
                # next available species ID.
                self.num_species[eid] += 1
                self.species_list[eid].append(self.atom_tag[atom])
                self.atom_species_id[atom] = self.num_species[eid]
            else:
                self.atom_species_id[atom] = found

    def count_element_atoms(self):
        """Get the number of atoms of each element in the system.

        Populates ``self.element_count`` as a 1-indexed list of length
        ``num_elements + 1`` (index 0 is the None placeholder).  Each entry
        ``element_count[e]`` holds the total number of atoms whose
        ``atom_element_id`` equals *e*.

        Must be called after ``create_element_list`` has established
        ``num_elements`` and ``atom_element_id``.
        """
        # Initialize the element count to zero for every element slot.
        self.element_count = [None] + [0] * self.num_elements

        for atom in range(1, self.num_atoms + 1):
            self.element_count[self.atom_element_id[atom]] += 1

    def map_element_number(self):
        """Map each atom's element name to its atomic number Z.

        Iterates over all atoms and looks up atomic_z[atom] by matching
        atom_element_name[atom] against the element database name table.
        The result is stored in self.atomic_z[1..num_atoms] (1-indexed;
        index 0 is a None placeholder).

        Prerequisite: the element database (self.element_data) must already
        be populated, i.e. setup_data_base_ref() must have been called first.

        Note: in the Perl original, getElementNumber walked the
        $elementNames_ref list directly; here we delegate to
        element_data.get_element_z() which performs the same lookup.
        """
        if len(self.atomic_z) < self.num_atoms + 1:
            self.atomic_z = [None] * (self.num_atoms + 1)
        for atom in range(1, self.num_atoms + 1):
            self.atomic_z[atom] = self.get_element_number(
                self.atom_element_name[atom])

    def compute_implicit_info(self):
        """Compute basis/potential summary data from element database.

        Scans all atoms and queries the element database (via self.atomic_z)
        to establish four global extrema used downstream by the OLCAO program:

        Attributes set:
          min_pot_alpha       Minimum Gaussian exponent found in any atom's
                              potential expansion.  Initialized to the sentinel
                              1000.0 (larger than any physical exponent) so the
                              first real value always wins.
          max_num_pot_alphas  Maximum number of Gaussian terms in any atom's
                              potential expansion.
          max_num_atom_alphas Maximum number of Gaussian terms in any atom's
                              wavefunction expansion.  get_num_terms_wf(z)
                              returns a list; only index [0] is the count of
                              wavefunction alphas.
          max_num_vale_states Maximum number of valence states for any single
                              atom, counting full angular-momentum degeneracy
                              (2l+1: s=1, p=3, d=5, f=7) and summing over all
                              three basis sets: MB (minimal), FB (full), EB
                              (extended), corresponding to basis indices 1, 2, 3.

        Called after all atoms and species are resolved so that self.atomic_z
        is fully populated.

        Note: the Perl original declared but never used local variables
        $element and $orbitalType — these have been omitted in Python.
        """
        ed = self.element_data
        # Sentinel: 1000.0 is larger than any real Gaussian exponent.
        self.min_pot_alpha = 1000.0
        self.max_num_pot_alphas = 0
        self.max_num_atom_alphas = 0
        self.max_num_vale_states = 0
        for atom in range(1, self.num_atoms + 1):
            z = self.atomic_z[atom]

            temp = ed.min_term_pot[z]
            if temp < self.min_pot_alpha:
                self.min_pot_alpha = temp

            temp = ed.num_terms_pot[z]
            if temp > self.max_num_pot_alphas:
                self.max_num_pot_alphas = temp

            # get_num_terms_wf returns a list; [0] is the wavefunction alpha count.
            temp_array = ed.get_num_terms_wf(z)
            if temp_array[0] > self.max_num_atom_alphas:
                self.max_num_atom_alphas = temp_array[0]

            # Sum valence states over MB+FB+EB (basis=1,2,3), weighted by
            # angular-momentum degeneracy 2l+1.
            temp_value = 0
            for basis in range(1, 4):  # MB+FB+EB
                temp_array = ed.get_vale_orbitals(basis, z)
                temp_value += 1 * temp_array[0]  # s  (2·0+1 = 1)
                temp_value += 3 * temp_array[1]  # p  (2·1+1 = 3)
                temp_value += 5 * temp_array[2]  # d  (2·2+1 = 5)
                temp_value += 7 * temp_array[3]  # f  (2·3+1 = 7)
            if temp_value > self.max_num_vale_states:
                self.max_num_vale_states = temp_value

    def compute_crystal_parameters(self):
        """Build an orthorhombic P1 cell that encloses all atoms.

        Used when reading non-crystalline formats (XYZ, PDB, HIN, etc.) that
        carry no lattice information.  The box is orthorhombic (all angles
        pi/2) but is deliberately assigned space group '1_a' / P1 so that it
        is treated as a general triclinic cell by the rest of the code — no
        symmetry operations will be applied.  Cell lengths are set to
        (max - min + buffer) along each Cartesian axis, where ``buffer`` is a
        user-configurable padding value (default: see ``set_buffer``).

        Call sequence (mirrors Perl):
          get_min_max_xyz → set space group → set angles → get_angle_sine
          → set mag → get_abc_vectors → make_inv_or_recip_lattice (inverse)
          → shift_xyz_center
        """
        # Find the maximum and minimum values of all atoms for each x,y,z axis.
        self.get_min_max_xyz(1, self.num_atoms)

        # We assume an orthorhombic box that contains all the atoms with a
        # buffer on all sides of the system.  Even though the box is
        # orthorhombic, we assign it space group 1_a (P1) so that no symmetry
        # operations are ever applied to it.
        self.space_group = '1_a'
        self.space_group_name = 'P1'
        self.space_group_num = 1
        self.space_group_sub_num = 1

        # Define the angles (all 90 degrees = orthorhombic).
        for axis in range(1, 4):
            self.angle[axis] = PI / 2.0
            self.angle_deg[axis] = 90.0

        # Compute the sine of each angle.
        self.get_angle_sine()

        # Define the cell dimensions: span of atom positions plus buffer.
        # The extra face_safety_margin (0.01 Angstroms per axis) prevents
        # atoms on the positive face of the bounding box from landing at
        # fractional coordinate exactly 1.0 after shift_xyz_center recentres
        # the atom cloud inside the cell.  Such an atom would otherwise be
        # folded back to fractional 0.0 by any subsequent apply_space_group
        # call (even for P1), which detaches the atom from its molecule and
        # corrupts later per-molecule bounding box computations.  The fold
        # tolerance inside applySpaceGroup is smallThresh*100 = 1e-6 in
        # fractional units, so the absolute safety margin must satisfy
        # (margin / 2) / mag > 1e-6; using 0.01 Angstroms keeps us well above
        # that threshold for cells up to several thousand Angstroms.
        face_safety_margin = 0.01
        for axis in range(1, 4):
            self.mag[axis] = (self.max_pos[axis] - self.min_pos[axis]
                              + self.buffer + face_safety_margin)

        # Define the real lattice parameters (a,b,c) in x,y,z vector form.
        self.get_abc_vectors()

        # Generate the inverse of the real-space lattice.
        # (Perl originally called makeLatticeInv; renamed makeInvOrRecipLattice.)
        self.make_inv_or_recip_lattice(make_recip=False)

        # Demand that there be no atoms with negative positions and that the
        # system be centered inside the cell.
        self.shift_xyz_center(1, self.num_atoms)

    def get_element_number(self, element_name):
        """Return the atomic number Z for the given element abbreviation.

        Perl name: getElementNumber

        Prerequisite: the link to the element data base must already be
        established (i.e. ``self.element_data`` populated, typically via
        ``setup_data_base_ref``).  The Perl comment reads: "This subroutine
        assums that the link to the element data base has been made."

        Mechanism: the Perl implementation performs a linear scan through
        ``elementNames_ref[1..N]`` (1-indexed), comparing each name against
        ``givenElement``, and returns the loop index on a match.  Because
        ``elementNames_ref`` is indexed from 1 and ordered by atomic number,
        that index is exactly Z.  Python delegates the same scan to
        ``element_data.get_element_z()``, which preserves the same semantics.

        Note: the Perl code had ``return $element; last;`` — the ``last``
        (break) is unreachable after ``return`` but was written to signal the
        intent to stop scanning once a match is found.

        Parameters
        ----------
        element_name : str
            Lower-case element abbreviation (e.g. 'si', 'o').

        Returns
        -------
        int
            Atomic number Z (same as the 1-based index in the ElementData
            tables), or None if not found.
        """
        return self.element_data.get_element_z(element_name)

    def assign_coval_radii(self):
        """Snapshot covalent radii from ElementData into a flat per-atom array.

        Copies ``element_data.coval_radii[Z]`` for each atom's atomic number Z
        into ``self.coval_radii[atom]``, producing a 1-indexed list of length
        ``num_atoms + 1`` (index 0 is the None placeholder).  Units: Angstroms.

        Bonding-factor note: ``ElementData`` exposes ``apply_bond_factor(f)``
        which scales every entry in its internal coval_radii table by *f*
        in-place.  If that was called before this method, the scaled values are
        snapshotted here.  The Perl bonding code notes: "the bonding factor was
        already included when we obtained the covalent radii."

        The bonding criterion used by ``obtain_atomic_interaction`` is:
            distance <= coval_radii[atom1] + coval_radii[ext_to_central[atom2]]
        so two atoms are considered bonded when their inter-atomic distance is
        within the sum of their (optionally scaled) covalent radii.

        Called by ``create_bonding_list``, ``create_pore_list``, and
        ``create_coordination_list`` immediately before
        ``obtain_atomic_interaction``.
        """
        self.coval_radii = [None] * (self.num_atoms + 1)
        for atom in range(1, self.num_atoms + 1):
            z = self.atomic_z[atom]
            self.coval_radii[atom] = self.element_data.coval_radii[z]

    def assign_color_vtk(self):
        """Snapshot per-atom RGBA colors from ElementData into a flat array.

        Copies ``element_data.color_vtk[Z]`` for each atom's atomic number Z
        into ``self.color_vtk[atom]``, producing a 1-indexed list of length
        ``num_atoms + 1`` (index 0 is the None placeholder).

        Each entry is a 4-component RGBA list [r, g, b, a].  The Perl
        implementation copies the four components explicitly as indices [0]..[3];
        Python uses ``list()`` to copy all components at once.

        The Perl calling context describes the purpose as: "Create a color map
        in case the model needs to be plotted."  The color data is used by
        downstream visualization code (e.g. VTK-based renderers) to assign
        element-specific colors to atoms.

        Called by ``create_bonding_list``, ``create_ext_bonding_list``, and
        ``create_q_list`` immediately before ``obtain_atomic_interaction``.
        """
        self.color_vtk = [None] * (self.num_atoms + 1)
        for atom in range(1, self.num_atoms + 1):
            z = self.atomic_z[atom]
            self.color_vtk[atom] = list(self.element_data.color_vtk[z])

    # ======================================================================
    # Section: Analysis
    # ======================================================================

    def create_min_dist_matrix(self, num_items1=None, items1_abc=None,
                               items1_xyz=None, num_items2=0,
                               items2_abc=None, items2_xyz=None):
        """Compute the minimum pairwise distance matrix.

        Two optional item groups can be provided and are joined internally
        before computing extended positions and distances.  When called with
        no arguments (the common case from apply_filter) all atoms in the
        model are used as the single group.

        Populates self.min_dist and self.self_min_dist via
        obtain_atomic_interaction (interaction flag 1).

        Note on Perl parameter ordering: the Perl signature interleaves the
        two groups differently — (numItems1, numItems2, items1ABC_ref,
        items2ABC_ref, items1XYZ_ref, items2XYZ_ref) — all ABC refs before
        all XYZ refs.  The Python signature uses the more natural per-group
        interleaving (num_items1, items1_abc, items1_xyz, ...).

        Parameters
        ----------
        num_items1 : int, optional
            Number of items in the first group.  Defaults to self.num_atoms.
        items1_abc : list, optional
            1-indexed fractional coordinates for group 1.
            Defaults to self.fract_abc.
        items1_xyz : list, optional
            1-indexed direct-Cartesian coordinates for group 1.
            Defaults to self.direct_xyz.
        num_items2 : int, optional
            Number of items in the second group (0 if omitted).
        items2_abc : list, optional
            1-indexed fractional coordinates for group 2.
        items2_xyz : list, optional
            1-indexed direct-Cartesian coordinates for group 2.
        """
        # Defaults: operate on all atoms.
        if num_items1 is None:
            num_items1 = self.num_atoms
            items1_abc = self.fract_abc
            items1_xyz = self.direct_xyz

        # Join the two lists into a single combined 1-indexed list.
        items_total_abc = [None]
        items_total_xyz = [None]
        num_items_total = 0

        for item in range(1, num_items1 + 1):
            num_items_total += 1
            items_total_abc.append(items1_abc[item])
            items_total_xyz.append(items1_xyz[item])

        if num_items2 and items2_abc and items2_xyz:
            for item in range(1, num_items2 + 1):
                num_items_total += 1
                items_total_abc.append(items2_abc[item])
                items_total_xyz.append(items2_xyz[item])

        # Determine the number of cells in each direction needed to account
        #   for all intercell interactions.
        self.compute_num_cells(num_items_total, items_total_abc)

        # Compute the position of each item in each of the non-origin cells
        #   where it is needed because it interacts with an item in the
        #   origin cell.
        self.compute_extended_pos(num_items_total, items_total_abc)

        # Compute and record the interatomic distances between all items.
        #   Flag 1 is used to record the min_dist matrix from the general
        #   purpose obtain_atomic_interaction method.
        self.obtain_atomic_interaction(
            1, num_items_total, self.num_items_ext,
            items_total_xyz, self.ext_direct_xyz_list)

    # ------------------------------------------------------------------
    # Section: Pair-filtering helpers (used by compute_num_cells /
    #          initialize_interaction_data)
    # ------------------------------------------------------------------

    def _border_coords(self):
        """Return the per-atom coordinate array matching border_coord_type."""
        ct = getattr(self, 'border_coord_type', 'xyz')
        if ct == 'xyz':
            return self.direct_xyz
        elif ct == 'abc':
            return self.direct_abc
        return self.fract_abc

    def item_out_of_bounds(self, item, coords=None):
        """Determine if this item should be excluded because it is outside
        the defined boundaries (if any were defeined).

        In the case of a calculation with a constraint on the positions of the
        items (e.g. RPDF), items other than those either inside or outside the
        defined borders are ignored, as requested via border_zone.

        border_zone == 1: exclude atoms *outside* the defined borders
            (keep only the atoms inside the box).
        border_zone == 2: exclude atoms *inside* the defined borders
            (keep only the atoms outside the box).
        border_zone == 0: no constraint; never skip (Python-only guard, not
            present in Perl — Perl callers check border_zone before calling).

        Parameters
        ----------
        item : int
            1-indexed item number.
        coords : list, optional
            1-indexed coordinate array to test against.  If None, the
            appropriate attribute (direct_xyz, direct_abc, or fract_abc) is
            chosen according to self.border_coord_type.

        Returns
        -------
        bool
            True if the item should be skipped, False if it should be kept.
        """
        if self.border_zone == 0:
            return False
        if coords is None:
            coords = self._border_coords()
        if item >= len(coords) or coords[item] is None:
            return False

        if self.border_zone == 1:   # Exclude atoms outside defined borders.
            # Assume that we will not skip this item.
            for axis in range(1, 4):
                v = coords[item][axis]
                if (v < self.border_low[axis] - EPSILON or
                        v > self.border_high[axis] + EPSILON):
                    return True
            return False
        else:                        # border_zone == 2: Exclude atoms inside defined borders.
            # Assume that we will skip this item.
            for axis in range(1, 4):
                v = coords[item][axis]
                if (v < self.border_low[axis] - EPSILON or
                        v > self.border_high[axis] + EPSILON):
                    return False
            return True

    def item_not_requested(self, item, tag_list, request_list):
        """Determine if this item should be excluded because it is not labeled
        with the requested tag (if such a request was made).

        Returns True (skip) when the item's tag is absent from request_list,
        and False (do not skip) when request_list is empty or the tag matches.

        In Perl the empty-list test is ``$#$requestList < 0`` (last index < 0),
        which is equivalent to ``len(request_list) <= 1`` here because the list
        is 1-indexed (index 0 is always a None placeholder).

        # In the case of a calculation with a constraint on the types of items
        # to consider, we here ignore items other than the requested ones.

        Parameters
        ----------
        item : int
            1-indexed item index.
        tag_list : list
            1-indexed list of tag strings (e.g. atom_element_name).
        request_list : list
            1-indexed (or empty) list of allowed tags.  If empty or contains
            only a None placeholder, every item is accepted.

        Returns
        -------
        bool
            True  → skip this item (not requested).
            False → include this item.
        """
        # Assume that we will skip this item.
        # If the request list has negative length (i.e. is empty) then every
        # item is not skipped.
        if not request_list or len(request_list) <= 1:
            return False
        tag = tag_list[item] if item < len(tag_list) else None
        # If the tag of the current item in the tag list is present in the
        # list of requested tags then we will not skip this item.
        for idx in range(1, len(request_list)):
            if tag == request_list[idx]:
                return False
        return True

    def pair_not_of_element(self, item1, item2):
        """Determine if this atom pair should be excluded because they are not
        of the correct element pair (if such a request was made).

        Parameters
        ----------
        item1, item2 : int
            Central-cell atom indices (1-indexed).

        Returns
        -------
        bool
            True if the pair should be skipped; False otherwise.
        """
        # Assume that we will not skip the item2.
        skip_item = False

        # In the case of an RPDF calculation with a constraint on the types of
        # elements to consider, we here ignore atom pairs that do not exactly
        # match the requested pair.
        if self.rpdf_select == 1:
            n1 = self.atom_element_name[item1]
            n2 = self.atom_element_name[item2]
            s1 = self.rpdf_select_atoms[1]
            s2 = self.rpdf_select_atoms[2]
            # Skip to the next j if we do not have a match between
            # (i name and atom1) + (j name and atom2) or
            # (i name and atom2) + (j name and atom1).
            if not ((n1 == s1 and n2 == s2) or (n1 == s2 and n2 == s1)):
                skip_item = True

        return skip_item

    def create_bonding_list(self):
        """Build the bonding list for the central cell using covalent radii.

        This is the primary bonding routine for most analysis tasks.  It
        orchestrates six steps:

        1. Determine how many periodic-image cells are needed in each direction
           to capture all interactions within limit_dist (compute_num_cells).
        2. Compute Cartesian positions of every atom in those image cells
           (compute_extended_pos → populates num_items_ext and
           ext_direct_xyz_list).
        3. Assign each atom's covalent radius from the element database
           (assign_coval_radii).
        4. Assign per-element VTK colours in case the model is visualised
           (assign_color_vtk).
        5. Walk all central-atom / extended-atom pairs; record those within
           the sum of covalent radii as bonds.  Flag 2 tells
           obtain_atomic_interaction to record extended bonding data in
           bonded_ext and bond_length_ext (obtain_atomic_interaction).
        6. Fold the extended bond partners back to their central-cell
           equivalents, populating num_bonds and bonded (map_ext_to_central).

        Contrast with create_ext_bonding_list, which skips step 6 and
        therefore leaves bonds expressed in terms of extended-cell indices
        rather than central-cell indices.

        Populates: bonded_ext, bond_length_ext, num_bonds, bonded.
        """
        # Determine the number of cells in each direction needed to account
        # for all intercell interactions.
        self.compute_num_cells()

        # Compute the position of each atom in each of the non-origin cells
        # where it is needed because it interacts with an atom in the origin
        # cell.
        self.compute_extended_pos()

        # Assign the covalent radius of each atom from the database.
        self.assign_coval_radii()

        # Create a color map in case the model needs to be plotted.
        self.assign_color_vtk()

        # Compute and record which atoms are bonded to which other atoms.
        # Flag #2 is used to record extended bonding information from the
        # general purpose obtain_atomic_interaction routine.
        self.obtain_atomic_interaction(
            2, self.num_atoms, self.num_items_ext,
            self.direct_xyz, self.ext_direct_xyz_list)

        # Map the extended atom positions back to the original central cell.
        self.map_ext_to_central()

    def compute_pore_map(self, resolution, num_mesh_points):
        """Compute a volumetric map of accessible pore space.

        Marks every grid point in ext_pore_map that falls within the covalent
        radius of any atom (including periodic images) with 1.  Grid points
        outside all atomic spheres remain 0.

        The parameters use Cartesian xyz axes, not fractional abc axes ("not
        abc yet" in the Perl).  The caller (e.g. the porosity script) converts
        abc-based mesh specifications to xyz before calling this routine.

        Parameters
        ----------
        resolution : list
            Step size in angstroms along [None, x, y, z] (1-indexed).
        num_mesh_points : list
            Number of mesh points along [None, x, y, z] (1-indexed).

        Algorithm
        ---------
        1. compute_num_cells — determine how many periodic images are needed in
           each +/− direction to capture all intercell interactions.
        2. compute_extended_pos — build the full extended-cell atom list
           (num_items_ext entries) with Cartesian positions.
        3. For each extended atom:
           a. Look up atomic Z via ext_to_central_item_map (the extended list
              may contain images; the central-cell copy carries the element).
           b. Retrieve the already-scaled covalent radius and square it.
           c. Find the nearest grid-index (index_num) along each axis by
              truncating position/resolution toward zero (int()).
           d. Compute mini_cube_size: the covalent radius divided by the grid
              spacing, rounded to the nearest integer.  This defines a cubic
              bounding box that *circumscribes* the atomic sphere — grid points
              outside this cube are guaranteed to be outside the sphere too.
           e. Triple loop (x, y, z) over the bounding cube.  Only consider
              mesh points inside the central system cell (1..num_mesh_points).
           f. Compute the squared Euclidean distance between the mesh point
              (in angstroms: index * resolution) and the atom position.
              The Perl code had a faster Inline-C helper ``CdistanceSqrd`` that
              was commented out, leaving a pure-Perl fallback; that fallback
              included a ``limitDist`` early-exit (if any single coordinate
              difference exceeded limit_dist, it set distanceSqrd = BIG_REAL
              and broke out of the axis loop).  Python uses the direct
              three-component squared distance without that optimisation.
           g. If dist_sqrd < covalent_radius², mark the grid point as 1.

        Notes
        -----
        - ext_pore_map is a 3-D array indexed [x][y][z], 1-based, pre-allocated
          by the caller before invoking this method.
        - The triple loop in the Perl source is written "without indentation"
          (all three for-loops at the same indent level) as a deliberate style
          choice to keep the inner body visually unindented.  Python uses normal
          nesting.
        """
        self.compute_num_cells()
        self.compute_extended_pos()

        ed = self.element_data

        for atom in range(1, self.num_items_ext + 1):
            # Map extended atom back to its central-cell copy to get element Z.
            central = self.ext_to_central_item_map[atom]
            z = self.atomic_z[central]
            curr_coval_radius = ed.coval_radii[z]
            curr_coval_radius_sqrd = curr_coval_radius * curr_coval_radius

            # Find the nearest grid index for this atom in each axis by
            # truncating position/resolution toward zero.
            index_num = [None, 0, 0, 0]
            # mini_cube_size circumscribes a sphere of size equal to the
            # covalent radius; rounded to the nearest grid interval.
            mini_cube_size = [None, 0, 0, 0]
            for axis in range(1, 4):
                index_num[axis] = int(
                    self.ext_direct_xyz_list[atom][axis] / resolution[axis])
                raw = curr_coval_radius / resolution[axis]
                if (raw - int(raw)) > 0.5:
                    mini_cube_size[axis] = int(raw) + 1
                else:
                    mini_cube_size[axis] = int(raw)

            # Triple loop over x,y,z of the mini-cube dimensions.
            # Only consider mesh points inside the central system cell.
            for x_point in range(index_num[1] - mini_cube_size[1],
                                 index_num[1] + mini_cube_size[1] + 1):
                if x_point < 1 or x_point > num_mesh_points[1]:
                    continue
                for y_point in range(index_num[2] - mini_cube_size[2],
                                     index_num[2] + mini_cube_size[2] + 1):
                    if y_point < 1 or y_point > num_mesh_points[2]:
                        continue
                    for z_point in range(index_num[3] - mini_cube_size[3],
                                         index_num[3] + mini_cube_size[3] + 1):
                        if z_point < 1 or z_point > num_mesh_points[3]:
                            continue

                        # Compute the squared distance between this mesh point
                        # (position = index * resolution, in angstroms) and the
                        # atom.  coord_diff per axis corresponds to the Perl
                        # @coordDiff array used by its CdistanceSqrd fallback.
                        dx = (x_point * resolution[1]
                              - self.ext_direct_xyz_list[atom][1])
                        dy = (y_point * resolution[2]
                              - self.ext_direct_xyz_list[atom][2])
                        dz = (z_point * resolution[3]
                              - self.ext_direct_xyz_list[atom][3])
                        dist_sqrd = dx*dx + dy*dy + dz*dz

                        if dist_sqrd < curr_coval_radius_sqrd:
                            self.ext_pore_map[x_point][y_point][z_point] = 1

    def create_ext_bonding_list(self):
        """Build the extended-cell bonding list for coordination/RPDF/Ylm analysis.

        This is the counterpart to create_bonding_list (flag 1).  Where
        create_bonding_list only considers the origin cell, this routine
        considers every periodic image that can interact with an atom in the
        origin cell, enabling analysis that requires neighbour shells beyond
        the immediate unit cell (e.g. RPDF, extended coordination, Ylm/BOO).

        Step 1 – determine the number of periodic images (cells) needed in
        each direction (+a/−a, +b/−b, +c/−c) to capture all intercell
        interactions up to the cutoff distance.  The results are stored in
        neg_bit and pos_bit (per axis).

        Step 2 – compute the Cartesian position of every atom in every
        non-origin periodic image that can interact with any origin-cell atom.
        Results populate ext_direct_xyz_list, ext_fract_abc_list,
        central2ext_item_map, ext2central_item_map, and num_atoms_ext.

        Step 3 – assign the covalent radius of each element from the database
        (prerequisite for the distance-based bond criterion).

        Step 4 – assign a VTK colour to each element in case the model is
        visualised later.

        Step 5 – find and record all bonded pairs using flag=2, which directs
        obtain_atomic_interaction to record extended (inter-image) bonding
        information rather than the origin-only bonding used by flag=1.
        """
        # Determine the number of cells in each direction needed to account
        # for all intercell interactions.
        self.compute_num_cells()

        # Compute the position of each atom in each non-origin cell where it
        # interacts with an atom in the origin cell.
        self.compute_extended_pos()

        # Assign the covalent radius of each atom from the database.
        self.assign_coval_radii()

        # Create a colour map in case the model needs to be plotted.
        self.assign_color_vtk()

        # Compute and record which atoms are bonded to which other atoms.
        # Flag 2 records extended bonding information from the general-purpose
        # obtain_atomic_interaction routine (vs. flag 1 for origin-cell only).
        self.obtain_atomic_interaction(2, self.num_atoms, self.num_atoms_ext,
                                       self.direct_xyz, self.ext_direct_xyz_list)

    def compute_rpdf(self):
        """Compute the radial pair distribution function g(r).

        Called by the ``rpdf`` script after the caller has optionally set:
          - ``set_limit_dist(d)``      — max distance (default 10.0 Å)
          - ``set_border(zone, ...)``  — fractional-ABC bounding box filter
          - ``set_rpdf_select_atoms(el1, el2)`` — element-pair filter

        Three-step pipeline (mirrors the Perl body exactly):

        1. **compute_num_cells** — determines per-atom neg_bit/pos_bit: how
           many periodic images in each ±a/b/c direction must be searched.
           Uses limit_dist / (mag[axis] * sin_angle[axis]) to extend the
           cutoff for oblique cells (hexagonal etc.).  Atoms outside the
           bounding box or not of the requested element type receive zero
           bits and contribute no extended images.

        2. **compute_extended_pos** — populates ext_direct_xyz_list with the
           Cartesian positions of every needed periodic image and builds the
           ext2central / central2ext index maps.  Returns num_atoms_ext.

        3. **obtain_atomic_interaction(flag=3)** — three internal phases:

           *Initialize (type 3)*:
             - skip_item1[i] set by itemOutOfBounds against borderCoords
               (central-cell bounding box filter; always False if no -limit).
             - skip_item2[j] set by itemOutOfBounds against borderCoordsExt
               (extended-image bounding box filter).
             - skip_pair[i][j] set by pairNotOfElement: True unless the
               element names of atom i and the central-cell image of atom j
               both match the two requested elements (or no filter was set).
               Flag value 1 signals compute_rpdf to filter by element pair
               rather than the covalent-radius bonding criterion used by flag
               2.  When set_rpdf_select_atoms is not called, all pairs pass.
             - rpdf[1..1000] initialised to 0.0.

           *Record loop*:
             For every pair (i, j) with distance < limit_dist and
             int(distance * 100) != 0 (i.e. non-zero, avoids self-image
             artefacts at r=0):
               rpdf[int(distance * 100)] += 1
             The factor of 100 maps Ångströms onto a 0.01 Å grid of 1000
             bins covering 0.01–10.00 Å.  A commented-out alternative in the
             Perl code accumulated raw distances and broadened the full list
             afterwards; that approach was abandoned in favour of the simpler
             bin-then-broaden scheme.

           *Finalize (type 3)*:
             - Apply Gaussian broadening with sigma=0.01 to the raw integer
               histogram stored in rpdf[1..1000], producing rpdfBroadened.
             - Divide each point by num_atoms * r² where r = point * 0.01 Å:
                 rpdf[point] = rpdfBroadened[point] / (num_atoms * r * r)
               The 1/r² factor corrects for the spherical-shell volume growth
               so the result is the true pair distribution function g(r).

        Results stored in self.rpdf[1..1000] (1-indexed, index 0 = None).
        The ``rpdf`` script then prints each point as:
          "%-6.3f %10.5f\\n" % (point/100.0, rpdf[point])
        """
        # Determine the number of cells in each direction needed to account
        # for all inter-cell interactions.
        self.compute_num_cells()

        # Compute the position of each atom in each non-origin cell where it
        # is needed because it interacts with an atom in the origin cell.
        self.num_atoms_ext = self.compute_extended_pos()

        # Compute and record the RPDF.  Flag 3 specifies RPDF mode for the
        # general-purpose obtain_atomic_interaction routine.
        self.obtain_atomic_interaction(3, self.num_atoms, self.num_atoms_ext,
                                       self.direct_xyz,
                                       self.ext_direct_xyz_list)

    def create_q_list(self):
        """Build the neighbor list needed to compute Bond Orientational Order (BOO).

        BOO parameters Q^n (also written Q_l) characterize the local symmetry of
        an atom's coordination shell.  This method mirrors Perl's createQList and
        sets up the extended-atom geometry needed before the actual Q^n accumulation
        and normalization steps (init_ylm_coeff / accumulate_q_bar /
        q_bar_per_bond / q_bar_normalize / q_bar_correlation).

        Pipeline (five steps, matching Perl inline comments verbatim):

        1. compute_num_cells
           Determine the number of cells in each direction needed to account for
           all intercell interactions.

        2. compute_extended_pos
           Compute the position of each atom in each of the non-origin cells where
           it is needed because it interacts with an atom in the origin cell.

        3. assign_coval_radii
           Assign the covalent radius of each atom from the database.  The radii
           are used by obtain_atomic_interaction to decide which atom pairs form
           bonds.

        4. assign_color_vtk
           Create a color map in case the model needs to be plotted.

        5. obtain_atomic_interaction (flag=4)
           Compute and record BOO-specific information.  Flag 4 distinguishes
           this call from the bonding-list (flag=1), extended-bonding (flag=2),
           RPDF (flag=3), atom-mesh-distance (flag=5), and bond-mesh-distance
           (flag=6) variants of the same general-purpose subroutine.

        Prerequisites: set_limit_dist, set_ylm_l must be called before this
        method to configure the interaction cutoff and the l quantum number for
        the spherical harmonic expansion.
        """
        # Determine the number of cells in each direction needed to account for
        # all intercell interactions.
        self.compute_num_cells()

        # Compute the position of each atom in each of the non-origin cells
        # where it is needed because it interacts with an atom in the origin cell.
        self.compute_extended_pos()

        # Assign the covalent radius of each atom from the database.
        self.assign_coval_radii()

        # Create a color map in case the model needs to be plotted.
        self.assign_color_vtk()

        # Compute and record the BOO.  Flag #4 is used to record BOO specific
        # information for this general purpose subroutine.
        self.obtain_atomic_interaction(
            4, self.num_atoms, self.num_items_ext,
            self.direct_xyz, self.ext_direct_xyz_list)

    def compute_atom_mesh_dist(self):
        """Compute distances from each atom (of the requested element) to scan mesh points.

        Three-step pipeline (mirrors Perl computeAtomMeshDist):

        1. compute_num_cells — determine the number of supercell replicas needed
           in each direction to account for all inter-cell interactions given the
           current limit distance.
        2. compute_extended_pos — compute the position of each atom in each of the
           non-origin cells where it is needed because it interacts with an atom in
           the origin cell.
        3. obtain_atomic_interaction(flag=5) — compute and record the distance
           between each atom (of the requested element) and the points of the
           defined scan line or mesh.  Flag 5 selects this action from the general-
           purpose obtain_atomic_interaction dispatcher.

        Prerequisites: set_limit_dist, set_scan_points (or set_scan_element) must
        be called before invoking this method.
        """
        # Determine the number of cells in each direction needed to account for
        # all intercell interactions.
        self.compute_num_cells(self.num_atoms, self.fract_abc)

        # Compute the position of each atom in each of the non-origin cells
        # where it is needed because it interacts with an atom in the origin cell.
        self.compute_extended_pos(self.num_atoms, self.fract_abc)

        # Compute and record the distance between each atom (of the requested
        # element) and the points of the defined scan line or mesh.  Flag 5
        # requests this type of action from the general-purpose subroutine.
        self.obtain_atomic_interaction(
            5, self.num_scan_points, self.num_items_ext,
            self.scan_points, self.ext_direct_xyz_list)

    def compute_bond_mesh_dist(self):
        """Compute the distance between each atom and the defined scan line or mesh.

        Identical pipeline to compute_atom_mesh_dist except that num_bonds_total
        (not num_atoms) is passed to compute_num_cells.  Using the bond count
        sizes the periodic-image search radius to cover all bond positions rather
        than just the raw atom count.

        Steps (matching Perl inline comments verbatim):
          1. Determine the number of cells in each direction needed to account
             for all intercell interactions.
          2. Compute the position of each atom in each of the non-origin cells
             where it is needed because it interacts with an atom in the origin
             cell.
          3. Compute and record the distance between each atom (of the requested
             element) and the points of the defined scan line or mesh.  Flag #5
             requests this type of action from the general-purpose subroutine.

        Prerequisites: set_limit_dist, set_scan_points (or set_scan_element)
        must be called before this method.
        """
        # Determine the number of cells in each direction needed to account for
        #   all intercell interactions.
        self.compute_num_cells(self.num_bonds_total, self.fract_abc)

        # Compute the position of each atom in each of the non-origin cells where
        #   it is needed because it interacts with an atom in the origin cell.
        self.num_atoms_ext = self.compute_extended_pos(self.num_atoms,
                                                       self.fract_abc)

        # Compute and record the distance between each atom (of the requested
        #   element) and the points of the defined scan line or mesh.  The #5 is
        #   the flag requesting this type of action from the general purpose
        #   subroutine.
        self.obtain_atomic_interaction(
            5, self.num_scan_points, self.num_atoms_ext,
            self.scan_points, self.ext_direct_xyz_list)

    def compute_num_cells(self, num_items=None, fract_abc=None):
        """Determine how many cell replications are needed (neg_bit, pos_bit).

        Computes per-item neg_bit[item][axis] and pos_bit[item][axis] — the
        number of periodic images needed in the negative / positive direction
        along each lattice axis so that all neighbours within limit_dist are
        covered.

        To make a significant efficiency gain we determine which periodic cells
        each item needs to be searched for in.  Consider an item at the origin
        of a large cell.  Other items in the origin cell will search for this
        item in the origin cell (x0) and at the +x1 cell, but not at the -x1
        cell.  Items near the origin item are still too far from the replication
        of the origin item in the -x1 cell to need to count it.  However, items
        in the origin cell but near the +x1 border will want to include the
        origin item in the +x1 replicated cell when they search for interacting
        items.  This means we need the position of the origin item in the origin
        cell and in the +x1 cell, but not in the -x1 cell.

        The difficult aspect is making this work for arbitrary lattice
        parameters.  Consider a hexagonal cell and an item near one border.  We
        cannot just add the limit distance in the a direction to see if the atom
        needs to be searched for in the -a1 cell (this test is counter-intuitive:
        if original position + limit_dist crosses into the +a1 cell then this
        item must be searched for in the -a1 cell).  The a,b,c position +
        limit_dist in the a direction may not penetrate into the next cell, but
        the a,b,c position + limit_dist perpendicular to the wall might.

        Hexagonal cell schematic (c is out of the screen)::

               b-----
              /    ./ 
             /     /
            c-----a

        The distance to the cell wall in the a direction is greater than the
        distance perpendicular to the wall.  If the limit distance is just the
        right length, an item would not be searched for when it should be.  To
        correct this, the limit distance is increased by a factor of 1/sin(angle)
        where the angle is alpha, beta, or gamma as appropriate when added to
        the a,b,c item positions.

        For each axis the two relevant angles must be considered and the minimum
        sine taken:
          - a-axis → beta (a,c) and gamma (a,b)
          - b-axis → alpha (b,c) and gamma (a,b)
          - c-axis → beta (a,c) and alpha (b,c)

        Parameters
        ----------
        num_items : int, optional
            Number of items to consider.  Defaults to self.num_atoms.
        fract_abc : list, optional
            1-indexed fractional coordinates.  Defaults to self.fract_abc.
        """
        if num_items is None:
            num_items = self.num_atoms
        if fract_abc is None:
            fract_abc = self.fract_abc

        # Get the minimum sine for angles involving each axis.
        # sine_angle[1]=alpha (b,c), sine_angle[2]=beta (a,c),
        # sine_angle[3]=gamma (a,b).
        sine_min = [None,
            # a-axis: beta, gamma  (a,c  a,b)
            min(self.sine_angle[2], self.sine_angle[3]),
            # b-axis: alpha, gamma (b,c  a,b)
            min(self.sine_angle[1], self.sine_angle[3]),
            # c-axis: beta, alpha  (a,c  b,c)
            min(self.sine_angle[2], self.sine_angle[1]),
        ]
        # Compute the maximum modification to the limit distance for each axis.
        ext_limit = [None] + [
            self.limit_dist / self.mag[ax] / sine_min[ax]
            for ax in range(1, 4)
        ]

        self.neg_bit = [None] * (num_items + 1)
        self.pos_bit = [None] * (num_items + 1)

        rpdf_sel = getattr(self, 'rpdf_select_atoms', [])

        for item in range(1, num_items + 1):
            self.neg_bit[item] = [None, 0, 0, 0]
            self.pos_bit[item] = [None, 0, 0, 0]

            # Check if this position is outside the defined bounding box.
            if self.item_out_of_bounds(item):
                continue

            # Check if this position is one of the requested tags.
            # Always true if no specific tags were requested.
            if self.item_not_requested(item, self.atom_element_name, rpdf_sel):
                continue

            # Compute the number of cells in each direction (a,b,c).
            for axis in range(1, 4):
                fa = fract_abc[item][axis]
                self.neg_bit[item][axis] = math.floor(fa + ext_limit[axis])
                self.pos_bit[item][axis] = math.floor(
                    1.0 - fa + ext_limit[axis])
            # Each bit gives the number of periodic cells in that direction
            # that need to be searched in *by* other items *for* this item.

    def compute_extended_pos(self, num_items=None, fract_abc=None):
        """Generate the extended (replicated) atom list.

        For each item in the central cell, loop through all the necessary
        replicated periodic-image cells (determined per-item by neg_bit and
        pos_bit from compute_num_cells) and store the direct-space x,y,z
        coordinates of that item in each image cell.  The result is a flat
        list of every (item, image-cell) combination — the "extended" system.

        The Perl implementation also maintained a 5-D positional array indexed
        as extDirectXYZItem[xyz][cellBitA+negA][cellBitB+negB][cellBitC+negC][item],
        which was used to recover the xyz position of item i in a specific
        image cell without scanning the flat list.  Python uses only the flat
        ext_direct_xyz_list / ext_fract_abc_list approach (sufficient for all
        current callers).

        Fractional coordinate of an image cell entry:
            tf[a] = fract_abc[item][a] + cellBitA   (similarly for b, c)
        where cellBitA ranges over -neg_bit[item][1] .. +pos_bit[item][1].
        The image corresponding to cellBit == (0,0,0) is the central-cell
        copy of the item.

        Cartesian conversion uses real_lattice[abc][xyz] (both 1-indexed):
            tx[xyz] = tf[1]*rl[1][xyz] + tf[2]*rl[2][xyz] + tf[3]*rl[3][xyz]

        Two index maps are maintained (both 1-indexed):

        central_to_ext_item_map[item] — maps a central-cell item number to its
            position in the flat extended list.  E.g. item 1 in the central
            cell is "item 1 duh", but item 1 in the extended list may be a
            periodic image in the -negBit direction; this map records where in
            the extended list the actual origin-cell copy lives.
            (Index = original item number, value = extended item number.)

        ext_to_central_item_map[ext] — maps an extended-list item back to its
            original central-cell item number.  E.g. extended item 1 is some
            periodic replica of some central-cell atom; this records which one.
            (Index = extended item number, value = original central-cell item.)

        Populates
        ---------
        self.ext_direct_xyz_list : list of [None, x, y, z], 1-indexed
        self.ext_fract_abc_list  : list of [None, a, b, c], 1-indexed
        self.central_to_ext_item_map : list, 1-indexed by central-cell item
        self.ext_to_central_item_map : list, 1-indexed by extended item
        self.num_items_ext  : total count of extended positions
        self.num_atoms_ext  : alias kept in sync for backward-compat callers

        Parameters
        ----------
        num_items : int, optional
            Defaults to self.num_atoms.
        fract_abc : list, optional
            1-indexed fractional coordinates.  Defaults to self.fract_abc.

        Returns
        -------
        int
            Number of extended positions generated.
        """
        if num_items is None:
            num_items = self.num_atoms
        if fract_abc is None:
            fract_abc = self.fract_abc

        self.num_items_ext = 0
        self.ext_direct_xyz_list = [None]
        self.ext_fract_abc_list  = [None]
        self.ext_to_central_item_map  = [None]
        self.central_to_ext_item_map  = [None] * (num_items + 1)

        rl = self.real_lattice  # rl[abc][xyz], both 1-indexed

        # Loop through all the necessary replicated cells for each item and
        # store the direct-space x,y,z coordinates of the item in that cell.
        for item in range(1, num_items + 1):
            neg = self.neg_bit[item]
            pos = self.pos_bit[item]
            fa  = fract_abc[item]

            for ca in range(-neg[1], pos[1] + 1):
                for cb in range(-neg[2], pos[2] + 1):
                    for cc in range(-neg[3], pos[3] + 1):
                        # Fractional coords shifted by the image-cell offset.
                        tf = [None,
                              fa[1] + ca,
                              fa[2] + cb,
                              fa[3] + cc]

                        # Convert fractional → direct Cartesian xyz.
                        tx = [None, 0.0, 0.0, 0.0]
                        for xyz in range(1, 4):
                            tx[xyz] = (tf[1] * rl[1][xyz] +
                                       tf[2] * rl[2][xyz] +
                                       tf[3] * rl[3][xyz])

                        # Add this item to the extended (flat) list.
                        self.num_items_ext += 1
                        ext = self.num_items_ext
                        self.ext_direct_xyz_list.append(
                            [None, tx[1], tx[2], tx[3]])
                        self.ext_fract_abc_list.append(
                            [None, tf[1], tf[2], tf[3]])

                        # ext→central: which central-cell atom does this
                        # extended entry represent?
                        self.ext_to_central_item_map.append(item)

                        # central→ext: the (0,0,0) image is the central-cell
                        # copy; record its flat-list index for this item.
                        if ca == 0 and cb == 0 and cc == 0:
                            self.central_to_ext_item_map[item] = ext

        # Keep num_atoms_ext in sync for callers that use it.
        self.num_atoms_ext = self.num_items_ext
        return self.num_items_ext

    def map_ext_to_central(self):
        """Translate bond partners from extended-cell indices to central-cell indices.

        After ``create_bonding_list`` runs, ``self.bonded_ext[atom][bond]``
        holds the index of each bond partner in the *extended* (replicated)
        cell list.  ``self.ext_to_central_item_map`` maps every extended-cell
        item number back to the original central-cell atom number it was
        replicated from — i.e. the index is the extended item number and the
        value is the original central-cell atom number.

        This routine populates ``self.bonded[atom][bond]`` using that map, so
        that all bond-partner references become central-cell indices.  Central-
        cell indexing is often better for sharing data with other programs or
        for output where periodic-image information is not needed (see the
        Perl comment preceding ``obtainAtomicInteraction``).

        Both ``bonded_ext`` and ``bonded`` use 1-based indexing for atoms and
        bonds; index 0 is a None placeholder.  This method rebuilds
        ``self.bonded`` from scratch on each call (Perl overwrote in-place;
        Python reinitialises to keep index-0 None semantics clean).
        """
        self.bonded = [None] + [[None] for _ in range(self.num_atoms)]
        for atom in range(1, self.num_atoms + 1):
            for bond in range(1, self.num_bonds[atom] + 1):
                # bonded_ext[atom][bond] is the extended-cell item number;
                # ext_to_central_item_map converts it to the central-cell atom.
                ext_idx = self.bonded_ext[atom][bond]
                self.bonded[atom].append(
                    self.ext_to_central_item_map[ext_idx])

    def obtain_atomic_interaction(self, interaction_type,
                                   num_items1, num_items2,
                                   item1_xyz, item2_xyz):
        """Compute one of several types of pair interactions within limit_dist.

        The purpose of this subroutine is to compute one of a number of
        possible types of interactions between a set of points within the
        original cell (possibly the atom positions) and another set of points
        both within the original cell and in neighboring replicated cells.

        Interaction types
        -----------------
        1 — Minimal distance matrix.
            Computes the minimum distance between item pairs, including the
            effects of periodic boundary conditions.  At least one item from
            each pair must be from the central cell.

        2 — Bonding (extended-cell link list).
            Computes which items are "linked" to which other items and how many
            links each item has.  The index numbers for items in the central
            cell are used to index the links between items.  For example, if
            bonded_ext[1][4] == 77 it means that the fourth link (bond) for
            item #1 in the central-cell list is linked (bonded) to item #77
            from the extended-cell list.  This is often used to identify which
            pairs of atoms are bonded.  The data can be reduced so that all
            indices refer only to the central cell; sometimes one listing is
            better than another for sharing data with other programs or people.

        3 — Radial pair distribution function (RPDF).
            Computes the RPDF for the given set of items with the accompanying
            tag and space restrictions.

        4 — Bond orientational order (BOO).
            Used for the bond orientational order calculation.
            Not fully usable at the moment.

        5 — Scan distances.
            Computes distances between points in a scan line and atomic
            positions.

        Parameters
        ----------
        interaction_type : int
            1=min-dist matrix, 2=bonding, 3=RPDF,
            4=bond-orientational order, 5=scan distances.
        num_items1 : int
            Number of central-cell items.
        num_items2 : int
            Number of extended items.
        item1_xyz : list
            1-indexed direct-xyz coordinates for central items.
        item2_xyz : list
            1-indexed direct-xyz coordinates for extended items.

        Notes
        -----
        Distance pruning uses a two-phase strategy inherited from the Perl
        original.  The Perl code once had an optional Inline-C version
        (``Cdistance``) for speed; the pure-Perl fallback checked each Cartesian
        axis individually against ``limit_dist`` before computing the full
        Euclidean distance.  Python retains that same two-phase logic: a cheap
        per-axis ``abs(d) > lim`` early-out, followed by the squared-distance
        check against ``limit_dist_sqrd``.  In the C version a sentinel value of
        BIG_REAL was returned for out-of-range pairs; the Python port simply uses
        ``continue`` instead.

        The progress indicator (dots and pipe characters printed to stdout) is
        intentional: the double loop over item1 × item2 can take a while for
        large systems, so the running count keeps the user from wondering whether
        the process has stalled.
        """
        self.initialize_interaction_data(
            interaction_type, num_items1, num_items2)

        lim   = self.limit_dist
        lim2  = self.limit_dist_sqrd

        for item1 in range(1, num_items1 + 1):
            # Print a running count of progress so the user does not get too
            # bored waiting.  '|' every 10th item, '.' otherwise; atom number
            # printed and newline every 50th item.
            if item1 % 10 == 0:
                print('|', end='', flush=True)
            else:
                print('.', end='', flush=True)
            if item1 % 50 == 0:
                print(f' {item1}')

            if self.skip_item1[item1]:
                continue

            xyz1 = item1_xyz[item1]

            for item2 in range(1, num_items2 + 1):
                if self.skip_item2[item2]:
                    continue
                if self.skip_pair[item1][item2]:
                    continue

                xyz2 = item2_xyz[item2]
                dx = xyz1[1] - xyz2[1]
                dy = xyz1[2] - xyz2[2]
                dz = xyz1[3] - xyz2[3]

                # Phase 1: cheap per-axis pruning — if any single component
                # already exceeds limit_dist the Euclidean distance must also
                # exceed it, so skip immediately without computing the sqrt.
                if abs(dx) > lim or abs(dy) > lim or abs(dz) > lim:
                    continue

                # Phase 2: full squared-distance check before taking the sqrt.
                dist2 = dx * dx + dy * dy + dz * dz
                if dist2 >= lim2:
                    continue
                distance = math.sqrt(dist2)

                self.record_interaction_data(
                    interaction_type, distance,
                    [None, dx, dy, dz], item1, item2)

        if num_items1 % 50 != 0:
            print()
        self.finalize_interaction_data(
            interaction_type, num_items1, num_items2)

    def initialize_interaction_data(self, interaction_type,
                                     num_items1, num_items2):
        """Allocate and zero-initialise interaction data arrays for one
        pass through obtain_atomic_interaction.

        Called once per analysis run to reset all arrays before the main
        pairwise distance loop.  The five interaction types each require
        a different set of skip flags and accumulator arrays:

          1 — Minimal distance matrix.
              Finds the minimal distances between central-cell item pairs,
              accounting for periodicity via the extended cell.  Watch
              carefully for the distinction between
              ext_to_central_item_map[item2] (the central-cell index of an
              extended item) and just item2 (its extended index).  We must
              search all extended items but record results for central-cell
              items only.

          2 — Extended bonding list (check for extended bonding).
              No skipping; every central/extended pair is examined.
              Initialises bond counts for each central-cell atom.

          3 — Radial pair distribution function (RPDF).
              Atoms outside a user-defined bounding box are skipped.
              Pairs not matching the requested elemental types are also
              skipped.  The rpdf accumulator (1000 points) is zeroed so
              that it can be accumulated across the pairwise loop.

          4 — Bond orientational order (BOO / Ylm / q_bar).
              No skipping; every pair is examined.
              q_bar[item][m] = [real, imag] is the accumulating sum of
              spherical-harmonic values Ylm for each bonded neighbour.
              (Perl note: a complex module was available but gave trouble,
              so real/imag parts are stored separately as a "lame method".)
              Ylm coefficients are initialised via init_ylm_coeff().

          5 — Minimum distance from extended atoms to grid (scan) points.
              Scan-point items (item1) are never skipped.  Extended items
              (item2) are skipped when their central-cell atom's element
              does not match scan_element.  ext_scan_dist is initialised
              to BIG_REAL so the minimum can be found by comparison.

        Parameters
        ----------
        interaction_type : int
            1–5 as described above; raises ValueError for any other value.
        num_items1 : int
            Number of central-cell items (or scan-grid points for type 5).
        num_items2 : int
            Number of extended-cell items.

        Arrays set (all 1-indexed; index 0 is None or False placeholder):
            skip_item1[item1], skip_item2[item2], skip_pair[item1][item2]
            Type 1: self_min_dist[item], min_dist[item1][item2]
            Type 2: num_bonds[item], bonded_ext[item], bond_length_ext[item]
            Type 3: rpdf[point]  (1..1000)
            Type 4: num_bonds[item], bonded_ext[item],
                    q_bar[item][m_idx]  (m_idx in 1..2*l+1)
            Type 5: ext_scan_dist[item1][item2]
        """
        # skip_item1/2: 1-indexed bool lists (index 0 unused/False)
        # skip_pair:    skip_pair[item1][item2], both 1-indexed

        if interaction_type == 1:   # minimal distance matrix
            # Reset the minimal distance matrix and self-min-distance array
            # (Perl used `undef @selfMinDist; undef @minDist` to discard old
            # data before rebuilding; here we simply overwrite the attributes.)
            self.self_min_dist = [None] + [BIG_REAL] * num_items1
            self.min_dist = [
                [None] + [BIG_REAL] * num_items1
                for _ in range(num_items1 + 1)
            ]
            self.min_dist[0] = None

            # Don't skip any atoms, but skip pairs where the extended item2
            # maps back to a central atom whose index is less than item1.
            # This avoids double-counting symmetric pairs.
            self.skip_item1 = [False] * (num_items1 + 1)
            self.skip_item2 = [False] * (num_items2 + 1)
            self.skip_pair  = [
                [False] * (num_items2 + 1) for _ in range(num_items1 + 1)
            ]
            for i1 in range(1, num_items1 + 1):
                for i2 in range(1, num_items2 + 1):
                    if self.ext_to_central_item_map[i2] < i1:
                        self.skip_pair[i1][i2] = True

        elif interaction_type == 2:  # check for extended bonding
            # Don't skip any atoms or pairs of atoms.
            self.skip_item1 = [False] * (num_items1 + 1)
            self.skip_item2 = [False] * (num_items2 + 1)
            self.skip_pair  = [
                [False] * (num_items2 + 1) for _ in range(num_items1 + 1)
            ]
            # Initialize the count of the number of bonds for each atom
            # in the central cell.
            self.num_bonds      = [None] + [0] * num_items1
            # 1-indexed inner list grows as bonds are recorded:
            # bonded_ext[item] = [None, ext_idx1, ext_idx2, ...]
            self.bonded_ext     = [None] + [
                [None] for _ in range(num_items1)
            ]
            self.bond_length_ext = [None] + [
                [None] for _ in range(num_items1)
            ]

        elif interaction_type == 3:  # rpdf
            ct   = getattr(self, 'border_coord_type', 'xyz')
            cext = 'ext_' + ('direct_xyz_list'
                             if ct == 'xyz' else 'fract_abc_list')
            bc    = None if self.border_zone == 0 else self._border_coords()
            bcext = None if self.border_zone == 0 else getattr(self, cext, None)

            # Skip certain atoms if they are outside the defined bounding box.
            self.skip_item1 = [
                False if bc is None else self.item_out_of_bounds(i, bc)
                for i in range(num_items1 + 1)
            ]
            self.skip_item2 = [
                False if bcext is None else self.item_out_of_bounds(i, bcext)
                for i in range(num_items2 + 1)
            ]
            # Skip atom pairs if they are not of the requested elemental types.
            self.skip_pair = [
                [False] * (num_items2 + 1) for _ in range(num_items1 + 1)
            ]
            for i1 in range(1, num_items1 + 1):
                for i2 in range(1, num_items2 + 1):
                    self.skip_pair[i1][i2] = self.pair_not_of_element(
                        i1, self.ext_to_central_item_map[i2])
            # Initialize the rpdf values to zero so they can be accumulated.
            self.rpdf = [None] + [0.0] * 1000

        elif interaction_type == 4:  # bond orientational order
            # Don't skip any atoms or pairs.
            self.skip_item1 = [False] * (num_items1 + 1)
            self.skip_item2 = [False] * (num_items2 + 1)
            self.skip_pair  = [
                [False] * (num_items2 + 1) for _ in range(num_items1 + 1)
            ]
            # Initialize the count of the number of bonds for each atom.
            self.num_bonds  = [None] + [0] * num_items1
            self.bonded_ext = [None] + [[None] for _ in range(num_items1)]
            m = 2 * self.ylm_l + 1
            # Initialize q_bar for each atom.  q_bar is a complex
            # multi-component vector that accumulates the sum of the spherical
            # harmonics Ylm for each bonded neighbour.  Real and imaginary
            # parts are stored separately (index 1=real, 2=imag) to match
            # Perl's ``@qBar[$atom][$k][1..2]`` convention; slot 0 is an
            # unused sentinel so the module-wide 1-indexed style is kept
            # even at the innermost level.
            # q_bar[item][m_idx] = [None, real, imag];  m_idx in 1..2*l+1.
            self.q_bar = [None] + [
                [None] + [[None, 0.0, 0.0] for _ in range(m)]
                for _ in range(num_items1)
            ]
            # Initialize the coefficients for the Ylm.
            self.init_ylm_coeff()

        elif interaction_type == 5:  # min dist of ext items to grid positions
            # Do not skip any scan points (item1 = grid positions).
            self.skip_item1 = [False] * (num_items1 + 1)
            self.skip_item2 = [False] * (num_items2 + 1)
            # Skip extended items that are not of the requested scan element.
            for i2 in range(1, num_items2 + 1):
                central = self.ext_to_central_item_map[i2]
                self.skip_item2[i2] = (
                    self.atom_element_name[central] != self.scan_element)
            # Do not skip any pairs.
            self.skip_pair = [
                [False] * (num_items2 + 1) for _ in range(num_items1 + 1)
            ]
            # Initialize the distance between each mesh point and each
            # extended item to BIG_REAL so the minimum can be found.
            self.ext_scan_dist = [
                [None] + [BIG_REAL] * num_items2
                for _ in range(num_items1 + 1)
            ]
            self.ext_scan_dist[0] = None

        else:
            raise ValueError(
                f'Unknown interaction type {interaction_type}.')

    def record_interaction_data(self, interaction_type, distance,
                                 diff, item1, item2):
        """Store one pair interaction result dispatched by interaction_type.

        Called once per (item1, item2) pair by obtain_atomic_interaction.
        The five interaction types correspond to the five analysis modes:

          1 — Minimal distance matrix.  Populate self_min_dist (self-pairs)
              or min_dist[item1][central2] (cross-pairs).
          2 — Bonding.  If distance <= sum of covalent radii, record bond:
              increment num_bonds, append item2 to bonded_ext, append distance
              to bond_length_ext.
          3 — Radial pair distribution function (rpdf).  Bin the distance at
              integer index int(distance * 100); skip zero (same site).
          4 — Bond orientational order.  Same bonding criterion as type 2; on
              a bond, compute spherical angles (theta, phi) from the diff
              vector and accumulate q_bar for item1.  Bond length is not
              stored here — it may not be needed for BOO.
          5 — Minimum distance of extended items to mesh/grid positions.
              Store the scalar distance in ext_scan_dist[item1][item2].

        Parameters
        ----------
        interaction_type : int
            Which analysis mode is active (1–5).
        distance : float
            Scalar distance between item1 and item2.
        diff : list
            1-indexed [None, dx, dy, dz] displacement vector (item1 − item2),
            used only for type 4 (spherical-angle computation).
        item1 : int
            Central-cell item index (1-indexed).
        item2 : int
            Extended-cell item index (1-indexed).

        Notes
        -----
        IMPORTANT / TRICKY: item1 indexes into the central-cell atom list and
        item2 indexes into the *extended* (periodic) atom list.  These two
        lists are NOT aligned — atom 1 in the central cell does NOT correspond
        to atom 1 in the extended list.  Use ext_to_central_item_map[item2]
        (c2 below) whenever you need the central-cell counterpart of item2.
        """
        if interaction_type == 1:  # Minimal distance matrix.
            c2 = self.ext_to_central_item_map[item2]
            if item1 == c2:
                # Self-interaction: record the smallest non-zero distance for
                # item1 against its own periodic images.
                if distance != 0 and distance < self.self_min_dist[item1]:
                    self.self_min_dist[item1] = distance
            else:
                # Cross-pair: record the value in the minimal distance matrix.
                if distance < self.min_dist[item1][c2]:
                    self.min_dist[item1][c2] = distance

        elif interaction_type == 2:  # Check for bonding.
            c2 = self.ext_to_central_item_map[item2]
            # If the distance is <= the sum of the covalent radii of the two
            # atoms then we say these atoms are bonded.  Note that the bonding
            # factor was already included when we obtained the covalent radii.
            if (distance != 0 and
                    distance <= self.coval_radii[item1] +
                                 self.coval_radii[c2]):
                # Increment the number of bonds for the central-cell atom.
                self.num_bonds[item1] += 1
                # Store the extended-cell item number for this bond.
                # (Perl used 1-indexed array assignment; Python uses append.)
                self.bonded_ext[item1].append(item2)
                # Record the distance between these two items (bond length).
                self.bond_length_ext[item1].append(distance)

        elif interaction_type == 3:  # Radial pair distribution function.
            # Bin the distance at graph scale (graph point = 100 * distance);
            # skip index 0 (same-site, zero distance).
            # Note: an alternative approach (commented out in Perl) would push
            # raw distances and apply Gaussian broadening after the array is
            # filled, but the histogram approach is used instead.
            idx = int(distance * 100)
            if idx != 0:
                self.rpdf[idx] += 1

        elif interaction_type == 4:  # Bond orientational order.
            c2 = self.ext_to_central_item_map[item2]
            # Same bonding criterion as type 2.
            if (distance != 0 and
                    distance <= self.coval_radii[item1] +
                                 self.coval_radii[c2]):
                # Increment the number of bonds for the central-cell atom.
                self.num_bonds[item1] += 1
                # Store the extended-cell item number for this bond.
                # (Bond length is not recorded here — it may not be needed
                # for a bond orientational order calculation.)
                self.bonded_ext[item1].append(item2)
                # Get the spherical angles theta and phi from the diff vector.
                # diff is already 1-indexed [None, dx, dy, dz], so pass it
                # directly to the 1-indexed spherical_angles helper.
                theta, phi = self.spherical_angles(diff)
                # Accumulate the values of q_bar for the central-cell atom.
                self.accumulate_q_bar(item1, theta, phi)

        elif interaction_type == 5:  # Min dist of extended items to grid positions.
            # Record the distance between the mesh point and the extended item.
            self.ext_scan_dist[item1][item2] = distance

        else:
            raise ValueError(
                f"Unknown interaction type {interaction_type}.  Aborting.")

    def finalize_interaction_data(self, interaction_type,
                                   num_items1, num_items2):
        """Post-process raw interaction data collected by obtain_atomic_interaction.

        Called once after the full atom-pair loop in obtain_atomic_interaction
        completes.  The required post-processing depends on which kind of
        interaction was accumulated:

        interaction_type == 1  (makeinput minimal distance matrix)
            record_interaction_data filled only the upper-triangle
            (item1 < item2) of self.min_dist.  This step:
              • sets the diagonal to 0 (self-distance)
              • mirrors the upper triangle into the lower triangle so the
                matrix is fully symmetric.

        interaction_type == 2  (check for bonding)
            Nothing to do; bonding flags were set directly during
            record_interaction_data.

        interaction_type == 3  (radial pair distribution function, RPDF)
            self.rpdf[1..1000] holds raw pair counts binned at 0.01 Å
            intervals.  This step:
              • applies Gaussian broadening (sigma = 0.01) via
                gaussian_broaden to smooth the delta-function spikes
              • divides each broadened bin by num_items1 * r^2 to
                normalise for the increasing shell volume at larger r
                (the index-to-distance conversion is r = point * 0.01,
                so dividing by r^2 is equivalent to dividing by
                index^2 / 100^2 in the Perl formulation).

        interaction_type == 4  (bond orientational order, BOO)
            Finalises the q_bar vectors accumulated per bond:
              • q_bar_per_bond   — divides each atom's accumulated q_bar
                                   by its bond count
              • q_bar_normalize  — normalises each q_bar vector to unit
                                   length
              • q_bar_correlation — accumulates the dot-product correlation
                                   factor across bonded pairs

        interaction_type == 5  (minimum distance of extended items to grid
                                positions)
            Nothing to do; distances were recorded directly.

        Parameters
        ----------
        interaction_type : int
            Selects which post-processing branch to execute (1–5).
        num_items1 : int
            Number of primary items (atoms or elements) in the interaction.
        num_items2 : int
            Number of secondary items (used by some branches; ignored by
            others).
        """
        if interaction_type == 1:
            # Complete the minimal distance matrix for the self reference.
            for item in range(1, num_items1 + 1):
                self.min_dist[item][item] = 0

            # Complete the symmetric portion of the minimal distance matrix.
            for i1 in range(1, num_items1):
                for i2 in range(i1 + 1, num_items1 + 1):
                    self.min_dist[i2][i1] = self.min_dist[i1][i2]

        elif interaction_type == 2:
            pass  # nothing to do

        elif interaction_type == 3:
            # Apply Gaussian broadening to the raw RPDF points.  self.rpdf
            # is already a 1-indexed ``[None, v1, ..., v1000]`` list so it
            # can be handed straight to the 1-indexed ``gaussian_broaden``
            # without any slice-and-shift bridging.
            sigma = 0.01
            broadened = self.gaussian_broaden(self.rpdf, sigma)

            # Divide by r^2 to account for the effect of more atoms being
            # present at each greater distance.  To divide by r^2 we divide
            # by index^2 / 100^2 because the index numbers are 100 times
            # greater than the distances (r = point * 0.01).
            for point in range(1, 1001):
                r = point * 0.01
                self.rpdf[point] = broadened[point] / (num_items1 * r * r)

        elif interaction_type == 4:
            # Divide the q_bar for each atom by the number of bonds for
            # that atom.
            self.q_bar_per_bond()

            # Normalise each q_bar vector.
            self.q_bar_normalize()

            # Accumulate the q_bar correlation factor.
            self.q_bar_correlation()

        elif interaction_type == 5:
            pass  # nothing to do

        else:
            raise ValueError(
                f"Unknown interaction type {interaction_type}.  Aborting."
            )

    def create_coordination_list(self):
        """Build per-atom elemental coordination strings.

        For each atom in the central cell, constructs a compact coordination
        string of the form "<count><elem><count><elem>..." where elements are
        sorted in ASCII order (e.g. "2N1Si" means 2 nitrogen + 1 silicon
        neighbors).  The result is stored in self.coordination[atom].

        Bonding prerequisite logic (from Perl):
          - If num_bonds is empty (len==1 due to None placeholder), no bond
            data exists yet → call create_bonding_list() to compute from
            scratch (extended positions + distance thresholds).
          - If num_bonds is populated but bonded (central-cell mapping) is
            empty, the extended calculation exists but has not been mapped
            back to the central cell → call map_ext_to_central().

        The border-box test (item_out_of_bounds) skips atoms that lie outside
        the region the user defined via set_border / set_buffer.  The Perl
        comment: "Compare the position of this atom with the box the user
        defined (or the default).  This also considers if the user wanted only
        atoms inside the box or only outside.  The return value is 1 if the
        atom is out of the desired region and it is 0 if the atom is inside
        the desired region."

        The Perl implementation collected element names into @bondedElem then
        called getUniqueArray(\\@bondedElem, 0, 1) with flag 0 = compare by
        ASCII characters and flag 1 = sort the unique elements by those ASCII
        characters.  The Python dict + sorted() reproduces this exactly.

        Populates: self.coordination  (1-indexed list, index 0 = None)
        """
        # If the length of the num_bonds array equals zero (len==1: only the
        # None placeholder), then we have not yet computed bonds → do that now.
        # If however num_bonds is populated but bonded is empty, the extended
        # calculation exists but must be mapped back to the central cell.
        if len(self.num_bonds) == 1:
            self.create_bonding_list()
        elif len(self.bonded) == 1:
            self.map_ext_to_central()

        self.coordination = [None] + [''] * self.num_atoms

        for atom in range(1, self.num_atoms + 1):
            # Compare atom position against the user-defined border box.
            # Skip (return value 1 = out of desired region).
            if self.item_out_of_bounds(atom):
                continue

            # Construct a dict of element names for atoms bonded to `atom`.
            # Perl built @bondedElem (1-indexed, elem 0 = "" sentinel) then
            # called getUniqueArray; here a plain dict achieves the same result.
            elem_counts = {}
            for bond in range(1, self.num_bonds[atom] + 1):
                elem = self.atom_element_name[self.bonded[atom][bond]]
                elem_counts[elem] = elem_counts.get(elem, 0) + 1

            # Construct the coordination string: for each unique element (sorted
            # by ASCII character value, matching Perl getUniqueArray flag 1),
            # append "<count><elemName>".
            coord_str = ''
            for elem in sorted(elem_counts):
                coord_str += str(elem_counts[elem]) + elem
            self.coordination[atom] = coord_str

    def create_coordination_summary(self):
        """Create a coordination summary for each element in the system.

        Populates self.coordination_summary[1..num_elements] with a
        human-readable string describing the distribution of coordination
        environments across all atoms of each element.

        Each summary string has the form:
            '[count1]coordStr1  [count2]coordStr2  ...'
        where coordStr is the coordination tag produced by
        create_coordination_list (e.g. 'O4' for oxygen coordinated to four
        neighbours) and count is the number of atoms of this element that
        carry that coordination.

        Prerequisites:
            create_coordination_list must have been called (or will be called
            automatically on first entry when self.coordination is empty).

        Perl analogue: createCoordinationSummary (StructureControl.pm ~7403).
        The Perl version used getUniqueArray to extract unique coordination
        strings and their occurrence counts from a per-element list; Python
        replaces that with an inline dict-based count, which is equivalent.
        """
        # Check to be sure that the coordination data has been obtained.
        # self.coordination is initialised as [None], so length == 1 means
        # it has not yet been filled by create_coordination_list.
        if len(self.coordination) == 1:
            self.create_coordination_list()

        # Initialize a per-element list of coordination strings.  Index 0 is
        # the None placeholder (1-indexed convention used throughout OLCAO).
        # Append the coordination of each atom onto the coordination summary
        # for the associated element; skip atoms with empty coordination
        # strings (atoms outside the border or with no bonds).
        coord_by_elem = [None] + [[] for _ in range(self.num_elements)]
        for atom in range(1, self.num_atoms + 1):
            coord = self.coordination[atom]
            if coord == '':
                continue
            coord_by_elem[self.atom_element_id[atom]].append(coord)

        # Get the unique coordinations for each element along with the number
        # of times that each coordination occurs, then construct a string for
        # this element that describes the different coordinations it has.
        # Initialize the coordination summary for this element to be empty
        # before accumulating.
        self.coordination_summary = [None] + [''] * self.num_elements
        for element in range(1, self.num_elements + 1):
            counts = {}
            for coord in coord_by_elem[element]:
                counts[coord] = counts.get(coord, 0) + 1
            summary = ''
            for coord in sorted(counts):
                summary += f'[{counts[coord]}]{coord}  '
            self.coordination_summary[element] = summary

    def compute_qn(self, bridges, ions):
        """Compute Q^n structural parameters for network glasses.

        Q^n (for a network-forming ion) is the count of bridging neighbours
        bonded to that ion.  A neighbour qualifies as a bridge only when it is
        in the `bridges` element list AND has more than one bond — i.e. it
        truly bridges between two polyhedra rather than being a terminal
        (non-bridging) species.

        The computed per-atom values are then aggregated into absolute counts,
        fractional occupancies, and an overall network connectivity scalar for
        the whole system, considering only the targeted ions.

        We assume that Q^n is meaningful for n in the range 0..8.  Overbonding
        to n=5 is expected to be rare, n=6 practically unheard of, and n=7,8
        virtually impossible — but values are tracked for all of them.

        Populated attributes
        --------------------
        atom_qn : list (1-indexed)
            Q^n value for each atom.  -1 for non-ion or out-of-bounds atoms.
        absolute_sys_qn : list (0-indexed, length 9)
            ``absolute_sys_qn[n]`` = number of ions with Q^n == n, for
            n = 0..8.  The index *is* the Q^n value n, following Perl's
            ``$absoluteSysQn[0..8]`` convention — a deliberate 0-indexed
            exception to the surrounding 1-indexed style because ``n``
            here is literally the array index, not an atom/axis label.
        fractional_sys_qn : list (0-indexed, length 9)
            ``fractional_sys_qn[n]`` = absolute_sys_qn[n] / num_qn_atoms.
        num_qn_atoms : int
            Total number of participating (non-skipped) ions.
        net_conn_qn : float
            Overall network connectivity: ``sum_n { n * fractional_sys_qn[n] }``.

        Parameters
        ----------
        bridges : list of str
            Lowercase element names that act as bridging atoms (e.g. ['o']).
        ions : list of str
            Lowercase element names of the network-forming ions (e.g. ['si']).

        Prerequisites
        -------------
        Call set_limit_dist and set_border before this method.  Bonding data
        are computed automatically if not already present.
        """
        # If the length of the num_bonds array equals zero, then we have not
        # yet computed the bonds between atoms and so we should do that.  If
        # however, the array length is finite, then we have to check to be
        # sure that the extended calculation is mapped back to the central cell.
        if len(self.num_bonds) <= 1:
            self.create_bonding_list()
        elif len(self.bonded) <= 1:
            self.map_ext_to_central()

        bridge_set = set(bridges)
        ion_set = set(ions)

        # Compute Q^n for each atom.
        self.atom_qn = [None] + [0] * self.num_atoms
        for atom in range(1, self.num_atoms + 1):
            # If the atom is not in the list of allowed ions, then skip it.
            if self.atom_element_name[atom] not in ion_set:
                self.atom_qn[atom] = -1
                continue

            # Compare the position of this atom with the box the user defined
            # (or the default).  This also considers if the user wanted only
            # atoms inside the box or only outside.  The return value is 1 if
            # the atom is out of the desired region and 0 if it is inside.
            if self.item_out_of_bounds(atom):
                self.atom_qn[atom] = -1
                continue

            # Perform a direct and simple observation of the number of bridging
            # atoms bonded to the current atom.  (Note: the Perl comment
            # reads "non-bridging" here, but the code — and the intent — is to
            # count bridging bonds, not non-bridging ones.)
            self.atom_qn[atom] = 0
            for bond in range(1, self.num_bonds[atom] + 1):
                # Get the atom id number of the bonded atom.
                bonded_atom = self.bonded[atom][bond]

                # Compare this to the list of bridge elements.
                if self.atom_element_name[bonded_atom] in bridge_set:
                    # Check if this bridge atom has more than one bond.
                    # (I.e., it bridges from the current atom to _any_ other
                    # atom(s).)  If it does, increment the Q^n count for atom.
                    if self.num_bonds[bonded_atom] > 1:
                        self.atom_qn[atom] += 1

        # Now that the Q^n for each atom has been determined, we can compute
        # the absolute and fractional Q^n for the whole system.  (Obviously,
        # we consider the targeted ions only.)

        # Allocate the Q^n histogram using Perl's 0-indexed convention:
        # absolute_sys_qn[n] = number of atoms with Q^n value n, for
        # n in 0..8.  Here n *is* the index — there is no atom/axis
        # sentinel — so using slot 0 as a live data entry matches Perl's
        # @absoluteSysQn[0..8] exactly.
        self.absolute_sys_qn = [0] * 9

        # Now, accumulate Q^n from each atom and count the total number of
        # participating ions.
        self.num_qn_atoms = 0
        for atom in range(1, self.num_atoms + 1):
            if self.atom_qn[atom] == -1:  # Ignore non-target atoms.
                continue
            self.num_qn_atoms += 1
            self.absolute_sys_qn[self.atom_qn[atom]] += 1

        # Once the absolute Q^n has been computed, we can easily obtain the
        # fractional Q^n, using the same 0-indexed layout.
        self.fractional_sys_qn = [0.0] * 9
        if self.num_qn_atoms > 0:
            for n in range(9):
                self.fractional_sys_qn[n] = (
                    self.absolute_sys_qn[n] / self.num_qn_atoms)

        # From the fractional Q^n, we can compute the overall network
        # connectivity.
        self.net_conn_qn = sum(
            n * self.fractional_sys_qn[n] for n in range(9))

    def compute_ring_distribution(self, min_ring_len, max_ring_len):
        """Identify rings in the bonding network and tabulate their sizes.

        Prerequisite: create_bonding_list() must have been called (called
        internally) and compute_min_network_dist_matrix() must be callable.

        Performs a modified depth-first traversal starting from each atom.
        For each target ring length from min_ring_len to max_ring_len, every
        atom is used as the root of a traversal that finds all shortest-path
        rings of that length.

        The traversal has two phases:
          Outward phase (level_delta = +1): push bonded atoms whose network
            distance from the root increases by one at each step, until we
            reach atoms that are ring_midpoint steps away.
          Inward phase (level_delta = -1): push bonded atoms whose network
            distance from the root decreases by one at each step, until we
            return to the root atom.

        The depth "searches" are really traversals so that we find every
        single possible ring.

        A ring is complete when the traversal stack yields the root atom again
        and the ring buffer contains ring_len + 1 atoms (head repeated at
        tail).  Completed rings are passed to save_ring(), which sorts the
        atom list and discards duplicates.

        Once all rings of a given length through an atom have been found, that
        atom is added to completed_vertices and skipped as a bond target for
        subsequent starting atoms.

        Populates self.ring_counts and self.rings.

        Parameters
        ----------
        min_ring_len : int
            Smallest ring size to search for.
        max_ring_len : int
            Largest ring size to search for.
        """
        self.create_bonding_list()
        self.compute_min_network_dist_matrix()

        # ring_counts[ring_len] = number of unique rings of that length.
        # rings[ring_len] is a 1-indexed list; rings[ring_len][i] is a sorted
        # list of atom indices for the i-th ring of that length.
        self.ring_counts = [None] + [0] * max_ring_len
        self.rings = [None] + [[] for _ in range(max_ring_len)]
        # Pad each rings[ring_len] with a None sentinel at index 0.
        for ring_len in range(1, max_ring_len + 1):
            self.rings[ring_len] = [None]

        # Backtrack stack — sentinel 0 at the bottom; tracks how many
        # branches were added at each level of the DFS so we know how
        # far to unwind curr_ring when backtracking.
        back_track_stack = [0]

        for ring_len in range(min_ring_len, max_ring_len + 1):
            # Define the midpoint of the ring we are seeking and also set an
            # even/odd flag.
            if ring_len % 2 == 0:
                even_ring = True
                ring_midpoint = ring_len // 2
            else:
                even_ring = False
                ring_midpoint = (ring_len - 1) // 2

            # Initialize the list of vertices for which all rings of the
            # current size have been found.
            completed_vertices = [0] * (self.num_atoms + 1)

            for atom in range(1, self.num_atoms + 1):
                # Perform a modified depth first search for every atom that
                # is ring_midpoint away from the current starting atom.  The
                # modification is that every step out to the distant atom must
                # have a distance from the starting vertex atom that is one
                # higher than the last.
                # Once an atom that is ring_midpoint away from the starting
                # atom has been found, we perform another depth first search
                # back to the starting atom with the modified requirement that
                # all steps must decrease the value of their network distance
                # from the starting atom by one (for each step) and we cannot
                # have any duplicate vertices in the ring list.
                # The depth "searches" are really traversals so that we find
                # every single possible ring.

                # Initialize curr_out_level to -1 so that when the first
                # atom is added to the ring it will be at level 0.
                curr_out_level = -1
                out_stack = []
                level_delta = 1  # +1 = going out; -1 = coming back
                curr_ring = []

                # Initialize the outgoing stack with the first (root) atom.
                out_stack.append(atom)

                while out_stack:
                    # Any atom that had been added to the stack will
                    # definitely be part of some ring.  Pop the top atom and
                    # add it to the current ring buffer.
                    curr_atom = out_stack.pop()
                    curr_ring.append(curr_atom)

                    # If the atom already exists in the ring, set a duplicate
                    # member flag.  We check all positions except the one we
                    # just appended (curr_ring[:-1]).
                    dup_member = curr_atom in curr_ring[:-1]

                    # ----------------------------------------------------------
                    # Complete-ring check: the ring buffer has ring_len+1
                    # atoms (head appears again at tail) and head == tail.
                    # ----------------------------------------------------------
                    if (len(curr_ring) == ring_len + 1 and
                            curr_ring[0] == curr_ring[ring_len]):
                        # Record that we found a new ring.
                        self.ring_counts[ring_len] += 1

                        # Save the current ring minus the tail (which
                        # duplicates the head).
                        curr_ring.pop()
                        self.save_ring(ring_len, curr_ring)

                        # The stack may have additional elements in it.  We
                        # need to backtrack the ring to the point where the
                        # divergence occurred.  Consider a stack that grew as
                        # A,B,G,C,D,E,F,A and a ring that grew as
                        # A,B,C,D,E,F,A.  The stack has already popped
                        # A,F,E,D,C and so contains only A,B,G.  Thus, we
                        # need to pop the ring all the way back to B and then
                        # add on G.
                        while (back_track_stack and
                               back_track_stack[-1] == 1):
                            back_track_stack.pop()
                            curr_ring.pop()
                        if len(back_track_stack) > 1:
                            back_track_stack[-1] -= 1

                        # Compute the new level delta.
                        level_delta = self.compute_ring_level_delta(
                            ring_len, ring_midpoint, len(curr_ring) - 1)

                        # The current out level is the level of the last atom
                        # in the ring (unless there are no atoms left).
                        if curr_ring:
                            curr_out_level = (
                                self.min_net_dist[atom][curr_ring[-1]])
                        continue

                    # ----------------------------------------------------------
                    # Update the current out level for the atom we just popped.
                    # ----------------------------------------------------------
                    # While going out, any atom that we pop will be at a
                    # "one higher" level than wherever we were, except in the
                    # case that we are in the middle of the turning point for
                    # an odd-length ring.
                    # This "if" checks if the current level (before going up)
                    # is not equal to the midpoint.  If it is equal, that
                    # would indicate that the previous step out was at the
                    # first part of an odd ring turning point.  In that case,
                    # we do not increase the level but we do trigger the flag
                    # to start reducing levels.
                    if curr_out_level != ring_midpoint:
                        curr_out_level += level_delta  # level of curr_atom
                    else:
                        level_delta = -1

                    # If the ring has even length and we are at the midpoint,
                    # then we also need to start reducing the levels.  Note
                    # that this "if" is different from the last one because
                    # this is asked *after* the current level has been
                    # incremented to the level of the current (just popped)
                    # atom.
                    if even_ring and curr_out_level == ring_midpoint:
                        level_delta = -1

                    # ----------------------------------------------------------
                    # Consider each atom bonded to the current atom.  As long
                    # as the bonded atom meets certain conditions, add it to
                    # the stack where it will eventually be part of some ring.
                    # Assume that we will not add any bonded atoms to the
                    # stack.
                    # ----------------------------------------------------------
                    curr_stack_adds = 0
                    for bond in range(1, self.num_bonds[curr_atom] + 1):
                        # For convenience, get the atom number of the bonded
                        # atom and its minimum network distance from the root.
                        bonded_atom = self.bonded[curr_atom][bond]
                        bonded_mnd = self.min_net_dist[atom][bonded_atom]

                        # In all cases, never add an atom that has already
                        # had all of its rings constructed.
                        if completed_vertices[bonded_atom]:
                            continue
                        # In all cases, never add bonded atoms if the current
                        # "root" atom is a duplicate member in the ring.
                        if dup_member:
                            continue

                        if level_delta == 1:
                            # In all cases, never add an atom that has a
                            # level greater than the ring midpoint.
                            if bonded_mnd > ring_midpoint:
                                continue
                            # In all cases, never add an atom that has a
                            # level less than the current level.
                            if bonded_mnd < curr_out_level:
                                continue
                            # For an even midpoint ring, never add an atom
                            # that has a level the same as the current level.
                            if even_ring and bonded_mnd == curr_out_level:
                                continue
                            # For an odd midpoint ring, never add an atom
                            # that has the same level as the current level
                            # and whose bonded level is less than the midpoint.
                            # (There is an exception when the bonded level
                            # equals the ring midpoint for odd rings.  E.g.,
                            # consider a ring with network distances 012210 —
                            # it has 5 atoms because the two 0-distance ends
                            # are the same atom.)
                            if (not even_ring and
                                    bonded_mnd == curr_out_level and
                                    bonded_mnd < ring_midpoint):
                                continue
                            # Now the only atoms left are "normal" next steps
                            # for the current ring or atoms at the midpoint
                            # of an odd ring that are at the same level.
                        else:  # level_delta == -1
                            # In all cases, never add an atom that has a
                            # level greater than or equal to the current level.
                            if bonded_mnd >= curr_out_level:
                                continue
                            # Now the only atoms left are "normal" next steps
                            # for the current ring.

                        out_stack.append(bonded_atom)
                        curr_stack_adds += 1

                    # ----------------------------------------------------------
                    # If no bonds were added then curr_atom is a dead end and
                    # should not be included in a possible ring.  However, if
                    # bonds were added to the stack, we need to track how many
                    # were added so that when we backtrack the ring we know
                    # how far to go.
                    # ----------------------------------------------------------
                    if curr_stack_adds == 0:
                        curr_ring.pop()
                        if not curr_ring:
                            continue
                        curr_out_level = (
                            self.min_net_dist[atom][curr_ring[-1]])

                        # The stack may have additional elements.  Backtrack
                        # the ring to the point where the divergence occurred.
                        # (Same A,B,G,C,D,E,F,A example as the complete-ring
                        # block above.)  If we pop all elements off the ring,
                        # then we are done.
                        while (back_track_stack and
                               back_track_stack[-1] == 1):
                            back_track_stack.pop()
                            curr_ring.pop()
                            if curr_ring:
                                curr_out_level = (
                                    self.min_net_dist[atom][curr_ring[-1]])

                        if not curr_ring:
                            continue
                        if len(back_track_stack) > 1:
                            back_track_stack[-1] -= 1

                        # Compute the new level delta.
                        level_delta = self.compute_ring_level_delta(
                            ring_len, ring_midpoint, len(curr_ring) - 1)
                    else:
                        back_track_stack.append(curr_stack_adds)

                # Once this atom is done, add it to the list of vertices that
                # should not be included in any future rings of this size.
                completed_vertices[atom] = 1

    def compute_ring_level_delta(self, target_ring_len, ring_midpoint,
                                  curr_ring_max_index):
        """Determine the level-delta direction (outward or inward) after
        backtracking to a ring of current depth curr_ring_max_index.

        Called after completing a ring or hitting a dead end, once curr_ring
        has been unwound to a new branch point.  Decides whether the traversal
        from that point should continue going outward (+1) or inward (-1).

        Parameters
        ----------
        target_ring_len : int
            Length of the ring being sought.
        ring_midpoint : int
            Midpoint level of the ring (ring_len//2 for even,
            (ring_len-1)//2 for odd).
        curr_ring_max_index : int
            Index of the last atom currently in curr_ring (i.e.
            len(curr_ring) - 1), used as a proxy for depth in the traversal.

        Returns
        -------
        int
            +1 (going outward) or -1 (coming back inward).
        """
        if target_ring_len % 2 == 0:
            # If the target ring length is even, we have an easy switch
            # between "going out" and "coming back" right at the midpoint.
            # If the level that the ring was popped back to is greater than
            # or equal to the ring midpoint, set level_delta = -1.  Else
            # (the ring was popped back to a level less than the midpoint),
            # set level_delta = 1.
            if curr_ring_max_index >= ring_midpoint:
                return -1
            else:
                return 1
        else:
            # If the ring length is odd:
            # If the level that the ring was popped back to is greater than
            # the ring midpoint, set level_delta = -1.  (I.e., we are on the
            # second of two equally distant atoms, or beyond.)  Else (the
            # ring was popped back to a level less than or equal to the
            # midpoint), set level_delta = 1.
            if curr_ring_max_index > ring_midpoint:
                return -1
            else:
                return 1

    def save_ring(self, ring_len, curr_ring):
        """Record a completed ring, discarding duplicates.

        Called by compute_ring_distribution immediately after incrementing
        ring_counts[ring_len].  Sorts curr_ring so that atom-by-atom
        comparison with previously saved rings is straightforward.

        If this is the first ring of the given length (ring_counts == 1), it
        is saved unconditionally.  Otherwise every previously saved ring of
        that length is compared against the sorted candidate.  If a duplicate
        is found, the erroneous increment to ring_counts is reversed (−1) and
        the ring is discarded.

        Parameters
        ----------
        ring_len : int
            Length of the ring (number of atoms).
        curr_ring : list of int
            Atom indices (1-indexed) forming the ring, in traversal order.
            The list should NOT include the repeated head atom at the tail.
        """
        # Sort the current ring atoms so that we can make an easy comparison
        # with other already saved rings.
        sorted_ring = sorted(curr_ring[:ring_len])

        # Just prior to calling this subroutine, ring_counts was incremented
        # by one.  So, if ring_counts equals one, it is the first ring to be
        # saved — we can just save the sorted ring and exit.
        if self.ring_counts[ring_len] == 1:
            self.rings[ring_len].append(sorted_ring)
            return

        # If we did not exit, determine if the ring we just counted should
        # really be saved — i.e., is it a duplicate of an already-saved ring?

        # We start by assuming the ring we want to add is different from every
        # other ring in the current list of rings.
        # rings[ring_len] is 1-indexed (index 0 is the None sentinel), so
        # iterate over indices 1 .. ring_counts[ring_len]-1.
        same = False
        for ring_idx in range(1, self.ring_counts[ring_len]):
            # For this specific ring, assume the candidate is the same as it.
            # Check for a match atom-by-atom.  Because the rings are sorted,
            # we can do a direct list comparison.
            if self.rings[ring_len][ring_idx] == sorted_ring:
                # If there was a match for this specific ring, we found a
                # global match and can stop looking.
                same = True
                break

        if not same:
            self.rings[ring_len].append(sorted_ring)
        else:
            # Revert the erroneous ring count increment.
            self.ring_counts[ring_len] -= 1

    def compute_bond_angles_ext(self):
        """Compute all bond angles in the extended cell.

        For each atom in the central cell, loops over all unique pairs of
        bonded extended-cell neighbours (bond1 < bond2) and computes the
        angle at the vertex (central) atom using the law of cosines:

            c^2 = a^2 + b^2 - 2ab*cos(gamma)

        where gamma is the angle we wish to discover, a is the distance from
        the vertex atom to bonded atom #1, b is the distance from the vertex
        atom to bonded atom #2, and c is the distance between the two bonded
        atoms.

        Populates:
          self.num_bond_angles[atom]          -- count of angles at each atom
          self.bond_angles_ext[atom][1..]     -- angles in degrees (1-indexed,
                                                 index 0 is None placeholder)

        Prereqs: create_bonding_list / create_ext_bonding_list must have been
        called so that bonded_ext, bond_length_ext, ext_direct_xyz_list, and
        num_bonds are populated.
        """
        # Initialize per-atom angle count and angle list (1-indexed; index 0
        # is the None placeholder matching the Perl array convention).
        self.num_bond_angles = [None] + [0      for _ in range(self.num_atoms)]
        self.bond_angles_ext = [None] + [[None] for _ in range(self.num_atoms)]

        for atom in range(1, self.num_atoms + 1):
            # Compare the position of this atom with the box the user defined
            # (or the default).  This also considers if the user wanted only
            # atoms inside the box or only outside.  The return value is 1 if
            # the atom is out of the desired region and it is 0 if the atom is
            # inside the desired region.
            if self.item_out_of_bounds(atom):
                continue

            # Initialize the count of the number of bond angles for this atom.
            # (Already set to 0 above; the comment is Perl-authoritative.)

            # Use the law of cosines to compute the angle between this vertex
            # atom and its bonded atom pairs:  c^2 = a^2 + b^2 - 2ab*cos(gamma).
            # Gamma = the angle we wish to discover, a = distance from vertex
            # atom to bonded atom #1, b = distance from vertex atom to bonded
            # atom #2 and c = distance between two bonded atoms.
            for bond1 in range(1, self.num_bonds[atom]):
                # Short name for the distance between vertex and this bonded atom.
                a = self.bond_length_ext[atom][bond1]
                ext1 = self.bonded_ext[atom][bond1]

                for bond2 in range(bond1 + 1, self.num_bonds[atom] + 1):
                    # Increment the count of the number of bond angles.
                    self.num_bond_angles[atom] += 1

                    # Short name for the distance between vertex and this bonded atom.
                    b = self.bond_length_ext[atom][bond2]
                    ext2 = self.bonded_ext[atom][bond2]

                    # Now we must compute the distance between bonded atom #1 and
                    # bonded atom #2.  First we get the difference vector between
                    # the a and b positions, then we find the magnitude of that
                    # distance.
                    c = math.sqrt(sum(
                        (self.ext_direct_xyz_list[ext1][ax] -
                         self.ext_direct_xyz_list[ext2][ax]) ** 2
                        for ax in range(1, 4)
                    ))

                    # Apply the law of cosines:  c^2 = a^2 + b^2 - 2ab cos(gamma)
                    # to solve for gamma (the angle between two bonds) where c is
                    # the distance between the two bonded atoms, a is the distance
                    # from the central (vertex) atom to one bonded atom and b is
                    # the distance from the vertex to the other.
                    cos_gamma = (a*a + b*b - c*c) / (2.0 * a * b)
                    # Guard against numerical noise that pushes cos_gamma below -1,
                    # which would make acos undefined.  When cos_gamma ≈ -1 the
                    # bond angle is exactly 180°.
                    if abs(cos_gamma + 1.0) < EPSILON:
                        gamma = math.pi
                    else:
                        gamma = math.acos(cos_gamma)

                    # Record for general use in degrees.
                    self.bond_angles_ext[atom].append(
                        gamma * 180.0 / math.pi)

    # (item_not_requested, pair_not_of_element, item_out_of_bounds are
    #  implemented above near create_bonding_list)

    def compute_min_network_dist_matrix(self):
        """Compute shortest-path distances through the bonding network.

        Consider each atom and construct a minimal network distance array for
        it.  In contrast to the regular minimal distance (which is Euclidean),
        this distance is the count of integer steps that must be taken on a
        bonding network to travel from one node to another.  If a node is
        unreachable from another node, then the distance is -1.

        The network lengths are determined by a breadth-first traversal of the
        bonding network.  The bonding network is defined in the context of
        periodic boundary conditions, but all atom numbers are from the central
        cell (i.e., between 1 and num_atoms inclusive).

        Populates self.min_net_dist, which is 1-indexed:
        min_net_dist[atom1][atom2] is the minimum number of bond hops from
        atom1 to atom2.  Unreachable pairs have distance -1.  Distance from an
        atom to itself is 0.
        """
        self.min_net_dist = [None]  # index 0 unused

        for atom in range(1, self.num_atoms + 1):
            # visited[i] = hop distance from `atom` to atom i; -1 = unvisited.
            # Distance to myself starts as zero.
            visited = [-1] * (self.num_atoms + 1)
            visited[atom] = 0

            queue = [atom]

            while queue:
                # Get the vertex that we will now use as the current starting
                # point for extending the traversal path.
                curr_atom = queue.pop(0)

                # Compute the depth of the path.
                curr_depth = visited[curr_atom] + 1

                # Record the depth for each bonded atom that has not already
                # been visited.
                for bond in range(1, self.num_bonds[curr_atom] + 1):
                    # For convenience, get the atom number of the bonded atom.
                    bonded_atom = self.bonded[curr_atom][bond]

                    # Determine if the atom at the other end of this bond has
                    # been visited yet.
                    if visited[bonded_atom] == -1:
                        # Since it has not been visited, record the path length.
                        visited[bonded_atom] = curr_depth

                        # Add the bonded atom to the queue so it can later be
                        # used as a starting point vertex for future depth steps.
                        queue.append(bonded_atom)

            # At this point, the breadth-first traversal should be complete and
            # we can add this path length list (visited array) to the minimal
            # network distance matrix.
            self.min_net_dist.append(visited)

    # ======================================================================
    # Section: Math and vector utilities
    # ======================================================================

    def get_plane_normal(self, p1, p2, p3):
        """Return a normal vector to the plane defined by three points.

        The normal is the cross product of two edge vectors:
          d1 = p2 - p1  (first difference vector)
          d2 = p3 - p1  (second difference vector)
        The result is NOT normalized to unit length.

        The Perl original also computed the triangle centroid (midpoint of p1,p2
        then 1/3 of the way from p3 toward that midpoint), but never used it —
        dead code that is correctly omitted here.

        Parameters
        ----------
        p1, p2, p3 : sequence of float
            Three non-collinear points, each 1-indexed
            ``[None, x, y, z]``.  Perl stored these as
            ``@p[1..3]``; Python preserves that convention with a
            ``None`` sentinel at slot 0.

        Returns
        -------
        list of float
            1-indexed normal vector ``[None, nx, ny, nz]`` (not unit
            length).
        """
        # Build the two edge vectors in 1-indexed form so they can be fed
        # directly into cross_product without offset bookkeeping.
        d1 = [None] + [p2[axis] - p1[axis] for axis in range(1, 4)]
        d2 = [None] + [p3[axis] - p1[axis] for axis in range(1, 4)]
        return self.cross_product(d1, d2)

    def dot_product(self, a, b):
        """Return the scalar dot product of two 3-vectors.

        Perl equivalent: dotProduct (line 8282), which operated on
        1-indexed array references ``$vector_ref->[1..3]``.  The Python
        port preserves that convention so the same mathematical
        expression reads identically in either language.

        Parameters
        ----------
        a, b : sequence of float
            Input vectors in 1-indexed form ``[None, x, y, z]``.

        Returns
        -------
        float
        """
        return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

    def cross_product(self, a, b):
        """Return the cross product a x b.

        Computes the standard right-hand-rule cross product of two
        3-vectors using the 1-indexed convention established throughout
        this module (``v[1]`` = x, ``v[2]`` = y, ``v[3]`` = z).

        Parameters
        ----------
        a, b : sequence of float
            1-indexed input vectors ``[None, x, y, z]``.

        Returns
        -------
        list of float
            1-indexed result ``[None, cx, cy, cz]``.  Slot 0 is a
            sentinel (``None``) matching Perl's ``$product[0] = 0.0``.

        Notes
        -----
        Perl equivalent: crossProduct (StructureControl.pm).
        """
        return [
            None,
            a[2]*b[3] - a[3]*b[2],   # cx = ay*bz - az*by
            a[3]*b[1] - a[1]*b[3],   # cy = az*bx - ax*bz
            a[1]*b[2] - a[2]*b[1],   # cz = ax*by - ay*bx
        ]

    def cross_product_mag(self, a, b):
        """Return the magnitude of the cross product |a x b|.

        Parameters
        ----------
        a, b : sequence of float
            1-indexed input vectors ``[None, x, y, z]``.

        Returns
        -------
        float
            Euclidean norm of the cross product vector.

        Notes
        -----
        Perl equivalent: crossProductMag (StructureControl.pm). Like
        the Perl original this routine sums components ``[1..3]`` of the
        cross-product result.
        """
        c = self.cross_product(a, b)
        return math.sqrt(c[1]**2 + c[2]**2 + c[3]**2)

    def normalized_cross_product(self, a, b):
        """Return the unit vector in the direction of a x b.

        Computes the cross product c = a x b via cross_product(), then
        divides each component by the magnitude |c| to produce a unit vector.
        This is the same operation as cross_product_mag() but returns the full
        normalized vector rather than just the scalar magnitude.

        Parameters
        ----------
        a, b : sequence of float
            1-indexed input vectors ``[None, x, y, z]``.

        Returns
        -------
        list of float
            1-indexed unit cross product vector
            ``[None, x, y, z]``.

        Notes
        -----
        Perl equivalent: normalizedCrossProduct.  The Perl original used
        1-indexed arrays with a 0.0 placeholder at index [0]; Python
        mirrors that layout with a ``None`` placeholder.
        """
        c = self.cross_product(a, b)
        # Magnitude of the cross product vector (sum over slots 1..3).
        mag = math.sqrt(c[1]**2 + c[2]**2 + c[3]**2)
        # Divide each component by the magnitude to normalize, preserving
        # the [None, x, y, z] sentinel layout.
        return [None] + [c[axis] / mag for axis in range(1, 4)]

    def get_vector_angle(self, a, b):
        """Return the angle in radians between vectors a and b.

        Computes the angle from:
            theta = arccos( (v1 . v2) / (|v1| |v2|) )

        The result is corrected for slight deviations from 90 degrees that may
        occur because some vector elements are expressed as large numbers rather
        than infinity.  If the computed angle is within EPSILON of pi/2 it is
        snapped to exactly pi/2.

        The cosine is clamped to [-1, 1] before ``acos`` to guard against
        floating-point rounding outside that domain (a Python-side
        robustness improvement over the Perl original).

        Parameters
        ----------
        a, b : sequence of float
            1-indexed input vectors ``[None, x, y, z]``.

        Returns
        -------
        float
            Angle in radians in [0, pi].
        """
        # Compute the magnitudes of the two vectors (sum over slots 1..3).
        mag_a = math.sqrt(a[1]**2 + a[2]**2 + a[3]**2)
        mag_b = math.sqrt(b[1]**2 + b[2]**2 + b[3]**2)

        # Compute the angle between the vectors from:
        #   theta = arccos(v1 . v2) / (|v1| |v2|).
        # Clamp to [-1, 1] to guard against floating-point rounding (Python addition).
        cos_t = (a[1]*b[1] + a[2]*b[2] + a[3]*b[3]) / (mag_a * mag_b)
        angle = math.acos(max(-1.0, min(1.0, cos_t)))

        # Correct for slight deviations from 90 degrees which may occur because
        #   we do not express some vector elements as infinity.  Instead they are
        #   expressed as large numbers.
        if abs(angle - PI / 2.0) < EPSILON:
            angle = PI / 2.0
        return angle

    # Perl's ``getUniqueArray`` has been removed from the Python port.
    # It was dead code in this module and everywhere in the converted
    # scripts: ``create_coordination_list`` and ``create_coordination_summary``
    # reimplement the unique-with-counts logic inline with a Python dict,
    # and no external script (makeinput, pdb2skl, condense, bond_analysis,
    # make_reactions, mod_struct) called it.  Preserving it only as a
    # lexicographic uniqueifier also silently dropped Perl's count-return
    # contract.  Any future caller that wants the sorted unique values of
    # a list should use ``sorted(set(values))`` directly; callers needing
    # per-element counts should use ``collections.Counter``.

    def spherical_angles(self, vector):
        """Convert a Cartesian vector to spherical angles (theta, phi).

        Called by the Ylm / BOO accumulation pipeline to convert a bond
        difference vector into the two spherical-coordinate angles needed
        by the real spherical harmonic Y_l^m(theta, phi).

        Perl counterpart: sphericalAngles (StructureControl.pm ~8077).
        The Perl version accepts a 1-indexed reference ``[$x,$y,$z]``
        and clamps in-place.  The Python port keeps the same 1-indexed
        input convention so the caller can pass a bond difference vector
        directly without any 0-to-1 offset translation, and builds a
        local clamped copy so the caller's original vector is left
        untouched (a Python-side non-mutation improvement).

        Parameters
        ----------
        vector : sequence of float
            1-indexed Cartesian difference vector ``[None, x, y, z]``.

        Returns
        -------
        tuple of float
            (theta, phi) in radians.
            theta — polar angle measured from the +z axis (0 .. pi).
            phi   — azimuthal angle measured from +x toward +y (-pi .. pi).

        Notes
        -----
        Epsilon clamping: components whose absolute value is below EPSILON
        are forced to exactly 0.0 to correct for numerical round-off before
        the atan2 calls.

        Undefined phi: when both x and y are zero the point lies on the
        z-axis and phi is geometrically undefined; following the Perl
        original we return phi = 0.0 in that case.
        """
        # Correct for numerical error (Perl: "foreach $axis (1..3) { if abs < epsilon }").
        # Build a fresh clamped copy in 1-indexed form so the caller's
        # vector is not mutated and the axis bindings read the same as
        # the surrounding 1-indexed code.
        clamped = [None] + [0.0 if abs(vector[axis]) < EPSILON else vector[axis]
                            for axis in range(1, 4)]
        x, y, z = clamped[1], clamped[2], clamped[3]

        # Compute theta: polar angle from +z.
        theta = math.atan2(math.sqrt(x**2 + y**2), z)

        # Compute phi: azimuthal angle from +x.
        # In the case that phi is undefined (point on z-axis) we set it to zero.
        phi = 0.0 if (x == 0.0 and y == 0.0) else math.atan2(y, x)

        return (theta, phi)

    def init_ylm_coeff(self):
        """Initialise the real spherical-harmonic normalization coefficients for
        angular-momentum order ylm_l (set by set_ylm_l).

        Populates self.ylm_coeff as a 1-indexed list (index 0 is None) with
        2*ylm_l + 1 entries, one per magnetic quantum number m ordered
        m = -l, -(l-1), ..., 0, ..., l-1, +l.

        These are the angular-normalization prefactors for the real solid
        spherical harmonics Y_l^m used in Bond-Order Parameter (BOO/Ylm)
        analysis (see accumulate_q_bar, q_bar_per_bond, q_bar_normalize,
        q_bar_correlation).  The actual Y_l^m value for a given bond direction
        is obtained by multiplying the coefficient by the appropriate angular
        function of (theta, phi).

        Supported orders:
          l = 0  ->  1 coefficient  (isotropic; Y_0^0 = (1/2) sqrt(1/pi))
          l = 1  ->  3 coefficients (p-like; Y_1^{-1}, Y_1^0, Y_1^{+1})
          l = 6  ->  13 coefficients (used for icosahedral / FCC order;
                     see Steinhardt et al., Phys. Rev. B 28, 784 (1983))

        Raises ValueError for unsupported ylm_l values.  The Perl source
        (initYlmCoeff, StructureControl.pm:8496) silently left ylm_coeff
        unchanged for unsupported values.
        """
        if self.ylm_l == 0:
            # l=0: single isotropic term Y_0^0
            self.ylm_coeff = [None,
                0.5 * math.sqrt(1.0 / PI),
            ]
        elif self.ylm_l == 1:
            # l=1: three p-like terms (m = -1, 0, +1)
            self.ylm_coeff = [None,
                 0.5 * math.sqrt(3.0 / 2.0 / PI),
                 0.5 * math.sqrt(3.0 / PI),
                -0.5 * math.sqrt(3.0 / 2.0 / PI),
            ]
        elif self.ylm_l == 6:
            # Initialize the values of the Ylm for order 6.
            self.ylm_coeff = [None,
                 1.0/64.0 * math.sqrt(3003.0 / PI),
                 3.0/32.0 * math.sqrt(1001.0 / PI),
                 3.0/32.0 * math.sqrt(91.0 / 2.0 / PI),
                 1.0/32.0 * math.sqrt(1365.0 / PI),
                 1.0/64.0 * math.sqrt(1365.0 / PI),
                 1.0/16.0 * math.sqrt(273.0 / 2.0 / PI),
                 1.0/32.0 * math.sqrt(13.0 / PI),
                -1.0/16.0 * math.sqrt(273.0 / 2.0 / PI),
                 1.0/64.0 * math.sqrt(1365.0 / PI),
                -1.0/32.0 * math.sqrt(1365.0 / PI),
                 3.0/32.0 * math.sqrt(91.0 / 2.0 / PI),
                -3.0/32.0 * math.sqrt(1001.0 / PI),
                 1.0/64.0 * math.sqrt(3003.0 / PI),
            ]
        else:
            raise ValueError(f'init_ylm_coeff: unsupported ylm_l={self.ylm_l}')

    def accumulate_q_bar(self, atom, theta, phi):
        """Accumulate one bond's spherical-harmonic contribution into q_bar.

        Called by obtain_atomic_interaction (flag=4, via create_q_list) once for
        every bond of every atom.  For each call, atom is the central atom and
        (theta, phi) are the polar/azimuthal angles of the bond vector in the
        lab frame.

        This method adds the real spherical harmonic Y_l^m(theta, phi) — scaled
        by the precomputed normalization coefficient ylm_coeff[m_index] — to the
        running sum in self.q_bar[atom].  After all bonds have been accumulated,
        q_bar_per_bond divides by the bond count to form the per-bond average
        q̄_lm, and q_bar_normalize/q_bar_correlation complete the BOO analysis.

        Array layout
        ------------
        q_bar[atom] is a 1-indexed list.  Each entry q_bar[atom][k] is a
        1-indexed real/imag pair ``[None, real, imag]``: slot 1 is the
        real component and slot 2 is the imaginary component, matching
        Perl's ``$qBar[$atom][$k][1..2]`` layout exactly.  Slot 0 is an
        unused ``None`` sentinel so the module-wide 1-indexed convention
        is honoured at every nesting level.
        Index k runs from 1 to 2*l+1; the corresponding magnetic quantum number
        is m = k - (l+1), i.e. k=1 → m=-l, k=l+1 → m=0, k=2l+1 → m=+l.

        Parameters
        ----------
        atom : int
            Central atom index (1-indexed).
        theta : float
            Polar angle in radians (from +z axis).
        phi : float
            Azimuthal angle in radians (from +x axis).

        Notes
        -----
        Perl debug artifacts (not replicated here):
          - A commented-out print at the top:
            #print STDOUT "atom=$atom theta=$theta   phi=$phi\\n"
          - An active debug print inside the l=0 branch:
            print STDOUT "$qBar[$atom][1][1]\\n"

        For l=6 at m=0 (k=7): Perl writes c[7]*poly[7] to BOTH the real and
        imaginary components (mathematically sin(0·phi)=0, so the imaginary
        part should be zero; this appears to be a Perl bug faithfully preserved
        here).
        """
        c  = self.ylm_coeff
        qb = self.q_bar[atom]
        ct = math.cos(theta)
        st = math.sin(theta)

        if self.ylm_l == 0:
            # l=0: Y_0^0 = c[1] (isotropic); imaginary component is always zero
            qb[1][1] += c[1]
            qb[1][2] += 0.0

        elif self.ylm_l == 1:
            # REAL parts (cos(m*phi) * angular factor)
            qb[1][1] += c[1] * math.cos(-1.0 * phi) * st
            qb[2][1] += c[2] * ct
            qb[3][1] += c[3] * math.cos( 1.0 * phi) * st
            # IMAGINARY parts (sin(m*phi) * angular factor; m=0 term is zero)
            qb[1][2] += c[1] * math.sin(-1.0 * phi) * st
            qb[2][2] += 0.0
            qb[3][2] += c[3] * math.sin( 1.0 * phi) * st

        elif self.ylm_l == 6:
            ct2 = ct**2; ct3 = ct**3; ct4 = ct**4; ct5 = ct**5; ct6 = ct**6
            st2 = st**2; st3 = st**3; st4 = st**4; st5 = st**5; st6 = st**6
            # Theta-only Associated Legendre polynomial factor for each m_index
            # k=1..13, corresponding to m = k-7 = -6..+6.  Matches the explicit
            # 13-line REAL block and 13-line IMAGINARY block in the Perl source.
            poly = [None,
                st6,                                      # k=1  m=-6
                st5 * ct,                                 # k=2  m=-5
                st4 * (11.0*ct2 - 1.0),                   # k=3  m=-4
                st3 * (11.0*ct3 - 3.0*ct),                # k=4  m=-3
                st2 * (33.0*ct4 - 18.0*ct2 + 1.0),        # k=5  m=-2
                st  * (33.0*ct5 - 30.0*ct3 + 5.0*ct),     # k=6  m=-1
                231.0*ct6 - 315.0*ct4 + 105.0*ct2 - 5.0,  # k=7  m= 0
                st  * (33.0*ct5 - 30.0*ct3 + 5.0*ct),     # k=8  m=+1
                st2 * (33.0*ct4 - 18.0*ct2 + 1.0),        # k=9  m=+2
                st3 * (11.0*ct3 - 3.0*ct),                # k=10 m=+3
                st4 * (11.0*ct2 - 1.0),                   # k=11 m=+4
                st5 * ct,                                  # k=12 m=+5
                st6,                                       # k=13 m=+6
            ]
            # Phi frequency: m = k - 7 (k=1→-6, k=7→0, k=13→+6)
            # REAL component (slot 1) uses cos(m*phi); IMAGINARY (slot 2)
            # uses sin(m*phi).  Exception at k=7 (m=0): Perl assigned
            # poly*coeff to both real and imaginary (sin(0)=0 so
            # imaginary should be 0; Perl bug preserved).
            for k in range(1, 14):
                m = float(k - 7)
                if k == 7:  # m=0: replicate Perl (both components get poly*coeff)
                    qb[k][1] += c[k] * poly[k]
                    qb[k][2] += c[k] * poly[k]
                else:
                    qb[k][1] += c[k] * math.cos(m * phi) * poly[k]  # REAL
                    qb[k][2] += c[k] * math.sin(m * phi) * poly[k]  # IMAGINARY

        else:
            raise ValueError(f'accumulate_q_bar: unsupported ylm_l={self.ylm_l}')

    def q_bar_per_bond(self):
        """Divide the q_bar for each atom by the number of bonds for that atom.

        This is step 1 of 3 in the bond-orientational-order (BOO) finalization
        sequence called from finalize_interaction_data (interaction_type == 4):
            1. q_bar_per_bond      — average over bonds (this method)
            2. q_bar_normalize     — normalize the averaged vector
            3. q_bar_correlation   — accumulate the correlation factor

        q_bar[atom][ylm_m] holds a 1-indexed ``[None, real, imag]`` pair
        where slot 1 is the real component and slot 2 is the imaginary
        component (matching Perl's ``@qBar[$atom][$Ylm_m][1..2]``).
        After this call each component is the per-bond average of the
        accumulated Ylm sums from accumulate_q_bar.

        Perl debug prints (commented out in the original) bracketed the inner
        loop to show q_bar values before and after the division; they are
        preserved here as a reminder of the data layout:
            # print(f"num_bonds[{atom}] = {self.num_bonds[atom]}")
            # print(self.q_bar[atom][1][1], self.q_bar[atom][2][1], ...)
        """
        for atom in range(1, self.num_atoms + 1):
            for ylm_m in range(1, 2 * self.ylm_l + 2):
                # Slots 1/2 hold the real/imag components; divide each
                # separately to average over the bond count.
                for part in range(1, 3):
                    self.q_bar[atom][ylm_m][part] /= self.num_bonds[atom]

    def q_bar_normalize(self):
        """Normalize each atom's q_bar complex vector in-place.

        q_bar[atom][ylm_m] holds a complex number stored as a 1-indexed
        ``[None, real, imag]`` pair (slot 1 = real, slot 2 = imaginary,
        slot 0 = sentinel).  The full vector for one atom has
        2*(2*l+1) real-valued components (real and imaginary parts for
        each of the 2l+1 magnetic quantum numbers m = -l … +l).

        The magnitude is the Euclidean norm of that real-valued vector:
            magnitude = sqrt( sum_{m=1}^{2l+1} ( re_m^2 + im_m^2 ) )
        Both the real and imaginary parts are then divided by this single
        scalar magnitude, yielding a unit vector in complex (2*(2l+1))-space.

        Note: the Perl source kept a commented-out alternative normalization
        that used abs() on each component before dividing — this was an
        experimental variant that was not adopted:
            # $qBar[$atom][$Ylm_m][1] = abs($qBar[$atom][$Ylm_m][1]) / $magnitude;
            # $qBar[$atom][$Ylm_m][2] = abs($qBar[$atom][$Ylm_m][2]) / $magnitude;
        The Perl source also kept several commented-out debug prints (showing
        magnitude per atom and the first three m components before and after
        normalization); they are omitted here.

        Prereq: accumulate_q_bar (or q_bar_per_bond) must have been called
        first to populate self.q_bar.
        """
        for atom in range(1, self.num_atoms + 1):
            # Initialize the magnitude for accumulation.
            magnitude = 0.0

            # Accumulate the sum of squares of the vector components
            # (slot 1 = real, slot 2 = imag).
            for ylm_m in range(1, 2 * self.ylm_l + 2):
                magnitude += (self.q_bar[atom][ylm_m][1] ** 2
                              + self.q_bar[atom][ylm_m][2] ** 2)

            # Obtain the magnitude of the vector in complex space.
            magnitude = math.sqrt(magnitude)

            # Divide each component by the magnitude to make it a normalized vector.
            for ylm_m in range(1, 2 * self.ylm_l + 2):
                for part in range(1, 3):
                    self.q_bar[atom][ylm_m][part] /= magnitude

    def q_bar_correlation(self):
        """Compute the bond orientational order parameter correlation for every atom.

        For each atom, accumulates the "dot product" of its normalised q_bar
        vector with the normalised q_bar of every bonded neighbour (Perl:
        "Accumulate the dot product of $atom's normalized qBar with the
        normalized qBar of each bonded atom"), then divides by the number
        of bonds (Perl: "Divide by the total number of bonded atoms").
        Results are stored in ``self.q_order[atom]`` (1-indexed; index 0 is
        None).

        Prerequisites
        -------------
        ``q_bar_normalize`` must have been called first so that each atom's
        ``q_bar`` vector is a unit vector in the 2*(2l+1)-dimensional real
        space of Re/Im components.

        Calls ``map_ext_to_central`` first so that ``self.bonded`` contains
        central-cell atom indices.

        Complex arithmetic note
        -----------------------
        The Perl expands the per-m contribution as four separate products::

            qBar[atom][m][1]  *  qBar[bonded][m][1]    # Re_a * Re_c
          + qBar[atom][m][1]  * -qBar[bonded][m][2]    # Re_a * (-Im_c)
          + qBar[atom][m][2]  *  qBar[bonded][m][1]    # Im_a * Re_c
          + qBar[atom][m][2]  *  qBar[bonded][m][2]    # Im_a * Im_c

        The Python port uses the same slot layout: slot 1 = real, slot 2
        = imaginary.  Expanded: a*c - a*d + b*c + b*d.  Note this differs
        from the standard complex inner product (which would be
        a*c + b*d); this specific formula is preserved verbatim from the
        Perl original.

        A commented-out debug print in the Perl (omitted here)::

            #print STDOUT "qOrder[$atom] = $qOrder[$atom]\\n";

        Output
        ------
        self.q_order : list
            1-indexed list of floats; ``q_order[atom]`` is the bond
            orientational order parameter for that atom.
        """
        self.map_ext_to_central()

        # Initialize the bond orientational order parameter for all atoms.
        self.q_order = [None] + [0.0] * self.num_atoms

        for atom in range(1, self.num_atoms + 1):
            # Accumulate the dot product of atom's normalized q_bar with the
            # normalized q_bar of each bonded atom.
            for bond in range(1, self.num_bonds[atom] + 1):
                bonded_atom = self.bonded[atom][bond]
                for ylm_m in range(1, 2 * self.ylm_l + 2):
                    a = self.q_bar[atom][ylm_m][1]        # Re(atom)
                    b = self.q_bar[atom][ylm_m][2]        # Im(atom)
                    c = self.q_bar[bonded_atom][ylm_m][1] # Re(bonded)
                    d = self.q_bar[bonded_atom][ylm_m][2] # Im(bonded)
                    self.q_order[atom] += a*c - a*d + b*c + b*d

            # Divide by the total number of bonded atoms.
            self.q_order[atom] /= self.num_bonds[atom]

    # Perl's ``stableSort`` helper has been removed from the Python port.
    # It was dead code: nothing inside structure_control.py, or anywhere in
    # the converted scripts (makeinput, pdb2skl, condense, bond_analysis,
    # make_reactions, mod_struct), called it.  The single in-file consumer
    # — ``sort_atoms`` — builds its own inline ``sorted(...)`` call, and
    # the 0-indexed signature of the previous port diverged from Perl's
    # 1-indexed in-place convention.  Rather than carry a silently
    # drifted public API, the method is deleted.  Any future caller that
    # genuinely needs a stable sort should call Python's built-in
    # ``sorted`` (which is stable by construction — Timsort).

    def gaussian_broaden(self, data, sigma):
        """Apply Gaussian broadening to a 1-indexed 1-D data array.

        Each element of ``data`` is treated as a delta-function spike
        centred at its index position.  A normalised Gaussian of width
        ``sigma`` is placed at each spike and accumulated into the
        output array, giving a smoothed spectrum.

        The broadening formula is::

            result[gv] += data[point] * exp(-expTerm) / (sigma * sqrt(pi))

        where ``expTerm = ((point - gv) / 100.0)^2 / sigma^2``.

        The /100.0 factor converts the integer index difference into
        Angstrom units (the caller's bin width is 0.01 Å — see
        ``compute_rpdf``).

        Contributions where ``expTerm >= 20`` are skipped as negligible
        (``exp(-20) ≈ 2e-9``).  ``1 / (sigma*sqrt(pi))`` is the
        normalisation constant so the integrated area is preserved.

        Perl note: the original ``gaussianBroaden`` accepted explicit
        ``$start`` / ``$end`` loop bounds and a separate ``$graph_ref``
        accumulator (written back in place).  The Python version drops
        the bound parameters and instead consumes the whole array, but
        preserves Perl's 1-indexed layout (slot 0 is an unused
        sentinel) so callers can hand it ``self.rpdf`` or any other
        per-bin list built with the module-wide convention directly.
        Callers that need a subset should copy the slice into a fresh
        1-indexed list before calling.

        Parameters
        ----------
        data : list of float
            1-indexed spike values ``[None, v1, v2, ..., vN]``.  Slot 0
            is ignored; slots 1..N carry the data.
        sigma : float
            Gaussian broadening width (standard deviation) in Angstroms.

        Returns
        -------
        list of float
            1-indexed broadened data ``[None, b1, b2, ..., bN]`` of the
            same length as ``data`` (slot 0 left as ``None``).
        """
        n = len(data) - 1  # number of live bins (slot 0 is sentinel)
        sigma_sqrt_pi = sigma * math.sqrt(PI)
        result = [None] + [0.0] * n
        for point in range(1, n + 1):
            curr = data[point]
            for gv in range(1, n + 1):
                # Determine the difference between the point in question and
                # each place on the graph scale (bin width = 0.01 Å).
                diff = (point - gv) / 100.0

                # Determine the exponential term for the broadening factor.
                # This is the "blah" part of gaussian = normalConst * e^(-blah).
                exp_term = (diff * diff) / (sigma * sigma)

                # If gv is far enough from point that the contribution of a
                # Gaussian centred at point would be negligible, don't bother
                # completing the broadening.  Note that 1/sigma_sqrt_pi is the
                # normalisation constant.
                if exp_term < 20.0:
                    result[gv] += curr * math.exp(-exp_term) / sigma_sqrt_pi
        return result

    # ======================================================================
    # Section: Private utility methods
    # ======================================================================

    @staticmethod
    def _prep_line(f, full_line=False):
        """Read one line from an open file and return a list of tokens.

        Perl origin: ``prepLine`` (StructureControl.pm).

        The Perl version accepted three arguments — a file handle, a
        pre-loaded line string, and a regex splitter — so it could work
        either from a file or from an already-read string, with an arbitrary
        delimiter.  In Python neither feature is necessary: callers that
        already hold a string simply call ``str.split()`` directly, and
        whitespace splitting covers every delimiter used in practice.  The
        signature is therefore reduced to a file object plus the ``full_line``
        flag described below.

        The Perl implementation also shifted off a leading empty token that
        ``split`` could produce when the line started with the delimiter (e.g.
        a leading space with a space-based splitter).  Python's ``str.split()``
        with no argument handles this automatically by consuming all leading
        and internal whitespace, so no explicit shift is needed.

        Parameters
        ----------
        f : file object
            Open text file positioned at the next line to read.
        full_line : bool, optional
            When True, return the entire stripped line as a single-element
            list rather than splitting on whitespace.  Used when callers need
            the raw content of a line (e.g. a title string that may contain
            internal spaces).

        Returns
        -------
        list of str
            Tokens from the line, or a one-element list containing the full
            stripped line when ``full_line`` is True.
        """
        line = f.readline().strip()
        if full_line:
            return [line]
        return line.split()

    @staticmethod
    def _find_number(values, target):
        """Search a list for a numeric target and return its 1-based index.

        Perl signature: findNumber(itemToFind, arrayRef, startIndex, endIndex)
        The Perl version accepted an array reference plus explicit startIndex /
        endIndex bounds, allowing a partial-range scan within a larger array.
        The Python version simplifies this: callers pass the desired sub-list
        directly and the entire passed list is searched.

        Parameters
        ----------
        values : list of float or int
            List to search.  Pass a slice if only a sub-range is needed.
        target : float or int
            Value to find (numeric equality, ==).

        Returns
        -------
        int
            1-indexed position of target in values, or 0 if not found.
            (Mirrors Perl's convention: 0 means "not found"; the returned index
            equals what Perl's $arrayIndex would have been for a 1-indexed array
            starting at startIndex=1.)
        """
        # Assume that the item will not be found in the array.
        # (Perl: $found = 0)

        # Search the array for the item.  Return immediately on first match
        # (mirrors Perl's `last` statement inside the foreach loop).
        for i, v in enumerate(values):
            if v == target:
                return i + 1
        return 0

    @staticmethod
    def _find_string(values, target):
        """Search a list for a string target and return its 1-indexed position.

        Perl equivalent: findString($itemToFind, $array, $startIndex, $endIndex)

        The Perl version accepted explicit start/end indices for a bounded search
        of an array reference.  Python simplifies this to a full-list scan, which
        covers all call-sites in the codebase (callers always pass the full range).

        Returns 0 if the target is not found — 0 is the conventional "not found"
        sentinel (Perl: "Assume that the item will not be found in the array.").
        The search stops at the first match (Perl: `last` after assignment).

        Parameters
        ----------
        values : list of str
            List to search (may be 1-indexed with None at [0]).
        target : str
            String to find.

        Returns
        -------
        int
            1-indexed position of target in values, or 0 if not found.
        """
        # Assume that the item will not be found in the array.
        for i, v in enumerate(values):
            if v == target:
                return i + 1  # first match; Perl: last (break) after recording index
        return 0

    @staticmethod
    def _int_remainder(a, b):
        """Return the integer remainder of ``a`` divided by ``b``.

        Perl name: ``intRemainder``

        The Perl source (line 8896) defines two parameters (``$num1``,
        ``$num2``) but its active return statement is::

            return ($num2 - int($num2));

        which computes the *fractional part of* ``$num2`` and completely
        ignores ``$num1`` — almost certainly a copy-paste bug.  The
        commented-out alternative in the Perl source::

            # return (($num1/$num2) - int($num1/$num2));

        computes the fractional part of the *quotient* ``num1/num2``, which
        is also not a standard integer remainder.

        This function has no callers in the current codebase.  The Python
        implementation uses ``a % b``, which matches the conventional
        meaning of "integer remainder" and the function's name.

        Parameters
        ----------
        a, b : int or float

        Returns
        -------
        int or float
            ``a % b``
        """
        return a % b

    @staticmethod
    def _min(a, b):
        """Return the smaller of two values.

        Parameters
        ----------
        a, b : comparable

        Returns
        -------
        Same type as inputs.
        """
        return a if a < b else b
