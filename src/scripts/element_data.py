#!/usr/bin/env python3

"""element_data.py -- Periodic table and OLCAO basis/potential parameters.

This module is the Python equivalent of the Perl ElementData.pm module.  It
provides a single class, ElementData, that reads the OLCAO elements database
file ($OLCAO_DATA/elements.dat) and exposes all of its data through a clean,
documented interface.

TYPICAL USAGE::

    from element_data import ElementData

    ed = ElementData()
    ed.init_element_data()          # reads $OLCAO_DATA/elements.dat

    z = ed.get_element_z('si')      # returns 14
    mass = ed.atomic_masses[14]     # 28.086 (direct attribute access)
    radii = ed.coval_radii          # full 1-indexed list

INDEXING CONVENTION:
    All element-indexed arrays are 1-indexed so that the array index equals
    the atomic number Z (H=1, He=2, ..., Lr=103).  Index 0 of every such
    list is None and must not be used.  This matches the convention used
    throughout the OLCAO Fortran source and the original Perl scripts, and
    it keeps loop indices directly interpretable as physical quantities.

    The quantum number l (qn_l) is 0-indexed everywhere: s=0, p=1, d=2, f=3.

MULTI-DIMENSIONAL ARRAY LAYOUTS:
    core_orbitals[element]          -> list [s, p, d, f] of ints
    vale_charge[element]            -> list [s, p, d, f] of floats
    core_charge[element]            -> list [s, p, d, f] of floats
    num_terms_wf[element]           -> list [s, p, d, f] of ints
    vale_orbitals[basis][element]   -> list [s, p, d, f] of ints
                                       basis: 1=MB, 2=FB, 3=EB
    orbital_terms[qn_l][element]    -> str (Gaussian selection mask)
    lj_pair_coeffs[element]         -> tuple (epsilon, sigma)
    color_vtk[element]              -> list [r, g, b, alpha]
"""

import os


class ElementData:
    """Stores periodic table data and OLCAO basis/potential parameters.

    All data are read from the plain-text file $OLCAO_DATA/elements.dat.
    The file is structured as tagged sections (e.g. NUM_ELEMENTS, MASS, ...)
    each followed by one value per element.

    After construction, call init_element_data() before accessing any
    attributes.  Calling init_element_data() more than once is safe -- the
    second and later calls return immediately without re-reading the file.

    Attributes
    ----------
    num_elements : int
        Total number of elements in the database (typically 103).
    max_qn_l : int
        Highest angular momentum quantum number present (typically 3 for f).
    element_names : list of str
        Abbreviated element names, lower-case (e.g. 'h', 'he', 'li').
        1-indexed; element_names[0] is None.
    element_full_names : list of str
        Full element names, lower-case (e.g. 'hydrogen').
        1-indexed; element_full_names[0] is None.
    atomic_masses : list of float
        Atomic mass in atomic mass units.  1-indexed.
    coval_radii : list of float
        Covalent radius in Angstroms.  1-indexed.  Can be scaled in place
        by apply_bond_factor().
    atomic_radii : list of float
        Atomic radius in Angstroms.  1-indexed.
    neut_scatt : list of float
        Neutron scattering factor (Neutron News, Vol. 3, No. 3, 1992).
        1-indexed.
    num_uj_electrons : list of int
        Number of electrons in the highest d or f orbital in the ground
        state.  Used for DFT+U/J calculations.  1-indexed.
    lj_pair_coeffs : list of tuple
        Lennard-Jones pair coefficients (epsilon, sigma) for LAMMPS.
        1-indexed; lj_pair_coeffs[element] = (epsilon, sigma).
    core_orbitals : list of list
        core_orbitals[element] = [n_s, n_p, n_d, n_f] counts of core
        orbitals.  1-indexed on element; qn_l is 0-indexed.
    vale_orbitals : list
        vale_orbitals[basis][element] = [n_s, n_p, n_d, n_f] counts of
        valence orbitals beyond the previous basis set.
        basis 1=MB (minimal beyond core), 2=FB (full beyond MB),
        3=EB (extended beyond FB).  Both basis and element are 1-indexed.
    core_charge : list of list
        core_charge[element][qn_l] = number of core electrons in orbital
        type qn_l.  1-indexed on element; qn_l is 0-indexed.
    vale_charge : list of list
        vale_charge[element][qn_l] = number of valence electrons in orbital
        type qn_l.  1-indexed on element; qn_l is 0-indexed.
    num_terms_wf : list of list
        num_terms_wf[element][qn_l] = number of Gaussian terms used to
        represent the radial wave function for orbital type qn_l.
        1-indexed on element; qn_l is 0-indexed.
    min_term_wf : list of float
        Minimum Gaussian exponent (alpha) for the radial wave function
        expansion.  1-indexed.
    max_term_wf : list of float
        Maximum Gaussian exponent (alpha) for the radial wave function
        expansion.  1-indexed.
    num_terms_pot : list of int
        Number of Gaussian terms used to represent the potential function.
        1-indexed.
    min_term_pot : list of float
        Minimum Gaussian exponent (alpha) for the potential function.
        1-indexed.
    max_term_pot : list of float
        Maximum Gaussian exponent (alpha) for the potential function.
        1-indexed.
    orbital_terms : list of list of str
        orbital_terms[qn_l][element] = compact boolean string identifying
        which Gaussians in the orbital expansion are active (e.g. '11010').
        qn_l is 0-indexed; element is 1-indexed (orbital_terms[qn_l][0]
        is None).
    color_vtk : list of list
        color_vtk[element] = [r, g, b, alpha] CPK color used in VTK
        visualisation.  All four values are on the 0-255 scale.
        All elements are fully opaque (alpha=255).  1-indexed.
    color_dx : list of float
        OpenDX color index (scale 1-100) per element.  1-indexed.
    grey_dx : list of float
        OpenDX greyscale index (scale 1-100) per element.  1-indexed.
    """

    def __init__(self):
        self._init_data = False

        # Scalars
        self.num_elements = 0
        self.max_qn_l = 0

        # 1-indexed lists (index 0 is a None placeholder; do not use it)
        self.element_names = [None]
        self.element_full_names = [None]
        self.atomic_masses = [None]
        self.coval_radii = [None]
        self.atomic_radii = [None]
        self.neut_scatt = [None]
        self.num_uj_electrons = [None]
        self.lj_pair_coeffs = [None]    # entries: (epsilon, sigma)
        self.core_orbitals = [None]     # entries: [s, p, d, f]
        self.core_charge = [None]       # entries: [s, p, d, f]
        self.vale_charge = [None]       # entries: [s, p, d, f]
        self.num_terms_wf = [None]      # entries: [s, p, d, f]
        self.min_term_wf = [None]
        self.max_term_wf = [None]
        self.num_terms_pot = [None]
        self.min_term_pot = [None]
        self.max_term_pot = [None]
        self.color_vtk = [None]         # entries: [r, g, b, alpha], all 0-255
        self.color_dx = [None]
        self.grey_dx = [None]

        # vale_orbitals[basis][element][qn_l]
        #   Outer index 0 is None (basis is 1-indexed: 1=MB, 2=FB, 3=EB).
        #   Middle index 0 is None (element is 1-indexed).
        #   Inner index is qn_l, 0-indexed (s=0, p=1, d=2, f=3).
        self.vale_orbitals = [None, None, None, None]

        # orbital_terms[qn_l][element]
        #   Outer index is qn_l, 0-indexed.  Length set in init_element_data()
        #   once max_qn_l is known.
        #   Inner index 0 is None (element is 1-indexed).
        self.orbital_terms = []

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _prep_line(f, full_line=False):
        """Read one line from open file f and return a list of tokens.

        Strips leading/trailing whitespace before splitting.  Blank lines
        (e.g. trailing newlines) are skipped automatically because strip()
        on an empty line returns '' and split() on '' returns [].

        Parameters
        ----------
        f : file object
            An open text file positioned at the next line to read.
        full_line : bool
            When True, return the entire stripped line as a single-element
            list rather than splitting on whitespace.  This is needed for
            the orbital_terms strings, which are compact boolean masks
            (e.g. '11010') that must be preserved verbatim.

        Returns
        -------
        list of str
        """
        line = f.readline().strip()
        if full_line:
            return [line]
        return line.split()

    @staticmethod
    def _expect_tag(values, expected, filename):
        """Raise ValueError when the first token does not match expected.

        This is called immediately after reading each section header line
        in elements.dat.  If the file is malformed or a future version of
        elements.dat reorders sections, this gives a clear error message.

        Parameters
        ----------
        values : list of str
            Tokens from the header line (output of _prep_line).
        expected : str
            The tag string that must appear as values[0].
        filename : str
            Path to the file being read; included in the error message.
        """
        if values[0] != expected:
            raise ValueError(
                f"Expected tag '{expected}' in {filename}, "
                f"got '{values[0]}'"
            )

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------

    def init_element_data(self):
        """Read all element data from $OLCAO_DATA/elements.dat.

        Populates every attribute listed in the class docstring.  The file
        is plain ASCII with tagged sections; each section header is a
        single keyword line followed by one data line per element.

        This method is idempotent -- calling it more than once is safe and
        the second call returns without re-reading the file.

        Raises
        ------
        KeyError
            If the OLCAO_DATA environment variable is not set.
        FileNotFoundError
            If elements.dat cannot be found at $OLCAO_DATA/elements.dat.
        ValueError
            If an unexpected section tag is encountered in the file.
        """
        if self._init_data:
            return
        self._init_data = True

        data_dir = os.environ.get('OLCAO_DATA', '')
        file_path = os.path.join(data_dir, 'elements.dat')

        with open(file_path, 'r') as f:
            fn = file_path   # shorthand used in _expect_tag calls

            # ----------------------------------------------------------
            # Header scalars
            # ----------------------------------------------------------

            self._expect_tag(self._prep_line(f), 'NUM_ELEMENTS', fn)
            self.num_elements = int(self._prep_line(f)[0])

            self._expect_tag(self._prep_line(f), 'MAX_QN_L', fn)
            self.max_qn_l = int(self._prep_line(f)[0])

            # Convenience shorthands used throughout this method.
            n = self.num_elements
            l = self.max_qn_l

            # ----------------------------------------------------------
            # Element identity
            # ----------------------------------------------------------

            # Abbreviated names: 'h', 'he', 'li', ...
            self._expect_tag(self._prep_line(f), 'ELEMENT_NAMES', fn)
            for _ in range(1, n + 1):
                self.element_names.append(self._prep_line(f)[0])

            # Full names: 'hydrogen', 'helium', ...
            self._expect_tag(self._prep_line(f), 'ELEMENT_FULL_NAMES', fn)
            for _ in range(1, n + 1):
                self.element_full_names.append(self._prep_line(f)[0])

            # ----------------------------------------------------------
            # Physical properties
            # ----------------------------------------------------------

            # Atomic masses in atomic mass units (amu)
            self._expect_tag(self._prep_line(f), 'MASS', fn)
            for _ in range(1, n + 1):
                self.atomic_masses.append(float(self._prep_line(f)[0]))

            # Covalent radii in Angstroms
            self._expect_tag(self._prep_line(f), 'COVALENT_RADII', fn)
            for _ in range(1, n + 1):
                self.coval_radii.append(float(self._prep_line(f)[0]))

            # Atomic radii in Angstroms
            self._expect_tag(self._prep_line(f), 'ATOMIC_RADII', fn)
            for _ in range(1, n + 1):
                self.atomic_radii.append(float(self._prep_line(f)[0]))

            # Neutron scattering factors
            #   Source: Neutron News, Vol. 3, No. 3, 1992, pp. 29-37
            self._expect_tag(self._prep_line(f), 'NEUT_SCATT', fn)
            for _ in range(1, n + 1):
                self.neut_scatt.append(float(self._prep_line(f)[0]))

            # Number of electrons in the highest d or f orbital (ground state)
            #   Used for DFT+U/J (Hubbard correction) calculations.
            self._expect_tag(self._prep_line(f), 'NUM_UJ_ELECTRONS', fn)
            for _ in range(1, n + 1):
                self.num_uj_electrons.append(int(self._prep_line(f)[0]))

            # Lennard-Jones pair coefficients for LAMMPS MD simulations.
            #   Each entry is stored as a tuple (epsilon, sigma).
            self._expect_tag(self._prep_line(f), 'LJ_PAIR_COEFFS', fn)
            for _ in range(1, n + 1):
                vals = self._prep_line(f)
                self.lj_pair_coeffs.append(
                    (float(vals[0]), float(vals[1]))
                )

            # ----------------------------------------------------------
            # Basis set descriptors
            # ----------------------------------------------------------

            # Core orbital counts: core_orbitals[element] = [n_s, n_p, n_d, n_f]
            self._expect_tag(self._prep_line(f), 'CORE_ORBITALS', fn)
            for _ in range(1, n + 1):
                vals = self._prep_line(f)
                self.core_orbitals.append(
                    [int(vals[qn_l]) for qn_l in range(l + 1)]
                )

            # Valence orbital counts beyond each successive basis level.
            #   Three passes: MB beyond core, FB beyond MB, EB beyond FB.
            #   vale_orbitals[basis][element][qn_l]
            expected_tags = ['MB_BEYOND_CORE', 'FB_BEYOND_MB', 'EB_BEYOND_FB']
            for basis in range(1, 4):
                self._expect_tag(
                    self._prep_line(f), expected_tags[basis - 1], fn
                )
                basis_list = [None]  # 1-indexed; slot 0 unused
                for _ in range(1, n + 1):
                    vals = self._prep_line(f)
                    basis_list.append(
                        [int(vals[qn_l]) for qn_l in range(l + 1)]
                    )
                self.vale_orbitals[basis] = basis_list

            # ----------------------------------------------------------
            # Charge distributions
            # ----------------------------------------------------------

            # Core charge per s,p,d,f orbital type
            self._expect_tag(self._prep_line(f), 'CORE_CHARGE', fn)
            for _ in range(1, n + 1):
                vals = self._prep_line(f)
                self.core_charge.append(
                    [float(vals[qn_l]) for qn_l in range(l + 1)]
                )

            # Valence charge per s,p,d,f orbital type
            self._expect_tag(self._prep_line(f), 'VALE_CHARGE', fn)
            for _ in range(1, n + 1):
                vals = self._prep_line(f)
                self.vale_charge.append(
                    [float(vals[qn_l]) for qn_l in range(l + 1)]
                )

            # ----------------------------------------------------------
            # Gaussian expansion parameters (radial wave functions)
            # ----------------------------------------------------------

            # Number of Gaussian terms per s,p,d,f orbital type
            self._expect_tag(self._prep_line(f), 'NUM_RWF_TERMS_SPDF', fn)
            for _ in range(1, n + 1):
                vals = self._prep_line(f)
                self.num_terms_wf.append(
                    [int(vals[qn_l]) for qn_l in range(l + 1)]
                )

            # Min/max Gaussian exponent alpha for the wave function expansion
            self._expect_tag(self._prep_line(f), 'MIN_RWF_ALPHA', fn)
            for _ in range(1, n + 1):
                self.min_term_wf.append(float(self._prep_line(f)[0]))

            self._expect_tag(self._prep_line(f), 'MAX_RWF_ALPHA', fn)
            for _ in range(1, n + 1):
                self.max_term_wf.append(float(self._prep_line(f)[0]))

            # ----------------------------------------------------------
            # Gaussian expansion parameters (potential functions)
            # ----------------------------------------------------------

            # Number of Gaussian terms for the potential function
            self._expect_tag(self._prep_line(f), 'NUM_POT_TERMS', fn)
            for _ in range(1, n + 1):
                self.num_terms_pot.append(int(self._prep_line(f)[0]))

            # Min/max Gaussian exponent alpha for the potential expansion
            self._expect_tag(self._prep_line(f), 'MIN_POT_ALPHA', fn)
            for _ in range(1, n + 1):
                self.min_term_pot.append(float(self._prep_line(f)[0]))

            self._expect_tag(self._prep_line(f), 'MAX_POT_ALPHA', fn)
            for _ in range(1, n + 1):
                self.max_term_pot.append(float(self._prep_line(f)[0]))

            # ----------------------------------------------------------
            # Orbital Gaussian selection masks
            # ----------------------------------------------------------

            # orbital_terms[qn_l][element] = compact boolean string.
            #   Each character in the string is '0' or '1' indicating
            #   whether the corresponding Gaussian term is included in the
            #   basis expansion for that orbital type.  For example,
            #   '11010' means terms 1, 2, and 4 are active out of 5 total.
            #
            #   The file is written element-major (all qn_l lines for
            #   element 1, then element 2, etc.), but the array is stored
            #   qn_l-major so that indexing is orbital_terms[qn_l][element].
            self._expect_tag(self._prep_line(f), 'GAUSS_TERMS_SPDF', fn)
            self.orbital_terms = [
                [None] * (n + 1) for _ in range(l + 1)
            ]
            for element in range(1, n + 1):
                for qn_l in range(l + 1):
                    self.orbital_terms[qn_l][element] = (
                        self._prep_line(f, full_line=True)[0]
                    )

            # ----------------------------------------------------------
            # Visualization colors
            # ----------------------------------------------------------

            # VTK colors using the CPK chemistry color scheme.
            #   color_vtk[element] = [r, g, b, alpha], all on the 0-255 scale.
            #   All elements are fully opaque (alpha=255).
            self._expect_tag(self._prep_line(f), 'VTK_COLOR', fn)
            for _ in range(1, n + 1):
                vals = self._prep_line(f)
                self.color_vtk.append([float(vals[i]) for i in range(4)])

            # OpenDX color index (scale 1-100)
            self._expect_tag(self._prep_line(f), 'ODX_COLOR', fn)
            for _ in range(1, n + 1):
                self.color_dx.append(float(self._prep_line(f)[0]))

            # OpenDX greyscale index (scale 1-100)
            self._expect_tag(self._prep_line(f), 'ODX_GREY', fn)
            for _ in range(1, n + 1):
                self.grey_dx.append(float(self._prep_line(f)[0]))

    # ------------------------------------------------------------------
    # Modifiers
    # ------------------------------------------------------------------

    def apply_bond_factor(self, bonding_factor):
        """Scale all covalent radii by a multiplicative factor.

        This is used to tighten or loosen the bond-detection criterion
        throughout OLCAO.  A factor greater than 1.0 allows longer bonds
        to be detected; a factor less than 1.0 restricts detection to
        shorter bonds.  Has no effect if bonding_factor is exactly 1.0.

        Parameters
        ----------
        bonding_factor : float
            Multiplicative scale factor applied to every covalent radius.
        """
        if bonding_factor != 1.0:
            for element in range(1, len(self.coval_radii)):
                self.coval_radii[element] *= bonding_factor

    # ------------------------------------------------------------------
    # Lookup methods
    # ------------------------------------------------------------------

    def get_element_z(self, element_name):
        """Return the atomic number Z for a given element abbreviation.

        The lookup is case-insensitive, so 'Si', 'SI', and 'si' all return
        14.  Returns None if the element abbreviation is not found.

        Parameters
        ----------
        element_name : str
            Element abbreviation (e.g. 'h', 'Si', 'FE').

        Returns
        -------
        int or None
            Atomic number (1-indexed) or None if not found.
        """
        name = element_name.lower()
        for element in range(1, self.num_elements + 1):
            if name == self.element_names[element]:
                return element
        return None

    # ------------------------------------------------------------------
    # Accessor methods
    # These mirror the Perl ElementData getter subroutines and provide a
    # stable named API for StructureControl and other callers.  For simple
    # scalar or list lookups, direct attribute access is equally valid.
    # ------------------------------------------------------------------

    def get_num_elements(self):
        """Return the total number of elements in the database."""
        return self.num_elements

    def get_max_qn_l(self):
        """Return the highest angular momentum quantum number (typically 3)."""
        return self.max_qn_l

    def get_core_orbitals(self, element):
        """Return [n_s, n_p, n_d, n_f] core orbital counts for element Z.

        Parameters
        ----------
        element : int
            Atomic number (1-indexed).
        """
        return self.core_orbitals[element]

    def get_vale_orbitals(self, basis, element):
        """Return [n_s, n_p, n_d, n_f] valence orbital counts for a basis.

        Parameters
        ----------
        basis : int
            1 = minimal beyond core (MB),
            2 = full beyond MB (FB),
            3 = extended beyond FB (EB).
        element : int
            Atomic number (1-indexed).
        """
        return self.vale_orbitals[basis][element]

    def get_orbital_terms(self, qn_l, element):
        """Return the Gaussian selection string for orbital qn_l of element.

        Parameters
        ----------
        qn_l : int
            Angular momentum quantum number, 0-indexed (s=0, p=1, d=2, f=3).
        element : int
            Atomic number (1-indexed).
        """
        return self.orbital_terms[qn_l][element]

    def get_num_terms_wf(self, element):
        """Return [n_s, n_p, n_d, n_f] Gaussian term counts for wave function.

        Parameters
        ----------
        element : int
            Atomic number (1-indexed).
        """
        return self.num_terms_wf[element]

    def get_vale_charge(self, element):
        """Return [q_s, q_p, q_d, q_f] valence electron counts for element.

        Parameters
        ----------
        element : int
            Atomic number (1-indexed).
        """
        return self.vale_charge[element]
