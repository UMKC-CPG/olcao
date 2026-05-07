"""Tests for element_data.py -- ElementData class.

These tests verify that ElementData correctly reads $OLCAO_DATA/elements.dat
and exposes the data through its attributes and accessor methods.

Run with:
    cd /path/to/olcao
    python src/tests/test_o.py src/tests/test_element_data.py -v

The tests are grouped into sections matching the sections of elements.dat so
it is easy to see which part of the reader is being exercised.  Known values
for a few well-understood elements (H, C, Si, Fe) are checked explicitly so
that any regression in the data file or reader is immediately obvious.
"""

import os
import sys
import pytest

# Allow imports from src/scripts without installing the package.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from element_data import ElementData


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def ed():
    """Return a fully initialised ElementData instance.

    scope='module' means the file is read only once for the entire test
    module, which keeps the test suite fast.
    """
    data_dir = os.environ.get('OLCAO_DATA', '')
    if not data_dir:
        pytest.skip('OLCAO_DATA environment variable is not set.')
    instance = ElementData()
    instance.init_element_data()
    return instance


# ---------------------------------------------------------------------------
# Idempotency
# ---------------------------------------------------------------------------

class TestIdempotency:
    """init_element_data() must be safe to call more than once."""

    def test_double_init_does_not_raise(self, ed):
        ed.init_element_data()  # second call on already-initialised instance

    def test_values_unchanged_after_double_init(self, ed):
        mass_before = ed.atomic_masses[1]
        ed.init_element_data()
        assert ed.atomic_masses[1] == mass_before


# ---------------------------------------------------------------------------
# Header scalars
# ---------------------------------------------------------------------------

class TestHeaderScalars:
    """NUM_ELEMENTS and MAX_QN_L sections."""

    def test_num_elements(self, ed):
        assert ed.num_elements == 103

    def test_max_qn_l(self, ed):
        # OLCAO uses s, p, d, f so max l = 3.
        assert ed.max_qn_l == 3


# ---------------------------------------------------------------------------
# Element identity
# ---------------------------------------------------------------------------

class TestElementNames:
    """ELEMENT_NAMES and ELEMENT_FULL_NAMES sections."""

    def test_index_zero_is_none(self, ed):
        # The 1-indexed convention requires index 0 to be unused (None).
        assert ed.element_names[0] is None
        assert ed.element_full_names[0] is None

    def test_list_length(self, ed):
        # Length should be num_elements + 1 (one None placeholder at 0).
        assert len(ed.element_names) == ed.num_elements + 1
        assert len(ed.element_full_names) == ed.num_elements + 1

    def test_first_element_hydrogen(self, ed):
        assert ed.element_names[1] == 'h'
        assert ed.element_full_names[1] == 'hydrogen'

    def test_last_element_lawrencium(self, ed):
        assert ed.element_names[103] == 'lr'
        assert ed.element_full_names[103] == 'lawrencium'

    def test_known_elements(self, ed):
        # A handful of spot checks across the periodic table.
        expected = {6: 'c', 8: 'o', 14: 'si', 26: 'fe', 79: 'au'}
        for z, name in expected.items():
            assert ed.element_names[z] == name, (
                f'element_names[{z}]: expected {name!r}, '
                f'got {ed.element_names[z]!r}'
            )

    def test_names_are_lowercase(self, ed):
        for z in range(1, ed.num_elements + 1):
            assert ed.element_names[z] == ed.element_names[z].lower(), (
                f'element_names[{z}] is not lower-case: '
                f'{ed.element_names[z]!r}'
            )


# ---------------------------------------------------------------------------
# Physical properties
# ---------------------------------------------------------------------------

class TestAtomicMasses:
    """MASS section."""

    def test_index_zero_is_none(self, ed):
        assert ed.atomic_masses[0] is None

    def test_hydrogen_mass(self, ed):
        assert ed.atomic_masses[1] == pytest.approx(1.0079, rel=1e-4)

    def test_carbon_mass(self, ed):
        assert ed.atomic_masses[6] == pytest.approx(12.011, rel=1e-4)

    def test_silicon_mass(self, ed):
        assert ed.atomic_masses[14] == pytest.approx(28.0855, rel=1e-4)

    def test_all_masses_positive(self, ed):
        for z in range(1, ed.num_elements + 1):
            assert ed.atomic_masses[z] > 0, f'atomic_masses[{z}] <= 0'


class TestCovalentRadii:
    """COVALENT_RADII section."""

    def test_index_zero_is_none(self, ed):
        assert ed.coval_radii[0] is None

    def test_hydrogen_coval_radius(self, ed):
        assert ed.coval_radii[1] == pytest.approx(0.32, rel=1e-4)

    def test_silicon_coval_radius(self, ed):
        assert ed.coval_radii[14] == pytest.approx(1.11, rel=1e-4)

    def test_all_radii_non_negative(self, ed):
        # Some elements (e.g. noble gases like Rn, Z=86) have no measured
        # covalent radius and are stored as 0.0 in elements.dat.
        for z in range(1, ed.num_elements + 1):
            assert ed.coval_radii[z] >= 0, f'coval_radii[{z}] < 0'


class TestAtomicRadii:
    """ATOMIC_RADII section."""

    def test_index_zero_is_none(self, ed):
        assert ed.atomic_radii[0] is None

    def test_all_radii_non_negative(self, ed):
        # Some heavy elements (e.g. Fr, Z=87) have no measured atomic radius
        # and are stored as 0.0 in elements.dat.
        for z in range(1, ed.num_elements + 1):
            assert ed.atomic_radii[z] >= 0, f'atomic_radii[{z}] < 0'


class TestNeutronScattering:
    """NEUT_SCATT section."""

    def test_index_zero_is_none(self, ed):
        assert ed.neut_scatt[0] is None

    def test_length(self, ed):
        assert len(ed.neut_scatt) == ed.num_elements + 1


class TestUJElectrons:
    """NUM_UJ_ELECTRONS section."""

    def test_index_zero_is_none(self, ed):
        assert ed.num_uj_electrons[0] is None

    def test_values_are_non_negative_integers(self, ed):
        for z in range(1, ed.num_elements + 1):
            val = ed.num_uj_electrons[z]
            assert isinstance(val, int), f'num_uj_electrons[{z}] not int'
            assert val >= 0, f'num_uj_electrons[{z}] < 0'


class TestLJPairCoeffs:
    """LJ_PAIR_COEFFS section."""

    def test_index_zero_is_none(self, ed):
        assert ed.lj_pair_coeffs[0] is None

    def test_entries_are_two_tuples(self, ed):
        for z in range(1, ed.num_elements + 1):
            entry = ed.lj_pair_coeffs[z]
            assert len(entry) == 2, (
                f'lj_pair_coeffs[{z}] should have 2 values, got {len(entry)}'
            )

    def test_silicon_coeffs(self, ed):
        eps, sig = ed.lj_pair_coeffs[14]
        assert eps == pytest.approx(0.03, rel=1e-4)
        assert sig == pytest.approx(2.5, rel=1e-4)


# ---------------------------------------------------------------------------
# Basis set descriptors
# ---------------------------------------------------------------------------

class TestCoreOrbitals:
    """CORE_ORBITALS section."""

    def test_index_zero_is_none(self, ed):
        assert ed.core_orbitals[0] is None

    def test_entry_length_equals_max_qn_l_plus_one(self, ed):
        # Each entry is a list [n_s, n_p, n_d, n_f].
        for z in range(1, ed.num_elements + 1):
            assert len(ed.core_orbitals[z]) == ed.max_qn_l + 1

    def test_hydrogen_has_no_core_orbitals(self, ed):
        # H: 1s valence only, no core.
        assert ed.core_orbitals[1] == [0, 0, 0, 0]

    def test_silicon_core_orbitals(self, ed):
        # Si: 1s2 2s2 2p6 core -> s=2, p=1 (shells), d=0, f=0.
        core = ed.get_core_orbitals(14)
        assert core[0] == 2   # s
        assert core[1] == 1   # p
        assert core[2] == 0   # d
        assert core[3] == 0   # f

    def test_all_values_non_negative(self, ed):
        for z in range(1, ed.num_elements + 1):
            for qn_l, val in enumerate(ed.core_orbitals[z]):
                assert val >= 0, (
                    f'core_orbitals[{z}][{qn_l}] is negative: {val}'
                )


class TestValeOrbitals:
    """MB_BEYOND_CORE / FB_BEYOND_MB / EB_BEYOND_FB sections."""

    def test_outer_index_zero_is_none(self, ed):
        assert ed.vale_orbitals[0] is None

    def test_inner_index_zero_is_none_for_each_basis(self, ed):
        for basis in range(1, 4):
            assert ed.vale_orbitals[basis][0] is None, (
                f'vale_orbitals[{basis}][0] should be None'
            )

    def test_entry_length_equals_max_qn_l_plus_one(self, ed):
        for basis in range(1, 4):
            for z in range(1, ed.num_elements + 1):
                entry = ed.vale_orbitals[basis][z]
                assert len(entry) == ed.max_qn_l + 1, (
                    f'vale_orbitals[{basis}][{z}] length wrong'
                )

    def test_accessor_matches_direct_access(self, ed):
        for basis in range(1, 4):
            for z in [1, 6, 14, 26]:
                assert ed.get_vale_orbitals(basis, z) == \
                       ed.vale_orbitals[basis][z]


# ---------------------------------------------------------------------------
# Charge distributions
# ---------------------------------------------------------------------------

class TestCoreCharge:
    """CORE_CHARGE section."""

    def test_index_zero_is_none(self, ed):
        assert ed.core_charge[0] is None

    def test_hydrogen_core_charge_all_zero(self, ed):
        assert all(v == 0.0 for v in ed.core_charge[1])

    def test_all_values_non_negative(self, ed):
        for z in range(1, ed.num_elements + 1):
            for qn_l, val in enumerate(ed.core_charge[z]):
                assert val >= 0, (
                    f'core_charge[{z}][{qn_l}] is negative: {val}'
                )


class TestValeCharge:
    """VALE_CHARGE section."""

    def test_index_zero_is_none(self, ed):
        assert ed.vale_charge[0] is None

    def test_hydrogen_vale_charge(self, ed):
        # H has one valence electron in the s orbital.
        assert ed.vale_charge[1][0] == pytest.approx(1.0, rel=1e-4)  # s
        assert ed.vale_charge[1][1] == pytest.approx(0.0, abs=1e-10)  # p

    def test_accessor_matches_direct_access(self, ed):
        for z in [1, 6, 14, 26]:
            assert ed.get_vale_charge(z) == ed.vale_charge[z]


# ---------------------------------------------------------------------------
# Gaussian expansion parameters
# ---------------------------------------------------------------------------

class TestNumTermsWF:
    """NUM_RWF_TERMS_SPDF section."""

    def test_index_zero_is_none(self, ed):
        assert ed.num_terms_wf[0] is None

    def test_all_values_non_negative(self, ed):
        # Light elements do not have all orbital types.  For example, H (Z=1)
        # has no d or f orbitals, so num_terms_wf[1][2] and [1][3] are 0.
        for z in range(1, ed.num_elements + 1):
            for qn_l, val in enumerate(ed.num_terms_wf[z]):
                assert val >= 0, (
                    f'num_terms_wf[{z}][{qn_l}] is negative: {val}'
                )

    def test_accessor_matches_direct_access(self, ed):
        for z in [1, 6, 14, 26]:
            assert ed.get_num_terms_wf(z) == ed.num_terms_wf[z]


class TestWFAlphaBounds:
    """MIN_RWF_ALPHA and MAX_RWF_ALPHA sections."""

    def test_index_zero_is_none(self, ed):
        assert ed.min_term_wf[0] is None
        assert ed.max_term_wf[0] is None

    def test_min_less_than_max(self, ed):
        for z in range(1, ed.num_elements + 1):
            assert ed.min_term_wf[z] < ed.max_term_wf[z], (
                f'min_term_wf[{z}] >= max_term_wf[{z}]'
            )

    def test_all_positive(self, ed):
        for z in range(1, ed.num_elements + 1):
            assert ed.min_term_wf[z] > 0, f'min_term_wf[{z}] <= 0'
            assert ed.max_term_wf[z] > 0, f'max_term_wf[{z}] <= 0'


class TestNumTermsPot:
    """NUM_POT_TERMS section."""

    def test_index_zero_is_none(self, ed):
        assert ed.num_terms_pot[0] is None

    def test_all_values_positive(self, ed):
        for z in range(1, ed.num_elements + 1):
            assert ed.num_terms_pot[z] > 0, (
                f'num_terms_pot[{z}] is not positive'
            )


class TestPotAlphaBounds:
    """MIN_POT_ALPHA and MAX_POT_ALPHA sections."""

    def test_index_zero_is_none(self, ed):
        assert ed.min_term_pot[0] is None
        assert ed.max_term_pot[0] is None

    def test_min_less_than_max(self, ed):
        for z in range(1, ed.num_elements + 1):
            assert ed.min_term_pot[z] < ed.max_term_pot[z], (
                f'min_term_pot[{z}] >= max_term_pot[{z}]'
            )


# ---------------------------------------------------------------------------
# Orbital Gaussian selection masks
# ---------------------------------------------------------------------------

class TestOrbitalTerms:
    """GAUSS_TERMS_SPDF section."""

    def test_outer_dimension_equals_max_qn_l_plus_one(self, ed):
        assert len(ed.orbital_terms) == ed.max_qn_l + 1

    def test_inner_index_zero_is_none(self, ed):
        for qn_l in range(ed.max_qn_l + 1):
            assert ed.orbital_terms[qn_l][0] is None

    def test_entries_are_strings(self, ed):
        for qn_l in range(ed.max_qn_l + 1):
            for z in range(1, ed.num_elements + 1):
                val = ed.orbital_terms[qn_l][z]
                assert isinstance(val, str), (
                    f'orbital_terms[{qn_l}][{z}] is not a str: {val!r}'
                )

    def test_entries_not_empty(self, ed):
        for qn_l in range(ed.max_qn_l + 1):
            for z in range(1, ed.num_elements + 1):
                assert ed.orbital_terms[qn_l][z] != '', (
                    f'orbital_terms[{qn_l}][{z}] is empty'
                )

    def test_accessor_matches_direct_access(self, ed):
        for qn_l in range(ed.max_qn_l + 1):
            for z in [1, 6, 14, 26]:
                assert ed.get_orbital_terms(qn_l, z) == \
                       ed.orbital_terms[qn_l][z]


# ---------------------------------------------------------------------------
# Visualization colors
# ---------------------------------------------------------------------------

class TestColorVTK:
    """VTK_COLOR section."""

    def test_index_zero_is_none(self, ed):
        assert ed.color_vtk[0] is None

    def test_entries_have_four_channels(self, ed):
        for z in range(1, ed.num_elements + 1):
            assert len(ed.color_vtk[z]) == 4, (
                f'color_vtk[{z}] does not have 4 channels'
            )

    def test_rgb_in_range_0_255(self, ed):
        for z in range(1, ed.num_elements + 1):
            for ch, val in enumerate(ed.color_vtk[z][:3]):
                assert 0 <= val <= 255, (
                    f'color_vtk[{z}][{ch}] out of 0-255 range: {val}'
                )

    def test_alpha_is_255(self, ed):
        # All elements should be fully opaque.
        for z in range(1, ed.num_elements + 1):
            assert ed.color_vtk[z][3] == 255, (
                f'color_vtk[{z}] alpha is not 255: {ed.color_vtk[z][3]}'
            )

    def test_known_colors(self, ed):
        # CPK color scheme spot checks.
        assert ed.color_vtk[1][:3] == [255.0, 255.0, 255.0]   # H: white
        assert ed.color_vtk[6][:3] == [144.0, 144.0, 144.0]   # C: grey


class TestColorDX:
    """ODX_COLOR section."""

    def test_index_zero_is_none(self, ed):
        assert ed.color_dx[0] is None

    def test_length(self, ed):
        assert len(ed.color_dx) == ed.num_elements + 1


class TestGreyDX:
    """ODX_GREY section."""

    def test_index_zero_is_none(self, ed):
        assert ed.grey_dx[0] is None

    def test_length(self, ed):
        assert len(ed.grey_dx) == ed.num_elements + 1


# ---------------------------------------------------------------------------
# Modifier: apply_bond_factor
# ---------------------------------------------------------------------------

class TestApplyBondFactor:
    """ElementData.apply_bond_factor() scales covalent radii in place."""

    def test_factor_one_is_no_op(self, ed):
        original = ed.coval_radii[14]
        ed.apply_bond_factor(1.0)
        assert ed.coval_radii[14] == original

    def test_scaling_and_inverse(self, ed):
        # Scale up then back down; result should round-trip to original.
        original = ed.coval_radii[6]
        ed.apply_bond_factor(1.5)
        assert ed.coval_radii[6] == pytest.approx(original * 1.5, rel=1e-10)
        ed.apply_bond_factor(1.0 / 1.5)
        assert ed.coval_radii[6] == pytest.approx(original, rel=1e-10)

    def test_all_elements_scaled(self, ed):
        # Collect originals, scale, verify every element changed.
        originals = ed.coval_radii[1:]          # drop the None at index 0
        ed.apply_bond_factor(2.0)
        for z in range(1, ed.num_elements + 1):
            assert ed.coval_radii[z] == pytest.approx(
                originals[z - 1] * 2.0, rel=1e-10
            )
        ed.apply_bond_factor(0.5)               # restore


# ---------------------------------------------------------------------------
# Lookup: get_element_z
# ---------------------------------------------------------------------------

class TestGetElementZ:
    """ElementData.get_element_z() returns atomic number by name."""

    def test_hydrogen_lower(self, ed):
        assert ed.get_element_z('h') == 1

    def test_hydrogen_upper(self, ed):
        assert ed.get_element_z('H') == 1

    def test_silicon_mixed_case(self, ed):
        assert ed.get_element_z('Si') == 14
        assert ed.get_element_z('SI') == 14
        assert ed.get_element_z('si') == 14

    def test_iron(self, ed):
        assert ed.get_element_z('fe') == 26

    def test_gold(self, ed):
        assert ed.get_element_z('au') == 79

    def test_lawrencium(self, ed):
        assert ed.get_element_z('lr') == 103

    def test_unknown_element_returns_none(self, ed):
        assert ed.get_element_z('xx') is None
        assert ed.get_element_z('') is None
