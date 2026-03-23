"""test_org_radwavefn.py -- Tests for Grasp2K wavefunction parser.

Uses Carbon MB test data from jobs/atoms/MB/C/.
"""

import os
import sys
import pytest

TESTS_DIR = os.path.dirname(__file__)
SCRIPTS_DIR = os.path.join(os.path.dirname(TESTS_DIR), 'scripts')
REPO_ROOT = os.path.dirname(os.path.dirname(TESTS_DIR))
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

CARBON_MB_DIR = os.path.join(REPO_ROOT, 'jobs', 'atoms', 'MB', 'C')

# Skip entire module if test data not available
pytestmark = pytest.mark.skipif(
    not os.path.isfile(os.path.join(CARBON_MB_DIR, 'rwfn.out')),
    reason='Carbon MB test data (jobs/atoms/MB/C/) not found'
)


@pytest.fixture(scope='module')
def carbon_isodata():
    from org_radwavefn import isodata_file
    return isodata_file(os.path.join(CARBON_MB_DIR, 'isodata'))


@pytest.fixture(scope='module')
def carbon_orbitals():
    from org_radwavefn import orbital_file
    return orbital_file(os.path.join(CARBON_MB_DIR, 'rwfn.out'))


@pytest.fixture(scope='module')
def carbon_system():
    from org_radwavefn import atomic_system
    return atomic_system(
        os.path.join(CARBON_MB_DIR, 'isodata'),
        os.path.join(CARBON_MB_DIR, 'rwfn.out')
    )


# ── isodata parsing ──────────────────────────────────────────────────

class TestIsodataParsing:
    def test_atomic_number(self, carbon_isodata):
        assert carbon_isodata.atomic_number == 6

    def test_mass_number(self, carbon_isodata):
        assert carbon_isodata.mass_number == 12

    def test_nuclear_mass(self, carbon_isodata):
        assert carbon_isodata.nucleus_mass == 12.0

    def test_nuclear_spin(self, carbon_isodata):
        assert carbon_isodata.nuclear_spin == 0.0

    def test_fermi_params_positive(self, carbon_isodata):
        assert carbon_isodata.fermi_dist_param_a > 0
        assert carbon_isodata.fermi_dist_param_c > 0


# ── rwfn.out orbital parsing ─────────────────────────────────────────

class TestOrbitalParsing:
    def test_orbital_count(self, carbon_orbitals):
        # Carbon MB: 1s_{-1}, 2s_{-1}, 2p_{1} (j=1/2), 2p_{-2} (j=3/2)
        assert carbon_orbitals.num_orbitals == 4

    def test_n_quantum_numbers(self, carbon_orbitals):
        assert carbon_orbitals.n_qnums == [1, 2, 2, 2]

    def test_k_quantum_numbers(self, carbon_orbitals):
        assert carbon_orbitals.k_qnums == [-1, -1, 1, -2]

    def test_max_l(self, carbon_orbitals):
        assert carbon_orbitals.max_l == 1  # p orbitals present

    def test_max_n(self, carbon_orbitals):
        assert carbon_orbitals.max_n == 2

    def test_orbital_list_format(self, carbon_orbitals):
        assert '1_-1' in carbon_orbitals.orbital_list
        assert '2_-1' in carbon_orbitals.orbital_list
        assert '2_1' in carbon_orbitals.orbital_list
        assert '2_-2' in carbon_orbitals.orbital_list

    def test_radial_grid_nonempty(self, carbon_orbitals):
        for i in range(carbon_orbitals.num_orbitals):
            assert len(carbon_orbitals.radial_comps[i]) > 0

    def test_large_components_match_grid_length(self, carbon_orbitals):
        for i in range(carbon_orbitals.num_orbitals):
            assert len(carbon_orbitals.large_comps[i]) == len(carbon_orbitals.radial_comps[i])

    def test_small_components_match_grid_length(self, carbon_orbitals):
        for i in range(carbon_orbitals.num_orbitals):
            assert len(carbon_orbitals.small_comps[i]) == len(carbon_orbitals.radial_comps[i])

    def test_eigenvalues_nonzero(self, carbon_orbitals):
        """All bound-state eigenvalues should be nonzero."""
        for ev in carbon_orbitals.energy_eigenvalues:
            assert ev != 0.0

    def test_interpolated_grid_is_sorted(self, carbon_orbitals):
        grid = carbon_orbitals.interpolated_radial_grid
        for i in range(len(grid) - 1):
            assert grid[i] <= grid[i + 1]


# ── single_component fitting mode ────────────────────────────────────

class TestSingleComponent:
    @pytest.fixture(autouse=True)
    def _setup(self):
        """Create a fresh orbital_file for each test (fitting mutates state)."""
        from org_radwavefn import orbital_file
        self.orb = orbital_file(os.path.join(CARBON_MB_DIR, 'rwfn.out'))

    def test_single_component_count(self):
        self.orb.single_component()
        # Carbon MB: 1s, 2s, 2p → 3 fitting functions
        assert len(self.orb.fitting_orbital_names) == 3

    def test_single_component_names(self):
        self.orb.single_component()
        assert '1s' in self.orb.fitting_orbital_names
        assert '2s' in self.orb.fitting_orbital_names
        assert '2p' in self.orb.fitting_orbital_names

    def test_single_component_lvals(self):
        self.orb.single_component()
        # 2 s-type (1s, 2s), 1 p-type (2p)
        assert self.orb.fitting_lvals[0] == 2  # s
        assert self.orb.fitting_lvals[1] == 1  # p

    def test_fitting_max_l(self):
        self.orb.single_component()
        assert self.orb.fitting_max_l == 1  # max l is p


# ── atomic_system convenience class ──────────────────────────────────

class TestAtomicSystem:
    def test_element_symbol(self, carbon_system):
        assert carbon_system.element == 'C'

    def test_has_orbital_info(self, carbon_system):
        assert carbon_system.orbital_info.num_orbitals == 4

    def test_has_atomic_info(self, carbon_system):
        assert carbon_system.atomic_info.atomic_number == 6
