"""test_dirac_helpers.py -- Tests for Dirac quantum number conversion helpers
and ElementData integration with the Grasp2K pipeline scripts.
"""

import os
import sys
import pytest

# Ensure src/scripts is importable
TESTS_DIR = os.path.dirname(__file__)
SCRIPTS_DIR = os.path.join(os.path.dirname(TESTS_DIR), 'scripts')
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)


# -- Dirac helper tests ---------------------------------------------------

class TestLktol:
    """lktol: kappa -> angular momentum l for the large Dirac component."""

    def test_kappa_neg1_gives_l0(self):
        from org_radwavefn import lktol
        assert lktol(-1) == 0  # s orbital

    def test_kappa_1_gives_l1(self):
        from org_radwavefn import lktol
        assert lktol(1) == 1  # p orbital

    def test_kappa_neg2_gives_l1(self):
        from org_radwavefn import lktol
        assert lktol(-2) == 1  # p orbital

    def test_kappa_2_gives_l2(self):
        from org_radwavefn import lktol
        assert lktol(2) == 2  # d orbital

    def test_kappa_neg3_gives_l2(self):
        from org_radwavefn import lktol
        assert lktol(-3) == 2  # d orbital

    def test_kappa_3_gives_l3(self):
        from org_radwavefn import lktol
        assert lktol(3) == 3  # f orbital


class TestSktol:
    """sktol: kappa -> angular momentum l for the small Dirac component."""

    def test_kappa_neg1_gives_l1(self):
        from org_radwavefn import sktol
        assert sktol(-1) == 1  # small component of s is p

    def test_kappa_1_gives_l0(self):
        from org_radwavefn import sktol
        assert sktol(1) == 0

    def test_kappa_neg2_gives_l2(self):
        from org_radwavefn import sktol
        assert sktol(-2) == 2

    def test_kappa_2_gives_l1(self):
        from org_radwavefn import sktol
        assert sktol(2) == 1


class TestLtok:
    """ltok: l -> both kappa values for the large component."""

    def test_l0_s_orbital(self):
        from org_radwavefn import ltok
        assert ltok(0) == (-1, -1)  # only one kappa for s

    def test_l1_p_orbital(self):
        from org_radwavefn import ltok
        assert ltok(1) == (1, -2)  # kappa=1 (j=1/2) and kappa=-2 (j=3/2)

    def test_l2_d_orbital(self):
        from org_radwavefn import ltok
        assert ltok(2) == (2, -3)

    def test_l3_f_orbital(self):
        from org_radwavefn import ltok
        assert ltok(3) == (3, -4)


class TestLtolsym:
    """ltolsym: l -> spectroscopic symbol."""

    def test_all_symbols(self):
        from org_radwavefn import ltolsym
        assert ltolsym(0) == 's'
        assert ltolsym(1) == 'p'
        assert ltolsym(2) == 'd'
        assert ltolsym(3) == 'f'
        assert ltolsym(4) == 'g'


# -- ElementData integration tests ----------------------------------------

@pytest.fixture(scope='module')
def element_data():
    """Initialize ElementData from $OLCAO_DATA/elements.dat."""
    data_dir = os.environ.get('OLCAO_DATA', '')
    if not data_dir:
        # Auto-detect from repo layout
        repo_root = os.path.dirname(os.path.dirname(TESTS_DIR))
        candidate = os.path.join(repo_root, 'share')
        if os.path.isdir(candidate):
            os.environ['OLCAO_DATA'] = candidate
            data_dir = candidate
        else:
            pytest.skip('OLCAO_DATA not set')

    from element_data import ElementData
    ed = ElementData()
    ed.init_element_data()
    return ed


class TestElementDataCoreVale:
    """Verify core/vale orbital counts for known elements."""

    def test_carbon_core(self, element_data):
        z = element_data.get_element_z('C')
        assert z == 6
        core = element_data.core_orbitals[z]
        assert core[0] == 1  # 1s core
        assert core[1] == 0  # no p core
        assert core[2] == 0
        assert core[3] == 0

    def test_silicon_core(self, element_data):
        z = element_data.get_element_z('Si')
        assert z == 14
        core = element_data.core_orbitals[z]
        assert core[0] == 2  # 1s, 2s core
        assert core[1] == 1  # 2p core
        assert core[2] == 0
        assert core[3] == 0

    def test_aluminum_core(self, element_data):
        z = element_data.get_element_z('Al')
        assert z == 13
        core = element_data.core_orbitals[z]
        assert core[0] == 2  # 1s, 2s core
        assert core[1] == 1  # 2p core

    def test_iron_core(self, element_data):
        z = element_data.get_element_z('Fe')
        assert z == 26
        core = element_data.core_orbitals[z]
        assert core[0] == 3  # 1s, 2s, 3s core
        assert core[1] == 2  # 2p, 3p core
        assert core[2] == 0  # 3d is valence
        assert core[3] == 0

    def test_carbon_mb_valence(self, element_data):
        z = element_data.get_element_z('C')
        vale_mb = element_data.vale_orbitals[1][z]
        assert vale_mb[0] == 1  # 2s valence
        assert vale_mb[1] == 1  # 2p valence

    def test_silicon_mb_valence(self, element_data):
        z = element_data.get_element_z('Si')
        vale_mb = element_data.vale_orbitals[1][z]
        assert vale_mb[0] == 1  # 3s valence
        assert vale_mb[1] == 1  # 3p valence

    def test_iron_mb_valence(self, element_data):
        z = element_data.get_element_z('Fe')
        vale_mb = element_data.vale_orbitals[1][z]
        assert vale_mb[0] == 1  # 4s
        assert vale_mb[1] == 1  # 4p
        assert vale_mb[2] == 1  # 3d
