"""Integration tests for StructureControl file I/O.

Tests in this module exercise reading structure files and verifying that the
parsed data (atom counts, element identities, lattice parameters) matches
known-good values.  Roundtrip tests write a structure and read it back to
confirm that no information is lost.

Run::

    python src/tests/test_o.py src/tests/test_structure_control_io.py -v

or select by marker::

    python src/tests/test_o.py -m integration -v
"""

import pytest

pytestmark = pytest.mark.integration


# ===========================================================================
# BN cubic (zincblende, SG 216, 2×2×2 prim)
# ===========================================================================
# Verified values from a live run:
#   16 atoms, 2 elements (b, n)
#   rhombohedral cell: a=b=c=5.1137 Å, alpha=beta=gamma=60°

class TestReadBNCubic:

    def test_atom_count(self, sc_bn_cubic):
        assert sc_bn_cubic.num_atoms == 16

    def test_element_count(self, sc_bn_cubic):
        assert sc_bn_cubic.num_elements == 2

    def test_element_names(self, sc_bn_cubic):
        names = set(
            sc_bn_cubic.atom_element_name[i]
            for i in range(1, sc_bn_cubic.num_atoms + 1)
        )
        assert names == {'b', 'n'}

    def test_equal_numbers_of_b_and_n(self, sc_bn_cubic):
        counts = {}
        for i in range(1, sc_bn_cubic.num_atoms + 1):
            name = sc_bn_cubic.atom_element_name[i]
            counts[name] = counts.get(name, 0) + 1
        assert counts['b'] == counts['n'] == 8

    def test_lattice_magnitudes(self, sc_bn_cubic):
        sc_bn_cubic.abc_alpha_beta_gamma()
        for axis in range(1, 4):
            assert sc_bn_cubic.mag[axis] == pytest.approx(5.1137, rel=1e-4)

    def test_lattice_angles_60_deg(self, sc_bn_cubic):
        sc_bn_cubic.abc_alpha_beta_gamma()
        for axis in range(1, 4):
            assert sc_bn_cubic.angle_deg[axis] == pytest.approx(60.0, abs=1e-3)

    def test_all_fract_abc_in_range(self, sc_bn_cubic):
        """All fractional coordinates must be in [0, 1)."""
        for atom in range(1, sc_bn_cubic.num_atoms + 1):
            for coord in range(1, 4):
                val = sc_bn_cubic.fract_abc[atom][coord]
                assert 0.0 <= val < 1.0 + 1e-8, (
                    f'fract_abc[{atom}][{coord}] = {val} out of [0, 1)')


# ===========================================================================
# Si diamond (SG 227_a, 1×1×1 full)
# ===========================================================================
# Verified: 8 atoms, 1 element (si), cubic a=5.4307 Å

class TestReadSiDiamond:

    def test_atom_count(self, sc_si_diamond):
        assert sc_si_diamond.num_atoms == 8

    def test_element_count(self, sc_si_diamond):
        assert sc_si_diamond.num_elements == 1

    def test_all_atoms_are_si(self, sc_si_diamond):
        for i in range(1, sc_si_diamond.num_atoms + 1):
            assert sc_si_diamond.atom_element_name[i] == 'si'

    def test_cubic_lattice(self, sc_si_diamond):
        sc_si_diamond.abc_alpha_beta_gamma()
        a = 5.43070
        for axis in range(1, 4):
            assert sc_si_diamond.mag[axis] == pytest.approx(a, rel=1e-4)

    def test_angles_are_90_deg(self, sc_si_diamond):
        sc_si_diamond.abc_alpha_beta_gamma()
        for axis in range(1, 4):
            assert sc_si_diamond.angle_deg[axis] == pytest.approx(90.0, abs=1e-4)


# ===========================================================================
# BeO hexagonal (SG 186, 1×1×1 full)
# ===========================================================================
# Verified: 4 atoms, 2 elements (be, o), hexagonal a=b=2.694, c=4.393, γ=120°

class TestReadBeOHex:

    def test_atom_count(self, sc_beo_hex):
        assert sc_beo_hex.num_atoms == 4

    def test_element_count(self, sc_beo_hex):
        assert sc_beo_hex.num_elements == 2

    def test_element_names(self, sc_beo_hex):
        names = set(
            sc_beo_hex.atom_element_name[i]
            for i in range(1, sc_beo_hex.num_atoms + 1)
        )
        assert names == {'be', 'o'}

    def test_a_equals_b(self, sc_beo_hex):
        sc_beo_hex.abc_alpha_beta_gamma()
        assert sc_beo_hex.mag[1] == pytest.approx(sc_beo_hex.mag[2], rel=1e-6)

    def test_a_magnitude(self, sc_beo_hex):
        sc_beo_hex.abc_alpha_beta_gamma()
        assert sc_beo_hex.mag[1] == pytest.approx(2.6940, rel=1e-4)

    def test_c_magnitude(self, sc_beo_hex):
        sc_beo_hex.abc_alpha_beta_gamma()
        assert sc_beo_hex.mag[3] == pytest.approx(4.3930, rel=1e-4)

    def test_gamma_120_deg(self, sc_beo_hex):
        sc_beo_hex.abc_alpha_beta_gamma()
        assert sc_beo_hex.angle_deg[3] == pytest.approx(120.0, abs=1e-3)

    def test_alpha_beta_90_deg(self, sc_beo_hex):
        sc_beo_hex.abc_alpha_beta_gamma()
        assert sc_beo_hex.angle_deg[1] == pytest.approx(90.0, abs=1e-3)
        assert sc_beo_hex.angle_deg[2] == pytest.approx(90.0, abs=1e-3)


# ===========================================================================
# C₂ molecule (P1 box)
# ===========================================================================
# Verified: 2 atoms, 1 element (c), cubic 5 Å box

class TestReadC2Molecule:

    def test_atom_count(self, sc_c2_molecule):
        assert sc_c2_molecule.num_atoms == 2

    def test_element_count(self, sc_c2_molecule):
        assert sc_c2_molecule.num_elements == 1

    def test_both_atoms_carbon(self, sc_c2_molecule):
        for i in range(1, 3):
            assert sc_c2_molecule.atom_element_name[i] == 'c'

    def test_cubic_box_5_angstrom(self, sc_c2_molecule):
        sc_c2_molecule.abc_alpha_beta_gamma()
        for axis in range(1, 4):
            assert sc_c2_molecule.mag[axis] == pytest.approx(5.0, rel=1e-5)

    def test_bond_length_1_5_angstrom(self, sc_c2_molecule):
        """The two C atoms are separated by 1.5 Å (= 3.25 - 1.75)."""
        # Cartesian positions from the skl file: c at (1.75, 2.5, 2.5)
        # and (3.25, 2.5, 2.5).  After reading, they live in fract_abc.
        sc_c2_molecule.abc_alpha_beta_gamma()
        import math
        x1 = sc_c2_molecule.direct_xyz[1][1]
        x2 = sc_c2_molecule.direct_xyz[2][1]
        y1 = sc_c2_molecule.direct_xyz[1][2]
        y2 = sc_c2_molecule.direct_xyz[2][2]
        z1 = sc_c2_molecule.direct_xyz[1][3]
        z2 = sc_c2_molecule.direct_xyz[2][3]
        dist = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
        assert dist == pytest.approx(1.5, rel=1e-4)


# ===========================================================================
# FeS₂ pyrite (SG 205, 2×2×1 prim)
# ===========================================================================
# Verified: 48 atoms, 2 elements (fe, s), orthorhombic-like cell

class TestReadFeS2:

    def test_atom_count(self, sc_fes2):
        assert sc_fes2.num_atoms == 48

    def test_element_count(self, sc_fes2):
        assert sc_fes2.num_elements == 2

    def test_element_names(self, sc_fes2):
        names = set(
            sc_fes2.atom_element_name[i]
            for i in range(1, sc_fes2.num_atoms + 1)
        )
        assert names == {'fe', 's'}

    def test_two_s_per_fe(self, sc_fes2):
        """FeS₂ stoichiometry: twice as many S as Fe."""
        counts = {}
        for i in range(1, sc_fes2.num_atoms + 1):
            name = sc_fes2.atom_element_name[i]
            counts[name] = counts.get(name, 0) + 1
        assert counts['s'] == 2 * counts['fe']


# ===========================================================================
# Roundtrip tests: read → print_olcao → read back
# ===========================================================================

class TestRoundtrip:
    """print_olcao followed by read_input_file must preserve key quantities."""

    def test_bn_atom_count(self, make_sc, sc_from_file, tmp_path):
        sc1 = make_sc('bn_cubic.skl')
        out = str(tmp_path / 'bn_rt.skl')
        sc1.print_olcao(filename=out, style='frac')
        sc2 = sc_from_file(out)
        assert sc2.num_atoms == sc1.num_atoms

    def test_bn_element_names(self, make_sc, sc_from_file, tmp_path):
        sc1 = make_sc('bn_cubic.skl')
        out = str(tmp_path / 'bn_rt.skl')
        sc1.print_olcao(filename=out, style='frac')
        sc2 = sc_from_file(out)
        names1 = set(sc1.atom_element_name[i] for i in range(1, sc1.num_atoms + 1))
        names2 = set(sc2.atom_element_name[i] for i in range(1, sc2.num_atoms + 1))
        assert names1 == names2

    def test_bn_lattice_magnitudes(self, make_sc, sc_from_file, tmp_path):
        sc1 = make_sc('bn_cubic.skl')
        sc1.abc_alpha_beta_gamma()
        mags_before = [sc1.mag[i] for i in range(1, 4)]

        out = str(tmp_path / 'bn_rt.skl')
        sc1.print_olcao(filename=out, style='frac')
        sc2 = sc_from_file(out)
        sc2.abc_alpha_beta_gamma()

        for i in range(3):
            assert sc2.mag[i + 1] == pytest.approx(mags_before[i], rel=1e-6)

    def test_c2_fractional_positions(self, make_sc, sc_from_file, tmp_path):
        """Fractional atom positions must survive a write/read cycle."""
        sc1 = make_sc('c2_molecule.skl')
        out = str(tmp_path / 'c2_rt.skl')
        sc1.print_olcao(filename=out, style='frac')
        sc2 = sc_from_file(out)

        for atom in range(1, sc1.num_atoms + 1):
            for coord in range(1, 4):
                assert sc2.fract_abc[atom][coord] == pytest.approx(
                    sc1.fract_abc[atom][coord], abs=1e-7
                )

    def test_beo_atom_count(self, make_sc, sc_from_file, tmp_path):
        sc1 = make_sc('beo_hexagonal.skl')
        out = str(tmp_path / 'beo_rt.skl')
        sc1.print_olcao(filename=out, style='frac')
        sc2 = sc_from_file(out)
        assert sc2.num_atoms == sc1.num_atoms
