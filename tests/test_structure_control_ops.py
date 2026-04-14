"""Tests for structural operations: supercell, bonding, coordinate transforms.

These tests verify that multi-step operations produce correct results for
well-understood reference structures.

Run::

    python src/tests/test_o.py src/tests/test_structure_control_ops.py -v

or::

    python src/tests/test_o.py -m integration -v
"""

import math
import pytest

pytestmark = pytest.mark.integration


# ===========================================================================
# Supercell expansion
# ===========================================================================

class TestSupercell:
    """apply_supercell (called internally by read_input_file via read_olcao_skl)
    must produce the correct atom count and lattice scaling."""

    def test_bn_2x2x2_atom_count(self, sc_bn_cubic):
        """BN primitive (2 atoms) + SG 216 → 8 atoms → 2×2×2 supercell = 16."""
        assert sc_bn_cubic.num_atoms == 16

    def test_si_1x1x1_atom_count(self, sc_si_diamond):
        """Si (1 atom) + SG 227_a full → 8 atoms in conventional cell."""
        assert sc_si_diamond.num_atoms == 8

    def test_beo_1x1x1_atom_count(self, sc_beo_hex):
        """BeO (2 atoms) + SG 186 full → 4 atoms."""
        assert sc_beo_hex.num_atoms == 4

    def test_fes2_2x2x1_atom_count(self, sc_fes2):
        """FeS₂ (2 atoms: Fe + S) + SG 205 prim → 12 atoms per cell
        → 2×2×1 supercell = 48 atoms."""
        assert sc_fes2.num_atoms == 48


# ===========================================================================
# Bonding
# ===========================================================================
# Note: create_bonding_list() requires set_limit_dist() to be called first.
# It uses limit_dist (via compute_num_cells) to determine how many periodic
# image cells to include.  A value of 6.0 Å covers all first-shell neighbours
# for the structures tested here.
#
# Important: the make_sc() fixture uses read_input_file() which calls
# map_element_number(), so atomic_z is populated before assign_coval_radii()
# is called inside create_bonding_list().

LIMIT_DIST = 6.0   # Å — sufficient for all first-neighbour bonds below


class TestBondingBNCubic:
    """BN zincblende has exactly 4 B-N nearest neighbours per atom."""

    @pytest.fixture
    def sc_with_bonds(self, make_sc):
        sc = make_sc('bn_cubic.skl')
        sc.set_limit_dist(LIMIT_DIST)
        sc.create_bonding_list()
        return sc

    def test_each_atom_has_4_bonds(self, sc_with_bonds):
        sc = sc_with_bonds
        for i in range(1, sc.num_atoms + 1):
            assert sc.num_bonds[i] == 4, (
                f'atom {i} ({sc.atom_element_name[i]}) has '
                f'{sc.num_bonds[i]} bonds, expected 4'
            )

    def test_bond_lengths_are_equal(self, sc_with_bonds):
        """All B-N bonds must have the same length: a√3/4 ≈ 1.5657 Å."""
        sc = sc_with_bonds
        expected = 3.615900000 * math.sqrt(3) / 4   # ≈ 1.5657 Å
        for atom in range(1, sc.num_atoms + 1):
            for bond in range(1, sc.num_bonds[atom] + 1):
                length = sc.bond_length_ext[atom][bond]
                assert length == pytest.approx(expected, rel=1e-3), (
                    f'bond length for atom {atom} bond {bond}: '
                    f'{length:.6f} Å, expected {expected:.6f} Å'
                )

    def test_b_bonds_to_n_only(self, sc_with_bonds):
        """In zincblende, each B bonds only to N and vice versa."""
        sc = sc_with_bonds
        for atom in range(1, sc.num_atoms + 1):
            my_elem = sc.atom_element_name[atom]
            expected_partner = 'n' if my_elem == 'b' else 'b'
            for bond in range(1, sc.num_bonds[atom] + 1):
                partner = sc.bonded[atom][bond]
                partner_elem = sc.atom_element_name[partner]
                assert partner_elem == expected_partner, (
                    f'atom {atom} ({my_elem}) bond {bond} goes to '
                    f'atom {partner} ({partner_elem}), expected {expected_partner}'
                )


class TestBondingC2Molecule:
    """C₂ dimer: the two atoms should be bonded to each other (d ≈ 1.5 Å)."""

    @pytest.fixture
    def sc_with_bonds(self, make_sc):
        sc = make_sc('c2_molecule.skl')
        sc.set_limit_dist(LIMIT_DIST)
        sc.create_bonding_list()
        return sc

    def test_each_atom_has_at_least_one_bond(self, sc_with_bonds):
        sc = sc_with_bonds
        for i in range(1, sc.num_atoms + 1):
            assert sc.num_bonds[i] >= 1, (
                f'C atom {i} has no bonds'
            )

    def test_bond_length_approx_1_5_angstrom(self, sc_with_bonds):
        sc = sc_with_bonds
        for bond in range(1, sc.num_bonds[1] + 1):
            length = sc.bond_length_ext[1][bond]
            assert length == pytest.approx(1.5, rel=0.01)


# ===========================================================================
# Coordinate transform consistency
# ===========================================================================

class TestCoordConsistency:
    """For every atom in each loaded structure, fract_abc and direct_xyz must
    be self-consistent (i.e. converting one gives the other)."""

    def _check_consistency(self, sc):
        sc.abc_alpha_beta_gamma()
        for atom in range(1, sc.num_atoms + 1):
            # fract_abc[atom] is itself a 1-indexed [None, a, b, c] row,
            # so pass it straight into the 1-indexed converter.
            xyz_from_frac = sc.fract_abc2direct_xyz(sc.fract_abc[atom])
            for axis in range(1, 4):
                assert xyz_from_frac[axis] == pytest.approx(
                    sc.direct_xyz[atom][axis], abs=1e-6), (
                    f'atom {atom} axis {axis}: '
                    f'from_frac={xyz_from_frac[axis]:.8f}, '
                    f'stored={sc.direct_xyz[atom][axis]:.8f}'
                )

    def test_bn_cubic(self, sc_bn_cubic):
        self._check_consistency(sc_bn_cubic)

    def test_si_diamond(self, sc_si_diamond):
        self._check_consistency(sc_si_diamond)

    def test_beo_hexagonal(self, sc_beo_hex):
        self._check_consistency(sc_beo_hex)

    def test_c2_molecule(self, sc_c2_molecule):
        self._check_consistency(sc_c2_molecule)


# ===========================================================================
# Element and species data
# ===========================================================================

class TestElementAndSpecies:
    """element_list, atom_element_id, and species_list must be consistent."""

    def test_bn_element_list_has_two_entries(self, sc_bn_cubic):
        assert sc_bn_cubic.num_elements == 2
        # element_list is 1-indexed; [0] is None.
        assert sc_bn_cubic.element_list[0] is None
        assert len(sc_bn_cubic.element_list) == sc_bn_cubic.num_elements + 1

    def test_bn_all_atom_ids_in_range(self, sc_bn_cubic):
        for i in range(1, sc_bn_cubic.num_atoms + 1):
            eid = sc_bn_cubic.atom_element_id[i]
            assert 1 <= eid <= sc_bn_cubic.num_elements

    def test_si_single_species(self, sc_si_diamond):
        """Si diamond has one element with one species."""
        assert sc_si_diamond.num_elements == 1
        eid = 1
        assert sc_si_diamond.num_species[eid] == 1

    def test_atomic_z_hydrogen_is_1(self, fresh_sc):
        """map_element_number maps 'h' → Z=1."""
        result = fresh_sc.get_element_number('h')
        assert result == 1

    def test_atomic_z_carbon_is_6(self, fresh_sc):
        result = fresh_sc.get_element_number('c')
        assert result == 6

    def test_atomic_z_silicon_is_14(self, fresh_sc):
        result = fresh_sc.get_element_number('si')
        assert result == 14
