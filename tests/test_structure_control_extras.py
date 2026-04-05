"""Additional tests for StructureControl methods not covered elsewhere.

Covers:
  - Vector math: normalized_cross_product, get_vector_angle, spherical_angles
  - Utilities: stable_sort, gaussian_broaden
  - Rotation: define_rot_matrix, rotate_one_point, rotate_arb_axis
  - Structure cutting: cut_block, cut_sphere
  - Write methods: print_vasp, print_pdb, print_cif, print_lmp, print_isaacs

Run::

    python src/tests/test_o.py src/tests/test_structure_control_extras.py -v

or by marker::

    python src/tests/test_o.py -m unit -v
    python src/tests/test_o.py -m integration -v
"""

import math
import os
import pytest


# ===========================================================================
# Vector math — unit tests
# ===========================================================================

class TestNormalizedCrossProduct:
    """normalized_cross_product(a, b) → unit vector parallel to a × b."""

    pytestmark = pytest.mark.unit

    def test_x_cross_y_is_z(self, fresh_sc):
        n = fresh_sc.normalized_cross_product([1, 0, 0], [0, 1, 0])
        assert n == pytest.approx([0, 0, 1], abs=1e-12)

    def test_result_is_unit_vector(self, fresh_sc):
        a = [1.0, 2.0, 3.0]
        b = [4.0, 5.0, 6.0]
        n = fresh_sc.normalized_cross_product(a, b)
        mag = math.sqrt(n[0]**2 + n[1]**2 + n[2]**2)
        assert mag == pytest.approx(1.0, rel=1e-10)

    def test_parallel_to_cross_product(self, fresh_sc):
        a = [2.0, 1.0, 0.0]
        b = [0.0, 3.0, 1.0]
        n = fresh_sc.normalized_cross_product(a, b)
        c = fresh_sc.cross_product(a, b)
        c_mag = math.sqrt(sum(x**2 for x in c))
        c_unit = [x / c_mag for x in c]
        assert n == pytest.approx(c_unit, rel=1e-10)

    def test_orthogonal_to_both_inputs(self, fresh_sc):
        a = [1.0, 0.0, 0.0]
        b = [0.0, 1.0, 1.0]
        n = fresh_sc.normalized_cross_product(a, b)
        assert fresh_sc.dot_product(n, a) == pytest.approx(0.0, abs=1e-12)
        assert fresh_sc.dot_product(n, b) == pytest.approx(0.0, abs=1e-12)

    def test_anticommutative_direction(self, fresh_sc):
        a = [1.0, 2.0, 3.0]
        b = [4.0, 5.0, 6.0]
        n_ab = fresh_sc.normalized_cross_product(a, b)
        n_ba = fresh_sc.normalized_cross_product(b, a)
        for i in range(3):
            assert n_ab[i] == pytest.approx(-n_ba[i], rel=1e-10)


class TestGetVectorAngle:
    """get_vector_angle(a, b) → angle in radians."""

    pytestmark = pytest.mark.unit

    def test_parallel_is_zero(self, fresh_sc):
        assert fresh_sc.get_vector_angle([1, 0, 0], [1, 0, 0]) == pytest.approx(0.0, abs=1e-12)

    def test_antiparallel_is_pi(self, fresh_sc):
        assert fresh_sc.get_vector_angle([1, 0, 0], [-1, 0, 0]) == pytest.approx(math.pi, rel=1e-10)

    def test_perpendicular_x_y_is_half_pi(self, fresh_sc):
        assert fresh_sc.get_vector_angle([1, 0, 0], [0, 1, 0]) == pytest.approx(math.pi / 2.0, rel=1e-10)

    def test_perpendicular_x_z_is_half_pi(self, fresh_sc):
        assert fresh_sc.get_vector_angle([1, 0, 0], [0, 0, 1]) == pytest.approx(math.pi / 2.0, rel=1e-10)

    def test_45_degrees(self, fresh_sc):
        # [1,1,0] and [0,1,0]: dot = 1, |a|=sqrt(2), |b|=1 → angle = pi/4
        angle = fresh_sc.get_vector_angle([1, 1, 0], [0, 1, 0])
        assert angle == pytest.approx(math.pi / 4.0, rel=1e-8)

    def test_scale_invariant(self, fresh_sc):
        """Scaling a vector should not change the angle."""
        angle1 = fresh_sc.get_vector_angle([1, 0, 0], [0, 1, 0])
        angle2 = fresh_sc.get_vector_angle([5, 0, 0], [0, 3, 0])
        assert angle1 == pytest.approx(angle2, rel=1e-10)

    def test_commutative(self, fresh_sc):
        a = [1.0, 2.0, 3.0]
        b = [4.0, 5.0, 6.0]
        assert fresh_sc.get_vector_angle(a, b) == pytest.approx(
            fresh_sc.get_vector_angle(b, a), rel=1e-10)

    def test_result_in_range_0_pi(self, fresh_sc):
        angle = fresh_sc.get_vector_angle([1, 2, 3], [4, -1, 2])
        assert 0.0 <= angle <= math.pi


class TestSphericalAngles:
    """spherical_angles(vector) → (theta, phi) in radians."""

    pytestmark = pytest.mark.unit

    def test_z_axis_theta_zero(self, fresh_sc):
        theta, phi = fresh_sc.spherical_angles([0, 0, 1])
        assert theta == pytest.approx(0.0, abs=1e-12)

    def test_z_axis_phi_zero(self, fresh_sc):
        theta, phi = fresh_sc.spherical_angles([0, 0, 1])
        assert phi == pytest.approx(0.0, abs=1e-12)

    def test_neg_z_axis_theta_pi(self, fresh_sc):
        theta, phi = fresh_sc.spherical_angles([0, 0, -1])
        assert theta == pytest.approx(math.pi, abs=1e-12)

    def test_x_axis_theta_half_pi(self, fresh_sc):
        theta, phi = fresh_sc.spherical_angles([1, 0, 0])
        assert theta == pytest.approx(math.pi / 2.0, abs=1e-12)

    def test_x_axis_phi_zero(self, fresh_sc):
        theta, phi = fresh_sc.spherical_angles([1, 0, 0])
        assert phi == pytest.approx(0.0, abs=1e-12)

    def test_y_axis_theta_half_pi(self, fresh_sc):
        theta, phi = fresh_sc.spherical_angles([0, 1, 0])
        assert theta == pytest.approx(math.pi / 2.0, abs=1e-12)

    def test_y_axis_phi_half_pi(self, fresh_sc):
        theta, phi = fresh_sc.spherical_angles([0, 1, 0])
        assert phi == pytest.approx(math.pi / 2.0, abs=1e-12)

    def test_neg_x_axis_phi_pi(self, fresh_sc):
        theta, phi = fresh_sc.spherical_angles([-1, 0, 0])
        assert abs(phi) == pytest.approx(math.pi, abs=1e-12)

    def test_scaling_does_not_change_angles(self, fresh_sc):
        """Multiplying the vector by a positive scalar should not change angles."""
        t1, p1 = fresh_sc.spherical_angles([1, 1, 1])
        t2, p2 = fresh_sc.spherical_angles([3, 3, 3])
        assert t1 == pytest.approx(t2, rel=1e-10)
        assert p1 == pytest.approx(p2, rel=1e-10)


# ===========================================================================
# Utility functions — unit tests
# ===========================================================================

class TestStableSort:
    """stable_sort(values, keys) → indices sorted by keys (stable on ties)."""

    pytestmark = pytest.mark.unit

    def test_already_sorted(self, fresh_sc):
        values = ['a', 'b', 'c']
        keys   = [1, 2, 3]
        result = fresh_sc.stable_sort(values, keys)
        assert result == [0, 1, 2]

    def test_reverse_order(self, fresh_sc):
        values = ['a', 'b', 'c']
        keys   = [3, 2, 1]
        result = fresh_sc.stable_sort(values, keys)
        assert result == [2, 1, 0]

    def test_general_permutation(self, fresh_sc):
        # keys: [3, 1, 4, 2] → sorted order of indices: 1, 3, 0, 2
        values = ['a', 'b', 'c', 'd']
        keys   = [3, 1, 4, 2]
        result = fresh_sc.stable_sort(values, keys)
        assert result == [1, 3, 0, 2]

    def test_stable_on_ties(self, fresh_sc):
        """Equal keys must preserve relative order of the tied items."""
        values = ['a', 'b', 'c']
        keys   = [1, 1, 2]
        result = fresh_sc.stable_sort(values, keys)
        # Both 'a' (idx 0) and 'b' (idx 1) have key=1; original order preserved.
        assert result == [0, 1, 2]

    def test_tie_at_end(self, fresh_sc):
        values = ['a', 'b', 'c', 'd']
        keys   = [2, 1, 3, 3]
        result = fresh_sc.stable_sort(values, keys)
        # key=1→idx 1, key=2→idx 0, key=3 tied→ 2 before 3 (stable)
        assert result == [1, 0, 2, 3]

    def test_single_element(self, fresh_sc):
        assert fresh_sc.stable_sort(['x'], [99]) == [0]

    def test_returns_all_indices(self, fresh_sc):
        values = list(range(10))
        keys   = list(reversed(range(10)))
        result = fresh_sc.stable_sort(values, keys)
        assert sorted(result) == list(range(10))


class TestGaussianBroaden:
    """gaussian_broaden(data, sigma) — Gaussian smoothing of a 1-D array."""

    pytestmark = pytest.mark.unit

    def test_all_zeros_stays_zero(self, fresh_sc):
        data = [0.0] * 50
        result = fresh_sc.gaussian_broaden(data, 0.05)
        assert all(v == pytest.approx(0.0) for v in result)

    def test_output_same_length(self, fresh_sc):
        data = [0.0] * 30
        data[15] = 1.0
        result = fresh_sc.gaussian_broaden(data, 0.05)
        assert len(result) == len(data)

    def test_peak_at_spike_location(self, fresh_sc):
        """Result maximum must coincide with the input spike."""
        n = 101
        data = [0.0] * n
        data[50] = 1.0
        result = fresh_sc.gaussian_broaden(data, 0.05)
        assert result[50] == max(result)

    def test_symmetric_around_spike(self, fresh_sc):
        """Broadened result must be symmetric around the spike index."""
        n = 101
        data = [0.0] * n
        data[50] = 1.0
        result = fresh_sc.gaussian_broaden(data, 0.05)
        for k in range(1, 15):
            assert result[50 - k] == pytest.approx(result[50 + k], rel=1e-10)

    def test_area_preservation(self, fresh_sc):
        """Total spectral weight must be conserved after broadening.

        The broadening formula is: result[gv] += data[point] * exp(-term) / (sigma*sqrt(pi)).
        Each data[point] acts as a Dirac delta with weight data[point].  The
        corresponding Gaussian integrates to data[point] over all x.  Therefore:
            sum(data)  ≈  sum(result) × bin_width
        where bin_width = 0.01 Å (implicit in the formula).
        """
        n = 201
        data = [0.0] * n
        data[100] = 1.0
        sigma = 0.05
        result = fresh_sc.gaussian_broaden(data, sigma)
        bin_width = 0.01  # Å (fixed by the formula in the implementation)
        area_in  = sum(data)              # Dirac weight of the spike
        area_out = sum(result) * bin_width  # integral of the broadened Gaussian
        assert area_out == pytest.approx(area_in, rel=0.02)

    def test_wider_sigma_lower_peak(self, fresh_sc):
        """A larger sigma should produce a shorter, wider peak."""
        n = 201
        data = [0.0] * n
        data[100] = 1.0
        peak_narrow = max(fresh_sc.gaussian_broaden(data, 0.02))
        peak_wide   = max(fresh_sc.gaussian_broaden(data, 0.10))
        assert peak_narrow > peak_wide


# ===========================================================================
# Rotation — unit tests (define_rot_matrix / rotate_one_point)
# ===========================================================================

class TestDefineRotMatrix:
    """define_rot_matrix(axis, angle_deg) stores a 3×3 rotation matrix."""

    pytestmark = pytest.mark.unit

    def test_identity_at_zero_degrees(self, fresh_sc):
        fresh_sc.define_rot_matrix([0, 0, 1], 0.0)
        for i in range(1, 4):
            for j in range(1, 4):
                expected = 1.0 if i == j else 0.0
                assert fresh_sc.rot_matrix[i][j] == pytest.approx(expected, abs=1e-12)

    def test_90_deg_about_z_maps_x_to_y(self, fresh_sc):
        """90° about z: (1,0,0) row-maps to (0,1,0)."""
        fresh_sc.define_rot_matrix([0, 0, 1], 90.0)
        # First row of rot_matrix should be approx [0, 1, 0]
        assert fresh_sc.rot_matrix[1][1] == pytest.approx(0.0, abs=1e-12)
        assert fresh_sc.rot_matrix[1][2] == pytest.approx(1.0, abs=1e-12)
        assert fresh_sc.rot_matrix[1][3] == pytest.approx(0.0, abs=1e-12)

    def test_180_deg_about_z_negates_x_and_y(self, fresh_sc):
        fresh_sc.define_rot_matrix([0, 0, 1], 180.0)
        # M should be diag(-1, -1, 1) approximately
        assert fresh_sc.rot_matrix[1][1] == pytest.approx(-1.0, abs=1e-12)
        assert fresh_sc.rot_matrix[2][2] == pytest.approx(-1.0, abs=1e-12)
        assert fresh_sc.rot_matrix[3][3] == pytest.approx(1.0, abs=1e-12)

    def test_matrix_is_orthogonal(self, fresh_sc):
        """R^T R = I for any rotation matrix.  Axis must be a unit vector."""
        s = math.sqrt(2)
        fresh_sc.define_rot_matrix([1/s, 1/s, 0], 73.0)
        M = fresh_sc.rot_matrix
        for i in range(1, 4):
            for j in range(1, 4):
                dot = sum(M[k][i] * M[k][j] for k in range(1, 4))
                expected = 1.0 if i == j else 0.0
                assert dot == pytest.approx(expected, abs=1e-12)

    def test_matrix_has_determinant_plus_one(self, fresh_sc):
        """det(R) = +1 for a proper rotation."""
        fresh_sc.define_rot_matrix([1, 0, 0], 45.0)
        M = fresh_sc.rot_matrix
        det = (M[1][1] * (M[2][2]*M[3][3] - M[2][3]*M[3][2])
             - M[1][2] * (M[2][1]*M[3][3] - M[2][3]*M[3][1])
             + M[1][3] * (M[2][1]*M[3][2] - M[2][2]*M[3][1]))
        assert det == pytest.approx(1.0, abs=1e-12)


class TestRotateOnePoint:
    """rotate_one_point(point) applies rot_matrix to a single point."""

    pytestmark = pytest.mark.unit

    def test_identity_leaves_point_unchanged(self, fresh_sc):
        fresh_sc.define_rot_matrix([0, 0, 1], 0.0)
        pt = [3.0, 4.0, 5.0]
        result = fresh_sc.rotate_one_point(pt)
        assert result == pytest.approx(pt, rel=1e-10)

    def test_90_deg_z_maps_x_hat_to_y_hat(self, fresh_sc):
        fresh_sc.define_rot_matrix([0, 0, 1], 90.0)
        result = fresh_sc.rotate_one_point([1.0, 0.0, 0.0])
        assert result == pytest.approx([0.0, 1.0, 0.0], abs=1e-12)

    def test_90_deg_z_maps_y_hat_to_neg_x_hat(self, fresh_sc):
        fresh_sc.define_rot_matrix([0, 0, 1], 90.0)
        result = fresh_sc.rotate_one_point([0.0, 1.0, 0.0])
        assert result == pytest.approx([-1.0, 0.0, 0.0], abs=1e-12)

    def test_90_deg_z_leaves_z_invariant(self, fresh_sc):
        fresh_sc.define_rot_matrix([0, 0, 1], 90.0)
        result = fresh_sc.rotate_one_point([0.0, 0.0, 7.0])
        assert result == pytest.approx([0.0, 0.0, 7.0], abs=1e-12)

    def test_180_deg_x_negates_y_and_z(self, fresh_sc):
        fresh_sc.define_rot_matrix([1, 0, 0], 180.0)
        result = fresh_sc.rotate_one_point([0.0, 2.0, 3.0])
        assert result == pytest.approx([0.0, -2.0, -3.0], abs=1e-12)

    def test_rotation_preserves_vector_length(self, fresh_sc):
        s = math.sqrt(3)
        fresh_sc.define_rot_matrix([1/s, 1/s, 1/s], 57.0)
        pt = [2.0, -1.0, 3.5]
        result = fresh_sc.rotate_one_point(pt)
        mag_before = math.sqrt(sum(x**2 for x in pt))
        mag_after  = math.sqrt(sum(x**2 for x in result))
        assert mag_after == pytest.approx(mag_before, rel=1e-10)

    def test_360_deg_returns_original_point(self, fresh_sc):
        fresh_sc.define_rot_matrix([0, 1, 0], 360.0)
        pt = [1.2, 3.4, 5.6]
        result = fresh_sc.rotate_one_point(pt)
        assert result == pytest.approx(pt, abs=1e-10)


# ===========================================================================
# Rotation — integration tests (rotate_arb_axis acts on a whole structure)
# ===========================================================================

class TestRotateArbAxis:
    """rotate_arb_axis rotates every atom; the structure remains physically valid."""

    pytestmark = pytest.mark.integration

    @pytest.fixture
    def sc_c2(self, make_sc):
        return make_sc('c2_molecule.skl')

    def test_atom_count_unchanged(self, sc_c2):
        n_before = sc_c2.num_atoms
        sc_c2.rotate_arb_axis([0, 0, 1], 45.0)
        assert sc_c2.num_atoms == n_before

    def test_element_names_unchanged(self, sc_c2):
        names_before = [sc_c2.atom_element_name[i] for i in range(1, sc_c2.num_atoms + 1)]
        sc_c2.rotate_arb_axis([1, 0, 0], 90.0)
        names_after  = [sc_c2.atom_element_name[i] for i in range(1, sc_c2.num_atoms + 1)]
        assert names_after == names_before

    def test_360_deg_recovers_positions(self, sc_c2):
        """A full-circle rotation must leave every atom in its original place."""
        xyz_before = [[sc_c2.direct_xyz[a][k] for k in range(1, 4)]
                      for a in range(1, sc_c2.num_atoms + 1)]
        sc_c2.rotate_arb_axis([0, 0, 1], 360.0)
        for i, atom in enumerate(range(1, sc_c2.num_atoms + 1)):
            for k in range(3):
                assert sc_c2.direct_xyz[atom][k + 1] == pytest.approx(
                    xyz_before[i][k], abs=1e-8)

    def test_inter_atom_distance_preserved(self, sc_c2):
        """Bond lengths are invariant under rotation."""
        def dist(sc, a, b):
            return math.sqrt(sum(
                (sc.direct_xyz[a][k] - sc.direct_xyz[b][k])**2
                for k in range(1, 4)))
        d_before = dist(sc_c2, 1, 2)
        sc_c2.rotate_arb_axis([1, 1, 0], 73.0)
        d_after = dist(sc_c2, 1, 2)
        assert d_after == pytest.approx(d_before, rel=1e-8)

    def test_lattice_magnitudes_unchanged_after_rotation(self, make_sc):
        """Rotating atoms does not change the lattice vectors themselves."""
        sc = make_sc('si_diamond.skl')
        sc.abc_alpha_beta_gamma()
        mags_before = [sc.mag[i] for i in range(1, 4)]
        sc.rotate_arb_axis([1, 0, 0], 37.0)
        sc.abc_alpha_beta_gamma()
        for i in range(3):
            assert sc.mag[i + 1] == pytest.approx(mags_before[i], rel=1e-10)


# ===========================================================================
# Structure cutting — integration tests
# ===========================================================================
# C₂ molecule fixture positions (from the skl file, cart style):
#   atom 1: direct_xyz ≈ (1.75, 2.5, 2.5)
#   atom 2: direct_xyz ≈ (3.25, 2.5, 2.5)
#   box: 5 × 5 × 5 Å

class TestCutBlock:
    """cut_block removes atoms inside or outside a rectangular block."""

    pytestmark = pytest.mark.integration

    @pytest.fixture
    def sc_c2(self, make_sc):
        return make_sc('c2_molecule.skl')

    # block_border format: [None, [None, lo, hi], [None, lo, hi], [None, lo, hi]]
    # x range [0, 2.0] captures atom 1 (x≈1.75) but not atom 2 (x≈3.25).

    def test_zone1_keep_inside_reduces_to_one_atom(self, sc_c2):
        """zone=1: keep atoms inside the block → only atom 1 survives."""
        border = [None, [None, 0.0, 2.0], [None, 0.0, 5.0], [None, 0.0, 5.0]]
        sc_c2.cut_block(zone=1, abcxyz_flag=1, block_border=border)
        assert sc_c2.num_atoms == 1

    def test_zone0_cut_inside_reduces_to_one_atom(self, make_sc):
        """zone=0: cut atoms inside the block → only atom 2 survives."""
        sc = make_sc('c2_molecule.skl')
        border = [None, [None, 0.0, 2.0], [None, 0.0, 5.0], [None, 0.0, 5.0]]
        sc.cut_block(zone=0, abcxyz_flag=1, block_border=border)
        assert sc.num_atoms == 1

    def test_zone1_full_box_keeps_all(self, make_sc):
        """A block that contains all atoms should leave the structure unchanged."""
        sc = make_sc('si_diamond.skl')
        sc.abc_alpha_beta_gamma()
        a = sc.mag[1] + 1.0  # slightly larger than the cell
        border = [None, [None, 0.0, a], [None, 0.0, a], [None, 0.0, a]]
        n_before = sc.num_atoms
        sc.cut_block(zone=1, abcxyz_flag=1, block_border=border)
        assert sc.num_atoms == n_before

    def test_zone0_empty_box_keeps_all(self, make_sc):
        """Cutting 'inside' an empty (zero-size) block removes nothing."""
        sc = make_sc('si_diamond.skl')
        n_before = sc.num_atoms
        # A block whose low >= high matches no atom
        border = [None, [None, 10.0, 10.0], [None, 10.0, 10.0], [None, 10.0, 10.0]]
        sc.cut_block(zone=0, abcxyz_flag=1, block_border=border)
        assert sc.num_atoms == n_before

    def test_surviving_atom_element_correct(self, sc_c2):
        """After cutting, the surviving atom must still be carbon."""
        border = [None, [None, 0.0, 2.0], [None, 0.0, 5.0], [None, 0.0, 5.0]]
        sc_c2.cut_block(zone=1, abcxyz_flag=1, block_border=border)
        assert sc_c2.atom_element_name[1] == 'c'


class TestCutSphere:
    """cut_sphere removes atoms inside or outside a sphere."""

    pytestmark = pytest.mark.integration

    @pytest.fixture
    def sc_c2(self, make_sc):
        return make_sc('c2_molecule.skl')

    # C₂: atoms at (1.75,2.5,2.5) and (3.25,2.5,2.5), separation 1.5 Å.
    # Centering on atom 1 (target_atom=1):
    #   distance from atom 1 to itself = 0
    #   distance from atom 1 to atom 2 = 1.5 Å

    def test_zone1_tight_sphere_keeps_one(self, sc_c2):
        """zone=1, radius=1.0 Å centered on atom 1 → atom 2 (1.5 Å away) is cut."""
        sc_c2.cut_sphere(zone=1, abcxyz_flag=1,
                         sphere_rad=1.0, target_atom=1,
                         sphere_loc=[None, 0, 0, 0])
        assert sc_c2.num_atoms == 1

    def test_zone1_wide_sphere_keeps_both(self, make_sc):
        """zone=1, radius=2.0 Å centered on atom 1 → both atoms inside → 2 atoms."""
        sc = make_sc('c2_molecule.skl')
        sc.cut_sphere(zone=1, abcxyz_flag=1,
                      sphere_rad=2.0, target_atom=1,
                      sphere_loc=[None, 0, 0, 0])
        assert sc.num_atoms == 2

    def test_zone0_tight_sphere_keeps_one(self, make_sc):
        """zone=0, radius=1.0 Å centered on atom 1 → atom 1 is cut, atom 2 survives."""
        sc = make_sc('c2_molecule.skl')
        sc.cut_sphere(zone=0, abcxyz_flag=1,
                      sphere_rad=1.0, target_atom=1,
                      sphere_loc=[None, 0, 0, 0])
        assert sc.num_atoms == 1

    def test_atom_count_nonnegative(self, make_sc):
        """Cutting more than everything should leave 0 atoms, not negative."""
        sc = make_sc('c2_molecule.skl')
        sc.cut_sphere(zone=0, abcxyz_flag=1,
                      sphere_rad=100.0, target_atom=0,
                      sphere_loc=[None, 2.5, 2.5, 2.5])
        assert sc.num_atoms >= 0

    def test_explicit_xyz_location(self, make_sc):
        """Pass sphere_loc explicitly (target_atom=0) with xyz coordinates."""
        sc = make_sc('c2_molecule.skl')
        # Center sphere on atom 1's known position; radius 0.5 → only atom 1 inside.
        sc.abc_alpha_beta_gamma()
        x1 = sc.direct_xyz[1][1]
        y1 = sc.direct_xyz[1][2]
        z1 = sc.direct_xyz[1][3]
        sc.cut_sphere(zone=1, abcxyz_flag=1,
                      sphere_rad=0.5, target_atom=0,
                      sphere_loc=[None, x1, y1, z1])
        assert sc.num_atoms == 1


# ===========================================================================
# Write methods — integration tests
# ===========================================================================

class TestPrintVasp:
    """print_vasp writes a VASP POSCAR file with the correct structure."""

    pytestmark = pytest.mark.integration

    @pytest.fixture
    def sc_si(self, make_sc):
        return make_sc('si_diamond.skl')

    def test_file_is_created(self, sc_si, tmp_path):
        out = str(tmp_path / 'POSCAR')
        sc_si.print_vasp(filename=out)
        assert os.path.isfile(out)

    def test_line_2_is_scaling_factor(self, sc_si, tmp_path):
        out = str(tmp_path / 'POSCAR')
        sc_si.print_vasp(filename=out)
        with open(out) as f:
            lines = f.readlines()
        assert lines[1].strip() == '1.000'

    def test_line_7_is_direct(self, sc_si, tmp_path):
        out = str(tmp_path / 'POSCAR')
        sc_si.print_vasp(filename=out)
        with open(out) as f:
            lines = f.readlines()
        assert lines[6].strip() == 'direct'

    def test_correct_number_of_coordinate_lines(self, sc_si, tmp_path):
        out = str(tmp_path / 'POSCAR')
        sc_si.print_vasp(filename=out)
        with open(out) as f:
            lines = [l for l in f.readlines() if l.strip()]
        # 1 title + 1 scale + 3 vectors + 1 counts + 1 keyword + N atom lines
        assert len(lines) == 7 + sc_si.num_atoms

    def test_atom_count_in_line_6(self, sc_si, tmp_path):
        """Line 6 must be the atom-count-per-element string (space-separated ints)."""
        out = str(tmp_path / 'POSCAR')
        sc_si.print_vasp(filename=out)
        with open(out) as f:
            lines = f.readlines()
        counts_line = lines[5].strip()
        total = sum(int(x) for x in counts_line.split())
        assert total == sc_si.num_atoms


class TestPrintPdb:
    """print_pdb writes a PDB file with CRYST1 and ATOM records."""

    pytestmark = pytest.mark.integration

    @pytest.fixture
    def sc_si(self, make_sc):
        return make_sc('si_diamond.skl')

    def test_file_is_created(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.pdb')
        sc_si.print_pdb(filename=out)
        assert os.path.isfile(out)

    def test_cryst1_record_present(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.pdb')
        sc_si.print_pdb(filename=out)
        with open(out) as f:
            content = f.read()
        assert 'CRYST1' in content

    def test_atom_records_count(self, sc_si, tmp_path):
        """Number of ATOM records must equal num_atoms."""
        out = str(tmp_path / 'si.pdb')
        sc_si.print_pdb(filename=out)
        with open(out) as f:
            lines = f.readlines()
        atom_lines = [l for l in lines if l.startswith('ATOM') or l.startswith('HETATM')]
        assert len(atom_lines) == sc_si.num_atoms

    def test_ter_record_at_end(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.pdb')
        sc_si.print_pdb(filename=out)
        with open(out) as f:
            content = f.read()
        assert 'TER' in content


class TestPrintCif:
    """print_cif writes a valid CIF file with cell parameters and atom sites."""

    pytestmark = pytest.mark.integration

    @pytest.fixture
    def sc_si(self, make_sc):
        return make_sc('si_diamond.skl')

    def test_file_is_created(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.cif')
        sc_si.print_cif(filename=out)
        assert os.path.isfile(out)

    def test_cell_length_a_present(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.cif')
        sc_si.print_cif(filename=out)
        with open(out) as f:
            content = f.read()
        assert '_cell_length_a' in content

    def test_atom_site_loop_present(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.cif')
        sc_si.print_cif(filename=out)
        with open(out) as f:
            content = f.read()
        assert '_atom_site' in content

    def test_data_block_header_present(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.cif')
        sc_si.print_cif(filename=out)
        with open(out) as f:
            first_line = f.readline().strip()
        assert first_line.startswith('data_')

    def test_number_of_atom_site_lines(self, sc_si, tmp_path):
        """There must be exactly num_atoms coordinate data lines in the atom loop."""
        out = str(tmp_path / 'si.cif')
        sc_si.print_cif(filename=out)
        with open(out) as f:
            lines = f.readlines()
        # After the _atom_site header lines, data lines start with element symbols
        # or sequential numbers.  Count non-blank, non-underscore, non-loop_ lines
        # that appear after the last '_atom_site' tag.
        last_site_tag = max(i for i, l in enumerate(lines) if '_atom_site' in l)
        data_lines = [l for l in lines[last_site_tag + 1:] if l.strip()]
        assert len(data_lines) == sc_si.num_atoms


class TestPrintLmp:
    """print_lmp writes a LAMMPS data file."""

    pytestmark = pytest.mark.integration

    @pytest.fixture
    def sc_si(self, make_sc):
        return make_sc('si_diamond.skl')

    def test_file_is_created(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.lmp')
        sc_si.print_lmp(filename=out)
        assert os.path.isfile(out)

    def test_atoms_count_line_present(self, sc_si, tmp_path):
        """LAMMPS header must contain 'N atoms'."""
        out = str(tmp_path / 'si.lmp')
        sc_si.print_lmp(filename=out)
        with open(out) as f:
            content = f.read()
        assert f'{sc_si.num_atoms} atoms' in content

    def test_atom_types_line_present(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.lmp')
        sc_si.print_lmp(filename=out)
        with open(out) as f:
            content = f.read()
        assert 'atom types' in content

    def test_atoms_section_present(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.lmp')
        sc_si.print_lmp(filename=out)
        with open(out) as f:
            content = f.read()
        assert 'Atoms' in content

    def test_masses_section_present(self, sc_si, tmp_path):
        out = str(tmp_path / 'si.lmp')
        sc_si.print_lmp(filename=out)
        with open(out) as f:
            content = f.read()
        assert 'Masses' in content


class TestPrintIsaacs:
    """print_isaacs writes a .ipf XML file and a companion .c3d file."""

    pytestmark = pytest.mark.integration

    @pytest.fixture
    def sc_si(self, make_sc):
        return make_sc('si_diamond.skl')

    def test_ipf_file_is_created(self, sc_si, tmp_path):
        stem = str(tmp_path / 'si')
        sc_si.print_isaacs(filename=stem)
        assert os.path.isfile(stem + '.ipf')

    def test_c3d_file_is_created(self, sc_si, tmp_path):
        stem = str(tmp_path / 'si')
        sc_si.print_isaacs(filename=stem)
        assert os.path.isfile(stem + '.c3d')

    def test_ipf_contains_atom_count(self, sc_si, tmp_path):
        stem = str(tmp_path / 'si')
        sc_si.print_isaacs(filename=stem)
        with open(stem + '.ipf') as f:
            content = f.read()
        assert f'<atoms>{sc_si.num_atoms}</atoms>' in content

    def test_ipf_is_xml(self, sc_si, tmp_path):
        stem = str(tmp_path / 'si')
        sc_si.print_isaacs(filename=stem)
        with open(stem + '.ipf') as f:
            first_line = f.readline()
        assert first_line.strip().startswith('<')

    def test_c3d_has_one_line_per_atom(self, sc_si, tmp_path):
        """The .c3d file has one coordinate line per atom (after the count header)."""
        stem = str(tmp_path / 'si')
        sc_si.print_isaacs(filename=stem)
        with open(stem + '.c3d') as f:
            data_lines = [l for l in f if l.strip()]
        # First line is the atom count; remaining lines are atom coordinates.
        assert len(data_lines) - 1 == sc_si.num_atoms
