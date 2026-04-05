"""Unit tests for StructureControl math and geometry utilities.

These tests exercise pure-computation methods that do not depend on reading
any structure file beyond the element database.  A StructureControl instance
is required because all methods are instance methods, but only the element
database (OLCAO_DATA/elements.dat) needs to be accessible.

Run::

    python src/tests/test_o.py src/tests/test_structure_control_unit.py -v

or select by marker::

    python src/tests/test_o.py -m unit -v
"""

import math
import pytest

pytestmark = pytest.mark.unit


# ===========================================================================
# Vector math
# ===========================================================================

class TestDotProduct:
    """dot_product(a, b) → scalar. Both vectors are 0-indexed [x, y, z]."""

    def test_orthogonal_x_y(self, fresh_sc):
        assert fresh_sc.dot_product([1, 0, 0], [0, 1, 0]) == pytest.approx(0.0)

    def test_orthogonal_x_z(self, fresh_sc):
        assert fresh_sc.dot_product([1, 0, 0], [0, 0, 1]) == pytest.approx(0.0)

    def test_self_dot_unit(self, fresh_sc):
        assert fresh_sc.dot_product([1, 0, 0], [1, 0, 0]) == pytest.approx(1.0)

    def test_self_dot_general(self, fresh_sc):
        # [1,2,3]·[1,2,3] = 1 + 4 + 9 = 14
        assert fresh_sc.dot_product([1, 2, 3], [1, 2, 3]) == pytest.approx(14.0)

    def test_antiparallel(self, fresh_sc):
        assert fresh_sc.dot_product([1, 0, 0], [-1, 0, 0]) == pytest.approx(-1.0)

    def test_general_vectors(self, fresh_sc):
        # [2, 3, 1]·[1, -1, 4] = 2 - 3 + 4 = 3
        assert fresh_sc.dot_product([2, 3, 1], [1, -1, 4]) == pytest.approx(3.0)

    def test_commutative(self, fresh_sc):
        a = [1.5, -2.0, 0.7]
        b = [3.1,  0.5, 4.2]
        assert fresh_sc.dot_product(a, b) == pytest.approx(fresh_sc.dot_product(b, a))

    def test_float_inputs(self, fresh_sc):
        assert fresh_sc.dot_product([0.5, 0.5, 0.0], [1.0, 1.0, 0.0]) == pytest.approx(1.0)


class TestCrossProduct:
    """cross_product(a, b) → [cx, cy, cz]. Both vectors are 0-indexed."""

    def test_x_cross_y_is_z(self, fresh_sc):
        assert fresh_sc.cross_product([1, 0, 0], [0, 1, 0]) == pytest.approx([0, 0, 1])

    def test_y_cross_x_is_neg_z(self, fresh_sc):
        assert fresh_sc.cross_product([0, 1, 0], [1, 0, 0]) == pytest.approx([0, 0, -1])

    def test_x_cross_z_is_neg_y(self, fresh_sc):
        assert fresh_sc.cross_product([1, 0, 0], [0, 0, 1]) == pytest.approx([0, -1, 0])

    def test_z_cross_x_is_y(self, fresh_sc):
        assert fresh_sc.cross_product([0, 0, 1], [1, 0, 0]) == pytest.approx([0, 1, 0])

    def test_self_cross_is_zero(self, fresh_sc):
        result = fresh_sc.cross_product([1, 2, 3], [1, 2, 3])
        assert result == pytest.approx([0, 0, 0], abs=1e-14)

    def test_anticommutative(self, fresh_sc):
        a = [1.0, 2.0, 3.0]
        b = [4.0, 5.0, 6.0]
        axb = fresh_sc.cross_product(a, b)
        bxa = fresh_sc.cross_product(b, a)
        for i in range(3):
            assert axb[i] == pytest.approx(-bxa[i])

    def test_general_result(self, fresh_sc):
        # [2,3,4] × [5,6,7]:
        # cx = 3*7 - 4*6 = 21 - 24 = -3
        # cy = 4*5 - 2*7 = 20 - 14 =  6
        # cz = 2*6 - 3*5 = 12 - 15 = -3
        assert fresh_sc.cross_product([2, 3, 4], [5, 6, 7]) == pytest.approx([-3, 6, -3])

    def test_result_orthogonal_to_inputs(self, fresh_sc):
        a = [1.0, 2.0, 3.0]
        b = [4.0, 5.0, 6.0]
        c = fresh_sc.cross_product(a, b)
        assert fresh_sc.dot_product(c, a) == pytest.approx(0.0, abs=1e-12)
        assert fresh_sc.dot_product(c, b) == pytest.approx(0.0, abs=1e-12)


class TestCrossProductMag:
    """cross_product_mag(a, b) → |a × b|."""

    def test_unit_vectors_perp(self, fresh_sc):
        assert fresh_sc.cross_product_mag([1, 0, 0], [0, 1, 0]) == pytest.approx(1.0)

    def test_parallel_vectors(self, fresh_sc):
        assert fresh_sc.cross_product_mag([1, 0, 0], [2, 0, 0]) == pytest.approx(0.0, abs=1e-14)

    def test_scaled_vectors(self, fresh_sc):
        # |[3,0,0] × [0,4,0]| = |[0,0,12]| = 12
        assert fresh_sc.cross_product_mag([3, 0, 0], [0, 4, 0]) == pytest.approx(12.0)

    def test_relation_to_cross_product(self, fresh_sc):
        a = [1.0, 2.0, 3.0]
        b = [4.0, 5.0, 6.0]
        c = fresh_sc.cross_product(a, b)
        expected_mag = math.sqrt(sum(x**2 for x in c))
        assert fresh_sc.cross_product_mag(a, b) == pytest.approx(expected_mag)


class TestGetPlaneNormal:
    """get_plane_normal(p1, p2, p3) — takes three 0-indexed [x,y,z] points."""

    def test_xy_plane_normal_is_z(self, fresh_sc):
        n = fresh_sc.get_plane_normal([0, 0, 0], [1, 0, 0], [0, 1, 0])
        # Normal must be parallel to [0,0,1]: first two components are zero.
        assert n[0] == pytest.approx(0.0, abs=1e-10)
        assert n[1] == pytest.approx(0.0, abs=1e-10)
        assert abs(n[2]) > 1e-10  # non-zero z component

    def test_xz_plane_normal_is_y(self, fresh_sc):
        n = fresh_sc.get_plane_normal([0, 0, 0], [1, 0, 0], [0, 0, 1])
        assert n[0] == pytest.approx(0.0, abs=1e-10)
        assert abs(n[1]) > 1e-10
        assert n[2] == pytest.approx(0.0, abs=1e-10)

    def test_normal_orthogonal_to_plane_vectors(self, fresh_sc):
        p1 = [0.0, 0.0, 0.0]
        p2 = [2.0, 1.0, 0.0]
        p3 = [0.0, 3.0, 1.0]
        n = fresh_sc.get_plane_normal(p1, p2, p3)
        d1 = [p2[i] - p1[i] for i in range(3)]
        d2 = [p3[i] - p1[i] for i in range(3)]
        assert fresh_sc.dot_product(n, d1) == pytest.approx(0.0, abs=1e-10)
        assert fresh_sc.dot_product(n, d2) == pytest.approx(0.0, abs=1e-10)


# ===========================================================================
# Lattice geometry
# ===========================================================================

class TestSetLatticeFromMagAngle:
    """set_lattice_from_mag_angle populates real_lattice and mag/angle arrays."""

    def test_cubic_mags(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle([5.0, 5.0, 5.0], [90.0, 90.0, 90.0])
        fresh_sc.abc_alpha_beta_gamma()
        for axis in range(1, 4):
            assert fresh_sc.mag[axis] == pytest.approx(5.0, rel=1e-10)

    def test_cubic_angles(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle([5.0, 5.0, 5.0], [90.0, 90.0, 90.0])
        fresh_sc.abc_alpha_beta_gamma()
        for axis in range(1, 4):
            assert fresh_sc.angle_deg[axis] == pytest.approx(90.0, abs=1e-10)

    def test_hexagonal_gamma_120(self, fresh_sc):
        # BeO: a=b=2.694, c=4.393, alpha=beta=90, gamma=120
        fresh_sc.set_lattice_from_mag_angle(
            [2.694, 2.694, 4.393], [90.0, 90.0, 120.0])
        fresh_sc.abc_alpha_beta_gamma()
        assert fresh_sc.angle_deg[3] == pytest.approx(120.0, abs=1e-8)

    def test_hexagonal_a_equals_b(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle(
            [2.694, 2.694, 4.393], [90.0, 90.0, 120.0])
        fresh_sc.abc_alpha_beta_gamma()
        assert fresh_sc.mag[1] == pytest.approx(fresh_sc.mag[2], rel=1e-10)

    def test_orthorhombic_distinct_mags(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle([3.0, 4.0, 5.0], [90.0, 90.0, 90.0])
        fresh_sc.abc_alpha_beta_gamma()
        assert fresh_sc.mag[1] == pytest.approx(3.0, rel=1e-10)
        assert fresh_sc.mag[2] == pytest.approx(4.0, rel=1e-10)
        assert fresh_sc.mag[3] == pytest.approx(5.0, rel=1e-10)


class TestAbcAlphaBetaGamma:
    """abc_alpha_beta_gamma computes mag and angle arrays from real_lattice."""

    def test_cubic_angles_are_90(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle([4.0, 4.0, 4.0], [90.0, 90.0, 90.0])
        fresh_sc.abc_alpha_beta_gamma()
        for axis in range(1, 4):
            assert fresh_sc.angle_deg[axis] == pytest.approx(90.0, abs=1e-8)

    def test_fcc_primitive_angles_are_60(self, fresh_sc):
        # FCC primitive cell: a=b=c, alpha=beta=gamma=60°.
        # Construct from the FCC primitive vectors for a=4.0 Å conventional:
        # a_prim = a/sqrt(2) = 2.8284
        a_prim = 4.0 / math.sqrt(2)
        fresh_sc.set_lattice_from_mag_angle(
            [a_prim, a_prim, a_prim], [60.0, 60.0, 60.0])
        fresh_sc.abc_alpha_beta_gamma()
        for axis in range(1, 4):
            assert fresh_sc.angle_deg[axis] == pytest.approx(60.0, abs=1e-6)

    def test_mags_stored_in_radians_too(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle([3.0, 4.0, 5.0], [90.0, 90.0, 90.0])
        fresh_sc.abc_alpha_beta_gamma()
        for axis in range(1, 4):
            assert fresh_sc.angle[axis] == pytest.approx(
                fresh_sc.angle_deg[axis] * math.pi / 180.0, rel=1e-10)


class TestMakeInvOrRecipLattice:
    """make_inv_or_recip_lattice: the real-space inverse satisfies M·M⁻¹ = I."""

    def test_cubic_inverse_is_diagonal(self, fresh_sc):
        # For a cubic cell with a=5, M = 5*I, M_inv = (1/5)*I.
        a = 5.0
        fresh_sc.set_lattice_from_mag_angle([a, a, a], [90.0, 90.0, 90.0])
        # real_lattice[abc][xyz], real_lattice_inv[xyz][abc]
        for xyz in range(1, 4):
            for abc in range(1, 4):
                expected = (1.0 / a) if xyz == abc else 0.0
                assert fresh_sc.real_lattice_inv[xyz][abc] == pytest.approx(
                    expected, abs=1e-10)

    def test_inverse_satisfies_identity(self, fresh_sc):
        # M * M_inv = I (in the 3×3 matrix sense).
        # real_lattice[i][k] * real_lattice_inv[k][j] = delta_ij
        fresh_sc.set_lattice_from_mag_angle([3.0, 4.0, 5.0], [90.0, 90.0, 90.0])
        for i in range(1, 4):
            for j in range(1, 4):
                dot = sum(
                    fresh_sc.real_lattice[i][k] * fresh_sc.real_lattice_inv[k][j]
                    for k in range(1, 4)
                )
                expected = 1.0 if i == j else 0.0
                assert dot == pytest.approx(expected, abs=1e-10)

    def test_hexagonal_inverse_satisfies_identity(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle(
            [2.694, 2.694, 4.393], [90.0, 90.0, 120.0])
        for i in range(1, 4):
            for j in range(1, 4):
                dot = sum(
                    fresh_sc.real_lattice[i][k] * fresh_sc.real_lattice_inv[k][j]
                    for k in range(1, 4)
                )
                expected = 1.0 if i == j else 0.0
                assert dot == pytest.approx(expected, abs=1e-10)


# ===========================================================================
# Coordinate transforms
# ===========================================================================

class TestFractABC2DirectXYZ:
    """fract_abc2direct_xyz(fract) — 0-indexed [a,b,c] → 0-indexed [x,y,z]."""

    def test_origin_maps_to_origin(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle([5.0, 5.0, 5.0], [90.0, 90.0, 90.0])
        xyz = fresh_sc.fract_abc2direct_xyz([0.0, 0.0, 0.0])
        assert xyz == pytest.approx([0.0, 0.0, 0.0], abs=1e-10)

    def test_cubic_half_half_half(self, fresh_sc):
        # In a cubic 5 Å cell, fract (0.5,0.5,0.5) → cart (2.5, 2.5, 2.5)
        fresh_sc.set_lattice_from_mag_angle([5.0, 5.0, 5.0], [90.0, 90.0, 90.0])
        xyz = fresh_sc.fract_abc2direct_xyz([0.5, 0.5, 0.5])
        assert xyz == pytest.approx([2.5, 2.5, 2.5], rel=1e-8)

    def test_cubic_corner(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle([3.0, 4.0, 5.0], [90.0, 90.0, 90.0])
        xyz = fresh_sc.fract_abc2direct_xyz([1.0, 1.0, 1.0])
        assert xyz == pytest.approx([3.0, 4.0, 5.0], rel=1e-8)


class TestDirectXYZ2FractABC:
    """direct_xyz2fract_abc(xyz) — 0-indexed [x,y,z] → fractional [a,b,c]."""

    def test_origin_maps_to_origin(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle([5.0, 5.0, 5.0], [90.0, 90.0, 90.0])
        frac = fresh_sc.direct_xyz2fract_abc([0.0, 0.0, 0.0])
        assert frac == pytest.approx([0.0, 0.0, 0.0], abs=1e-10)

    def test_cubic_corner(self, fresh_sc):
        fresh_sc.set_lattice_from_mag_angle([3.0, 4.0, 5.0], [90.0, 90.0, 90.0])
        frac = fresh_sc.direct_xyz2fract_abc([3.0, 4.0, 5.0])
        assert frac == pytest.approx([1.0, 1.0, 1.0], rel=1e-8)

    def test_roundtrip_cubic(self, fresh_sc):
        """xyz → fract → xyz should recover the original coordinates."""
        fresh_sc.set_lattice_from_mag_angle([4.0, 5.0, 6.0], [90.0, 90.0, 90.0])
        xyz_orig = [1.2, 3.4, 5.1]
        frac = fresh_sc.direct_xyz2fract_abc(xyz_orig)
        xyz_back = fresh_sc.fract_abc2direct_xyz(frac)
        assert xyz_back == pytest.approx(xyz_orig, rel=1e-8)

    def test_roundtrip_hexagonal(self, fresh_sc):
        """Roundtrip should work for non-orthogonal cells too."""
        fresh_sc.set_lattice_from_mag_angle(
            [2.694, 2.694, 4.393], [90.0, 90.0, 120.0])
        xyz_orig = [0.5, 1.2, 2.1]
        frac = fresh_sc.direct_xyz2fract_abc(xyz_orig)
        xyz_back = fresh_sc.fract_abc2direct_xyz(frac)
        assert xyz_back == pytest.approx(xyz_orig, rel=1e-7)
