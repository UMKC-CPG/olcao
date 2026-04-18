"""test_angle_utils.py -- Unit tests for the geometry-derived
angle clustering and force-constant helpers in angle_utils.py.

Covers PSEUDOCODE 10a (`cluster_angles`) and 10b (`get_angle_k`), the two
pure-computation helpers shared by `condense.create_lammps_files` (PSEUDOCODE
10c) and `make_reactions.py` (PSEUDOCODE 10d).  These tests are pure
unit tests: no fixture files, no element database, no
StructureControl import -- they depend only on the helpers themselves and on a
tiny stub bond-parameter table for the `get_angle_k` formula check.

Why a separate test file (instead of folding into one of the
`test_structure_control_*` modules): `angle_utils.py` is a distinct module in
`src/scripts/` and is not part of the StructureControl class, so it deserves its
own test surface. The conftest.py `SCRIPTS_DIR` insertion lets us import it
directly with `from angle_utils import ...`.
"""

import math

import pytest

from angle_utils import (
    AngleType, cluster_angles, get_angle_k)


# All tests in this file are pure-computation unit tests -- no file I/O, no
# element database, no fortran binaries.
pytestmark = pytest.mark.unit


# ---------------------------------------------------------
# cluster_angles -- PSEUDOCODE 10a
# ---------------------------------------------------------

def test_cluster_single_observation():
    """One observation collapses to one cluster of size 1.

    Verifies the trivial base case: a single observation must produce a single
    AngleType whose theta_0 is the observed theta, whose obs_count is 1, and
    whose base_tag is forwarded from the input verbatim.
    """
    obs = [(6, 8, 1, 109.5, "C O H")]
    types, type_map = cluster_angles(obs, tolerance=5.0)

    assert len(types) == 1
    assert type_map == [1]
    t = types[0]
    assert t.Z1 == 6 and t.Zv == 8 and t.Z2 == 1
    assert t.theta_0 == pytest.approx(109.5)
    assert t.obs_count == 1
    assert t.base_tag == "C O H"


def test_cluster_three_close_observations():
    """Three observations within tolerance merge into one
    cluster; theta_0 is the arithmetic mean of all merged
    observations and obs_count reflects the merge."""
    obs = [
        (6, 8, 1, 109.0, "C O H a"),
        (6, 8, 1, 109.5, "C O H b"),
        (6, 8, 1, 110.0, "C O H c"),
    ]
    types, type_map = cluster_angles(obs, tolerance=5.0)

    assert len(types) == 1
    assert type_map == [1, 1, 1]
    assert types[0].theta_0 == pytest.approx(109.5)
    assert types[0].obs_count == 3


def test_cluster_two_far_observations():
    """Two observations whose separation exceeds the
    tolerance form two separate clusters.  Type ids are assigned in finalize
    order, which matches sorted-theta
    order for non-merging inputs."""
    obs = [(6, 8, 1, 109.5, "tag_a"), (6, 8, 1, 120.0, "tag_b"),]
    types, type_map = cluster_angles(obs, tolerance=5.0)

    assert len(types) == 2
    # 109.5 finalizes first (smaller theta), so it gets type id 1; 120.0
    # finalizes after, so id 2.
    assert type_map == [1, 2]
    assert types[0].theta_0 == pytest.approx(109.5)
    assert types[1].theta_0 == pytest.approx(120.0)
    assert types[0].obs_count == 1
    assert types[1].obs_count == 1


def test_cluster_drift_chain_breaks():
    """A long chain of small steps cannot silently sweep
    into a single cluster.

    Each candidate is checked against the running mean (`within_tol`) and,
    defensively, against the 2*tolerance spread cap (`within_cap`).  The
    mathematical consequence is that a uniformly-spaced chain longer than a few
    times tolerance must break into more than one cluster.  This test verifies
    the chain-breaking behaviour without locking in the exact split point, which
    is an arithmetic detail of the running-mean walk.
    """
    obs = [(6, 8, 1, 100.0 + step, f"step{step}")
           for step in range(15)]
    types, type_map = cluster_angles(obs, tolerance=5.0)

    # The chain MUST break -- a 15-step linear chain with tolerance 5.0 cannot
    # sit inside any single cluster.
    assert len(types) > 1
    # Every observation must be assigned a non-zero type id (the "unassigned"
    # sentinel never appears for non-empty input).
    assert all(tid != 0 for tid in type_map)
    # Every emitted cluster must report at least one observation -- no empty
    # clusters slip through.
    assert all(t.obs_count >= 1 for t in types)


def test_cluster_multi_triplet_independence():
    """Observations from different element triplets are
    binned independently before clustering, so triplets
    cannot influence one another's type membership."""
    obs = [
        (6, 8, 1, 109.5, "C O H 1"),    # triplet A
        (6, 6, 6, 120.0, "C C C 1"),    # triplet B
        (6, 8, 1, 109.7, "C O H 2"),    # triplet A
        (6, 6, 6, 119.8, "C C C 2"),    # triplet B
        (1, 6, 1, 104.5, "H C H 1"),    # triplet C
    ]
    types, type_map = cluster_angles(obs, tolerance=5.0)

    # Three distinct triplets, each producing one cluster because each triplet's
    # members fit within tolerance.
    assert len(types) == 3
    ## Same-triplet pairs share type ids.
    assert type_map[0] == type_map[2]      # C-O-H pair shares
    assert type_map[1] == type_map[3]      # C-C-C pair shares

    ## Different-triplet members never share a type.
    assert type_map[4] != type_map[0]
    assert type_map[4] != type_map[1]
    # Cluster sizes match the per-triplet input counts.
    counts = sorted(t.obs_count for t in types)
    assert counts == [1, 2, 2]


def test_cluster_base_tag_from_first_sorted_member():
    """The representative base_tag is taken from the FIRST
    observation in sorted-theta order, NOT from the input
    order and NOT from lexicographic ordering of the tag
    strings themselves."""
    obs = [
        (6, 8, 1, 110.0, "second"),
        (6, 8, 1, 109.0, "first"),    # smallest theta
        (6, 8, 1, 109.5, "middle"),
    ]
    types, _ = cluster_angles(obs, tolerance=5.0)

    assert len(types) == 1
    assert types[0].base_tag == "first"


def test_cluster_empty_input_returns_empty_results():
    """An empty observation list produces empty outputs,
    not a crash and not a placeholder type."""
    types, type_map = cluster_angles([], tolerance=5.0)
    assert types == []
    assert type_map == []


# ---------------------------------------------------------
# get_angle_k -- PSEUDOCODE 10b
# ---------------------------------------------------------

def test_get_angle_k_geometric_mean_formula():
    """K_theta = stiffness * sqrt(K_arm1 * K_arm2) * scale.

    Drives the helper with a stub `get_bond_params` that returns known (K, r0)
    pairs so the formula shape can be verified independently of any production
    UFF lookup. Uses K_arm1 = 400, K_arm2 = 900 because their product (360000)
    has an exact integer square root (600), making the expected K_theta easy to
    read at a glance.
    """
    # Stub bond table: (z_a, z_b) -> (K, r0).  Only the bond stiffness K is
    # consumed by get_angle_k.
    bonds = {
        (6, 8): (400.0, 1.43),   # C-O arm
        (8, 6): (400.0, 1.43),   # symmetric entry
        (8, 1): (900.0, 0.96),   # O-H arm
        (1, 8): (900.0, 0.96),   # symmetric entry
    }
    def stub_bonds(za, zb):
        return bonds[(za, zb)]

    K_theta = get_angle_k(z1=6, zv=8, z2=1, stiffness=0.15, scale=1.0,
        get_bond_params=stub_bonds)

    ## K_theta = 0.15 * sqrt(400 * 900) * 1.0
    ##        = 0.15 * sqrt(360000)
    ##        = 0.15 * 600.0
    ##        = 90.0
    expected = 0.15 * math.sqrt(400.0 * 900.0) * 1.0
    assert K_theta == pytest.approx(expected)
    assert K_theta == pytest.approx(90.0)


def test_get_angle_k_scale_rescales_uniformly():
    """A change in `scale` multiplies K_theta proportionally
    while leaving the geometric-mean shape intact."""
    bonds = {
        (6, 8): (400.0, 0.0), (8, 6): (400.0, 0.0),
        (8, 1): (900.0, 0.0), (1, 8): (900.0, 0.0),
    }
    def stub_bonds(za, zb):
        return bonds[(za, zb)]

    K_one = get_angle_k(6, 8, 1, 0.15, 1.0, stub_bonds)
    K_two = get_angle_k(6, 8, 1, 0.15, 2.5, stub_bonds)
    assert K_two == pytest.approx(2.5 * K_one)


def test_get_angle_k_stiffness_rescales_uniformly():
    """A change in `stiffness` multiplies K_theta
    proportionally, mirroring the role of `scale`."""
    bonds = {
        (6, 8): (400.0, 0.0), (8, 6): (400.0, 0.0),
        (8, 1): (900.0, 0.0), (1, 8): (900.0, 0.0),
    }
    def stub_bonds(za, zb):
        return bonds[(za, zb)]

    K_low  = get_angle_k(6, 8, 1, 0.10, 1.0, stub_bonds)
    K_high = get_angle_k(6, 8, 1, 0.40, 1.0, stub_bonds)
    assert K_high == pytest.approx(4.0 * K_low)
