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
    AngleType, cluster_angles, get_angle_k,
    LocalRecord, FinalAngleType, cross_source_cluster)


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


# ---------------------------------------------------------
# cross_source_cluster -- PSEUDOCODE 10e
# ---------------------------------------------------------

def _record(source, local_id, z1, zv, z2, theta, obs_count,
            base_tag="tag"):
    """Construct one LocalRecord with sensible defaults.

    Most cross_source_cluster tests only vary the fields that
    matter for the scenario under test (theta, obs_count, source,
    local_id), so this factory keeps each test case tight and
    readable instead of re-specifying every slot at the call
    site.
    """
    return LocalRecord(
        z1=z1, zv=zv, z2=z2,
        theta_0_local=theta, obs_count=obs_count,
        base_tag=base_tag,
        source=source, local_type_id=local_id)


def test_cross_source_single_record():
    """One LocalRecord produces one FinalAngleType whose
    canonical theta_0 and obs_count total are forwarded
    verbatim, and a remap with exactly one entry."""
    rec = _record("lammps.dat", 1, 6, 8, 1, 109.5, 4,
                  base_tag="C 1 m O 1 m H 1 m")
    final_types, remap = cross_source_cluster(
        [rec], tolerance=5.0)

    assert len(final_types) == 1
    ft = final_types[0]
    assert ft.z1 == 6 and ft.zv == 8 and ft.z2 == 1
    assert ft.theta_0_final == pytest.approx(109.5)
    assert ft.obs_count_total == 4
    assert ft.base_tag == "C 1 m O 1 m H 1 m"
    assert remap == {("lammps.dat", 1): 1}


def test_cross_source_identical_theta_across_sources():
    """Two records from different sources with the same
    (triplet, theta) merge into one final type.  Their
    obs_count totals add; their (source, local_type_id) keys
    both map to the same final id."""
    recs = [
        _record("lammps.dat", 1, 6, 8, 1, 109.5, 3),
        _record("preRxn.data", 1, 6, 8, 1, 109.5, 2),
    ]
    final_types, remap = cross_source_cluster(recs, tolerance=5.0)

    assert len(final_types) == 1
    assert final_types[0].theta_0_final == pytest.approx(109.5)
    assert final_types[0].obs_count_total == 5
    assert remap[("lammps.dat", 1)] == 1
    assert remap[("preRxn.data", 1)] == 1


def test_cross_source_far_apart_stay_split():
    """Two records beyond tolerance remain as separate
    final types and remap to distinct ids."""
    recs = [
        _record("lammps.dat", 1, 6, 8, 1, 109.5, 3),
        _record("preRxn.data", 1, 6, 8, 1, 120.0, 3),
    ]
    final_types, remap = cross_source_cluster(recs, tolerance=5.0)

    assert len(final_types) == 2
    # 109.5 finalizes first (sorted-ascending), so it gets id 1.
    assert final_types[0].theta_0_final == pytest.approx(109.5)
    assert final_types[1].theta_0_final == pytest.approx(120.0)
    assert remap[("lammps.dat", 1)] == 1
    assert remap[("preRxn.data", 1)] == 2


def test_cross_source_obs_count_weighted_mean():
    """A cluster anchored by many observations pulls the
    final theta_0 more strongly than a sparse one.

    Uses theta_a=108.0 with obs_count=9 and theta_b=112.0 with
    obs_count=1.  The unweighted mean would be 110.0, but the
    obs_count-weighted mean is (9*108 + 1*112) / 10 = 108.4, so
    the sparse record pulls the final theta_0 less than one
    degree despite being 4 degrees from the dense anchor.
    """
    recs = [
        _record("lammps.dat", 1, 6, 8, 1, 108.0, 9),
        _record("preRxn.data", 1, 6, 8, 1, 112.0, 1),
    ]
    final_types, remap = cross_source_cluster(recs, tolerance=5.0)

    assert len(final_types) == 1
    assert final_types[0].theta_0_final == pytest.approx(108.4)
    assert final_types[0].obs_count_total == 10
    assert remap[("lammps.dat", 1)] == 1
    assert remap[("preRxn.data", 1)] == 1


def test_cross_source_spread_cap_breaks_drift_chain():
    """A chain of closely-spaced records cannot silently
    merge across a wide span.  With tolerance=2.5 and a
    spread_cap of 5.0, a 105..113 chain must break somewhere."""
    recs = [
        _record("s1", i + 1, 6, 8, 1, theta, 1)
        for i, theta in enumerate([105.0, 107.0, 109.0,
                                    111.0, 113.0])
    ]
    final_types, remap = cross_source_cluster(recs, tolerance=2.5)

    # The chain MUST break -- 8 degrees of span cannot sit in
    # one 5-degree spread cap.
    assert len(final_types) > 1
    # Every record must receive a remap entry.
    assert len(remap) == 5
    # Final type ids are 1-indexed and cover the len(final_types)
    # range without gaps.
    assert set(remap.values()) == set(
        range(1, len(final_types) + 1))


def test_cross_source_multi_triplet_independence():
    """Records from different (z1, zv, z2) triplets are
    binned independently, so triplets cannot influence one
    another's final-type membership."""
    recs = [
        _record("s1", 1, 6, 8, 1, 109.5, 1),   # triplet A
        _record("s1", 2, 6, 6, 6, 120.0, 1),   # triplet B
        _record("s2", 1, 6, 8, 1, 109.7, 1),   # triplet A
        _record("s2", 2, 6, 6, 6, 119.8, 1),   # triplet B
        _record("s3", 1, 1, 6, 1, 104.5, 1),   # triplet C
    ]
    final_types, remap = cross_source_cluster(recs, tolerance=5.0)

    # Three distinct triplets, each producing one final type.
    assert len(final_types) == 3
    # Same-triplet pairs from different sources map to the same
    # final id.
    assert remap[("s1", 1)] == remap[("s2", 1)]   # C-O-H pair
    assert remap[("s1", 2)] == remap[("s2", 2)]   # C-C-C pair
    # Different triplets never share a final id.
    assert remap[("s3", 1)] != remap[("s1", 1)]
    assert remap[("s3", 1)] != remap[("s1", 2)]


def test_cross_source_base_tag_from_first_sorted_member():
    """The representative base_tag of each final cluster is
    taken from the member with the smallest theta_0_local
    (the first after sort), not from input order."""
    recs = [
        _record("s1", 1, 6, 8, 1, 110.0, 1, base_tag="second"),
        _record("s2", 1, 6, 8, 1, 109.0, 1, base_tag="first"),
        _record("s3", 1, 6, 8, 1, 109.5, 1, base_tag="middle"),
    ]
    final_types, _ = cross_source_cluster(recs, tolerance=5.0)

    assert len(final_types) == 1
    assert final_types[0].base_tag == "first"


def test_cross_source_empty_input_returns_empty_results():
    """An empty local_records list produces empty outputs,
    not a crash and not a placeholder type."""
    final_types, remap = cross_source_cluster([], tolerance=5.0)
    assert final_types == []
    assert remap == {}


def test_cross_source_associativity_with_local_premerge():
    """Collapsing duplicate observations at the producer does
    not change what cross_source_cluster computes.

    Path A: three separate records with obs_count=1 each, all
    at theta=109.5 from one source.
    Path B: one record with obs_count=3 at theta=109.5 from
    that source.

    Both paths must produce the same final theta_0 and the same
    obs_count_total.  This is the associativity property that
    underpins make_reactions.py's tolerance=0 policy (DESIGN
    4.8.8 item 4a).
    """
    # Path A: three separate records.
    recs_a = [
        _record("s1", 1, 6, 8, 1, 109.5, 1),
        _record("s1", 2, 6, 8, 1, 109.5, 1),
        _record("s1", 3, 6, 8, 1, 109.5, 1),
    ]
    final_a, _ = cross_source_cluster(recs_a, tolerance=5.0)

    # Path B: one merged record.
    recs_b = [_record("s1", 1, 6, 8, 1, 109.5, 3)]
    final_b, _ = cross_source_cluster(recs_b, tolerance=5.0)

    assert len(final_a) == len(final_b) == 1
    assert final_a[0].theta_0_final == pytest.approx(
        final_b[0].theta_0_final)
    assert final_a[0].obs_count_total == (
        final_b[0].obs_count_total)
