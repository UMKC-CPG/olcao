"""angle_utils.py -- Geometry-derived angle clustering and the
matching UFF-derived angular force constants.

This module supplies the two pure-computation helpers shared by
`condense.create_lammps_files` and `make_reactions.py` to replace the legacy
`angles.dat` lookup with a procedure that discovers angle types directly from
the observed bond angles in each structure or reaction template.  See DESIGN
sections 4.8.3 and 4.8.4 for the physical and algorithmic rationale, and
PSEUDOCODE sections 10a and 10b for the line-by-line specification that this
module implements.

Why this module lives here (and not inside `condense.py` or
`structure_control.py`)
-----------------------------------------------------------
  * The helpers are shared by two top-level scripts
    (`condense.py` and `make_reactions.py`).  Putting them inside either one
    would force the other to import the larger module just to reach the helpers.
  * The work is force-field data preparation, not atomic-
    structure manipulation, so `structure_control.py` is the wrong domain
    regardless of size considerations.
  * `get_angle_k` needs a bond-stiffness lookup, but to keep
    `angle_utils.py` dependency-free the lookup is passed in as a callable
    parameter rather than imported from `condense.BondData`.  The producers each
    bind their own `get_bond_params` instance at the call site.

Public surface
--------------
AngleType
    NamedTuple holding the per-cluster output record produced by
    `cluster_angles`: `(Z1, Zv, Z2, theta_0, obs_count, base_tag)`.

LocalRecord
    NamedTuple holding one (source, local_type_id) entry that a producer feeds
    into `cross_source_cluster`: `(z1, zv, z2, theta_0_local, obs_count,
    base_tag, source, local_type_id)`.  See PSEUDOCODE 10e.

FinalAngleType
    NamedTuple holding one global cluster produced by `cross_source_cluster`:
    `(z1, zv, z2, theta_0_final, obs_count_total, base_tag)`.

cluster_angles(observations, tolerance)
    Greedy single-pass local clustering of observed bond angles grouped by
    element triplet (PSEUDOCODE 10a).  Returns the discovered angle types and
    a per-observation type-id map.  Used by PSEUDOCODE 10c
    (`create_lammps_files`) and 10d (`make_reactions.py`).

get_angle_k(z1, zv, z2, stiffness, scale, get_bond_params)
    Harmonic angular spring constant for one angle type, derived from the UFF
    bond stiffnesses of its two arms (PSEUDOCODE 10b).

cross_source_cluster(local_records, tolerance)
    Cross-source unification of the local types emitted by every producer
    (PSEUDOCODE 10e).  Same greedy-merge structure as `cluster_angles`, but
    the running mean is obs_count-weighted so a densely-observed local cluster
    anchors the final theta_0 more strongly than a sparse one.  Returns the
    global final types and a remap from (source, local_type_id) to global id.
    Used by `normalize_types` in `condense.py`.
"""

import math
from collections import namedtuple


# ---------------------------------------------------------
# Output record
# ---------------------------------------------------------
# Slot ordering deliberately matches PSEUDOCODE 10e so that `local_records` from
# individual producers and `final_types` from the cross-source merge can be
# compared and rewritten by index without per-source juggling.  Per-field
# documentation (preserved as a structured listing -- see CLAUDE.md Structured
# Comment Blocks):
##
##   Z1, Zv, Z2 : atomic numbers (canonicalized so Z1 <= Z2).
##   theta_0    : arithmetic mean (degrees) of all theta_obs values
##                merged into this cluster.
##   obs_count  : number of raw observations that fed this cluster --
##                the weight used by 10e for the obs_count-weighted
##                cross-source merge.
##   base_tag   : the FIRST member's base_tag, kept as the cluster's
##                representative tag prefix.  10c and 10d append the
##                rest-angle and the type id to this prefix; 10e/10f
##                carry it across sources unchanged.
AngleType = namedtuple('AngleType',
    ['Z1', 'Zv', 'Z2', 'theta_0', 'obs_count', 'base_tag'])


# LocalRecord is the per-source input consumed by
# `cross_source_cluster` (PSEUDOCODE 10e).  One LocalRecord is
# emitted per (source, local_type_id) pair after a producer
# finishes its own local clustering pass (PSEUDOCODE 10c for
# `create_lammps_files`, PSEUDOCODE 10d for `make_reactions.py`):
##
##   z1, zv, z2      : canonical atomic numbers (z1 <= z2).
##   theta_0_local   : rest angle of this local cluster in
##                     degrees, as written to disk by the
##                     producer.
##   obs_count       : number of raw observations feeding this
##                     local cluster.  For reaction templates
##                     the producer uses tolerance=0, so this
##                     counts bit-identical duplicates only;
##                     for lammps.dat the producer uses
##                     angle_cluster_tolerance, so this can be
##                     larger.  See DESIGN 4.8.10.
##   base_tag        : representative element/species/molecule
##                     tag prefix captured by `cluster_angles`.
##   source          : stable source identifier -- "lammps.dat"
##                     or the template file's basename.
##   local_type_id   : 1-indexed id of this local cluster within
##                     its source.  (source, local_type_id) is
##                     the key used in the remap returned by
##                     `cross_source_cluster`.
LocalRecord = namedtuple('LocalRecord',
    ['z1', 'zv', 'z2', 'theta_0_local', 'obs_count',
     'base_tag', 'source', 'local_type_id'])


# FinalAngleType is one entry in the output list of
# `cross_source_cluster`.  Each final type represents one
# chemically-distinct angle environment unified across every
# source:
##
##   z1, zv, z2       : canonical atomic numbers.
##   theta_0_final    : cross-source obs_count-weighted running
##                      mean of all contributing
##                      theta_0_local values.
##   obs_count_total  : sum of obs_count across every
##                      contributing record.
##   base_tag         : base_tag from the first contributing
##                      record, carried verbatim so that
##                      species_id and molecule_id metadata
##                      (which Z alone cannot reconstruct) stay
##                      attached to the final unified type.
FinalAngleType = namedtuple('FinalAngleType',
    ['z1', 'zv', 'z2', 'theta_0_final', 'obs_count_total',
     'base_tag'])


# ---------------------------------------------------------
# Local clustering helper -- PSEUDOCODE 10a
# ---------------------------------------------------------

def cluster_angles(observations, tolerance):
    """Cluster observed bond angles by element triplet.

    Implements PSEUDOCODE 10a.  Walks the observation list, bins each entry by
    (Z1, Zv, Z2), sorts each bin by theta, and greedily merges sorted neighbors
    into the running cluster while BOTH within_tol and within_cap hold.

      - within_tol: |theta_obs - cluster_mean| <= tolerance
      - within_cap: (cluster_max - cluster_min) <= 2.0 * tolerance

    The spread cap (within_cap) prevents a long chain of tightly clustered
    observations from silently sweeping values from opposite ends of a wide
    distribution into one cluster.  The same cap is applied in PSEUDOCODE 10e
    for the cross-source reclustering inside `normalize_types`, so the
    user-tunable `angle_cluster_tolerance` parameter has uniform semantics
    across the local and cross-source phases.

    Parameters
    ----------
    observations : sequence
        Each entry is read by index.  Expected slots, in order:

          - obs[0]: Z1 -- atomic number of one arm end.
          - obs[1]: Zv -- atomic number of the vertex atom.
          - obs[2]: Z2 -- atomic number of the other arm end.  The
            producer must canonicalize so that Z1 <= Z2.
          - obs[3]: theta_obs -- observed angle in degrees.
          - obs[4]: base_tag -- producer's tag prefix (element
            names, species ids, and molecule ids); 10c and 10d append the
            rest-angle and the type id later.

        Additional trailing slots (a1, atom, a2, ...) are ignored here;
        producers carry them through for use in their own later phases.

    tolerance : float
        Maximum deviation (degrees) between an observation and its cluster's
        running mean for the observation to merge into that cluster.  Defaults
        to 5.0 upstream, settable in condense.in via the keyword
        `angle_cluster_tolerance`.

    Returns
    -------
    angle_types : list of AngleType
        Discovered angle types in finalization order.  Each record carries `(Z1,
        Zv, Z2, theta_0, obs_count, base_tag)`; see AngleType's docstring for
        the slot semantics.

    angle_type_map : list of int
        Length == len(observations).  Entry i holds the 1-based id of the
        angle_type into which observation i was merged.  Type id 0 is reserved
        as the "unassigned" sentinel and is never written here for a non-empty
        `observations` list.
    """
    # Phase 1 -- Bin observations by element triplet.  Each bin entry retains
    # the observed theta, the original observation index (so we can write
    # angle_type_map at the end), and the producer's representative base_tag for
    # that observation (which becomes the cluster's representative base_tag when
    # the entry is the first member added).
    groups = {}
    for idx, obs in enumerate(observations):
        key = (obs[0], obs[1], obs[2])
        groups.setdefault(key, []).append(
            (obs[3], idx, obs[4]))

    # Allocate the type list and the per-observation id map. angle_type_map is
    # indexed by the original observation index and holds 1-based type ids; type
    # id 0 is the "unassigned" sentinel, which only appears when the observation
    # list is empty.
    angle_types = []
    angle_type_map = [0] * len(observations)
    spread_cap = 2.0 * tolerance

    # Phase 2 -- Greedy clustering inside each triplet bin.  Sort by theta, then
    # walk left-to-right, merging the candidate into the running cluster while
    # BOTH within_tol (running-mean check) and within_cap (spread check) hold.
    # When either fails, finalize the current cluster as a new type and start a
    # fresh cluster at the candidate.
    for key, entries in groups.items():
        entries.sort(key=lambda entry: entry[0])

        # Initialize the first cluster with entry 0.  The cluster_rep_base_tag
        # captures the FIRST member's base_tag and is propagated to the emitted
        # AngleType record so 10c and 10d can build tag tails from it and 10e
        # can carry it across sources.
        first_theta, first_idx, first_tag = entries[0]
        cluster_rep_base_tag = first_tag
        cluster_sum     = first_theta
        cluster_count   = 1
        cluster_min     = first_theta
        cluster_max     = first_theta
        cluster_members = [first_idx]

        for i in range(1, len(entries)):
            candidate_theta = entries[i][0]
            candidate_idx   = entries[i][1]
            candidate_tag   = entries[i][2]
            cluster_mean    = cluster_sum / cluster_count
            new_max = max(cluster_max, candidate_theta)
            new_min = min(cluster_min, candidate_theta)

            within_tol = abs(candidate_theta - cluster_mean) <= tolerance
            within_cap = (
                (new_max - new_min) <= spread_cap)

            if within_tol and within_cap:
                # Merge the candidate into the running cluster.  The running
                # mean is recomputed next iteration from cluster_sum and
                # cluster_count, so we only need to update those running totals
                # plus the spread bookkeeping.
                cluster_sum    += candidate_theta
                cluster_count  += 1
                cluster_min     = new_min
                cluster_max     = new_max
                cluster_members.append(candidate_idx)
            else:
                # Finalize the current cluster as an AngleType.  cluster_count
                # becomes the obs_count (slot 5); cluster_rep_base_tag becomes
                # the representative base_tag (slot 6).
                theta_0 = cluster_sum / cluster_count
                type_id = len(angle_types) + 1
                angle_types.append(AngleType(
                    Z1=key[0], Zv=key[1], Z2=key[2],
                    theta_0=theta_0,
                    obs_count=cluster_count,
                    base_tag=cluster_rep_base_tag))
                for member_idx in cluster_members:
                    angle_type_map[member_idx] = type_id

                # Start a fresh cluster at the candidate.
                cluster_rep_base_tag = candidate_tag
                cluster_sum     = candidate_theta
                cluster_count   = 1
                cluster_min     = candidate_theta
                cluster_max     = candidate_theta
                cluster_members = [candidate_idx]

        # Finalize the last cluster of this triplet bin. Slot ordering matches
        # the in-loop finalize above so 10e and 10f can rely on a single layout.
        theta_0 = cluster_sum / cluster_count
        type_id = len(angle_types) + 1
        angle_types.append(AngleType(
            Z1=key[0], Zv=key[1], Z2=key[2],
            theta_0=theta_0,
            obs_count=cluster_count,
            base_tag=cluster_rep_base_tag))
        for member_idx in cluster_members:
            angle_type_map[member_idx] = type_id

    return angle_types, angle_type_map


# ---------------------------------------------------------
# Angular force constant -- PSEUDOCODE 10b
# ---------------------------------------------------------

def get_angle_k(z1, zv, z2, stiffness, scale,
                get_bond_params):
    """Harmonic angular spring constant from UFF bond data.

    Implements PSEUDOCODE 10b.  Computes the force constant K_theta of a LAMMPS
    `angle_style harmonic` term from the two UFF bond stiffnesses of the angle's
    arms using the formula shown below (a geometric mean of the two arm
    stiffnesses).

        K_theta = stiffness * sqrt(K_arm1 * K_arm2) * scale

    The geometric mean of the two arm stiffnesses is a simple, physically
    motivated mixing rule that scales correctly when the arms differ in bond
    strength (e.g.  a C-H arm paired with a C=O arm).  The dimensionless
    `stiffness` coefficient calibrates the geometric mean against representative
    angular force constants from the literature; `scale` is a global
    user-tunable multiplier that lets the user soften or stiffen all angles
    uniformly without re-fitting the calibration.

    Parameters
    ----------
    z1, zv, z2 : int
        Atomic numbers of one arm end, the vertex, and the other arm end.  The
        two element-pair lookups go through (z1, zv) and (zv, z2); the relative
        order of z1 vs z2 does not matter to this helper since the geometric
        mean is symmetric in K_arm1, K_arm2.

    stiffness : float
        Dimensionless calibration constant -- the `angle_stiffness_coeff` from
        condense.in (defaults to 0.15 upstream).

    scale : float
        Global multiplier on all angle force constants -- the
        `angle_parameter_scale` from condense.in (defaults to 1.0 upstream).

    get_bond_params : callable
        Bond-parameter lookup of the form `get_bond_params(za, zb) -> (K, r0)`.
        Only the first return value (the bond stiffness K) is used here.  The
        lookup is passed in to keep this module dependency-free; in production
        it is bound to `condense.BondData.get_bond_params`.

    Returns
    -------
    float
        K_theta in the same units as the bond stiffnesses produced by
        `get_bond_params`, scaled by `stiffness` and `scale`.
    """
    K_arm1, _ = get_bond_params(z1, zv)
    K_arm2, _ = get_bond_params(zv, z2)
    return stiffness * math.sqrt(K_arm1 * K_arm2) * scale


# ---------------------------------------------------------
# Cross-source clustering -- PSEUDOCODE 10e
# ---------------------------------------------------------

def cross_source_cluster(local_records, tolerance):
    """Unify local angle types across every source into one
    global cluster list.

    Implements PSEUDOCODE 10e.  Collects the local types emitted
    by each producer (PSEUDOCODE 10c for lammps.dat, 10d for
    every reaction template), groups them by canonical (z1, zv,
    z2) triplet, and greedy-merges within each triplet group
    using an obs_count-weighted running mean.  A spread cap of
    ``2 * tolerance`` is applied on top of the running-mean
    tolerance check to stop a long chain of closely-spaced
    records from silently sweeping a wide distribution into one
    final cluster -- the same cap `cluster_angles` applies at
    the local step (PSEUDOCODE 10a), so local and cross-source
    clustering are semantically consistent.

    The obs_count weighting reflects the physical intuition that
    a cluster anchored by 200 observations is a better estimator
    of the underlying angle than one anchored by 3.  It is also
    associative with respect to the identity-only pre-merge
    performed by `make_reactions.py` at tolerance=0: collapsing
    duplicate observations at the producer does not change what
    `cross_source_cluster` computes here (see DESIGN 4.8.8
    item 4a for the proof).

    Parameters
    ----------
    local_records : iterable of LocalRecord (or compatible)
        One entry per (source, local_type_id) pair from every
        source.  Each entry must expose the slots listed in the
        LocalRecord NamedTuple docstring above.

    tolerance : float
        The user's `angle_cluster_tolerance` (default 5.0
        degrees).  Candidates merge when both
        |theta - running_mean| <= tolerance and the resulting
        cluster span (max - min) stays within ``2 * tolerance``.

    Returns
    -------
    final_types : list of FinalAngleType
        Discovered global angle types in finalization order.
        Each final type carries the cross-source canonical
        theta_0 and the contributing obs_count total.

    remap : dict
        ``{(source, local_type_id): final_type_id}`` covering
        every input LocalRecord.  `final_type_id` is 1-indexed
        and matches the (index + 1) position of the entry in
        `final_types`.  Downstream callers use this dict to
        rewrite per-angle type ids in every source file and to
        build the unified Angle Coeffs section of lammps.dat.
    """
    # Phase 1 -- bin local records by canonical triplet.  The
    # triplet is an (int, int, int) tuple suitable for dict keys
    # and mirrors the binning pattern cluster_angles uses at the
    # local step.
    groups = {}
    for rec in local_records:
        key = (rec.z1, rec.zv, rec.z2)
        groups.setdefault(key, []).append(rec)

    final_types = []
    remap = {}
    spread_cap = 2.0 * tolerance

    # Phase 2 -- greedy merge within each triplet group, sorted
    # by theta_0_local.  The running mean is obs_count-weighted;
    # finalize the current cluster whenever either the tolerance
    # or the spread cap would be violated by the next record.
    for key, entries in groups.items():
        entries.sort(key=lambda rec: rec.theta_0_local)

        first = entries[0]
        cluster_weighted_sum = (
            first.theta_0_local * first.obs_count)
        cluster_weight = first.obs_count
        cluster_min = first.theta_0_local
        cluster_max = first.theta_0_local
        members = [first]

        for i in range(1, len(entries)):
            candidate = entries[i]
            candidate_theta = candidate.theta_0_local
            candidate_weight = candidate.obs_count
            running_mean = cluster_weighted_sum / cluster_weight
            new_min = min(cluster_min, candidate_theta)
            new_max = max(cluster_max, candidate_theta)

            within_tolerance = (
                abs(candidate_theta - running_mean) <= tolerance)
            within_spread_cap = (
                (new_max - new_min) <= spread_cap)

            if within_tolerance and within_spread_cap:
                # Merge the candidate into the running cluster.
                # The weighted sum grows by theta * obs_count and
                # the weight grows by obs_count, so the updated
                # running_mean is the obs_count-weighted mean of
                # every theta_0_local added so far.
                cluster_weighted_sum += (
                    candidate_theta * candidate_weight)
                cluster_weight += candidate_weight
                cluster_min = new_min
                cluster_max = new_max
                members.append(candidate)
            else:
                _finalize_cross_source_cluster(
                    key, members,
                    cluster_weighted_sum, cluster_weight,
                    final_types, remap)
                # Reset the running totals to start a fresh
                # cluster anchored at the candidate.
                cluster_weighted_sum = (
                    candidate_theta * candidate_weight)
                cluster_weight = candidate_weight
                cluster_min = candidate_theta
                cluster_max = candidate_theta
                members = [candidate]

        _finalize_cross_source_cluster(
            key, members,
            cluster_weighted_sum, cluster_weight,
            final_types, remap)

    return final_types, remap


def _finalize_cross_source_cluster(
        key, members,
        cluster_weighted_sum, cluster_weight,
        final_types, remap):
    """Append one finalized cluster to final_types and register
    its (source, local_type_id) -> final_type_id entries in the
    remap dictionary.

    The function is the "finalize" helper referenced by
    PSEUDOCODE 10e.  It is deliberately local-scope (leading
    underscore) because its callers are closely bound to the
    running-cluster loop in `cross_source_cluster`; no other
    module should need it.

    Parameters
    ----------
    key : tuple of (int, int, int)
        The (z1, zv, zv, z2) triplet shared by every member of
        this cluster.  Carried through to the FinalAngleType
        record so the downstream K recomputation in PSEUDOCODE
        10f Phase C can look up arm bond stiffnesses.

    members : list of LocalRecord
        Every record that merged into this cluster, in the
        order they were added (sorted-by-theta order within
        the triplet group).

    cluster_weighted_sum, cluster_weight : float
        Running totals of ``theta * obs_count`` and
        ``obs_count`` accumulated by the outer loop.  Their
        ratio is the cross-source canonical theta_0.

    final_types : list of FinalAngleType
        Accumulator appended to in place.

    remap : dict
        Accumulator written to in place.  Each member's
        (source, local_type_id) maps to the new final_type_id.
    """
    theta_0_final = cluster_weighted_sum / cluster_weight
    obs_count_total = cluster_weight
    final_type_id = len(final_types) + 1

    # Take the first member's base_tag as the representative
    # prefix for the final type.  The base_tag carries
    # species_id and molecule_id metadata that the (z1, zv, z2)
    # triplet alone cannot reconstruct, so the representative
    # prefix is propagated verbatim rather than rebuilt.
    representative_base_tag = members[0].base_tag

    final_types.append(FinalAngleType(
        z1=key[0], zv=key[1], z2=key[2],
        theta_0_final=theta_0_final,
        obs_count_total=obs_count_total,
        base_tag=representative_base_tag))

    for member in members:
        remap[(member.source, member.local_type_id)] = (
            final_type_id)
