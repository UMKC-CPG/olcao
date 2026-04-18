# Pseudocode

> **Document hierarchy:** VISION -> ARCHITECTURE -> DESIGN
> -> **PSEUDOCODE** -> Code. For the design rationale behind
> these algorithms, see `DESIGN.md`.

---

## 1. Generate Tetrahedra (DESIGN 1.2)

```
function generateTetrahedra(nA, nB, nC):
    numTetrahedra = 6 * nA * nB * nC
    allocate tetrahedra(4, numTetrahedra)
    t = 0
    for a = 1 to nA:
        for b = 1 to nB:
            for c = 1 to nC:
                # 8 corners with periodic wrapping
                M1 = idx(a,     b,     c    )
                M2 = idx(a+1,   b,     c    )
                M3 = idx(a,     b+1,   c    )
                M4 = idx(a,     b,     c+1  )
                M5 = idx(a+1,   b+1,   c    )
                M6 = idx(a+1,   b,     c+1  )
                M7 = idx(a,     b+1,   c+1  )
                M8 = idx(a+1,   b+1,   c+1  )

                # 6 tetrahedra sharing diagonal M1-M8
                tetrahedra(:, t+1) = [M1, M2, M5, M8]
                tetrahedra(:, t+2) = [M1, M3, M5, M8]
                tetrahedra(:, t+3) = [M1, M3, M7, M8]
                tetrahedra(:, t+4) = [M1, M4, M7, M8]
                tetrahedra(:, t+5) = [M1, M4, M6, M8]
                tetrahedra(:, t+6) = [M1, M2, M6, M8]
                t = t + 6

function idx(a, b, c):
    # Periodic wrapping, 1-based indexing
    return getIndexFromIndices(
        mod(a-1, nA) + 1,
        mod(b-1, nB) + 1,
        mod(c-1, nC) + 1)
```

---

## 2. LAT TDOS (DESIGN 1.3)

The TDOS at each energy grid point is the sum of per-
corner DOS weights from `bloechlCornerDOSWt` (section
2a), summed over all bands and tetrahedra. The per-corner
weights are also used individually by the PDOS (section
8.3).

```
function computeTDOS_LAT(eigenValues, tetrahedra,
        numTetrahedra, tetraVol,
        energyGrid, numEnergyPoints,
        numStates, numSpins,
        fullKPToIBZKPMap):
    allocate tdos(numEnergyPoints, numSpins) = 0.0

    for spin = 1 to numSpins:
        for n = 1 to numStates:
            for T = 1 to numTetrahedra:
                # Map full-mesh corners to IBZ
                # eigenvalues.
                for c = 1 to 4:
                    kFull = tetrahedra(c, T)
                    kIBZ = fullKPToIBZKPMap(kFull)
                    eps(c) = eigenValues(
                        n, kIBZ, spin)

                # Sort eigenvalues ascending.
                sortedEps = sort(eps)

                for iE = 1 to numEnergyPoints:
                    E = energyGrid(iE)
                    if E < sortedEps(1) or
                            E >= sortedEps(4):
                        cycle

                    # Per-corner DOS weights. The
                    # TDOS uses only their sum.
                    cornerDOSWt_LAT(1:4) =
                        bloechlCornerDOSWt(
                            E, sortedEps)

                    tdos(iE, spin) +=
                        sum(cornerDOSWt_LAT)
                        * tetraVol / spin
                        / hartree

    # Diagnostic: integrated area should equal
    # the number of spin states in the energy
    # range. Use deltaDOS * hartree because
    # deltaDOS is in Hartree but the TDOS is
    # in states/eV.
    integratedArea = trapezoid(tdos)
        * deltaDOS * hartree

    return tdos
```

---

## 2a. Bloechl Corner DOS Weights (DESIGN 1.3)

`bloechlCornerDOSWt` computes the per-corner DOS
density weights `cornerDOSWt_LAT(1:4)` for one
tetrahedron at energy E. These are the energy
derivatives of the cumulative corner integration
weights `cornerIntgWt_LAT` (section 3a):

  cornerDOSWt_LAT(c) = d/dE [ cornerIntgWt_LAT(c) ]

Their sum equals the total per-tetrahedron DOS
(the `bloechlDOS` value from the original TDOS
implementation). This identity provides a built-in
self-consistency check.

The derivation follows from the product rule applied
to the cumulative weight expressions in section 3a.
For each case, we reuse the same intermediate
variables (t_j, s_j, a, b, c, d, v_I, v_II, v_III)
and also compute the total DOS `gTotal`, then
decompose it across the four corners.

```
function bloechlCornerDOSWt(E, eps):
    # eps = [e1, e2, e3, e4] sorted ascending.
    # Returns cornerDOSWt_LAT(1:4): the per-corner
    # DOS density weights at energy E.
    #
    # cornerDOSWt_LAT(c) is the spectral density
    # (units: 1/energy) attributed to sorted
    # corner c. Their sum equals the total DOS
    # per unit BZ volume for this tetrahedron.

    e1 = eps(1);  e2 = eps(2)
    e3 = eps(3);  e4 = eps(4)
    tol = 1.0e-12

    # Case 0: outside eigenvalue range.
    if E < e1 or E >= e4:
        return [0, 0, 0, 0]

    # ------------------------------------------
    # Case 1: e1 <= E < e2
    # ------------------------------------------
    # From section 3a, the cumulative weights are:
    #   w(j) = f * t_j / 4     for j = 2,3,4
    #   w(1) = f - w(2) - w(3) - w(4)
    # where f = t2*t3*t4, t_j = (E-e1)/(e_j-e1).
    #
    # Applying the product rule:
    #   d(f*t_j)/dE = df/dE * t_j + f * dt_j/dE
    #               = gTotal * t_j + f/(e_j - e1)
    # where gTotal = df/dE = 3*(E-e1)^2 / denom
    #   is the total DOS for this case.
    #
    if E < e2:
        denom = (e2-e1) * (e3-e1) * (e4-e1)
        if abs(denom) < tol:
            return [0, 0, 0, 0]
        t2 = (E - e1) / (e2 - e1)
        t3 = (E - e1) / (e3 - e1)
        t4 = (E - e1) / (e4 - e1)
        f = t2 * t3 * t4
        gTotal = 3.0 * (E - e1)**2 / denom

        g(2) = (gTotal*t2 + f/(e2-e1)) / 4
        g(3) = (gTotal*t3 + f/(e3-e1)) / 4
        g(4) = (gTotal*t4 + f/(e4-e1)) / 4
        g(1) = gTotal - g(2) - g(3) - g(4)
        return g

    # ------------------------------------------
    # Case 2: e2 <= E < e3 (middle range)
    # ------------------------------------------
    # From section 3a, the cumulative weights use
    # three sub-tetrahedra volumes v_I, v_II,
    # v_III with intersection parameters a, b,
    # c, d. Their derivatives are:
    #
    #   da/dE = 1/e31,  db/dE = 1/e41
    #   dc/dE = 1/e32,  dd/dE = 1/e42
    #
    #   dv_I/dE   = da*b + a*db
    #             = b/e31 + a/e41
    #   dv_II/dE  = da*d*(1-b) + a*dd*(1-b)
    #             + a*d*(-db)
    #             = d*(1-b)/e31 + a*(1-b)/e42
    #             - a*d/e41
    #   dv_III/dE = (-da)*c*d + (1-a)*dc*d
    #            + (1-a)*c*dd
    #             = -c*d/e31 + (1-a)*d/e32
    #             + (1-a)*c/e42
    #
    # Each cumulative weight w(j) is a linear
    # combination of v_I, v_II, v_III with
    # coefficients that also depend on a, b, c, d.
    # Applying the product rule to each term
    # and collecting:
    #
    if E < e3:
        e31 = e3-e1;  e41 = e4-e1
        e32 = e3-e2;  e42 = e4-e2
        if e31*e41 < tol or e32*e42 < tol:
            return [0, 0, 0, 0]

        a = (E-e1) / e31
        b = (E-e1) / e41
        c_var = (E-e2) / e32
        d_var = (E-e2) / e42

        # Sub-tetrahedra volumes.
        v_I   = a * b
        v_II  = a * d_var * (1 - b)
        v_III = (1 - a) * c_var * d_var

        # Volume derivatives.
        dv_I   = b/e31 + a/e41
        dv_II  = d_var*(1-b)/e31
                 + a*(1-b)/e42
                 - a*d_var/e41
        dv_III = -c_var*d_var/e31
                 + (1-a)*d_var/e32
                 + (1-a)*c_var/e42

        # Parameter derivatives.
        da = 1/e31;  db = 1/e41
        dc = 1/e32;  dd = 1/e42

        # Corner 1: w(1) = [v_I*(3-a-b)
        #   + v_II*(2-a-b) + v_III*(1-a)] / 4
        g(1) = (dv_I*(3-a-b)
                + v_I*(-da - db)
                + dv_II*(2-a-b)
                + v_II*(-da - db)
                + dv_III*(1-a)
                + v_III*(-da)) / 4

        # Corner 2: w(2) = [v_I + v_II*(2-d)
        #   + v_III*(3-c-d)] / 4
        g(2) = (dv_I
                + dv_II*(2-d_var)
                + v_II*(-dd)
                + dv_III*(3-c_var-d_var)
                + v_III*(-dc - dd)) / 4

        # Corner 3: w(3) = [v_I*a + v_II*a
        #   + v_III*(a+c)] / 4
        g(3) = (dv_I*a + v_I*da
                + dv_II*a + v_II*da
                + dv_III*(a+c_var)
                + v_III*(da + dc)) / 4

        # Corner 4: w(4) = [v_I*b
        #   + v_II*(b+d) + v_III*d] / 4
        g(4) = (dv_I*b + v_I*db
                + dv_II*(b+d_var)
                + v_II*(db + dd)
                + dv_III*d_var
                + v_III*dd) / 4

        return g

    # ------------------------------------------
    # Case 3: e3 <= E < e4
    # ------------------------------------------
    # From section 3a, the cumulative weights use
    # the unoccupied sub-tet fraction f_un and
    # parameters s_j = (e4-E)/(e4-e_j). Their
    # derivatives are:
    #   ds_j/dE = -1/(e4-e_j)
    #   df_un/dE = -gTotal  (where gTotal is the
    #     total DOS for this case)
    #
    # Applying the product rule to
    #   w(j) = 1/4 - f_un*s_j/4 for j=1,2,3:
    #   dw(j)/dE = -(df_un*s_j + f_un*ds_j)/4
    #            = (gTotal*s_j + f_un/(e4-e_j))/4
    #
    denom = (e4-e1) * (e4-e2) * (e4-e3)
    if abs(denom) < tol:
        return [0, 0, 0, 0]
    s1 = (e4 - E) / (e4 - e1)
    s2 = (e4 - E) / (e4 - e2)
    s3 = (e4 - E) / (e4 - e3)
    f_un = s1 * s2 * s3
    gTotal = 3.0 * (e4 - E)**2 / denom

    g(1) = (gTotal*s1 + f_un/(e4-e1)) / 4
    g(2) = (gTotal*s2 + f_un/(e4-e2)) / 4
    g(3) = (gTotal*s3 + f_un/(e4-e3)) / 4
    g(4) = gTotal - g(1) - g(2) - g(3)
    return g
```

---

## 3. electronPopulation_LAT (DESIGN 1.5)

```
function computeElectronPopulation_LAT(
        eigenValues, tetrahedra,
        numTetrahedra, tetraVol,
        eFermi, numStates,
        numKPoints, numSpins):
    # Computes electronPopulation_LAT(n, k, spin):
    # the LAT analog of electronPopulation for
    # integrated properties (effective charge,
    # bond order). Each entry gives the fractional
    # electron occupation of state (n, k) as
    # determined by tetrahedron integration.
    allocate electronPopulation_LAT(
        numStates, numKPoints, numSpins) = 0.0

    for spin = 1 to numSpins:
        for n = 1 to numStates:
            for T = 1 to numTetrahedra:
                corners(1:4) = tetrahedra(1:4, T)
                eps_raw(1:4) =
                    eigenValues(n, corners, spin)

                # Sort and track permutation
                sigma = argsort(eps_raw)
                eps(1:4) = eps_raw(sigma)

                # Corner integration weights for
                # occupied-state integration
                # (Bloechl eqs. 18-21, evaluated
                # at the Fermi energy)
                cornerIntgWt_LAT(1:4) =
                    bloechlCornerWeights(
                        eFermi, eps)

                for i = 1 to 4:
                    ki = corners(sigma(i))
                    electronPopulation_LAT(
                        n, ki, spin) +=
                        cornerIntgWt_LAT(i)
                            * tetraVol

    return electronPopulation_LAT
```

---

## 3a. Bloechl Corner Integration Weights (DESIGN 1.5)

### Motivation

Section 2 (LAT TDOS) computes a single number g(E)
for each tetrahedron: the total DOS contribution at
energy E. Section 3 (electronPopulation_LAT) calls
`bloechlCornerWeights(E, eps)` to decompose the
tetrahedron's occupation into four separate corner
weights. This section derives those weights from first
principles and presents the pseudocode.

**Why corners need separate weights.** Each corner of
a tetrahedron corresponds to a different k-point with
different eigenvector projections (Mulliken populations,
orbital character). The total DOS does not need this
decomposition because it depends only on eigenvalues.
Partial properties (effective charge, bond order, PDOS)
require knowing how much of each corner's projection
to include. The corner weights provide exactly this
decomposition.

### Definitions

Consider a tetrahedron with four corners having sorted
eigenvalues e1 <= e2 <= e3 <= e4. Within the
tetrahedron, the eigenvalue is linearly interpolated
via barycentric coordinates:

  epsilon(r) = lambda_1 * e1 + lambda_2 * e2
             + lambda_3 * e3 + lambda_4 * e4

where lambda_i >= 0 and sum(lambda_i) = 1.

The **corner integration weight** w_j(E) is defined as
the integral of the j-th barycentric coordinate over
the occupied region {epsilon <= E}, normalized by the
tetrahedron volume V_T:

  w_j(E) = (1/V_T) * integral_{epsilon<=E} lambda_j dV

Since sum(lambda_j) = 1, the four weights sum to the
total occupied fraction:

  f(E) = sum_j w_j(E) = (1/V_T) * Vol({epsilon <= E})

### Key property: vertex averaging

The integral of any linear function L(r) over a
tetrahedron equals the volume times the average of L
at the four vertices:

  integral_T L dV = V_T * [L(v1)+L(v2)+L(v3)+L(v4)]/4

Since barycentric coordinates are linear functions,
this means: if we decompose the occupied region into
sub-tetrahedra, each corner weight contribution from
a sub-tetrahedron S is:

  w_j(S) = (V_S / V_T)
         * [sum of lambda_j at S's 4 vertices] / 4

### Case 0: trivial bounds

  E < e1:   w_j = 0 for all j  (empty region)
  E >= e4:  w_j = 1/4 for all j  (full tetrahedron)

### Case 1: e1 <= E < e2

Only corner 1 has eigenvalue below E. The occupied
region is a small tetrahedron with apex at corner 1,
cut by the iso-energy surface epsilon = E. The surface
intersects the three edges from corner 1 at:

  edge 1->j at parameter t_j = (E-e1)/(ej-e1)
                                for j = 2, 3, 4

The sub-tetrahedron has four vertices with barycentric
coordinates (lambda_1, lambda_2, lambda_3, lambda_4):

  corner 1:          (1,     0,   0,   0)
  edge 1->2 at t_2:  (1-t_2, t_2, 0,   0)
  edge 1->3 at t_3:  (1-t_3, 0,   t_3, 0)
  edge 1->4 at t_4:  (1-t_4, 0,   0,   t_4)

Volume ratio: f = t_2 * t_3 * t_4

Applying vertex averaging (summing lambda_j across
the 4 vertices, multiplying by f/4):

  w_2 = f * t_2 / 4
  w_3 = f * t_3 / 4
  w_4 = f * t_4 / 4
  w_1 = f - w_2 - w_3 - w_4

Verification: sum(w_j) = f.

### Case 4: e3 <= E < e4 (complement of Case 1)

The *unoccupied* region is a small tetrahedron near
corner 4. Define:

  s_j = (e4 - E) / (e4 - ej)    for j = 1, 2, 3

Unoccupied sub-tetrahedron vertices:

  corner 4:          (0,   0,   0,   1)
  edge 4->1 at s_1:  (s_1, 0,   0,   1-s_1)
  edge 4->2 at s_2:  (0,   s_2, 0,   1-s_2)
  edge 4->3 at s_3:  (0,   0,   s_3, 1-s_3)

Unoccupied fraction: f_unocc = s_1 * s_2 * s_3
Occupied fraction:   f = 1 - f_unocc

The occupied weights are the whole-tetrahedron weights
(1/4 each) minus the unoccupied contributions:

  w_1 = 1/4 - f_unocc * s_1 / 4
  w_2 = 1/4 - f_unocc * s_2 / 4
  w_3 = 1/4 - f_unocc * s_3 / 4
  w_4 = f - w_1 - w_2 - w_3

Verification: sum(w_j) = f.

### Case 2: e2 <= E < e3 (middle range)

Corners 1 and 2 lie below E; corners 3 and 4 lie
above. The iso-energy surface cuts four edges:

  edge 1->3 at  a = (E-e1)/(e3-e1)     point A
  edge 1->4 at  b = (E-e1)/(e4-e1)     point B
  edge 2->3 at  c = (E-e2)/(e3-e2)     point C
  edge 2->4 at  d = (E-e2)/(e4-e2)     point D

The occupied region is a pentahedron with vertices
{corner 1, corner 2, A, B, C, D}. We decompose it
into three sub-tetrahedra:

  T_I   = (corner 1, corner 2, A, B)
  T_II  = (corner 2, A, B, D)
  T_III = (corner 2, A, C, D)

The volume ratios follow from the determinant of
the 4x4 barycentric coordinate matrix for each
sub-tetrahedron:

  v_I   = a * b
  v_II  = a * d * (1 - b)
  v_III = (1 - a) * c * d

Occupied fraction: f = v_I + v_II + v_III

The barycentric coordinates at each vertex:

  T_I:
    corner 1   (1,     0,     0,   0)
    corner 2   (0,     1,     0,   0)
    A          (1-a,   0,     a,   0)
    B          (1-b,   0,     0,   b)

  T_II:
    corner 2   (0,     1,     0,   0)
    A          (1-a,   0,     a,   0)
    B          (1-b,   0,     0,   b)
    D          (0,     1-d,   0,   d)

  T_III:
    corner 2   (0,     1,     0,   0)
    A          (1-a,   0,     a,   0)
    C          (0,     1-c,   c,   0)
    D          (0,     1-d,   0,   d)

Summing lambda_j over the four vertices of each
sub-tetrahedron, multiplying by v_k/4, and summing
over sub-tetrahedra gives the corner weights:

  w_1 = [v_I*(3-a-b) + v_II*(2-a-b)
         + v_III*(1-a)] / 4

  w_2 = [v_I + v_II*(2-d)
         + v_III*(3-c-d)] / 4

  w_3 = [v_I*a + v_II*a
         + v_III*(a+c)] / 4

  w_4 = [v_I*b + v_II*(b+d)
         + v_III*d] / 4

Verification: for each sub-tetrahedron, the sum of
all four lambda_j at any vertex is 1, so the sum over
all four vertices is 4. Therefore:
  sum(w_j) = (4*v_I + 4*v_II + 4*v_III) / 4
           = v_I + v_II + v_III = f.

### Continuity between cases

The formulas are continuous at the case boundaries:

- At E = e2: Case 2 reduces to Case 1 because
  c = d = 0, so v_II = v_III = 0 and
  f = v_I = a*b = (e2-e1)^2 / [(e3-e1)(e4-e1)].
  All four corner weights match.

- At E = e3: Case 2 reduces to Case 4 because
  a = 1, so v_III = 0 and the occupied fraction
  equals 1 - (e4-e3)^2 / [(e4-e1)(e4-e2)].
  All four corner weights match.

### Derivative consistency with TDOS

The energy derivative of sum(w_j) must equal the
total per-tetrahedron DOS. This relationship is now
built into `bloechlCornerDOSWt` (section 2a): the
four `cornerDOSWt_LAT` values are defined so that
`sum(cornerDOSWt_LAT) = gTotal`. This was verified
numerically for the middle range: with e1=0, e2=1,
e3=3, e4=5, both the derivative of
f = v_I + v_II + v_III and the TDOS formula give
g(2) = 51/120.

### Pseudocode

```
function bloechlCornerWeights(E, eps):
    # eps = [e1, e2, e3, e4] sorted ascending.
    # Returns w(1:4): the integrated corner
    # weights at energy E for one tetrahedron.
    #
    # w(i) is the fraction of the tetrahedron's
    # occupation attributed to sorted corner i.
    # sum(w) = f(E), the occupied volume fraction.

    e1 = eps(1);  e2 = eps(2)
    e3 = eps(3);  e4 = eps(4)
    tol = 1.0e-12

    # Case 0: trivial bounds
    if E < e1:
        return [0, 0, 0, 0]
    if E >= e4:
        return [0.25, 0.25, 0.25, 0.25]

    # Case 1: e1 <= E < e2
    if E < e2:
        denom = (e2-e1) * (e3-e1) * (e4-e1)
        if abs(denom) < tol:
            return [0, 0, 0, 0]
        t2 = (E - e1) / (e2 - e1)
        t3 = (E - e1) / (e3 - e1)
        t4 = (E - e1) / (e4 - e1)
        f = t2 * t3 * t4
        w(2) = f * t2 / 4
        w(3) = f * t3 / 4
        w(4) = f * t4 / 4
        w(1) = f - w(2) - w(3) - w(4)
        return w

    # Case 2: e2 <= E < e3
    if E < e3:
        e31 = e3-e1;  e41 = e4-e1
        e32 = e3-e2;  e42 = e4-e2
        if e31*e41 < tol or e32*e42 < tol:
            return [0, 0, 0, 0]
        a = (E-e1) / e31
        b = (E-e1) / e41
        c = (E-e2) / e32
        d = (E-e2) / e42

        v_I   = a * b
        v_II  = a * d * (1 - b)
        v_III = (1 - a) * c * d

        w(1) = (v_I*(3-a-b) + v_II*(2-a-b)
                + v_III*(1-a)) / 4
        w(2) = (v_I + v_II*(2-d)
                + v_III*(3-c-d)) / 4
        w(3) = (v_I*a + v_II*a
                + v_III*(a+c)) / 4
        w(4) = (v_I*b + v_II*(b+d)
                + v_III*d) / 4
        return w

    # Case 3: e3 <= E < e4
    denom = (e4-e1) * (e4-e2) * (e4-e3)
    if abs(denom) < tol:
        return [0.25, 0.25, 0.25, 0.25]
    s1 = (e4 - E) / (e4 - e1)
    s2 = (e4 - E) / (e4 - e2)
    s3 = (e4 - E) / (e4 - e3)
    f_un = s1 * s2 * s3
    w(1) = 0.25 - f_un * s1 / 4
    w(2) = 0.25 - f_un * s2 / 4
    w(3) = 0.25 - f_un * s3 / 4
    w(4) = (1 - f_un) - w(1) - w(2) - w(3)
    return w
```

---

## 4. Build Atom Permutation Table (DESIGN 2.4)

The atom permutation table records, for each point group
operation R and each atom A, which atom B = R(A) the
operation maps A to.  This is the single piece of
infrastructure needed for correct IBZ unfolding of all
shell-summed quantities (Q*, bond order, PDOS modes 0-2).

The algorithm works in fractional (abc) coordinates
because the point group operations are stored in that
basis.  Cartesian atom positions are converted to
fractional using invRealVectors (= recipVectors / 2*pi).

```
function buildAtomPerm(numPointOps, abcPointOps,
                       numAtomSites, atomSites,
                       invRealVectors):
    # Returns atomPerm(numPointOps, numAtomSites)
    #   where atomPerm(R, A) = B means operation R
    #   maps atom A to atom B.

    allocate atomPerm(numPointOps, numAtomSites)

    # Convert all atom positions from Cartesian (xyz)
    # to fractional (abc) coordinates.
    allocate abcAtomPos(3, numAtomSites)
    for A = 1 to numAtomSites:
        for i = 1 to 3:
            abcAtomPos(i, A) =
                sum(invRealVectors(i,:)
                    * atomSites(A)%cartPos(:))

    # For each operation and atom, apply R to the
    # fractional position and find the matching atom.
    for R = 1 to numPointOps:
        for A = 1 to numAtomSites:

            # Apply the point group rotation to atom
            # A's fractional position.
            for i = 1 to 3:
                rotPos(i) =
                    sum(abcPointOps(i,:,R)
                        * abcAtomPos(:,A))

            # Wrap the rotated position into [0,1).
            for i = 1 to 3:
                rotPos(i) = modulo(rotPos(i), 1.0)

            # Search for the atom at the rotated
            # position. R preserves species, so only
            # atoms of the same type can match.
            atomPerm(R, A) = -1  # sentinel
            for B = 1 to numAtomSites:
                if atomType(B) != atomType(A):
                    cycle

                # Compute the difference, wrapped
                # into [-0.5, 0.5) on each axis.
                for i = 1 to 3:
                    diff(i) = rotPos(i)
                             - abcAtomPos(i, B)
                    diff(i) = diff(i)
                             - nint(diff(i))

                if all(|diff(:)| < threshold):
                    atomPerm(R, A) = B
                    exit  # found the match

            # Safety check: every atom must have a
            # match. If not, the point group or the
            # atom positions are inconsistent.
            if atomPerm(R, A) == -1:
                error("No match for atom", A,
                      "under operation", R)

    deallocate abcAtomPos
    return atomPerm
```

---

## 4a. Build Inverse Atom Permutation (DESIGN 1.4, 2.4)

The inverse atom permutation invAtomPerm(R, B) gives
the atom A such that atomPerm(R, A) = B, i.e.,
A = R^{-1}(B). It is used during LAT PDOS tetrahedron
corner assembly to map channel indices from full-mesh
k-points back to their IBZ representatives (see
section 8). Built in O_AtomicSites alongside atomPerm.

```
function buildInvAtomPerm(numPointOps,
                          numAtomSites,
                          atomPerm):
    allocate invAtomPerm(numPointOps,
                         numAtomSites)

    for R = 1 to numPointOps:
        for A = 1 to numAtomSites:
            B = atomPerm(R, A)
            invAtomPerm(R, B) = A

    return invAtomPerm
```

---

## 5. Save fullKPToIBZOpMap (DESIGN 2.4)

This augments the existing IBZ folding loop in
`initializeKPointMesh`.  The current code saves
`fullKPToIBZKPMap(k_full) = k_IBZ` (which IBZ point
does this full-mesh point fold onto).  We additionally
save `fullKPToIBZOpMap(k_full) = R` (which point group
operation maps the IBZ representative to k_full, in
the forward direction: R(k_IBZ) = k_full).

The change is a single integer store at the point where
a match is found, plus identity for the IBZ
representative itself.

```
# Inside initializeKPointMesh, after allocating
# fullKPToIBZOpMap(numFullMeshKP):

# When a new IBZ representative i is found:
    fullKPToIBZOpMap(i) = identityOpIndex

# When mesh point j matches IBZ point i under
# operation m (the existing isMatch==1 branch):
    fullKPToIBZOpMap(j) = m
```

**Identity operation index.**  The space group
database guarantees that the identity is always the
first point group operation (verified across all 759
space group files in share/spaceDB).  Therefore
`identityOpIndex = 1` -- no runtime search is needed.

---

## 6. Corrected Effective Charge (DESIGN 2.4)

The current code accumulates Q* directly into
`atomCharge(A)` using only the IBZ k-point's Mulliken
projection.  The fix loops over the star of each IBZ
k-point and distributes the projection into the
permuted atom index.

The star of IBZ k-point k_IBZ is the set of full-mesh
k-points that fold onto it.  This set is not stored
explicitly -- it is traversed by scanning
fullKPToIBZKPMap.

The key change is in the innermost accumulation.
The outer loop structure (spin h, kpoint i, band j,
basis function l, atom k) remains the same.  The
Mulliken projection `oneValeRealAccum` for atom k at
IBZ k-point i is computed exactly as before.  The
difference is where it is accumulated.

```
# Current code (incorrect with IBZ):
#   atomCharge(k, h) += oneValeRealAccum
#                       * statePopulation

# Corrected code:
#
# After the band loop (j) completes for kpoint i,
# we have accumulated per-atom projections for
# this kpoint:
#   ibzAtomProj(k) = sum over bands of
#       oneValeRealAccum(k,j) * statePopulation(j)
#
# Count the star size (number of full-mesh points
# that fold to this IBZ kpoint):
#   starSize = count(fullKPToIBZKPMap(:) == i)
#
# Distribute across the star:
#   for each full-mesh kpoint f where
#           fullKPToIBZKPMap(f) == i:
#       R = fullKPToIBZOpMap(f)
#       for A = 1 to numAtomSites:
#           atomCharge(atomPerm(R, A), h) +=
#               ibzAtomProj(A) / starSize

# The normalizer is starSize, not numFullMeshKP.
# statePopulation already encodes the full BZ-
# integration weight for the star of kpoint i
# (kPointWeight for Gaussian, or the summed
# tetrahedron corner weights for LAT -- both
# proportional to starSize).  ibzAtomProj is
# therefore the total contribution from kpoint i.
# Dividing by starSize distributes this total
# equally among the star members.  The sum across
# the star recovers ibzAtomProj, preserving the
# total charge.
```

**Alternative (equivalent, avoids scanning):**  The
same result can be achieved within the existing
kpoint loop by changing only the accumulation target.
Instead of accumulating at index k, accumulate at
atomPerm(R, k) for each operation R in the star.
But the star is implicit -- it requires collecting
all full-mesh k-points for IBZ index i.

The cleanest implementation collects the per-atom
projection for one IBZ k-point, then distributes it
in a separate inner loop over the star.

---

## 7. Corrected Bond Order (DESIGN 2.4)

The same star-distribution pattern applies to bond
order.  The Mulliken overlap between atoms A and B
at IBZ k-point i is computed as before.  The
correction distributes it into the permuted atom
pair.

```
# Current code (incorrect with IBZ):
#   bondOrder(A_bonded, k) +=
#       oneValeRealAccum * statePopulation

# Corrected code:
#
# After computing bondOrderRaw(A, B) for IBZ
# kpoint i (accumulated over bands j):
#
#   starSize = count(fullKPToIBZKPMap(:) == i)
#
#   for each full-mesh kpoint f where
#           fullKPToIBZKPMap(f) == i:
#       R = fullKPToIBZOpMap(f)
#       for each bonded pair (A, B):
#           A_rot = atomPerm(R, A)
#           B_rot = atomPerm(R, B)
#           bondOrder(A_rot, B_rot) +=
#               bondOrderRaw(A, B) / starSize
```

**Integration with existing loop structure.**  The
current `computeBond` accumulates bond order inside
the band loop (j) interleaved with charge
accumulation.  The IBZ correction requires a second
pass over the star after all bands are processed for
a given IBZ k-point.  This suggests restructuring:

  1. For each IBZ kpoint i:
     a. Read eigenvectors and overlap (unchanged)
     b. Loop over bands j, accumulate raw per-atom
        charge ibzAtomProj(A) and raw per-pair bond
        order ibzBondRaw(A, B) at IBZ indices
     c. Distribute ibzAtomProj and ibzBondRaw across
        the star of i using atomPerm

Step (c) is the only new code.  Steps (a) and (b)
are the existing computation, with the accumulation
target changed from the final arrays to temporary
per-IBZ-kpoint buffers.

**Weight convention.**  The raw projection is
weighted by statePopulation (from either
electronPopulation_LAT or electronPopulation).
Both already encode the full star-weighted BZ-
integration weight for each IBZ kpoint:

- Gaussian: statePopulation includes
  kPointWeight(i) = 2 * starSize(i) /
  numFullMeshKP, so ibzAtomProj and ibzBondRaw
  are proportional to starSize.
- LAT: statePopulation includes
  electronPopulation_LAT(j,i,h), which sums
  tetrahedron corner weights from all full-mesh
  points in the star, again proportional to
  starSize.

The star distribution divides by starSize to
extract the per-full-mesh-point contribution, then
deposits it at the permuted atom (or atom pair).
The sum across starSize members recovers the
original ibzAtomProj (or ibzBondRaw), so the
overall charge and bond order totals are preserved.

---

## 8. LAT PDOS (DESIGN 1.4)

The LAT PDOS requires Mulliken projections at all four
corners of each tetrahedron simultaneously. Since
eigenvectors exist only at IBZ k-points, a two-pass
design is required: first compute and store projections
at IBZ k-points, then integrate over tetrahedra with
on-the-fly IBZ unfolding of the channel index.

### 8.1 Channel Permutation Table

For efficiency the channel permutation is precomputed
as a lookup table channelPermTable(R, alpha) so the
inner loop avoids repeated decode/encode. Mode 0
needs no permutation. Mode 1 uses invAtomPerm
directly. Mode 2 remaps the atom index while
preserving the l-shell offset within the atom.

```
function buildChannelPermTable(
        detailCodePDOS, numPointOps,
        cumulDOSTotal, cumulNumDOS,
        numAtomSites, invAtomPerm):

    allocate channelPermTable(numPointOps,
                              cumulDOSTotal)

    if detailCodePDOS == 0:
        # Per-type, per-l: identity (type-level
        # sums are invariant under R).
        for R = 1 to numPointOps:
            for alpha = 1 to cumulDOSTotal:
                channelPermTable(R, alpha) = alpha
        return channelPermTable

    if detailCodePDOS == 1:
        # Per-atom total: channel = atom index.
        for R = 1 to numPointOps:
            for A = 1 to numAtomSites:
                channelPermTable(R, A) =
                    invAtomPerm(R, A)
        return channelPermTable

    if detailCodePDOS == 2:
        # Per-atom, per-l: remap atom index,
        # preserve l-shell offset.
        for R = 1 to numPointOps:
            for A = 1 to numAtomSites:
                permA = invAtomPerm(R, A)
                baseOld = cumulNumDOS(A)
                baseNew = cumulNumDOS(permA)
                nOrbitals = cumulNumDOS(A+1)
                          - cumulNumDOS(A)
                for off = 1 to nOrbitals:
                    channelPermTable(R,
                        baseOld + off) =
                        baseNew + off
        return channelPermTable
```

### 8.2 Pass 1: Compute Projections

Stream through IBZ k-points, read eigenvectors and
overlap from HDF5, compute Mulliken projections, and
store into projArray(channel, band, kIBZ). The
Mulliken computation is identical to the existing code
in computeDOS (waveFnSqrd, oneValeRealAccum).

```
function computeProjections(inSCF, h,
        numKPoints, numStates, numAtomSites,
        numAtomStates, pdosIndex, valeDim,
        cumulDOSTotal, spin):
    allocate projArray(cumulDOSTotal,
                       numStates, numKPoints)
    projArray = 0.0

    for i = 1 to numKPoints:
        # Read eigenvectors + overlap for this
        # IBZ kpoint and spin orientation.
        readData(h, i, numStates, overlapCode=1)

        for j = 1 to numStates:
            valeDimIndex = 0
            for k = 1 to numAtomSites:
                for l = 1 to numAtomStates(k):
                    valeDimIndex += 1

                    # Compute Mulliken projection
                    # (existing waveFnSqrd *
                    # valeValeOL dot product).
                    oneValeRealAccum =
                        mullikenProjection(
                            valeDimIndex, j)

                    # Accumulate into the channel
                    # determined by pdosIndex.
                    ch = pdosIndex(valeDimIndex)
                    projArray(ch, j, i) +=
                        oneValeRealAccum
                        / real(spin)

    return projArray
```

### 8.3 Pass 2: Tetrahedron Integration

Loop over bands and tetrahedra. For each tetrahedron,
sort corner eigenvalues with tracked permutation,
compute `bloechlCornerDOSWt` at each energy point,
and accumulate weighted projections into pdosComplete.
The channel permutation table handles IBZ unfolding of
the projection index.

Note: this uses `bloechlCornerDOSWt` (section 2a),
which returns per-corner DOS density weights (units:
1/energy). The cumulative corner weights from
`bloechlCornerWeights` (section 3a) are NOT used here
-- those are for integrated properties only (section
3, `electronPopulation_LAT`).

```
function integratePDOS_LAT(projArray,
        channelPermTable,
        eigenValues, tetrahedra,
        numTetrahedra, tetraVol,
        fullKPToIBZKPMap, fullKPToIBZOpMap,
        energyScale, numEnergyPoints,
        numStates, cumulDOSTotal, spin):
    allocate pdosComplete(cumulDOSTotal,
                          numEnergyPoints)
    pdosComplete = 0.0

    for n = 1 to numStates:
        for T = 1 to numTetrahedra:
            # Look up corner info from the full
            # mesh, mapping to IBZ eigenvalues.
            for c = 1 to 4:
                kFull(c) = tetrahedra(c, T)
                kIBZ(c) =
                    fullKPToIBZKPMap(kFull(c))
                opIdx(c) =
                    fullKPToIBZOpMap(kFull(c))
                eps(c) =
                    eigenValues(n, kIBZ(c), h)

            # Sort eigenvalues ascending, tracking
            # permutation: sigma(i) = original
            # corner index in sorted position i.
            sigma = argsort(eps)
            sortedEps = eps(sigma)

            for iE = 1 to numEnergyPoints:
                E = energyScale(iE)

                # Skip if outside eigenvalue range.
                if E < sortedEps(1) or
                        E >= sortedEps(4):
                    cycle

                # Per-corner DOS density weights
                # for the sorted eigenvalues.
                cornerDOSWt_LAT(1:4) =
                    bloechlCornerDOSWt(
                        E, sortedEps)

                # Accumulate weighted projections
                # into pdosComplete. Each sorted
                # corner c maps back to original
                # corner sigma(c), whose IBZ kpoint
                # and operation index determine the
                # projection lookup.
                for c = 1 to 4:
                    orig = sigma(c)
                    R = opIdx(orig)
                    kIc = kIBZ(orig)

                    for alpha = 1 to cumulDOSTotal:
                        permA =
                            channelPermTable(
                                R, alpha)
                        pdosComplete(alpha, iE) +=
                            cornerDOSWt_LAT(c)
                            * tetraVol / hartree
                            * projArray(
                                permA, n, kIc)

    return pdosComplete
```

### 8.4 Normalization

For the LAT path, the corner DOS weights from
`bloechlCornerDOSWt` provide exact BZ integration
(no broadening artifacts). The electronFactor ratio
(currentPopulation / totalElectronsComputed) should
be ≈ 1.0. Compute and log it as a diagnostic but do
not apply it to pdosComplete. A ratio significantly
different from 1.0 signals an integration bug.

The "Spin States Calculated" diagnostic (integral of
totalSystemDos over the energy grid) must use
`deltaDOS * hartree` in the trapezoidal rule because
deltaDOS is stored in Hartree while the DOS is in
states/eV. This applies to both the Gaussian and LAT
paths.

---

## 9. UFF Bond Parameter Computation (DESIGN 4.2)

Given two atomic numbers, compute the UFF equilibrium
bond length and harmonic force constant.  The per-element
parameters (covalent radius r_i, effective charge Zstar_i,
GMP electronegativity chi_i) are read once from
`bond_parameters.dat` and stored in arrays indexed by
atomic number Z.

The prefactor 332.06 = 664.12 / 2 converts from the UFF
spring constant convention E = (1/2) k (r-r0)^2 to the
LAMMPS `bond_style harmonic` convention E = K (r-r0)^2.

```
# -------------------------------------------------
# Data structures (populated once by init_bond_data
# from bond_parameters.dat):
#
#   num_uff_elements : int
#       Number of elements in the table
#       (= maximum Z covered).
#   uff_r(Z)     : covalent radius (Angstroms)
#   uff_Zstar(Z) : effective charge
#   uff_chi(Z)   : GMP electronegativity (eV)
#
# These arrays are indexed by atomic number Z
# (1-based: uff_r(1) = hydrogen, etc.).
#
# The reader uses the Z column on each data line
# as the array index (not the sequential row
# number).  This makes the file order-independent.
# -------------------------------------------------

UFF_K_PREFACTOR = 332.06

function get_bond_params(z1, z2):
    # Compute UFF equilibrium bond length and
    # LAMMPS harmonic force constant for the
    # element pair (z1, z2).
    #
    # Inputs:
    #   z1, z2 : atomic numbers (order irrelevant;
    #            the formula is symmetric)
    #
    # Returns:
    #   K_ij : force constant (kcal/mol/A^2)
    #   r_ij : equilibrium bond length (Angstroms)
    #
    # Requires uff_r, uff_Zstar, uff_chi arrays
    # to be initialized from bond_parameters.dat.

    # --- Validate element coverage ---
    if z1 < 1 or z1 > num_uff_elements:
        error("Element Z =", z1,
              "not in bond_parameters.dat")
    if z2 < 1 or z2 > num_uff_elements:
        error("Element Z =", z2,
              "not in bond_parameters.dat")

    # --- Look up per-element parameters ---
    r1    = uff_r(z1)
    r2    = uff_r(z2)
    Zs1   = uff_Zstar(z1)
    Zs2   = uff_Zstar(z2)
    chi1  = uff_chi(z1)
    chi2  = uff_chi(z2)

    # --- Electronegativity correction ---
    # r_EN shortens the bond between elements of
    # unequal electronegativity.  For homonuclear
    # bonds (chi1 == chi2), r_EN = 0 and the bond
    # length is simply r1 + r2.
    denom_EN = chi1 * r1 + chi2 * r2
    if denom_EN > 0:
        r_EN = r1 * r2
               * (sqrt(chi1) - sqrt(chi2))**2
               / denom_EN
    else:
        r_EN = 0.0

    # --- Equilibrium bond length ---
    r_ij = r1 + r2 - r_EN

    # --- Force constant ---
    # Guard against zero or near-zero bond length
    # (should not occur for physical elements, but
    # protects against corrupt data).
    if r_ij <= 0:
        error("Non-positive bond length for Z =",
              z1, z2, "; check bond_parameters.dat")

    K_ij = UFF_K_PREFACTOR * Zs1 * Zs2
           / r_ij**3

    return K_ij, r_ij
```

**Usage in create_lammps_files and normalize_types.**
Both output paths contain a linear scan over
`hooke_bond_coeffs` to match element pairs to force
constants.  Both are replaced by a direct call to
`get_bond_params`:

```
# Current code (linear scan, in both paths):
#   for hb = 1 to num_hooke_bonds:
#       if atom1_z == hbc[hb][1]
#              and atom2_z == hbc[hb][2]:
#           k  = hbc[hb][3]
#           r0 = hbc[hb][4]
#
# New code (direct computation, both paths):
    K_ij, r_ij = bond_data.get_bond_params(
                     atom1_z, atom2_z)
    K_ij = K_ij * self.bond_parameter_scale
```

The `bond_parameter_scale` multiplier (default 1.0,
defined in `condenserc.py`, overridable in condense.in)
is applied after the UFF computation in **both**
output paths.  It scales only K_ij, not r_ij.

---

## 10. Angle Clustering and Force Constants (DESIGN 4.8)

Replace the `angles.dat` database lookup with a two-phase
procedure: (a) cluster observed bond angles by element
triplet to discover angle types, then (b) compute a force
constant for each type from UFF bond stiffnesses.

### 10a. Cluster Observed Angles by Triplet (DESIGN 4.8.3)

For each molecule, the bond analysis produces a list of
angles with atom indices and observed angle values.  The
`create_lammps_files` method already iterates over these
and extracts the element triplet (Z1, Zv, Z2).  The new
clustering step replaces the `angles.dat` lookup.

```
# -------------------------------------------------
# Data structures:
#
#   angle_cluster_tolerance : float
#       Maximum deviation (degrees) between an
#       observed angle and a cluster's running mean
#       for the angle to join that cluster.
#       Default: 5.0.  Read from condense.in.
#
#   spread_cap : float
#       Maximum allowed total span (max - min) of
#       theta values within any one cluster.  Set
#       to 2.0 * angle_cluster_tolerance.  Prevents
#       a long chain of closely-spaced observations
#       from silently sweeping values from opposite
#       ends of a wide distribution into a single
#       cluster.  The same cap is applied in 10e
#       for cross-source clustering.
#
#   Input: a list of angle observations, each
#       being (Z1, Zv, Z2, theta_obs, base_tag)
#       where Z1 <= Z2 (canonicalized).  The
#       base_tag is the producer's tag prefix
#       (element names, species ids, molecule
#       ids) for this specific atom triple.
#
#   Output:
#       angle_types : list of
#           (Z1, Zv, Z2, theta_0, obs_count,
#            base_tag)
#       angle_type_map : maps each observation
#           index to its angle_type index
#
#       obs_count is the number of observations
#       merged into the cluster.  The
#       representative base_tag is taken from
#       the first observation in the cluster.
#       The slot ordering (obs_count at slot 5,
#       base_tag at slot 6) matches 10e's
#       local_records and final_types tuples,
#       so 10c/10d/10e/10f use consistent slot
#       indices throughout.
# -------------------------------------------------

function cluster_angles(observations, tolerance):
    # Group observations by element triplet.
    # Each entry carries the observed theta, the
    # original observation index (for the
    # angle_type_map), and the producer's
    # representative base_tag for that observation.
    groups = {}
    for each (idx, obs) in enumerate(observations):
        key = (obs.Z1, obs.Zv, obs.Z2)
        groups[key].append(
            (obs.theta, idx, obs.base_tag))

    angle_types = []
    angle_type_map = array(len(observations))
    spread_cap = 2.0 * tolerance

    for each key in groups:
        # Sort angles within this triplet group.
        entries = groups[key]
        sort entries by theta ascending

        # Greedy clustering: walk the sorted list
        # and merge the candidate into the current
        # cluster while BOTH of these hold:
        #   (a) |theta - running_mean| <= tolerance
        #   (b) resulting (max - min) <= spread_cap
        # If either fails, finalize the current
        # cluster as a type and start a new cluster
        # at the candidate.  cluster_rep_base_tag is
        # captured from the first observation in
        # the cluster and propagated to the emitted
        # angle_type record so 10c/10d can build
        # tag tails and 10e can carry it across
        # sources.
        cluster_rep_base_tag = entries[0].base_tag
        cluster_sum   = entries[0].theta
        cluster_count = 1
        cluster_min   = entries[0].theta
        cluster_max   = entries[0].theta
        cluster_members = [entries[0].idx]

        for i = 1 to len(entries) - 1:
            cluster_mean = cluster_sum / cluster_count
            candidate_theta = entries[i].theta
            new_max = max(cluster_max, candidate_theta)
            new_min = min(cluster_min, candidate_theta)

            within_tol =
                |candidate_theta - cluster_mean|
                    <= tolerance
            within_cap =
                (new_max - new_min) <= spread_cap

            if within_tol and within_cap:
                # Merge into current cluster.
                cluster_sum += candidate_theta
                cluster_count += 1
                cluster_min = new_min
                cluster_max = new_max
                cluster_members.append(entries[i].idx)
            else:
                # Finalize current cluster as a type.
                # cluster_count is the number of raw
                # observations that fed this local
                # cluster (slot 5 = obs_count); the
                # first member's base_tag is the
                # representative prefix (slot 6).
                theta_0 = cluster_sum / cluster_count
                type_id = len(angle_types) + 1
                angle_types.append(
                    (key.Z1, key.Zv, key.Z2, theta_0,
                     cluster_count,
                     cluster_rep_base_tag))
                for m in cluster_members:
                    angle_type_map[m] = type_id

                # Start a new cluster at entry i.
                cluster_rep_base_tag =
                    entries[i].base_tag
                cluster_sum = candidate_theta
                cluster_count = 1
                cluster_min = candidate_theta
                cluster_max = candidate_theta
                cluster_members = [entries[i].idx]

        # Finalize the last cluster.  Slot
        # ordering matches the finalize above:
        # slot 5 = obs_count, slot 6 = base_tag.
        theta_0 = cluster_sum / cluster_count
        type_id = len(angle_types) + 1
        angle_types.append(
            (key.Z1, key.Zv, key.Z2, theta_0,
             cluster_count,
             cluster_rep_base_tag))
        for m in cluster_members:
            angle_type_map[m] = type_id

    return angle_types, angle_type_map
```

### 10b. Angle Force Constant (DESIGN 4.8.4)

Compute the harmonic angular spring constant for a
given angle type from the UFF bond stiffnesses of its
two arms.  This reuses `get_bond_params()` (section 9).

```
# -------------------------------------------------
# Data structures:
#
#   angle_stiffness_coeff : float
#       Dimensionless calibration constant that
#       converts the geometric mean of bond
#       stiffnesses into an angular stiffness.
#       Default: 0.15.  Read from condense.in.
#
#   angle_parameter_scale : float
#       Global multiplier on all angle force
#       constants.  Default: 1.0.
#       Read from condense.in.
# -------------------------------------------------

function get_angle_k(z1, zv, z2,
                     angle_stiffness_coeff,
                     angle_parameter_scale):
    # Compute the bond force constants for the
    # two arms of the angle: (z1, zv) and (zv, z2).
    K_arm1, _ = get_bond_params(z1, zv)
    K_arm2, _ = get_bond_params(zv, z2)

    # Geometric mean of arm stiffnesses, scaled
    # by the calibration constant and the global
    # user scale factor.
    K_angle = angle_stiffness_coeff
              * sqrt(K_arm1 * K_arm2)
              * angle_parameter_scale

    return K_angle
```

### 10c. Integration into create_lammps_files (DESIGN 4.8.8)

The existing angle loop in `create_lammps_files` extracts
(Z1, Zv, Z2) triplets and observed angles, then searches
`angles.dat` for a match.  The replacement runs the same
collect-cluster-emit structure as 10d, scoped to the
single lammps.dat file produced by `create_lammps_files`.
Both producers (10c and 10d) invoke the identical
`cluster_angles` helper from 10a, so local clustering
semantics are byte-identical.  Any residual theta_0
differences between 10c and 10d outputs are resolved by
10e during cross-source clustering inside
`normalize_types`.

```
# Phase 1: Collect all angle observations.  Replaces
# the angles.dat lookup loop.  The base_tag is built
# exactly as the current code builds it -- element
# names, species ids, and molecule ids -- with no
# rest-angle or type-id suffix appended yet (Phase 3
# adds those).
observations = []
for each atom with bond angles:
    for each angle_idx of atom:
        # Extract end atoms a1, a2 and vertex atom.
        z1 = element_z(a1)
        zv = element_z(atom)
        z2 = element_z(a2)
        if z1 > z2:
            swap z1, z2
            swap a1, a2
        theta_obs = bond_angles_ext[atom][angle_idx]
        base_tag = tag_string_for(a1, atom, a2)
        observations.append(
            (z1, zv, z2, theta_obs, base_tag,
             a1, atom, a2))

# Phase 2: Cluster locally using the shared helper
# from 10a.  Identical call signature and tolerance
# value as 10d.
angle_types, angle_type_map =
    cluster_angles(observations,
                   self.angle_cluster_tolerance)

# Phase 3: Build local type records with the
# cluster-mean theta_0 carried in the tag tail.
# These tables are local to this lammps.dat;
# normalize_types may merge them with types from
# reaction templates via 10e and rewrite both the
# tags and the per-angle ids in 10f.
num_local_angle_types = len(angle_types)
local_angle_tags   =
    [None] * (num_local_angle_types + 1)
local_angle_coeffs =
    [None] * (num_local_angle_types + 1)

for t = 1 to num_local_angle_types:
    atype = angle_types[t - 1]
    K = get_angle_k(
        atype.Z1, atype.Zv, atype.Z2,
        self.angle_stiffness_coeff,
        self.angle_parameter_scale)
    local_angle_coeffs[t] =
        [None, K, atype.theta_0]
    local_angle_tags[t] = (
        f"{atype.base_tag} "
        f"{atype.theta_0:.4f} {t}")

# Phase 4: Record per-atom angle entries with the
# local type ids.  normalize_types walks these and
# remaps the ids in 10f.  angle_bonded_atoms and
# ordered_angle_type follow the existing flat
# per-angle layout that the LAMMPS writer expects.
for i, obs in enumerate(observations):
    local_type_id = angle_type_map[i]
    angle_bonded_atoms.append(
        [None, obs.a1, obs.atom, obs.a2])
    ordered_angle_type.append(local_type_id)

# Export to normalize_types:
#   source tag = "lammps.dat"
#   local_angle_tags, local_angle_coeffs
#   angle_bonded_atoms, ordered_angle_type
#   per-local-type obs_count (slot 5 of each
#       entry in angle_types)
```

### 10d. Integration into make_reactions.py (DESIGN 4.8.8)

Mirrors 10c for the template-emission side of the
pipeline.  The existing Python port's angle construction
loop (around line 2518) iterates over angles in a
reaction template and searches `hooke_angle_coeffs` for a
matching row to build the tag tail.  The replacement runs
the same collect-cluster-emit structure as 10c, scoped to
one reaction template at a time.  Both producers invoke
the identical `cluster_angles` helper from 10a, so local
clustering semantics are byte-identical across sources
and any residual theta_0 differences between producers
are resolved by 10e downstream.

```
# Phase 1: Collect all angle observations for this
# reaction template.  Replaces the hooke_angle_coeffs
# scan near line 2518.  The base_tag is built exactly
# as today -- element names, species ids, and
# molecule ids, with no rest-angle or type-id suffix
# appended yet (Phase 3 adds those).
observations = []
for each vertex atom v in the template:
    for each (a1, a2) angle arm pair through v:
        z1 = element_z(a1)
        zv = element_z(v)
        z2 = element_z(a2)
        if z1 > z2:
            swap z1, z2
            swap a1, a2
        theta_obs = bond_angle(a1, v, a2)
        base_tag = tag_string_for(a1, v, a2)
        observations.append(
            (z1, zv, z2, theta_obs, base_tag,
             a1, v, a2))

# Phase 2: Cluster locally using the shared helper
# from 10a.  Same call signature as 10c, same
# angle_cluster_tolerance value.
angle_types, angle_type_map =
    cluster_angles(observations,
                   self.angle_cluster_tolerance)

# Phase 3: Build local type records with the
# cluster-mean theta_0 carried in the tag tail.
# These values are local to this template; the
# cross-source step in 10e may merge them with
# types from lammps.dat or from other templates.
num_local_angle_types = len(angle_types)
local_angle_tags   =
    [None] * (num_local_angle_types + 1)
local_angle_coeffs =
    [None] * (num_local_angle_types + 1)

for t = 1 to num_local_angle_types:
    atype = angle_types[t - 1]
    K = get_angle_k(
        atype.Z1, atype.Zv, atype.Z2,
        self.angle_stiffness_coeff,
        self.angle_parameter_scale)
    local_angle_coeffs[t] =
        [None, K, atype.theta_0]
    local_angle_tags[t] = (
        f"{atype.base_tag} "
        f"{atype.theta_0:.4f} {t}")

# Phase 4: Record per-atom angle entries with the
# local type ids.  normalize_types walks these
# and remaps the ids in 10f.
for i, obs in enumerate(observations):
    local_type_id = angle_type_map[i]
    angle_bonded[obs.v].append(
        [None, obs.a1, obs.a2])
    angle_tag_id[obs.v].append(local_type_id)

# Export to normalize_types:
#   source tag = "template:{name}"
#   local_angle_tags, local_angle_coeffs
#   angle_bonded, angle_tag_id
#   per-local-type obs_count (slot 5 of each
#       entry in angle_types)
```

### 10e. Cross-Source Angle Clustering (DESIGN 4.8.8 item 4a)

The first phase of angle handling inside
`normalize_types`.  Takes the per-source local cluster
centers emitted by 10c (lammps.dat) and 10d (each
reaction template) and merges any whose theta_0 values
represent the same physical angle.  This is what makes
bond/react type IDs consistent across sources.  The
algorithm is greedy merge with the same
`2 * tolerance` spread cap that 10a applies locally --
so local and cross-source clustering are semantically
consistent -- and adds observation-count weighting on
top, so a cluster anchored by many observations pulls
the final mean more strongly than a sparse one.

```
# -------------------------------------------------
# Data structures:
#
#   local_records : list, one entry per
#       (source, local_type_id) pair:
#           (z1, zv, z2,
#            theta_0_local,
#            obs_count,          # raw observations
#                                # feeding this
#                                # local cluster
#            base_tag,           # representative
#                                # tag prefix
#            source,             # "lammps.dat" or
#                                # "template:<name>"
#            local_type_id)
#
#   tolerance : float
#       angle_cluster_tolerance (default 5.0).
#
#   spread_cap : float
#       Max allowed total span (max-min) of
#       theta_0_local values within one final
#       cluster.  Default: 2.0 * tolerance.
#       Prevents greedy chaining from sweeping a
#       broad distribution into a single cluster.
#
#   Output:
#       final_types : list of
#           (z1, zv, z2,
#            theta_0_final,
#            obs_count_total,
#            representative_base_tag)
#       remap : dict
#           (source, local_type_id)
#               -> final_type_id
# -------------------------------------------------

function cross_source_cluster(local_records,
                              tolerance):
    # Group local records by canonical triplet.
    groups = {}
    for each rec in local_records:
        key = (rec.z1, rec.zv, rec.z2)
        groups[key].append(rec)

    final_types = []
    remap = {}
    spread_cap = 2.0 * tolerance

    for each key in groups:
        entries = groups[key]
        sort entries by theta_0_local ascending

        # Greedy merge, weighted by obs_count.
        cluster_w_sum = (entries[0].theta_0_local
                         * entries[0].obs_count)
        cluster_w     = entries[0].obs_count
        cluster_min   = entries[0].theta_0_local
        cluster_max   = entries[0].theta_0_local
        members       = [entries[0]]

        for i = 1 to len(entries) - 1:
            running_mean = cluster_w_sum / cluster_w
            candidate_theta = entries[i].theta_0_local
            new_max = max(cluster_max, candidate_theta)
            new_min = min(cluster_min, candidate_theta)

            within_tol =
                |candidate_theta - running_mean|
                    <= tolerance
            within_cap =
                (new_max - new_min) <= spread_cap

            if within_tol and within_cap:
                cluster_w_sum +=
                    candidate_theta * entries[i].obs_count
                cluster_w += entries[i].obs_count
                cluster_min = new_min
                cluster_max = new_max
                members.append(entries[i])
            else:
                finalize(key, members,
                         cluster_w_sum,
                         cluster_w,
                         final_types, remap)
                # Start a new cluster at entry i.
                cluster_w_sum = (candidate_theta
                    * entries[i].obs_count)
                cluster_w = entries[i].obs_count
                cluster_min = candidate_theta
                cluster_max = candidate_theta
                members = [entries[i]]

        finalize(key, members,
                 cluster_w_sum, cluster_w,
                 final_types, remap)

    return final_types, remap

function finalize(key, members,
                  cluster_w_sum, cluster_w,
                  final_types, remap):
    theta_0_final   = cluster_w_sum / cluster_w
    obs_count_total = cluster_w
    final_id = len(final_types) + 1
    # Take the first member's base_tag as the
    # representative prefix for the final type.
    # base_tag carries species_id / molecule_id
    # metadata that Z alone cannot reconstruct.
    representative_base_tag = members[0].base_tag
    final_types.append(
        (key.z1, key.zv, key.z2,
         theta_0_final, obs_count_total,
         representative_base_tag))
    for m in members:
        remap[(m.source, m.local_type_id)] =
            final_id
```

Decision notes embedded in the algorithm above:
- **Weighting.**  The running mean is `obs_count`-
  weighted, so a local cluster built from 200
  observations anchors the final theta_0 more
  strongly than one from 3.  This matches the
  physical intuition that the larger sample is a
  better estimator.
- **Spread cap.**  Greedy merge alone can chain
  across a wide distribution (e.g., observations at
  105, 107, 109, 111, 113 with tolerance 2.5 all
  collapse into one cluster spanning 8 degrees).
  The `spread_cap = 2 * tolerance` rule forces a
  split once the total span would exceed that
  bound, producing tighter clusters at distribution
  boundaries.
- **Canonical base_tag.**  Representative tag is
  taken from the first-encountered member rather
  than reconstructed from Z values, because the
  tag prefix carries species_id and molecule_id
  fields that Z alone does not encode.

### 10f. Tag Rewrite and Type-ID Remap (DESIGN 4.8.8 item 4c)

The second phase of angle handling inside
`normalize_types`, executed once 10e has produced
`(final_types, remap)`.  Every angle reference in every
source file is rewritten: the per-angle type id is
remapped to the global id, and the tag tail is replaced
with the final canonical theta_0 so any downstream tool
that inspects the tag sees a consistent value.  The
rewrite is deterministic given the cluster map, so
repeated runs on identical inputs produce byte-identical
output.

```
function rewrite_angles(sources, remap,
                        final_types,
                        angle_stiffness_coeff,
                        angle_parameter_scale):
    # Phase A: rewrite per-angle type ids in every
    # source.  lammps.dat carries an Angles section
    # with explicit type ids; each reaction
    # template carries a per-atom angle_tag_id
    # array.
    for each src in sources:
        if src is lammps.dat:
            for each angle entry in src.Angles:
                old_id = entry.type_id
                entry.type_id =
                    remap[(src.source_tag, old_id)]
        else:  # reaction template
            for each vertex atom v in src:
                for i in range(
                        len(src.angle_tag_id[v])):
                    old_id = src.angle_tag_id[v][i]
                    src.angle_tag_id[v][i] =
                        remap[(src.source_tag,
                               old_id)]

    # Phase B: build the unified global
    # unique_angle_tags table from final_types.
    # Each entry is
    #   "{representative_base_tag} {theta_0} {t}"
    # carrying the final canonical theta_0 and the
    # global type id.
    unique_angle_tags =
        [None] * (len(final_types) + 1)
    for t = 1 to len(final_types):
        ft = final_types[t - 1]
        unique_angle_tags[t] = (
            f"{ft.representative_base_tag} "
            f"{ft.theta_0_final:.4f} {t}")

    # Phase C: build the unified global
    # unique_angle_coeffs table via get_angle_k.
    # K_angle depends only on the triplet, so
    # recomputation here yields the same value
    # that any producer's local 10c/10d phase
    # computed -- cross-source merging does not
    # alter K_angle, only theta_0.
    unique_angle_coeffs =
        [None] * (len(final_types) + 1)
    for t = 1 to len(final_types):
        ft = final_types[t - 1]
        K = get_angle_k(
            ft.z1, ft.zv, ft.z2,
            angle_stiffness_coeff,
            angle_parameter_scale)
        unique_angle_coeffs[t] =
            [None, K, ft.theta_0_final]

    # Phase D: emit the cluster-map diagnostic.
    # For each final cluster, write:
    #   - global id
    #   - canonical theta_0_final
    #   - (z1, zv, z2)
    #   - every contributing
    #       (source, local_type_id,
    #        theta_0_local, obs_count) tuple
    # See DESIGN 4.8.8 item 4d.  This file is the
    # primary debuggability payback for routing
    # all clustering through normalize_types.
    write_cluster_map(final_types, remap,
                      local_records)

    return unique_angle_tags, unique_angle_coeffs
```

`normalize_types()`'s angle handling is therefore:
1. Gather `local_records` from every source.
2. `cross_source_cluster(local_records, tolerance)`
   -> `(final_types, remap)`.  (10e)
3. `rewrite_angles(sources, remap, final_types,
   angle_stiffness_coeff, angle_parameter_scale)`
   -> `(unique_angle_tags, unique_angle_coeffs)`.
   (10f)

No other changes are required inside
`normalize_types()`.  Bond handling (section 9) and
other type tables are unchanged.
