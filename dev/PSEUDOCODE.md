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
