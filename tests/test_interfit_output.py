"""Tests for interFit.py OLCAO output writer (Phase 4).

These tests validate the write_olcao_atom_type, write_olcao_pot_type, and
write_olcao_snippet methods by constructing mock orbital_fitting objects
with known fitting results and verifying the output format against the
atomicTypes.f90 and potTypes.f90 readers.

No actual Gaussian fitting is performed — the tests inject synthetic
coefficients and alpha values.
"""

import io
import math
import re
import pytest
import numpy as np


# ---------------------------------------------------------------------------
# Helpers to build mock orbital_fitting objects without reading Grasp2K files
# ---------------------------------------------------------------------------

class _MockAtomicInfo:
    def __init__(self, Z):
        self.atomic_number = Z


class _MockOrbitalInfo:
    def __init__(self, names, l_list, lvals, max_l):
        self.fitting_orbital_names = names
        self.fitting_l_list = l_list
        self.fitting_lvals = lvals
        self.fitting_max_l = max_l


class _MockGrasp:
    def __init__(self, Z, names, l_list, lvals, max_l):
        self.atomic_info = _MockAtomicInfo(Z)
        self.orbital_info = _MockOrbitalInfo(names, l_list, lvals, max_l)


def _make_mock_fitter(Z, names, l_list, lvals, max_l, num_terms, alphas,
                      coefficients):
    """Build a mock orbital_fitting with just enough state for the writer.

    Avoids the real __init__ (which reads isodata/rwfn.out) by constructing
    from scratch.
    """
    from element_data import ElementData

    class _Fitter:
        pass

    # Give it the writer methods from orbital_fitting
    import interFit
    _Fitter._fortran_e = interFit.fortran_e
    _Fitter._parse_orbital_info = interFit.orbital_fitting._parse_orbital_info
    _Fitter._classify_orbitals = interFit.orbital_fitting._classify_orbitals
    _Fitter.write_olcao_atom_type = interFit.orbital_fitting.write_olcao_atom_type
    _Fitter.write_olcao_pot_type = interFit.orbital_fitting.write_olcao_pot_type
    _Fitter.write_olcao_snippet = interFit.orbital_fitting.write_olcao_snippet

    obj = _Fitter()
    obj.grasp_description = _MockGrasp(Z, names, l_list, lvals, max_l)
    obj.best_num_terms = num_terms
    obj.best_alphas = alphas
    obj.best_coefficients = np.array(coefficients)
    obj.best_total_RMSE = 1e-4
    obj.best_orbital_RMSEs = [1e-4] * len(names)

    obj.element_data = ElementData()
    obj.element_data.init_element_data()
    return obj


# ---------------------------------------------------------------------------
# Fixture: Carbon (Z=6) MB single-component
# ---------------------------------------------------------------------------

@pytest.fixture
def carbon_mb():
    """Mock fitter for Carbon Z=6 with MB single-component orbitals.

    Orbitals: 1s (core), 2s (vale MB), 2p (vale MB).
    The fitting stores them highest-l-first, highest-n-first: [2p, 2s, 1s].
    num_terms: [12, 12, 0, 0, 0]  (s/p share 12 terms, no d/f/g).
    """
    Z = 6
    # Fitting order from single_component: highest l first, highest n first
    names = ['2p', '2s', '1s']
    l_list = [1, 0, 0]
    lvals = [2, 1, 0, 0, 0]  # 2 s-type, 1 p-type
    max_l = 1
    num_terms = [12, 12, 0, 0, 0]

    alphas = [0.12 * (50000.0 / 0.12) ** (i / 11.0) for i in range(12)]

    # 12 alphas x 3 orbitals — synthetic coefficients
    rng = np.random.default_rng(42)
    coefficients = rng.standard_normal((12, 3))

    return _make_mock_fitter(Z, names, l_list, lvals, max_l, num_terms,
                             alphas, coefficients)


# ---------------------------------------------------------------------------
# Fixture: Aluminium (Z=13) FB single-component
# ---------------------------------------------------------------------------

@pytest.fixture
def aluminium_fb():
    """Mock fitter for Al Z=13 with FB single-component orbitals.

    Core: 1s, 2s, 2p (core_orbitals = [2, 1, 0, 0]).
    Vale MB: 3s, 3p.
    Vale FB: 4s, 4p, 3d.
    num_terms: [20, 20, 12, 0, 0].
    """
    Z = 13
    # single_component order: highest l first, highest n first
    # d: 3d; p: 4p, 3p, 2p; s: 4s, 3s, 2s, 1s
    names = ['3d', '4p', '3p', '2p', '4s', '3s', '2s', '1s']
    l_list = [2, 1, 1, 1, 0, 0, 0, 0]
    lvals = [4, 3, 1, 0, 0]  # 4 s-type, 3 p-type, 1 d-type
    max_l = 2
    num_terms = [20, 20, 12, 0, 0]

    alphas = [0.12 * (100000.0 / 0.12) ** (i / 19.0) for i in range(20)]

    # 20 alphas x 8 orbitals — synthetic coefficients
    rng = np.random.default_rng(13)
    coefficients = rng.standard_normal((20, 8))

    return _make_mock_fitter(Z, names, l_list, lvals, max_l, num_terms,
                             alphas, coefficients)


# ---------------------------------------------------------------------------
# Tests: Fortran E format
# ---------------------------------------------------------------------------

class TestFortranEFormat:

    def test_positive(self):
        import interFit
        result = interFit.fortran_e(0.12)
        assert result.strip() == '0.12000000E+00'
        assert len(result) == 18

    def test_negative(self):
        import interFit
        result = interFit.fortran_e(-9.1196922e-06)
        assert result.strip().startswith('-0.')
        assert 'E-05' in result
        assert len(result) == 18

    def test_zero(self):
        import interFit
        result = interFit.fortran_e(0.0)
        assert result.strip() == '0.00000000E+00'

    def test_large(self):
        import interFit
        result = interFit.fortran_e(1e6)
        assert 'E+07' in result
        assert len(result) == 18

    def test_round_trip(self):
        """Fortran-formatted value should parse back to ~same float."""
        import interFit
        val = 3.14159265
        s = interFit.fortran_e(val).strip()
        parsed = float(s)
        assert abs(parsed - val) / val < 1e-7


# ---------------------------------------------------------------------------
# Tests: Orbital classification
# ---------------------------------------------------------------------------

class TestOrbitalClassification:

    def test_carbon_core_vale_split(self, carbon_mb):
        orbitals = carbon_mb._parse_orbital_info()
        core, vale = carbon_mb._classify_orbitals(orbitals, 6)
        # Carbon core: [1, 0, 0, 0] → only 1s is core
        assert len(core) == 1
        assert core[0]['name'] == '1s'
        assert core[0]['basis_code'] == 1
        # Valence: 2s, 2p
        assert len(vale) == 2
        vale_names = [o['name'] for o in vale]
        assert '2s' in vale_names
        assert '2p' in vale_names

    def test_aluminium_core_vale_split(self, aluminium_fb):
        orbitals = aluminium_fb._parse_orbital_info()
        core, vale = aluminium_fb._classify_orbitals(orbitals, 13)
        # Al core: [2, 1, 0, 0] → 1s, 2s, 2p
        assert len(core) == 3
        core_names = [o['name'] for o in core]
        assert '1s' in core_names
        assert '2s' in core_names
        assert '2p' in core_names
        # Valence: 3s(MB), 3p(MB), 4s(FB), 4p(FB), 3d(FB)
        assert len(vale) == 5
        # Check basis codes
        for o in vale:
            if o['name'] in ('3s', '3p'):
                assert o['basis_code'] == 1, f"{o['name']} should be MB"
            elif o['name'] in ('4s', '4p', '3d'):
                assert o['basis_code'] == 2, f"{o['name']} should be FB"

    def test_orbital_sort_order(self, aluminium_fb):
        orbitals = aluminium_fb._parse_orbital_info()
        # Should be sorted: l ascending, n ascending within l
        pairs = [(o['l'], o['n']) for o in orbitals]
        assert pairs == sorted(pairs)


# ---------------------------------------------------------------------------
# Tests: Atom type output format
# ---------------------------------------------------------------------------

class TestAtomTypeOutput:

    def test_tags_present(self, carbon_mb):
        buf = io.StringIO()
        carbon_mb.write_olcao_atom_type(buf)
        text = buf.getvalue()

        for tag in ['NUM_ATOM_TYPES', 'ATOM_TYPE_ID__SEQUENTIAL_NUMBER',
                     'ATOM_TYPE_LABEL', 'NUM_ALPHA_S_P_D_F', 'ALPHAS',
                     'NUM_CORE_RADIAL_FNS', 'NUM_VALE_RADIAL_FNS']:
            assert tag in text, f"Missing tag: {tag}"

    def test_label_format(self, carbon_mb):
        buf = io.StringIO()
        carbon_mb.write_olcao_atom_type(buf)
        text = buf.getvalue()
        # Label should be element_name + "1_1"
        assert 'c1_1' in text

    def test_alpha_counts_monotonic(self, aluminium_fb):
        """Alpha counts must be monotonically decreasing (s >= p >= d >= f)."""
        buf = io.StringIO()
        aluminium_fb.write_olcao_atom_type(buf)
        lines = buf.getvalue().splitlines()
        idx = lines.index('NUM_ALPHA_S_P_D_F')
        counts = list(map(int, lines[idx + 1].split()))
        assert len(counts) == 4
        for i in range(len(counts) - 1):
            assert counts[i] >= counts[i + 1], \
                f"Alpha counts not monotonically decreasing: {counts}"

    def test_alpha_count_matches_coefficients(self, aluminium_fb):
        """Each orbital's coefficient count must match its l-channel's alpha count."""
        buf = io.StringIO()
        aluminium_fb.write_olcao_atom_type(buf)
        lines = buf.getvalue().splitlines()

        # Parse alpha counts
        idx = lines.index('NUM_ALPHA_S_P_D_F')
        alpha_counts = list(map(int, lines[idx + 1].split()))

        # Parse total alphas
        idx = lines.index('ALPHAS')
        n_alphas = alpha_counts[0]
        alpha_lines = math.ceil(n_alphas / 4)

        # Find NL_RADIAL_FUNCTIONS blocks and verify coefficient counts
        nl_indices = [i for i, line in enumerate(lines)
                      if line.strip() == 'NL_RADIAL_FUNCTIONS']

        for nl_idx in nl_indices:
            i = nl_idx + 1
            while i < len(lines):
                line = lines[i].strip()
                # Check if we hit a new tag (all-caps with underscores)
                if re.match(r'^[A-Z_]+$', line):
                    break
                vals = line.split()
                if len(vals) == 2:
                    # numComponents, basisCode line
                    i += 1
                    # QN line: n, l, 2j, numStates, componentIndex
                    qn = list(map(int, lines[i].split()))
                    l = qn[1]
                    expected_coeffs = alpha_counts[l]
                    coeff_lines = math.ceil(expected_coeffs / 4)
                    # Read coefficient lines and count values
                    actual_count = 0
                    for cl in range(coeff_lines):
                        i += 1
                        actual_count += len(lines[i].split())
                    assert actual_count == expected_coeffs, \
                        f"Orbital n={qn[0]} l={l}: expected {expected_coeffs} " \
                        f"coefficients, got {actual_count}"
                i += 1

    def test_core_vale_counts(self, carbon_mb):
        buf = io.StringIO()
        carbon_mb.write_olcao_atom_type(buf)
        lines = buf.getvalue().splitlines()

        idx = lines.index('NUM_CORE_RADIAL_FNS')
        core_counts = list(map(int, lines[idx + 1].split()))
        assert core_counts == [1, 1, 1]  # 1s only, same for MB/FB/EB

        idx = lines.index('NUM_VALE_RADIAL_FNS')
        vale_counts = list(map(int, lines[idx + 1].split()))
        assert vale_counts == [2, 2, 2]  # 2s, 2p — all MB

    def test_al_core_vale_counts(self, aluminium_fb):
        buf = io.StringIO()
        aluminium_fb.write_olcao_atom_type(buf)
        lines = buf.getvalue().splitlines()

        idx = lines.index('NUM_CORE_RADIAL_FNS')
        core_counts = list(map(int, lines[idx + 1].split()))
        assert core_counts == [3, 3, 3]  # 1s, 2s, 2p

        idx = lines.index('NUM_VALE_RADIAL_FNS')
        vale_counts = list(map(int, lines[idx + 1].split()))
        # MB: 3s, 3p = 2; FB: +4s, 4p, 3d = 5; EB: same = 5
        assert vale_counts == [2, 5, 5]

    def test_quantum_numbers(self, carbon_mb):
        """Verify QN encoding: 2j=2l, numStates=(2l+1)*2, componentIdx=1."""
        buf = io.StringIO()
        carbon_mb.write_olcao_atom_type(buf)
        lines = buf.getvalue().splitlines()

        nl_indices = [i for i, line in enumerate(lines)
                      if line.strip() == 'NL_RADIAL_FUNCTIONS']

        for nl_idx in nl_indices:
            i = nl_idx + 1
            while i < len(lines):
                line = lines[i].strip()
                if re.match(r'^[A-Z_]+$', line):
                    break
                vals = line.split()
                if len(vals) == 2:
                    i += 1
                    qn = list(map(int, lines[i].split()))
                    n, l, two_j, num_states, comp_idx = qn
                    assert two_j == 2 * l, f"Expected 2j={2*l}, got {two_j}"
                    assert num_states == (2 * l + 1) * 2, \
                        f"Expected numStates={(2*l+1)*2}, got {num_states}"
                    assert comp_idx == 1
                    assert n >= 1
                i += 1


# ---------------------------------------------------------------------------
# Tests: Potential type output format
# ---------------------------------------------------------------------------

class TestPotTypeOutput:

    def test_tags_present(self, carbon_mb):
        buf = io.StringIO()
        carbon_mb.write_olcao_pot_type(buf)
        text = buf.getvalue()
        for tag in ['NUM_POTENTIAL_TYPES',
                     'POTENTIAL_TYPE_ID__SEQUENTIAL_NUMBER',
                     'POTENTIAL_TYPE_LABEL', 'NUCLEAR_CHARGE__ALPHA',
                     'COVALENT_RADIUS', 'NUM_ALPHAS', 'ALPHAS']:
            assert tag in text, f"Missing tag: {tag}"

    def test_nuclear_charge(self, carbon_mb):
        buf = io.StringIO()
        carbon_mb.write_olcao_pot_type(buf)
        lines = buf.getvalue().splitlines()
        idx = lines.index('NUCLEAR_CHARGE__ALPHA')
        vals = lines[idx + 1].split()
        assert float(vals[0]) == 6.0
        assert float(vals[1]) == 20.0  # standard nuclear alpha

    def test_pot_label(self, carbon_mb):
        buf = io.StringIO()
        carbon_mb.write_olcao_pot_type(buf)
        text = buf.getvalue()
        assert 'c1_1' in text

    def test_pot_alphas_are_min_max(self, carbon_mb):
        """ALPHAS line in potential section should be min and max values."""
        buf = io.StringIO()
        carbon_mb.write_olcao_pot_type(buf)
        lines = buf.getvalue().splitlines()
        idx = lines.index('ALPHAS')
        vals = lines[idx + 1].split()
        assert len(vals) == 2
        min_a, max_a = float(vals[0]), float(vals[1])
        assert min_a > 0
        assert max_a > min_a


# ---------------------------------------------------------------------------
# Tests: Combined snippet
# ---------------------------------------------------------------------------

class TestSnippetOutput:

    def test_snippet_contains_both_sections(self, carbon_mb, tmp_path):
        outfile = str(tmp_path / 'waveFn.dat')
        carbon_mb.write_olcao_snippet(outfile)

        with open(outfile) as f:
            text = f.read()

        assert 'NUM_ATOM_TYPES' in text
        assert 'NUM_POTENTIAL_TYPES' in text

    def test_snippet_raises_without_fitting(self):
        """Should raise RuntimeError if best_coefficients is None."""
        import interFit

        class _Empty:
            best_coefficients = None
        _Empty.write_olcao_snippet = interFit.orbital_fitting.write_olcao_snippet

        obj = _Empty()
        with pytest.raises(RuntimeError, match="No fitting results"):
            obj.write_olcao_snippet('test.dat')


# ---------------------------------------------------------------------------
# Tests: Reference format comparison (Al bulk)
# ---------------------------------------------------------------------------

class TestReferenceFormat:

    def test_al_alpha_field_widths(self, aluminium_fb):
        """Each alpha field should be 18 chars wide (Fortran E18.8)."""
        buf = io.StringIO()
        aluminium_fb.write_olcao_atom_type(buf)
        lines = buf.getvalue().splitlines()
        idx = lines.index('ALPHAS')
        # Check first alpha line
        alpha_line = lines[idx + 1]
        # 4 fields of 18 chars = 72 chars
        assert len(alpha_line) == 72, \
            f"Expected 72-char alpha line, got {len(alpha_line)}: '{alpha_line}'"

    def test_al_integer_field_widths(self, aluminium_fb):
        """Integer fields should be 8 chars wide."""
        buf = io.StringIO()
        aluminium_fb.write_olcao_atom_type(buf)
        lines = buf.getvalue().splitlines()
        idx = lines.index('NUM_ALPHA_S_P_D_F')
        # 4 fields of 8 chars = 32 chars
        assert len(lines[idx + 1]) == 32
