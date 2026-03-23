"""conftest.py -- Shared fixtures and path helpers for the OLCAO test suite.

All test files in src/tests/ automatically receive the fixtures defined here.

Required environment variable
------------------------------
OLCAO_DATA
    Path to the OLCAO share directory (must contain elements.dat, spaceDB/,
    atomicBDB/, atomicPDB/).  Any test that instantiates StructureControl or
    ElementData will be *skipped* automatically when this variable is unset.

Quick-start
-----------
    # Activate the project virtualenv, then from the repository root:
    python src/tests/test_o.py                   # run everything
    python src/tests/test_o.py -m unit           # fast, pure-computation tests only
    python src/tests/test_o.py -m integration    # file I/O and pipeline tests
    python src/tests/test_o.py -v                # verbose output

Fixture overview
----------------
olcao_data_dir      session  $OLCAO_DATA string; skips session if unset
fresh_sc            function freshly initialised StructureControl (no structure)
make_sc             function factory: make_sc('bn_cubic.skl') → loaded SC
sc_from_file        function factory: sc_from_file('/abs/path.skl') → loaded SC
sc_bn_cubic         session  BN zincblende (SG 216), 2×2×2 prim → 16 atoms
sc_si_diamond       session  Si diamond (SG 227_a), 1×1×1 full → 8 atoms
sc_beo_hex          session  BeO wurtzite (SG 186), 1×1×1 full → 4 atoms
sc_c2_molecule      session  C₂ in a 5 Å box, P1 → 2 atoms
sc_fes2             session  FeS₂ pyrite (SG 205), 2×2×1 prim → 48 atoms

The session-scoped fixtures are *read-only*.  Tests that modify the structure
(supercell expansion, rotation, etc.) must use make_sc() or fresh_sc so that
each test starts from a clean state.
"""

import os
import sys
import pytest

# ---------------------------------------------------------------------------
# Filesystem paths
# ---------------------------------------------------------------------------
TESTS_DIR      = os.path.dirname(__file__)
SCRIPTS_DIR    = os.path.join(os.path.dirname(TESTS_DIR), 'scripts')
FIXTURES_DIR   = os.path.join(TESTS_DIR, 'fixtures')
STRUCTURES_DIR = os.path.join(FIXTURES_DIR, 'structures')
REFERENCE_DIR  = os.path.join(FIXTURES_DIR, 'reference')

# Make src/scripts importable without installing the package.
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)


# ---------------------------------------------------------------------------
# Environment check
# ---------------------------------------------------------------------------

@pytest.fixture(scope='session')
def olcao_data_dir():
    """Return $OLCAO_DATA; skip the entire session if the variable is unset."""
    data_dir = os.environ.get('OLCAO_DATA', '')
    if not data_dir:
        pytest.skip(
            'OLCAO_DATA environment variable is not set.  '
            'Export it before running the tests:\n'
            '  export OLCAO_DATA=/home/rulisp/olcao/share'
        )
    return data_dir


# ---------------------------------------------------------------------------
# StructureControl import (done once per session)
# ---------------------------------------------------------------------------

@pytest.fixture(scope='session')
def _sc_import(olcao_data_dir):
    """Import and return the StructureControl class (once per session).

    Importing triggers ElementData initialisation which reads elements.dat,
    so it is done once and cached for the whole session.
    """
    from structure_control import StructureControl
    return StructureControl


# ---------------------------------------------------------------------------
# Function-scoped SC factories (each test gets a fresh instance)
# ---------------------------------------------------------------------------

@pytest.fixture
def fresh_sc(_sc_import):
    """Return a freshly initialised StructureControl with no structure loaded.

    Use this for tests that build up state programmatically (e.g. calling
    set_lattice_from_mag_angle) rather than reading a file.
    """
    return _sc_import()


@pytest.fixture
def make_sc(_sc_import):
    """Factory: read a .skl from fixtures/structures/ and return a loaded SC.

    Usage::

        def test_something(make_sc):
            sc = make_sc('bn_cubic.skl')
            assert sc.num_atoms == 16

    Each call produces a *fresh* StructureControl instance so that mutations
    in one test never affect another.

    Uses read_input_file() (not read_olcao_skl() directly) so that
    map_element_number() and compute_implicit_info() are called automatically,
    making atomic_z and basis-set metadata available immediately.
    """
    def _factory(filename):
        sc = _sc_import()
        sc.read_input_file(os.path.join(STRUCTURES_DIR, filename))
        return sc
    return _factory


@pytest.fixture
def sc_from_file(_sc_import):
    """Factory: read a .skl from an *absolute* path and return a loaded SC.

    Intended for roundtrip tests where the output file lives in tmp_path::

        def test_roundtrip(make_sc, sc_from_file, tmp_path):
            sc1 = make_sc('bn_cubic.skl')
            out = str(tmp_path / 'out.skl')
            sc1.print_olcao(filename=out, style='frac')
            sc2 = sc_from_file(out)
            assert sc2.num_atoms == sc1.num_atoms
    """
    def _factory(path):
        sc = _sc_import()
        sc.read_input_file(path)
        return sc
    return _factory


# ---------------------------------------------------------------------------
# Session-scoped read-only structure fixtures
# ---------------------------------------------------------------------------
# These are loaded once and shared across all tests in the session.
# Do NOT modify them in place; use make_sc() if you need a mutable copy.

@pytest.fixture(scope='session')
def sc_bn_cubic(_sc_import):
    """BN zincblende, SG 216, 2×2×2 prim → 16 atoms, rhombohedral cell."""
    sc = _sc_import()
    sc.read_input_file(os.path.join(STRUCTURES_DIR, 'bn_cubic.skl'))
    return sc


@pytest.fixture(scope='session')
def sc_si_diamond(_sc_import):
    """Si diamond cubic, SG 227_a, 1×1×1 full → 8 atoms, cubic cell."""
    sc = _sc_import()
    sc.read_input_file(os.path.join(STRUCTURES_DIR, 'si_diamond.skl'))
    return sc


@pytest.fixture(scope='session')
def sc_beo_hex(_sc_import):
    """BeO wurtzite, SG 186, 1×1×1 full → 4 atoms, hexagonal cell."""
    sc = _sc_import()
    sc.read_input_file(os.path.join(STRUCTURES_DIR, 'beo_hexagonal.skl'))
    return sc


@pytest.fixture(scope='session')
def sc_c2_molecule(_sc_import):
    """C₂ dimer in a 5 Å cubic box, P1, no symmetry → 2 atoms."""
    sc = _sc_import()
    sc.read_input_file(os.path.join(STRUCTURES_DIR, 'c2_molecule.skl'))
    return sc


@pytest.fixture(scope='session')
def sc_fes2(_sc_import):
    """FeS₂ pyrite, SG 205, 2×2×1 prim → 48 atoms."""
    sc = _sc_import()
    sc.read_input_file(os.path.join(STRUCTURES_DIR, 'fes2.skl'))
    return sc
