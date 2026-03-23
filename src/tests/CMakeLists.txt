# CMake file for OLCAO test suite.

# Document the current location
message("At the test suite.")

set (TEST_SCRIPTS
conftest.py
test_dirac_helpers.py
test_element_data.py
test_interfit_output.py
test_o.py
test_org_radwavefn.py
test_structure_control_extras.py
test_structure_control_io.py
test_structure_control_ops.py
test_structure_control_unit.py)

set (TEST_FIXTURES
fixtures/structures/beo_hexagonal.skl
fixtures/structures/bn_cubic.skl
fixtures/structures/c2_molecule.skl
fixtures/structures/fes2.skl
fixtures/structures/si_diamond.skl)

install (PROGRAMS ${TEST_SCRIPTS} DESTINATION tests)
install (FILES ${TEST_FIXTURES} DESTINATION tests/fixtures/structures)
