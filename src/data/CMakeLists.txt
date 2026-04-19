# CMake file for OLCAO databases.

# Document the current location
message("At the set of databases.")

# Note: angles.dat was retired in TODO C41 (DESIGN 4.8.8 items 1
# and 6).  Angle types are now discovered at run time by
# geometry-derived clustering in condense.py (PSEUDOCODE 10c/10e),
# so the database file is no longer installed.  The file may
# remain in the repository for historical reference but is not
# read by any active code path.
set (DATABASES
atomicBDB.tgz atomicPDB.tgz
bond_parameters.dat
elements.dat
isotopes.dat
curve_colors.dat curve_marks.dat curve_styles.dat precursorDB.tgz
spaceDB.tgz sybdDB.tgz
xc_code.dat)

install (PROGRAMS ${DATABASES} DESTINATION share)
