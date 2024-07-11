# CMake file for OLCAO databases.

# Document the current location
message("At the set of databases.")

set (DATABASES
angles.dat atomicBDB.tgz atomicPDB.tgz
bonds.dat
elements.dat
isotopes.dat
curve_colors.dat curve_marks.dat curve_styles.dat precursorDB.tgz
spaceDB.tgz sybdDB.tgz
xc_code.dat)

install (PROGRAMS ${DATABASES} DESTINATION share)
