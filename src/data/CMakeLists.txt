# CMake file for OLCAO databases.

# Document the current location
message("At the set of databases.")

set (DATABASES
angles.dat atomicBDB.tgz atomicPDB.tgz
bonds.dat
elements.dat
isotopes.dat
precursorDB.tgz
spaceDB.tgz sybdDB.tgz
xc_code.dat)

install (PROGRAMS ${DATABASES} DESTINATION share)
