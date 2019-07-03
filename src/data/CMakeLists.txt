# CMake file for OLCAO databases.

# Document the current location
message("At the set of databases.")

set (DATABASES
atomicBDB.tgz atomicPDB.tgz
precursorDB.tgz
spaceDB.tgz sybdDB.tgz
angles.dat
bonds.dat
elements.dat
xc_code.dat)

install (PROGRAMS ${DATABASES} DESTINATION share)
