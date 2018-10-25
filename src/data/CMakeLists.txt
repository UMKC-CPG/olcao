# CMake file for OLCAO databases.

# Document the current location
message("At the set of databases.")

set (DATABASES atomicBDB.tgz atomicPDB.tgz spaceDB.tgz sybdDB.tgz elements.dat
   bonds.dat angles.dat xc_code.dat)

install (PROGRAMS ${DATABASES} DESTINATION share)
