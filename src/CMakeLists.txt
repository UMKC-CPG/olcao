# Source level CMake file for OLCAO

# Document that we are at the source level.
message("At the source level.")

# List all subdirectories to be made.
add_subdirectory(applySpaceGroup)
add_subdirectory(atomSCF)
add_subdirectory(contract)
add_subdirectory(gaussFit)
add_subdirectory(makeKPoints)
add_subdirectory(scripts)
add_subdirectory(data)
add_subdirectory(olcao)
