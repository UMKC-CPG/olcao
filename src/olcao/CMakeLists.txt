# OLCAO level CMake file for OLCAO

# Document that we are at the OLCAO level.
message("At the OLCAO level.")

# List all subdirectories to be made.
add_subdirectory(real)
add_subdirectory(complex)
add_subdirectory(auxiliary)
