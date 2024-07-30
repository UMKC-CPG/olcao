# uOLCAO level CMake file for unified OLCAO

# Document that we are at the uOLCAO level.
message("At the uOLCAO level.")

# List all subdirectories to be made.
add_subdirectory(real)
add_subdirectory(complex)
add_subdirectory(auxiliary)
