# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library for the SKIRT core layer
# ------------------------------------------------------------------

# set the target name
set(TARGET skirtcore)

# list the source files in this directory
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

# create the library target
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# 1. Add the CUDA library dependency
target_link_libraries(${TARGET} PRIVATE skirtcuda)
target_include_directories(${TARGET} PRIVATE ../cuda)

# 2. Add SMILE dependencies (Modernized to target_include)
target_link_libraries(${TARGET} schema fundamentals)
target_include_directories(${TARGET} PRIVATE ../../SMILE/schema ../../SMILE/fundamentals)

# 3. add SKIRT library dependencies
target_link_libraries(${TARGET} utils)
target_include_directories(${TARGET} PRIVATE ../utils)

# adjust C++ compiler flags to our needs
include("../../SMILE/build/CompilerFlags.cmake")
