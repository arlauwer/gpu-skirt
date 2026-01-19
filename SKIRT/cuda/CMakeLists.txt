# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library for the cuda kernels
# ------------------------------------------------------------------

# set the target name
set(TARGET skirtcuda)

# list CUDA source files
file(GLOB SOURCES "*.cu")
file(GLOB HEADERS "*.cuh")

# Create a static library for CUDA kernels
add_library(${TARGET} STATIC ${SOURCES} ${HEADERS})

# Link against CUDA runtime
find_package(CUDAToolkit REQUIRED)
target_link_libraries(${TARGET} PUBLIC CUDA::cudart)

# Ensure the core can find these headers
include("../../SMILE/build/CompilerFlags.cmake")