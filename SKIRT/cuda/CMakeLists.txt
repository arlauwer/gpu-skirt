# //////////////////////////////////////////////////////////////////
# ///     The SKIRT project -- advanced radiative transfer       ///
# ///       © Astronomical Observatory, Ghent University         ///
# //////////////////////////////////////////////////////////////////

# ------------------------------------------------------------------
# Builds a library for the cuda kernels
# ------------------------------------------------------------------

# set the target name
set(CUDA_TARGET skirtcuda)

# list CUDA source files
file(GLOB CUDA_SOURCES "*.cu")
file(GLOB CUDA_HEADERS "*.cuh")

# Create a static library for CUDA kernels
add_library(${CUDA_TARGET} STATIC ${CUDA_SOURCES} ${CUDA_HEADERS})

# Link against CUDA runtime
find_package(CUDAToolkit REQUIRED)
target_link_libraries(${CUDA_TARGET} PUBLIC CUDA::cudart)

# Ensure the core can find these headers
target_include_directories(${CUDA_TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
