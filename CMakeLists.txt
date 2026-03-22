cmake_minimum_required(VERSION 3.10)
project(MNA_Solver VERSION 2.5.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find math library (libm)
find_library(MATH_LIBRARY m REQUIRED)

# =============================================================================
# Library Sources
# =============================================================================
set(MNA_SOURCES
    src/matrix.c
    src/solver/core.c
    src/solver/dc.c
    src/solver/ac.c
    src/solver/transient.c
    elements/passive.c
    elements/sources.c
    elements/nonlinear/nonlinear.c
    elements/nonlinear/magnetization.c
    elements/transformer.c
    elements/transformer_sat.c
)

set(MNA_HEADERS
    include/types.h
    include/matrix.h
    src/solver/core.h
    src/solver/dc.h
    src/solver/ac.h
    src/solver/transient.h
    elements/passive.h
    elements/sources.h
    elements/nonlinear/nonlinear.h
    elements/nonlinear/magnetization.h
    elements/transformer.h
    elements/transformer_sat.h
    mna_solver.h
)

# =============================================================================
# Static Library
# =============================================================================
add_library(mna_solver_static STATIC ${MNA_SOURCES})
target_include_directories(mna_solver_static PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/elements
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(mna_solver_static PUBLIC ${MATH_LIBRARY})
set_target_properties(mna_solver_static PROPERTIES
    OUTPUT_NAME mna_solver
    PREFIX lib
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# =============================================================================
# Shared Library
# =============================================================================
add_library(mna_solver_shared SHARED ${MNA_SOURCES})
target_include_directories(mna_solver_shared PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/elements
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(mna_solver_shared PUBLIC ${MATH_LIBRARY})
set_target_properties(mna_solver_shared PROPERTIES
    OUTPUT_NAME mna_solver
    VERSION ${PROJECT_VERSION}
    SOVERSION 2
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# =============================================================================
# CLI Application
# =============================================================================
set(CLI_SOURCES src/cli/netlist.c src/cli/cspice_cli.c)

add_executable(cspice ${CLI_SOURCES})
target_include_directories(cspice PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/elements
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(cspice PRIVATE mna_solver_static ${MATH_LIBRARY})
set_target_properties(cspice PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# =============================================================================
# Test Programs
# =============================================================================
add_executable(test_transformer test_transformer_csv.c)
target_include_directories(test_transformer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/elements
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(test_transformer PRIVATE mna_solver_shared ${MATH_LIBRARY})
set_target_properties(test_transformer PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# =============================================================================
# Installation
# =============================================================================
install(TARGETS mna_solver_static mna_solver_shared cspice
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(FILES ${MNA_HEADERS} DESTINATION include/mna_solver)
install(FILES mna_solver.h DESTINATION include)

# =============================================================================
# Build Summary
# =============================================================================
message(STATUS "")
message(STATUS "=== MNA Solver Configuration ===")
message(STATUS "  Version:      ${PROJECT_VERSION}")
message(STATUS "  Build Type:   ${CMAKE_BUILD_TYPE}")
message(STATUS "  C Standard:   C${CMAKE_C_STANDARD}")
message(STATUS "  Math Lib:     ${MATH_LIBRARY}")
message(STATUS "")
