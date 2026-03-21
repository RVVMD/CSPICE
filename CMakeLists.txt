cmake_minimum_required(VERSION 3.10)
project(MNA_Solver VERSION 2.5.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_library(MATH_LIBRARY m)

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
)

set(CLI_SOURCES
    src/cli/netlist.c
    src/cli/cspice_cli.c
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
    mna_solver.h
    src/cli/netlist.h
)

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

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/TR5_sim.c")
    add_executable(tr5_sim TR5_sim.c)
    target_include_directories(tr5_sim PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/elements
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    )
    target_link_libraries(tr5_sim PRIVATE mna_solver_static ${MATH_LIBRARY})
    set_target_properties(tr5_sim PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
endif()

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

install(TARGETS mna_solver_static mna_solver_shared
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

install(FILES ${MNA_HEADERS}
    DESTINATION include/mna_solver
)

message(STATUS "MNA Solver Version: ${PROJECT_VERSION}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Math Library: ${MATH_LIBRARY}")
