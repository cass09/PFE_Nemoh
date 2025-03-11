cmake_minimum_required(VERSION 3.0.2)

project(Nemoh Fortran)

message("EXECUTABLES_OUTPUT_DIR: " ${EXECUTABLES_OUTPUT_DIR})
if(DEFINED EXECUTABLES_OUTPUT_DIR)
    if(IS_ABSOLUTE ${EXECUTABLES_OUTPUT_DIR})
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLES_OUTPUT_DIR})
    else()
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/${EXECUTABLES_OUTPUT_DIR})
    endif()
else()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
endif()
message("CMAKE_RUNTIME_OUTPUT_DIRECTORY: " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if(NOT DEFINED CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "Release")
endif()

# TODO: remove the following block, which breaks normal use of compiler options
# Files that need pre-processing should have a capitalized extension (e.g. .F90 in place of .f90)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "-cpp")
    # Other flags that may be useful: -r8 (force double-precision on all real numbers)
    add_definitions(-D_INTEL_FORTRAN_COMPILE_RULE_)
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-cpp -DGNUFORT -ffree-line-length-none -static-libgfortran")
    # Other flags that may be useful: -fdefault-real-8 (force double-precision on all real numbers)
    add_definitions(-D_GNU_FORTRAN_COMPILE_RULE_)
endif()

# This is required to get _WIN32 defined because the previous block breaks the default behaviour
if(WIN32)
    add_compile_definitions(_WIN32)
endif()

set(COMMON_SRC_DIR ${PROJECT_SOURCE_DIR}/Common)
set(CORE_SRC_DIR ${PROJECT_SOURCE_DIR}/Solver/Core)
set(EXTERNAL_SRC_DIR ${PROJECT_SOURCE_DIR}/external)

add_subdirectory(Mesh)
add_subdirectory(preProcessor)
add_subdirectory(postProcessor)
add_subdirectory(Solver)

add_subdirectory(QTF)

include(CTest)

add_subdirectory(TestCases)
