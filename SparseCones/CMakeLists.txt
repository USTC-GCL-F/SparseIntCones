cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0074 NEW)
# set the project name
project(SparseCones)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Simplifier)

# search all source files in current directory
file(GLOB SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/Mesh/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/ShortestPath/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

if (GCLF_CMAKE_DEBUG)
  message(STATUS "SparseCones source files: ${SRCS}")
endif()

# OpenMP
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    find_package(OpenMP)
    if(NOT OPENMP_FOUND)
        message(FATAL_ERROR "OpenMP Not Found!")
    endif()
endif()

# add the library
add_library(${PROJECT_NAME} STATIC
  ${SRCS}
)

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${OPENMESH_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  "C:/Program Files/Mosek/10.0/tools/platform/win64x86/h"
)

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC
  GCLF-Simplification
  ${OPENMESH_LIBRARIES}
  "C:/Program Files/Mosek/10.0/tools/platform/win64x86/bin/mosek64_10_0.lib"
  "C:/Program Files/Mosek/10.0/tools/platform/win64x86/bin/fusion64_10_0.lib"
)

target_compile_definitions(${PROJECT_NAME} PUBLIC
)


if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  target_compile_definitions(${PROJECT_NAME} PUBLIC
    _USE_MATH_DEFINES
  )
endif()

# Compiler-specific options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # -frounding-math grant IEEE 754 compliance
  target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)

  # -g generate debug information
  target_compile_options(${PROJECT_NAME} PUBLIC -g)

# use intrinsic functions (CHECK WHAT TO DO FOR GCC !!!!!!!!)
# target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  # grant IEEE 754 compliance
  target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")

  # generate debug information
  target_compile_options(${PROJECT_NAME} PUBLIC "/Zi")

  # use intrinsic functions
  target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")

  # turn off annoying warnings
  target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
endif()