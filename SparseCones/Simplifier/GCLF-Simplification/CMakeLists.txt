cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0074 NEW)
# set the project name
project(GCLF-Simplification)

# search all source files in current directory
file(GLOB SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/Simplifier/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/Simplifier/Geom/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/Simplifier/SimplifyStages/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/Simplifier/SpaceSearch/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/Simplifier/Topo/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/Mesh/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/Utils/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/DE/*"
  "${CMAKE_CURRENT_SOURCE_DIR}/interface.h"
)

if (GCLF_CMAKE_DEBUG)
  message(STATUS "GCLF-Simplification source files: ${SRCS}")
endif()

# OpenMP
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    find_package(OpenMP)
    if(NOT OPENMP_FOUND)
        message(FATAL_ERROR "OpenMP Not Found!")
    endif()
endif()

# CGAL
set(CGAL_Boost_USE_STATIC_LIBS ON)
find_package(CGAL)
if(CGAL_FOUND)
    MESSAGE(STATUS "CGAL INCLUDE DIRS:${CGAL_INCLUDE_DIRS}")
else()
    MESSAGE(FATAL "Can't find CGAL.")
endif()

# Boost
# set(Boost_DEBUG ON)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost
    REQUIRED
    COMPONENTS filesystem json
)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    MESSAGE(STATUS "-------------------------------------------")
    MESSAGE(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}.")
    MESSAGE(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}.")
    MESSAGE(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}.")
    MESSAGE(STATUS "Boost_LIB_VERSION = ${Boost_LIB_VERSION}.")
    MESSAGE(STATUS "-------------------------------------------")
endif()

# add the executable
add_library(${PROJECT_NAME} STATIC
  ${SRCS}
)

# Public include directory
target_include_directories(${PROJECT_NAME} PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${OPENMESH_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC
  GCLF-Utils
  GCLF-Geometry
  ${OPENMESH_LIBRARIES}
  ${Boost_LIBRARIES}
)

target_compile_definitions(${PROJECT_NAME} PUBLIC
  #USE_TREE_SEARCH
)

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