# ========================================================================
# GCLF dependencies
# ========================================================================
if (WIN32)
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
      -D_CRT_SECURE_NO_WARNINGS
  )
endif ()
add_definitions(-DGCLF_DEBUG)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/GCLF-ThirdPartyLibs")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/GCLF-Utils")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/GCLF-Geometry")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/GCLF-Simplification")