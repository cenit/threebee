cmake_minimum_required (VERSION 2.8.7)
project(threebee)
enable_language(Fortran)

# make sure that the default is a RELEASE
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
  add_definitions(-DNO_OMP)
endif()

add_executable (poincare
${CMAKE_SOURCE_DIR}/src/poincare.f90
)

add_executable (threebee_legacy
${CMAKE_SOURCE_DIR}/src/three_bodies_legacy.f
)

add_executable (threebee_new
${CMAKE_SOURCE_DIR}/src/three_bodies_new.F90
)

add_executable (threebee_sho
${CMAKE_SOURCE_DIR}/src/three_bodies_sho.f90
)

install (TARGETS poincare	threebee_legacy threebee_new threebee_sho DESTINATION ${CMAKE_SOURCE_DIR}/bin/ )

