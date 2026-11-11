cmake_minimum_required(VERSION 3.15)
project(aldmc VERSION 1.0.0 LANGUAGES C CXX)

# aldmc: 3D voxel test-particle Monte Carlo feature solver with a C ABI (see README.md).
#   cmake -S lib/aldmc -B lib/aldmc/build -DCMAKE_BUILD_TYPE=Release
#   cmake --build lib/aldmc/build
#   ctest --test-dir lib/aldmc/build --output-on-failure
# Targets: aldmc (shared: libaldmc.dylib / .so / aldmc.dll) and aldmc_static (libaldmc.a, position independent).
# Parallelism is std::thread (no OpenMP), so Apple clang builds it as is.

option(ALDMC_BUILD_TESTS "Build the ctest suite" ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 99)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

# No fused multiply-add contraction, so a given seed gives the same trajectories on every platform the compiler
# supports (the random streams are the library's own xoshiro256**, also platform independent).
set(ALDMC_FP_FLAGS "")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  set(ALDMC_FP_FLAGS -ffp-contract=off)
elseif(MSVC)
  set(ALDMC_FP_FLAGS /fp:precise)
endif()

find_package(Threads REQUIRED)

set(ALDMC_SOURCES src/aldmc.cpp)
set(ALDMC_HEADER include/aldmc.h)

add_library(aldmc SHARED ${ALDMC_SOURCES} ${ALDMC_HEADER})
add_library(aldmc_static STATIC ${ALDMC_SOURCES} ${ALDMC_HEADER})

foreach(tgt aldmc aldmc_static)
  target_include_directories(${tgt} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
  target_compile_options(${tgt} PRIVATE ${ALDMC_FP_FLAGS})
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(${tgt} PRIVATE -Wall -Wextra)
  endif()
  target_link_libraries(${tgt} PUBLIC Threads::Threads)
  set_target_properties(${tgt} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER ${ALDMC_HEADER})
endforeach()

target_compile_definitions(aldmc PRIVATE ALDMC_BUILD_DLL)
set_target_properties(aldmc PROPERTIES OUTPUT_NAME aldmc)
if(WIN32)
  set_target_properties(aldmc_static PROPERTIES OUTPUT_NAME aldmc_static)
else()
  set_target_properties(aldmc_static PROPERTIES OUTPUT_NAME aldmc)
endif()

include(GNUInstallDirs)
install(TARGETS aldmc aldmc_static
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ALDMC_BUILD_TESTS)
  enable_testing()
  add_executable(test_aldmc tests/test_aldmc.cpp)
  target_link_libraries(test_aldmc PRIVATE aldmc_static)
  target_compile_options(test_aldmc PRIVATE ${ALDMC_FP_FLAGS})
  foreach(name
      version_and_defaults
      open_box_loses_all
      flat_langmuir_analytic
      flat_langmuir_deposit_scheme
      periodic_slit_conserves
      recombination_loses_all
      determinism_and_thread_invariance
      area_weights_cylinder
      hole_saturates
      error_codes)
    add_test(NAME aldmc.${name} COMMAND test_aldmc ${name})
  endforeach()
endif()

# the tests call exp/sqrt/floor: on glibc those live in libm and must be linked explicitly
find_library(ALDMC_LIBM m)
if(ALDMC_LIBM)
  if(TARGET test_aldmc)
    target_link_libraries(test_aldmc PRIVATE ${ALDMC_LIBM})
  endif()
  target_link_libraries(aldmc PRIVATE ${ALDMC_LIBM})
endif()
