cmake_minimum_required(VERSION 3.15)
project(aldkin VERSION 1.0.0 LANGUAGES C CXX)

# aldkin: multiscale ALD wall kinetics with a C ABI (see README.md).
#   cmake -S lib/aldkin -B lib/aldkin/build -DCMAKE_BUILD_TYPE=Release
#   cmake --build lib/aldkin/build
#   ctest --test-dir lib/aldkin/build --output-on-failure
# Targets: aldkin (shared: libaldkin.dylib / .so / aldkin.dll) and aldkin_static (libaldkin.a, position independent,
# the one the OpenFOAM adapter links). Options below switch the tests and the example off.

option(ALDKIN_BUILD_TESTS "Build the ctest suite" ON)
option(ALDKIN_BUILD_EXAMPLES "Build examples/plain_c_example.c" ON)

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

# Same floating-point contract as OpenFOAM's wmake rules (no fused multiply-add contraction), so the library's
# numbers are the boundary conditions' numbers bit for bit on every platform where wmake builds.
set(ALDKIN_FP_FLAGS "")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  set(ALDKIN_FP_FLAGS -ffp-contract=off)
elseif(MSVC)
  set(ALDKIN_FP_FLAGS /fp:precise)
endif()

set(ALDKIN_SOURCES src/aldkin.cpp)
set(ALDKIN_HEADER include/aldkin.h)

add_library(aldkin SHARED ${ALDKIN_SOURCES} ${ALDKIN_HEADER})
add_library(aldkin_static STATIC ${ALDKIN_SOURCES} ${ALDKIN_HEADER})

foreach(tgt aldkin aldkin_static)
  target_include_directories(${tgt} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
  target_compile_options(${tgt} PRIVATE ${ALDKIN_FP_FLAGS})
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(${tgt} PRIVATE -Wall -Wextra)
  endif()
  set_target_properties(${tgt} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER ${ALDKIN_HEADER})
endforeach()

target_compile_definitions(aldkin PRIVATE ALDKIN_BUILD_DLL)
set_target_properties(aldkin PROPERTIES OUTPUT_NAME aldkin)
if(WIN32)
  set_target_properties(aldkin_static PROPERTIES OUTPUT_NAME aldkin_static)
else()
  set_target_properties(aldkin_static PROPERTIES OUTPUT_NAME aldkin)
endif()

include(GNUInstallDirs)
install(TARGETS aldkin aldkin_static
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ALDKIN_BUILD_TESTS)
  enable_testing()
  add_executable(test_aldkin tests/test_aldkin.cpp)
  target_link_libraries(test_aldkin PRIVATE aldkin_static)
  target_compile_options(test_aldkin PRIVATE ${ALDKIN_FP_FLAGS})
  foreach(name
      version_and_defaults
      flat_langmuir_analytic
      flat_langmuir_euler_convergence
      flat_langmuir_desorption
      hole_constant_density_t90
      hole_euler_literal_step
      profiles_interpolate
      zero_density_and_zero_holes_are_safe
      clone_and_restart
      error_codes)
    add_test(NAME aldkin.${name} COMMAND test_aldkin ${name})
  endforeach()
endif()

if(ALDKIN_BUILD_EXAMPLES)
  add_executable(plain_c_example examples/plain_c_example.c)
  set_target_properties(plain_c_example PROPERTIES LINKER_LANGUAGE C)
  target_link_libraries(plain_c_example PRIVATE aldkin)
  if(ALDKIN_BUILD_TESTS)
    # the example is also a test: C linkage against the shared library and a run through the long case's drive
    add_test(NAME aldkin.plain_c_example
      COMMAND plain_c_example ${CMAKE_CURRENT_SOURCE_DIR}/examples/centre_face_drive.csv)
  endif()
endif()

# the C example and the tests call floor/exp/sqrt: on glibc those live in libm and must be linked explicitly
find_library(ALDKIN_LIBM m)
if(ALDKIN_LIBM)
  foreach(t IN ITEMS plain_c_example test_aldkin)
    if(TARGET ${t})
      target_link_libraries(${t} PRIVATE ${ALDKIN_LIBM})
    endif()
  endforeach()
  target_link_libraries(aldkin PRIVATE ${ALDKIN_LIBM})
endif()
