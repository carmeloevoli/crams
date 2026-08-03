cmake_minimum_required(VERSION 3.14)

file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION VERSION)
set(LANGS CXX)
project(CRAMS VERSION ${VERSION} LANGUAGES ${LANGS})

set(CMAKE_VERBOSE_MAKEFILE off)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -march=native -Wall -Wextra -Wpedantic")
set(CRAMS_EXTRA_INCLUDES)
set(CRAMS_EXTRA_LIBRARIES)

include_directories( include/ )
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/plog/include/plog/Log.h")
  message(FATAL_ERROR "plog not found. Run: cd external && sh get_plog.sh")
endif()
include_directories( external/plog/include/ )

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_HAS_LOCAL_CHANGES)
set(git_revision_cpp "${CMAKE_CURRENT_BINARY_DIR}/git_revision.cpp")
configure_file("src/utils/git_revision.cpp.in" "${git_revision_cpp}" @ONLY)

set( CRAMS_SRCS
     src/utils/csvreader.cpp
     src/utils/utilities.cpp
     src/physics/grammage.cpp
     src/inelastic.cpp
     src/fragmentation.cpp
     src/core/input.cpp
     src/physics/losses.cpp
     src/core/output.cpp
     src/particle.cpp
     src/particlelist.cpp
     src/particlesolver.cpp
     src/physics/primary.cpp
     src/secondary.cpp
     src/runner.cpp
     "${git_revision_cpp}"
)

# GSL (external)
find_package(GSL REQUIRED)
if(GSL_FOUND)
     set(HAVE_LIBGSL 1)
     list(APPEND CRAMS_EXTRA_INCLUDES ${GSL_INCLUDE_DIRS})
     list(APPEND CRAMS_EXTRA_LIBRARIES ${GSL_LIBRARIES})
else()
     message(FATAL_ERROR "Require GSL. Set GSL_HOME")
endif()

# copy data files in build
file(GLOB GRYPHON_DATA "data/*")
file(COPY ${GRYPHON_DATA} DESTINATION ${PROJECT_BINARY_DIR}/data)

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# make library
include_directories(include ${CRAMS_EXTRA_INCLUDES})
add_library(crams_core ${CRAMS_SRCS})

# make executables
add_executable(crams apps/main.cpp)
target_link_libraries (crams crams_core ${CRAMS_EXTRA_LIBRARIES})

enable_testing()

add_executable(test_cgs tests/test_cgs.cpp)
target_link_libraries(test_cgs ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_cgs COMMAND test_cgs)

add_executable(test_utilities tests/test_utilities.cpp)
target_link_libraries(test_utilities crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_utilities COMMAND test_utilities)

add_executable(test_pid tests/test_pid.cpp)
target_link_libraries(test_pid ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_pid COMMAND test_pid)

add_executable(test_numeric tests/test_numeric.cpp)
target_link_libraries(test_numeric ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_numeric COMMAND test_numeric)

add_executable(test_csvreader tests/test_csvreader.cpp)
target_link_libraries(test_csvreader crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_csvreader COMMAND test_csvreader)

add_executable(test_grid tests/test_grid.cpp)
target_link_libraries(test_grid ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_grid COMMAND test_grid)

add_executable(test_input tests/test_input.cpp)
target_link_libraries(test_input crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_input COMMAND test_input)

add_executable(test_example_ini tests/test_example_ini.cpp)
target_compile_definitions(test_example_ini PRIVATE EXAMPLE_INI="${CMAKE_SOURCE_DIR}/examples/crams.ini")
target_link_libraries(test_example_ini crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_example_ini COMMAND test_example_ini)

add_executable(test_particlelist tests/test_particlelist.cpp)
target_link_libraries(test_particlelist crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_particlelist COMMAND test_particlelist)

add_executable(test_particle tests/test_particle.cpp)
target_link_libraries(test_particle crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_particle COMMAND test_particle)

add_executable(test_grammage tests/test_grammage.cpp)
target_link_libraries(test_grammage crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_grammage COMMAND test_grammage)

add_executable(test_primary tests/test_primary.cpp)
target_link_libraries(test_primary crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_primary COMMAND test_primary)

add_executable(test_secondary tests/test_secondary.cpp)
target_link_libraries(test_secondary crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_secondary COMMAND test_secondary)

add_executable(test_losses tests/test_losses.cpp)
target_link_libraries(test_losses crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_losses COMMAND test_losses)

add_executable(test_inelastic tests/test_inelastic.cpp)
target_link_libraries(test_inelastic crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_inelastic COMMAND test_inelastic)

add_executable(test_models tests/test_models.cpp)
target_link_libraries(test_models crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_models COMMAND test_models)

add_test(NAME test_data_tables
         COMMAND ${CMAKE_COMMAND}
                 -DCRAMS_DATA_DIR=${CMAKE_SOURCE_DIR}/data
                 -P ${CMAKE_SOURCE_DIR}/tests/check_data_tables.cmake)

add_executable(test_output tests/test_output.cpp)
target_link_libraries(test_output crams_core ${CRAMS_EXTRA_LIBRARIES})
add_test(NAME test_output COMMAND test_output)

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... =   ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "CRAMS_VERSION ......... = ${VERSION}")
message (STATUS "CRAMS_INCLUDES ........ = ${CRAMS_EXTRA_INCLUDES}")
message (STATUS "CRAMS_LIBRARIES ....... = ${CRAMS_EXTRA_LIBRARIES}")
