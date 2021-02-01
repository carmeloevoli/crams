cmake_minimum_required(VERSION 3.2)

file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION VERSION)
set(LANGS CXX)
project(CRAMS VERSION ${VERSION} LANGUAGES ${LANGS})

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "/usr/local/Cellar/gcc/9.3.0_1/bin/g++-9" )
set( CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -std=c++11 -fopenmp -DTHREADS=4 -DEVOLI2019")
set( CRAMS_EXTRA_INCLUDES )
set( CRAMS_EXTRA_LIBRARIES )

include_directories( "${CMAKE_SOURCE_DIR}" )
include_directories( include/ )
include_directories( external/plog/include/ )

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_HAS_LOCAL_CHANGES)
set(git_revision_cpp "${CMAKE_CURRENT_BINARY_DIR}/git_revision.cpp")
configure_file("src/git_revision.cpp.in" "${git_revision_cpp}" @ONLY)

set( CRAMS_SRCS 
     src/xsecs/Evoli2019.cpp 
     src/xsecs/Tripathi99.cpp
#     src/chi2.cpp
     src/csvreader.cpp
     src/grammage.cpp
     src/inelastic.cpp
     src/input.cpp
     src/losses.cpp
     src/output.cpp
     src/particle.cpp
     src/particlelist.cpp
     src/particlesolver.cpp
     src/primary.cpp
     src/secondary.cpp
     src/utilities.cpp
     "${git_revision_cpp}"
)

# GSL (external)
find_package(GSL)
if(GSL_FOUND)
     set(HAVE_LIBGSL 1)
     list(APPEND CRAMS_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
     list(APPEND CRAMS_EXTRA_LIBRARIES ${GSL_LIBRARIES})
else(GSL_FOUND)
     message(FATAL_ERROR "Require GSL. Set GSL_HOME")
endif(GSL_FOUND)

# copy data files in build
file(GLOB GRYPHON_DATA "data/*")
file(COPY ${GRYPHON_DATA} DESTINATION ${PROJECT_BINARY_DIR}/data)

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# make library
include_directories (include ${CRAMS_EXTRA_INCLUDES})
add_library( CRAMS_LIB ${CRAMS_SRCS} ${CRAMS_EXTRA_INCLUDES} )

# make executables
add_executable(crams src/main.cpp)
target_link_libraries (crams CRAMS_LIB ${CRAMS_EXTRA_LIBRARIES})

add_executable(sandbox tests/sandbox.cpp)
target_link_libraries (sandbox CRAMS_LIB ${CRAMS_EXTRA_LIBRARIES})

add_executable(inelastic tests/inelasticPlot.cpp)
target_link_libraries (inelastic CRAMS_LIB ${CRAMS_EXTRA_LIBRARIES})

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "CRAMS_INCLUDES ........ = ${CRAMS_EXTRA_INCLUDES}")
message (STATUS "CRAMS_LIBRARIES ....... = ${CRAMS_EXTRA_LIBRARIES}")