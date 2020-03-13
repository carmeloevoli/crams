cmake_minimum_required(VERSION 3.2)

file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION VERSION)
set(LANGS CXX)
project(CRAMS VERSION ${VERSION} LANGUAGES ${LANGS})

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "/usr/local/Cellar/gcc/9.2.0_1/bin/g++-9" )
set( CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -std=c++11 -fopenmp -DTHREADS=4 -DDEBUG -DGALPROP2017")
set( CRAMS_EXTRA_INCLUDES )
set( CRAMS_EXTRA_LIBRARIES )

include_directories( "${CMAKE_SOURCE_DIR}" )
include_directories( include/ )

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_HAS_LOCAL_CHANGES)
set(git_revision_cpp "${CMAKE_CURRENT_BINARY_DIR}/git_revision.cpp")
configure_file("src/git_revision.cpp.in" "${git_revision_cpp}" @ONLY)

set( CRAMS_SRCS 
     src/chi2.cpp
	 src/grammage.cpp
	 src/inelastic.cpp
	 src/losses.cpp
     src/output.cpp
	 src/params.cpp
     src/particle.cpp
	 src/primary.cpp
	 src/secondary.cpp
	 src/spallation.cpp
	 src/Tripathi99.cpp
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

# MAKE EXECUTABLE
include_directories (include ${CRAMS_EXTRA_INCLUDES})

add_executable(CRAMS src/main.cpp ${CRAMS_SRCS} )
target_link_libraries (CRAMS ${CRAMS_EXTRA_LIBRARIES})

add_executable(sandbox src/sandbox.cpp ${CRAMS_SRCS} )
target_link_libraries (sandbox ${CRAMS_EXTRA_LIBRARIES})

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
