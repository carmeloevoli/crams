project(CRAMS CXX)

cmake_minimum_required(VERSION 2.6)

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_CXX_FLAGS "-O3 -std=c++11" )
set( CRAMS_EXTRA_INCLUDES )
set( CRAMS_EXTRA_LIBRARIES )

include_directories( "${CMAKE_SOURCE_DIR}" )
include_directories( include/ )

set( CRAMS_SRCS 
     src/main.cpp
     src/output.cpp
     src/utilities.cpp
     src/particle.cpp
	 src/losses.cpp
	 src/inelastic.cpp
	 src/grammage.cpp
	 src/params.cpp
	 src/primary.cpp
	 src/secondary.cpp
	 src/spallation.cpp
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

add_executable(CRAMS ${CRAMS_SRCS} )
target_link_libraries (CRAMS ${CRAMS_EXTRA_LIBRARIES})

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
