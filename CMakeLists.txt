cmake_minimum_required(VERSION 3.0.0)

# project name and it is C++ only
project(scATACseq CXX)


# libraries
## boost library
set( BOOST_ROOT "" CACHE PATH "Boost library path" )
set(Boost_NO_SYSTEM_PATHS on CACHE BOOL "Do not search system for Boost")
set(BOOST_INCLUDEDIR "<path_to_project_root>/lib/boost/include/")
set(BOOST_LIBRARYDIR "<path_to_project_root>/lib/")
find_package(Boost 1.70 COMPONENTS program_options REQUIRED)

## UnitTest++ library
find_library(UNITTEST_LIB
             NAMES "UnitTest++"
             PATHS "<path_to_project_root>/lib/UnitTest++/lib")
find_path(UNITTEST_INCLUDE
          NAMES "UnitTest++.h"
          PATHS "<path_to_project_root>/lib/UnitTest++/include")
include_directories(${UNITTEST_INCLUDE})

## zlib (for seqan Bam I/O)
find_package(ZLIB REQUIRED)

## SeqAn
find_package (SeqAn REQUIRED)

## threads
find_package(Threads REQUIRED)

add_subdirectory(src)

