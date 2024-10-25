cmake_minimum_required(VERSION 3.15)

project(odt LANGUAGES CXX C)

################### PROJECT SETUP ##########################

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)      # Release or Debug
endif()

# ---- set C++ language standard ----

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

#################### build options

option(BUILD_DOCS "Build ODT documentation" OFF)

set(CHEMISTRY               "SIMPLEDLR"             CACHE STRING "user mechanism source to add")
set(PARALLEL                "no"                    CACHE STRING "build with mpi")
set(SILENT                  "yes"                   CACHE STRING "code output config")

#################### preprocessor directives

if(PARALLEL STREQUAL "yes")
    message("MPI is on")
    add_definitions(-DDOMPI)
endif()

if(NOT CHEMISTRY STREQUAL "CANTERARR" )
    add_definitions(-DPROBLEMSPECIFICRR)
endif()

if(SILENT STREQUAL "yes" )
    message("Suppressing some runtime output")
    add_definitions(-DSILENT)
endif()

#################### project targets

add_subdirectory(src)

################### OPTIONAL PROJECT TARGETS #################

#-------- DOCUMENTATION --------

if(BUILD_DOCS)
    add_subdirectory(docs)
endif()
