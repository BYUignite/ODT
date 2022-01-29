################### PREAMBLE ###############################

cmake_minimum_required(VERSION 3.15)

project(ODT VERSION 2.0.0
        LANGUAGES CXX C)

################### PROJECT-WIDE SETUP ######################

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)       # default to Release build
endif()

set(CMAKE_CXX_STANDARD 11)              # set C++ language standard
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)          # tells find_package() where to look for Find[pkg].cmake files
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")                   # installation location base directory
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/run")     # runtime output location

# ---- build options ----
option(ODT_BUILD_TESTS      "Build ODT tests"       OFF)
option(ODT_BUILD_DOCS       "Build ODT docs"        OFF)

# ---- user options ----
set(ODT_CHEMISTRY   "SIMPLEDLR"     CACHE STRING    "User chemical mechanism")
set(ODT_PARALLEL    "OFF"           CACHE STRING    "Build with MPI")
set(ODT_SILENT      "OFF"           CACHE STRING    "Suppresses some runtime output")

################### DEPENDENCIES #########################

add_subdirectory(ext)

################### MAIN PROJECT TARGETS #####################

add_subdirectory(source)

################### OPTIONAL PROJECT TARGETS #####################

# ---- testing ----
if (ODT_BUILD_TESTS)
    add_subdirectory(tests)
endif()

# ---- documentation ----
if (ODT_BUILD_DOCS)
    add_subdirectory(docs)
endif()

#################### PREPROCESSOR DIRECTIVES ####################

if(ODT_PARALLEL STREQUAL "yes")
    message("MPI is on")
    add_definitions(-DDOMPI)
endif()

if(NOT ODT_CHEMISTRY STREQUAL "CANTERARR" )
    add_definitions(-DPROBLEMSPECIFICRR)
endif()

if(ODT_SILENT STREQUAL "yes" )
    message("ODT will suppress some runtime output")
    add_definitions(-DSILENT)
endif()

add_definitions(-D_USESTDVECTOR_)
add_definitions(-DDOCANTERA)