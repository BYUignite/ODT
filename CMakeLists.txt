################### PREAMBLE ###############################

cmake_minimum_required(VERSION 3.15)

project(ODT VERSION 2.0.0
        LANGUAGES CXX C)

################### PROJECT-WIDE SETUP ######################

include(FetchContent)
include(ExternalProject)
include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)       # default to Release build
endif()

set(CMAKE_CXX_STANDARD 11)              # set C++ language standard
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)          # tells find_package() where to look for Find[pkg].cmake files
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")                   # installation location base directory

# ---- build options ----
option(ODT_BUILD_TESTS      "Build ODT tests"       OFF)
option(ODT_BUILD_DOCS       "Build ODT docs"        OFF)

# ---- user options ----
set(ODT_CHEMISTRY       "SIMPLEDLR"     CACHE STRING    "User chemical mechanism")
option(ODT_PARALLEL     "Build with MPI support"                OFF)
option(ODT_VERBOSE      "Displays additional runtime output"    OFF)

################### DEPENDENCIES #########################

# -------- Boost + Threads --------
find_package(Boost REQUIRED COMPONENTS system)
find_package(Threads REQUIRED)

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