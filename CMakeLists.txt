################### CMAKE PROJECT ODT #########################

cmake_minimum_required(VERSION 3.15)

project(ODT
        VERSION 2.0
        LANGUAGES CXX C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)      # Release or Debug
endif()

################### BUILD & RUN DIRECTORIES #########################

set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/build")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/run")

###################### PROJECT OPTIONS & FLAGS #########################

set(CHEMISTRY               "SIMPLEDLR"              CACHE STRING "user chemical mechanism to add")
set(PARALLEL                "OFF"                    CACHE STRING "build with mpi")
set(SILENT                  "OFF"                    CACHE STRING "suppresses some runtime output")

# ---- options to build examples/tests/documentation ----
option(BUILD_EXAMPLES   "Build ODT examples"    OFF)
option(BUILD_TESTS      "Build ODT tests"       OFF)
option(BUILD_DOCS       "Build ODT docs"        OFF)

# ---- C++ standard & compiler flags -------------------
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} -ggdb3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")

#################### TESTS ###############################

if (BUILD_TESTS)
    include(Ctest)
    add_subdirectory(test)
endif()

#################### EXAMPLES ###############################

if (BUILD_EXAMPLES)
#    add_subdirectory(examples)
endif()

################### PACKAGES & LIBRARIES #########################

include(FetchContent)

find_package(Boost REQUIRED COMPONENTS system)
find_package(Threads REQUIRED)

#--------------------- YAML ------------------------
message(STATUS "Looking for yaml-cpp")
FetchContent_Declare(yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG yaml-cpp-0.6.3)
FetchContent_GetProperties(yaml-cpp)
if(NOT yaml-cpp_POPULATED)
    FetchContent_MakeAvailable(yaml-cpp)
endif()

#--------------------- SUNDIALS ------------------------
message(STATUS "Looking for sundials")
FetchContent_Declare(sundials
        GIT_REPOSITORY https://github.com/LLNL/sundials.git
        GIT_TAG v5.7.0)
FetchContent_GetProperties(sundials)
if(NOT sundials_POPULATED)
    FetchContent_Populate(sundials)
    set(BUILD_ARKODE             OFF     CACHE INTERNAL "")
    set(BUILD_CVODE              OFF     CACHE INTERNAL "")
    set(BUILD_CVODES             ON      CACHE INTERNAL "")
    set(BUILD_IDA                OFF     CACHE INTERNAL "")
    set(BUILD_IDAS               OFF     CACHE INTERNAL "")
    set(BUILD_KINSOL             OFF     CACHE INTERNAL "")
    set(BUILD_NVECTOR_SERIAL     ON      CACHE INTERNAL "")
    set(BUILD_NVECTOR_MANYVECTOR ON      CACHE INTERNAL "")
    set(EXAMPLES_ENABLE          OFF     CACHE INTERNAL "")
    set(EXAMPLES_INSTALL         OFF     CACHE INTERNAL "")
    set(EXAMPLES_ENABLE_CXX      OFF     CACHE INTERNAL "")
    set(EXAMPLES_ENABLE_C        OFF     CACHE INTERNAL "")
    set(EXAMPLES_ENABLE_F77      OFF     CACHE INTERNAL "")
    set(EXAMPLES_ENABLE_F90      OFF     CACHE INTERNAL "")
    set(EXAMPLES_ENABLE_F2003    OFF     CACHE INTERNAL "")
    set(ENABLE_LAPACK            ON      CACHE INTERNAL "")

    # Bring the populated content into the build
    add_subdirectory(${sundials_SOURCE_DIR} ${sundials_BINARY_DIR})
endif()

#--------------------- CANTERA ------------------------
message(STATUS "Looking for Cantera")
find_package(Cantera REQUIRED)

############################### DOCUMENTATION ###########################

if (BUILD_DOCS)
    find_package(Doxygen)
    if(Doxygen_FOUND)
        add_subdirectory(docs ${CMAKE_CURRENT_BINARY_DIR}/docs)
    else()
        message(STATUS "Doxygen not found---documentation will not be built")
    endif()
endif()

#################### ODT SOURCE #######################

add_subdirectory(source ${CMAKE_CURRENT_BINARY_DIR}/source)

#################### PREPROCESSOR DIRECTIVES ####################
# TODO do we need this with CMake? Is there a better way to do it?

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

add_definitions(-D_USESTDVECTOR_)