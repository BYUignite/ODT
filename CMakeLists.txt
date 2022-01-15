################### PREAMBLE ###############################

cmake_minimum_required(VERSION 3.15)

project(ODT VERSION 2.0.0
        LANGUAGES CXX C)

################### PROJECT-WIDE SETUP #######################

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)       # default to Release build
endif()

set(CMAKE_CXX_STANDARD 11)              # set C++ language standard
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)          # tells find_package() where to look for Find[pkg].cmake files
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")                   # installation location base directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/run")     # runtime output location

# ---- build options ----
option(BUILD_EXAMPLES   "Build ODT examples"    OFF)
option(BUILD_TESTS      "Build ODT tests"       OFF)
option(BUILD_DOCS       "Build ODT docs"        OFF)

# ---- user options ----
set(CHEMISTRY   "SIMPLEDLR"     CACHE STRING    "User chemical mechanism")
set(PARALLEL    "OFF"           CACHE STRING    "Build with MPI")
set(SILENT      "OFF"           CACHE STRING    "Suppresses some runtime output")

################### DEPENDENCIES #########################

include(FetchContent)

# ---- Dependency options ----
option(ODT_SYSTEM_YAML      "Use system-wide yaml installation"         ON)
option(ODT_SYSTEM_SUNDIALS  "Use system-wide Sundials installation"     ON)
option(ODT_SYSTEM_CANTERA   "Use system-wide Cantera installation"      ON)

# ---- DEBUGGING ONLY ----
#set(FETCHCONTENT_SOURCE_DIR_RADLIB  "/home/vbstephens/codes/radlib"  CACHE STRING "local radlib development location")
#set(FETCHCONTENT_SOURCE_DIR_SOOTLIB "/home/vbstephens/codes/sootlib" CACHE STRING "local sootlib development location")

# -------- Boost + Threads --------
find_package(Boost REQUIRED COMPONENTS system)
find_package(Threads REQUIRED)

# ---- FetchContent declarations ----
FetchContent_Declare(yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG 9a3624205e8774953ef18f57067b3426c1c5ada6)           # yaml-cpp-0.6.3
FetchContent_Declare(Sundials
        GIT_REPOSITORY https://github.com/LLNL/sundials.git
        GIT_TAG 887af4374af2271db9310d31eaa9b5aeff49e829)           # Sundials-5.3.0
FetchContent_Declare(radlib
        GIT_REPOSITORY https://github.com/BYUignite/radlib.git
        GIT_TAG da2420bc46a94e0c96f487735fda949bca96d07c)           # radlib-2.0.0
FetchContent_Declare(sootlib
        GIT_REPOSITORY https://github.com/BYUignite/sootlib.git
        GIT_TAG 5f5eeb6fdb1ef66e7d82fd136026bb20f06208e5)           # sootlib-0.1.0

# -------- YAML --------
message(STATUS "Looking for yaml-cpp")

if (ODT_SYSTEM_YAML)
    find_package(yaml-cpp)
    if (yaml-cpp_FOUND)
        message(STATUS "Found---using system installation of yaml-cpp")
    endif ()
endif ()

if (NOT ODT_SYSTEM_YAML OR NOT yaml-cpp_FOUND)
    FetchContent_MakeAvailable(yaml-cpp)
    if(yaml-cpp_FOUND)
        message(STATUS "Found---using yaml-cpp fetched from repository")
    else()
        message(FATAL_ERROR "yaml-cpp not found")
    endif()
endif ()

# -------- SUNDIALS --------
message(STATUS "Looking for Sundials")

if (ODT_SYSTEM_SUNDIALS)
    set(SUNDIALS_DIR "/usr/local/lib")
    find_package(SUNDIALS REQUIRED COMPONENTS cvode nvecserial sunmatrixdense)      # throws an error if required components not found
    if(SUNDIALS_FOUND)
        message(STATUS "Found---using system installation of Sundials")
    endif()
endif ()

if (NOT ODT_SYSTEM_SUNDIALS OR NOT SUNDIALS_FOUND)
    FetchContent_GetProperties(SUNDIALS)
    if (NOT SUNDIALS_POPULATED)
        FetchContent_Populate(SUNDIALS)                                     # populate content from repository
        set(BUILD_ARKODE             OFF     CACHE INTERNAL "")             # set options to only build the parts we need
        set(BUILD_CVODE              ON      CACHE INTERNAL "")
        set(BUILD_CVODES             OFF     CACHE INTERNAL "")
        set(BUILD_IDA                OFF     CACHE INTERNAL "")
        set(BUILD_IDAS               OFF     CACHE INTERNAL "")
        set(BUILD_KINSOL             OFF     CACHE INTERNAL "")
        set(BUILD_NVECTOR_SERIAL     ON      CACHE INTERNAL "")
        set(BUILD_NVECTOR_MANYVECTOR OFF     CACHE INTERNAL "")
        set(EXAMPLES_ENABLE          OFF     CACHE INTERNAL "")
        set(ENABLE_LAPACK            ON      CACHE INTERNAL "")             # required lapack parts will be installed automatically
        set(CMAKE_C_FLAGS            "-fPIC" CACHE INTERNAL "")             # required flag for building shared library
        add_subdirectory(${SUNDIALS_SOURCE_DIR} ${SUNDIALS_BINARY_DIR})     # Bring the populated content into the build
    endif ()
    if (SUNDIALS_FOUND)
        message(STATUS "Found---using Sundials fetched from repository")
    else ()
        message(FATAL_ERROR "Sundials not found")
    endif ()
endif ()

# -------- CANTERA --------
message(STATUS "Looking for Cantera")

if (ODT_SYSTEM_CANTERA)
    find_package(Cantera REQUIRED)      # throws an error if required components not found
    if (Cantera_FOUND)
        message(STATUS "Found---using system installation of Cantera")
    endif ()
endif ()

if (NOT ODT_SYSTEM_CANTERA OR NOT Cantera_FOUND)
    # TODO add Cantera via ExternalProject_Add commands here
    if (Cantera_FOUND)
        message(STATUS "Found---using Cantera fetched from repository")
    else ()
        message(FATAL_ERROR "Cantera not found")
    endif ()
endif ()

# -------- RADLIB --------
FetchContent_MakeAvailable(radlib)

# -------- SOOTLIB --------
FetchContent_MakeAvailable(sootlib)

################### MAIN PROJECT TARGETS #####################

add_subdirectory(source)

################### OPTIONAL PROJECT TARGETS #####################

# ---- examples ----
if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# ---- testing ----
if (BUILD_TESTS)
    include(Ctest)
    add_subdirectory(tests)
endif()

# ---- documentation ----
if (BUILD_DOCS)
    add_subdirectory(docs)
endif()

#################### PREPROCESSOR DIRECTIVES ####################

if(PARALLEL STREQUAL "yes")
    message("MPI is on")
    add_definitions(-DDOMPI)
endif()

if(NOT CHEMISTRY STREQUAL "CANTERARR" )
    add_definitions(-DPROBLEMSPECIFICRR)
endif()

if(SILENT STREQUAL "yes" )
    message("ODT will suppress some runtime output")
    add_definitions(-DSILENT)
endif()

add_definitions(-D_USESTDVECTOR_)
add_definitions(-DDOCANTERA)