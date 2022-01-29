################### DEPENDENCIES #########################

include(FetchContent)
include(ExternalProject)
include(GNUInstallDirs)

# ---- Dependency options ----
option(ODT_SYSTEM_YAML      "Use system-wide yaml installation"         OFF)
option(ODT_SYSTEM_SUNDIALS  "Use system-wide Sundials installation"     OFF)
option(ODT_SYSTEM_CANTERA   "Use system-wide Cantera installation"      OFF)
#set(ODT_CANTERA_INCLUDE_LOCATION    /usr/include    CACHE STRING "")
#set(ODT_CANTERA_LIBRARY_LOCATION    /usr/lib        CACHE STRING "")

# ---- set RPATH (runtime library locations) ----
file(RELATIVE_PATH reldir
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH $ORIGIN $ORIGIN/${reldir})

# ---- FetchContent declarations ----
FetchContent_Declare(yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG 9a3624205e8774953ef18f57067b3426c1c5ada6)           # yaml-cpp-0.6.3
FetchContent_Declare(sundials
        GIT_REPOSITORY https://github.com/LLNL/sundials.git
        GIT_TAG 887af4374af2271db9310d31eaa9b5aeff49e829)           # Sundials-5.3.0
FetchContent_Declare(radlib
        GIT_REPOSITORY https://github.com/BYUignite/radlib.git
        GIT_TAG da2420bc46a94e0c96f487735fda949bca96d07c)           # radlib-2.0.0
FetchContent_Declare(sootlib
        GIT_REPOSITORY https://github.com/BYUignite/sootlib.git
        GIT_TAG 5f5eeb6fdb1ef66e7d82fd136026bb20f06208e5)           # sootlib-0.1.0

# -------- Boost + Threads --------
find_package(Boost REQUIRED COMPONENTS system)
find_package(Threads REQUIRED)

# -------- RADLIB --------
message(STATUS "Fetching RadLib")
FetchContent_MakeAvailable(radlib)

# -------- SOOTLIB --------
message(STATUS "Fetching SootLib")
FetchContent_MakeAvailable(sootlib)

# -------- YAML --------
message(STATUS "Looking for yaml-cpp")
list(APPEND CMAKE_MESSAGE_INDENT "    ")

if (ODT_SYSTEM_YAML)
    find_package(yaml-cpp)
    if (yaml-cpp_FOUND)
        message(STATUS "Found---using system installation of yaml-cpp")
    endif ()
endif ()

if (NOT ODT_SYSTEM_YAML OR NOT yaml-cpp_FOUND)
    FetchContent_MakeAvailable(yaml-cpp)
    list(APPEND yaml-cpp_INCLUDE_DIR ${FETCHCONTENT_BASE_DIR}/yaml-cpp-build/include)
    message(STATUS "Found---using yaml fetched from repository")
endif ()

list(POP_BACK CMAKE_MESSAGE_INDENT)

# -------- SUNDIALS --------
message(STATUS "Looking for Sundials")
list(APPEND CMAKE_MESSAGE_INDENT "    ")

if (ODT_SYSTEM_SUNDIALS)
    #    set(SUNDIALS_DIR "/usr/local/lib")
    find_package(sundials REQUIRED COMPONENTS cvodes nvecserial sunmatrixdense)      # throws an error if required components not found
    if(sundials_FOUND)
        message(STATUS "Found---using system installation of Sundials")
    endif()
endif ()

if (NOT ODT_SYSTEM_SUNDIALS OR NOT sundials_FOUND)
    FetchContent_GetProperties(sundials)
    if (NOT sundials_POPULATED)
        FetchContent_Populate(sundials)                                     # populate content from repository
        set(BUILD_ARKODE             OFF     CACHE INTERNAL "")             # set options to only build the parts we need
        set(BUILD_CVODE              OFF     CACHE INTERNAL "")
        set(BUILD_CVODES             ON      CACHE INTERNAL "")
        set(BUILD_IDA                OFF     CACHE INTERNAL "")
        set(BUILD_IDAS               OFF     CACHE INTERNAL "")
        set(BUILD_KINSOL             OFF     CACHE INTERNAL "")
        set(BUILD_NVECTOR_SERIAL     ON      CACHE INTERNAL "")
        set(BUILD_NVECTOR_MANYVECTOR OFF     CACHE INTERNAL "")
        set(EXAMPLES_ENABLE_C        OFF     CACHE INTERNAL "")
        set(ENABLE_LAPACK            ON      CACHE INTERNAL "")             # required lapack parts will be installed automatically
        set(CMAKE_C_FLAGS            "-fPIC" CACHE INTERNAL "")             # required flag for building shared libraries

        add_subdirectory(${sundials_SOURCE_DIR} ${sundials_BINARY_DIR})                     # Bring the populated content into the build
        list(APPEND sundials_INCLUDE_DIR ${FETCHCONTENT_BASE_DIR}/sundials-build/include)   # tell CMake where to find sundials headers
        list(APPEND sundials_INCLUDE_DIR ${FETCHCONTENT_BASE_DIR}/sundials-src/include)

        message(STATUS "Found---using Sundials fetched from repository")
    endif ()
endif ()

list(POP_BACK CMAKE_MESSAGE_INDENT)

# -------- CANTERA --------
message(STATUS "Looking for Cantera")
list(APPEND CMAKE_MESSAGE_INDENT "    ")

#if (ODT_SYSTEM_CANTERA)
#    find_package(Cantera REQUIRED)      # throws an error if required components not found
#    if (Cantera_FOUND)
#        message(STATUS "Found---using system installation of Cantera")
#    endif ()
#endif ()

if(NOT cantera_POPULATED)
    FetchContent_Populate(cantera)
endif()
set(CANTERA_PATH ${cantera_SOURCE_DIR})
message(STATUS "CANTERA SOURCE DIR = ${cantera_INSTALL_DIR}")

include(ExternalProject)
ExternalProject_Add(cantera
        PREFIX      ${CANTERA_PATH}
        SOURCE_DIR  ${CANTERA_PATH}/src
        BINARY_DIR  ${CANTERA_PATH}
        INSTALL_DIR ${CANTERA_PATH}/install
        CONFIGURE_COMMAND ""
        BUILD_COMMAND scons build python_package=none matlab_toolbox=n f90_interface=n
        )

#ExternalProject_Add(Cantera
#        PREFIX              ${FETCHCONTENT_BASE_DIR}/cantera
#        GIT_REPOSITORY      https://github.com/Cantera/cantera
#        GIT_TAG             b0bace78223959cd3e5a15317734cacff7b0b0a2        # cantera 2.5.1
#        BINARY_DIR          ${FETCHCONTENT_BASE_DIR}/cantera/src/Cantera
#        CONFIGURE_COMMAND   ""
#        BUILD_COMMAND       scons build python_package=none matlab_toolbox=n f90_interface=n prefix=${PROJECT_SOURCE_DIR}
#        INSTALL_COMMAND     scons install
#        )

#list(APPEND Cantera_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
#list(APPEND Cantera_LIBRARY_DIR ${PROJECT_SOURCE_DIR}/lib)

#add_library(cantera SHARED IMPORTED GLOBAL)
#
#set_target_properties(cantera PROPERTIES
#        IMPORTED_LOCATION ${CANTERA_PATH}/lib/${CANTERA_PATH_SUFFIX}/cantera${CMAKE_LINK_LIBRARY_SUFFIX}
#        IMPORTED_IMPLIB ${CANTERA_PATH}/lib/${CANTERA_PATH_SUFFIX}/cantera${CMAKE_IMPORT_LIBRARY_SUFFIX}
#        INTERFACE_INCLUDE_DIRECTORIES ${CANTERA_PATH}/include
#)



list(POP_BACK CMAKE_MESSAGE_INDENT)