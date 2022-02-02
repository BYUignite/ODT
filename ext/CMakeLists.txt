################### DEPENDENCIES #########################




# ---- Dependency options ----
option(ODT_SYSTEM_YAML      "Use system-wide yaml installation"         OFF)
option(ODT_SYSTEM_SUNDIALS  "Use system-wide Sundials installation"     OFF)
option(ODT_SYSTEM_CANTERA   "Use system-wide Cantera installation"      OFF)
#set(ODT_CANTERA_INCLUDE_LOCATION    /usr/include    CACHE STRING "")
#set(ODT_CANTERA_LIBRARY_LOCATION    /usr/lib        CACHE STRING "")



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
FetchContent_Declare(cantera
        GIT_REPOSITORY https://github.com/Cantera/cantera
        GIT_TAG b0bace78223959cd3e5a15317734cacff7b0b0a2)           # cantera-2.5.1



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
        message(STATUS "Using system installation of yaml-cpp")
    endif ()
endif ()

if (NOT ODT_SYSTEM_YAML OR NOT yaml-cpp_FOUND)
    FetchContent_MakeAvailable(yaml-cpp)
    list(APPEND yaml-cpp_INCLUDE_DIR ${FETCHCONTENT_BASE_DIR}/yaml-cpp-build/include)
    message(STATUS "Using yaml fetched from repository")
endif ()

list(POP_BACK CMAKE_MESSAGE_INDENT)

# -------- SUNDIALS --------
message(STATUS "Looking for Sundials")
list(APPEND CMAKE_MESSAGE_INDENT "    ")

#if (ODT_SYSTEM_SUNDIALS)
#    set(SUNDIALS_DIR "/usr/local/lib")      # sundials default install location
#    find_package(sundials REQUIRED COMPONENTS cvodes nvecserial sunmatrixdense)      # throws an error if required components not found
#    if(sundials_FOUND)
#        message(STATUS "Using system installation of Sundials")
#    endif()
#endif ()

if (NOT ODT_SYSTEM_SUNDIALS OR NOT sundials_FOUND)
    message(STATUS "Using Sundials fetched from repository")
    FetchContent_GetProperties(sundials)
    if (NOT sundials_POPULATED)
        FetchContent_Populate(sundials)                                     # populate content from repository
        set(BUILD_ARKODE             OFF     CACHE INTERNAL "")             # set options to only build the parts we need
        set(BUILD_CVODE              OFF     CACHE INTERNAL "")
        set(BUILD_CVODES             ON      CACHE INTERNAL "")
        set(BUILD_IDA                OFF     CACHE INTERNAL "")
        set(BUILD_IDAS               OFF     CACHE INTERNAL "")
        set(BUILD_KINSOL             OFF     CACHE INTERNAL "")
        set(EXAMPLES_ENABLE_C        OFF     CACHE INTERNAL "")
        set(EXAMPLES_ENABLE_F77      OFF     CACHE INTERNAL "")
        set(EXAMPLES_ENABLE_F90      OFF     CACHE INTERNAL "")
        set(EXAMPLES_ENABLE_F2003    OFF     CACHE INTERNAL "")
        set(EXAMPLES_INSTALL         OFF     CACHE INTERNAL "")
        set(LAPACK_ENABLE            ON      CACHE INTERNAL "")             # required lapack parts will be installed automatically
        set(CMAKE_C_FLAGS            "-fPIC" CACHE INTERNAL "")             # required flag for building shared libraries
        set(BUILD_SHARED_LIBS        ON      CACHE INTERNAL "")
        set(BUILD_STATIC_LIBS        OFF     CACHE INTERNAL "")
        add_subdirectory(${sundials_SOURCE_DIR} ${sundials_BINARY_DIR})     # Bring the populated content into the build
    endif ()
endif ()

list(POP_BACK CMAKE_MESSAGE_INDENT)

# -------- CANTERA --------
message(STATUS "Looking for Cantera")
list(APPEND CMAKE_MESSAGE_INDENT "    ")

#if (ODT_SYSTEM_CANTERA)
#    find_package(cantera REQUIRED)      # throws an error if required package not found
#    if (cantera_FOUND)
#        message(STATUS "Using system installation of Cantera")
#    endif ()
#endif ()

if (NOT ODT_SYSTEM_CANTERA OR NOT cantera_FOUND)
    FetchContent_GetProperties(cantera)
    if (NOT cantera_POPULATED)
        FetchContent_Populate(cantera)
    endif()
    ExternalProject_Add(cantera
        SOURCE_DIR  ${cantera_SOURCE_DIR}
        BINARY_DIR  ${cantera_SOURCE_DIR}
        INSTALL_DIR ${CMAKE_CURRENT_LIST_DIR}/cantera
        CONFIGURE_COMMAND ""
        BUILD_COMMAND scons build python_package=n matlab_toolbox=n f90_interface=n googletest=none prefix=${CMAKE_CURRENT_LIST_DIR}/cantera
        INSTALL_COMMAND scons install
    )
    message(STATUS "Using Cantera fetched from repository")
endif()

add_library(cantera-ext SHARED IMPORTED GLOBAL)

set_target_properties(cantera-ext PROPERTIES
#        IMPORTED_LOCATION ${CMAKE_CURRENT_LIST_DIR}/cantera/lib/libcantera${CMAKE_LINK_LIBRARY_SUFFIX}
        IMPORTED_LOCATION ${CMAKE_CURRENT_LIST_DIR}/cantera/lib/libcantera.a
        IMPORTED_IMPLIB ${CMAKE_CURRENT_LIST_DIR}/cantera/lib/libcantera${CMAKE_IMPORT_LIBRARY_SUFFIX}
        INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_LIST_DIR}/cantera/include
)

list(POP_BACK CMAKE_MESSAGE_INDENT)