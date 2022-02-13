################### DEPENDENCIES #########################

# ---- ODT dependency options ----
option(ODT_SYSTEM_YAML          "Use system-wide yaml installation"         OFF)
option(ODT_SYSTEM_SUNDIALS      "Use system-wide Sundials installation"     OFF)
option(ODT_SYSTEM_CANTERA       "Use system-wide Cantera installation"      OFF)
option(YAML_CPP_BUILD_TESTS     "Build yaml-cpp tests"                      OFF)
option(YAML_CPP_BUILD_CONTRIB   "Build yaml-cpp contrib stuff"              OFF)

# ---- FetchContent declarations ----
FetchContent_Declare(yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG 9a3624205e8774953ef18f57067b3426c1c5ada6)           # yaml-cpp-0.6.3
FetchContent_Declare(sundials
        GIT_REPOSITORY https://github.com/LLNL/sundials.git
        GIT_TAG 887af4374af2271db9310d31eaa9b5aeff49e829)           # Sundials-5.3.0
FetchContent_Declare(radlib
        GIT_REPOSITORY https://github.com/BYUignite/radlib.git
        GIT_TAG bff973e22befdc1db16d55a2f04e675e35c4aea3)           # radlib-2.0.0
FetchContent_Declare(sootlib
        GIT_REPOSITORY https://github.com/BYUignite/sootlib.git
        GIT_TAG a9a96b68e676c65f0d05a17ee0f0f6c5dbb2a207)           # sootlib-0.1.0
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

# -------- CANTERA --------
message(STATUS "Looking for Cantera")
list(APPEND CMAKE_MESSAGE_INDENT "    ")

if (ODT_SYSTEM_CANTERA)
    find_package(cantera REQUIRED)      # throws an error if required package not found
    if (cantera_FOUND)
        message(STATUS "Using system installation of Cantera")
    endif ()
endif ()

if (NOT ODT_SYSTEM_CANTERA OR NOT cantera_FOUND)
    FetchContent_GetProperties(cantera)
    if (NOT cantera_POPULATED)
        FetchContent_Populate(cantera)
    endif()
    ExternalProject_Add(cantera
        SOURCE_DIR  ${cantera_SOURCE_DIR}
        BINARY_DIR  ${cantera_SOURCE_DIR}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND scons build python_package=n matlab_toolbox=n f90_interface=n googletest=none prefix=${FETCHCONTENT_BASE_DIR}/cantera-install system_sundials=n
        INSTALL_COMMAND scons install
    )
    message(STATUS "Using Cantera fetched from repository")
endif()

file(MAKE_DIRECTORY ${FETCHCONTENT_BASE_DIR}/cantera-install/include)

add_library(cantera-ext SHARED IMPORTED GLOBAL)
set_target_properties(cantera-ext PROPERTIES
        IMPORTED_LOCATION ${FETCHCONTENT_BASE_DIR}/cantera-install/lib/libcantera.a
        IMPORTED_IMPLIB ${FETCHCONTENT_BASE_DIR}/cantera-install/lib/libcantera.a
        INTERFACE_INCLUDE_DIRECTORIES ${FETCHCONTENT_BASE_DIR}/cantera-install/include
)

list(POP_BACK CMAKE_MESSAGE_INDENT)

# -------- SUNDIALS --------
message(STATUS "Looking for Sundials")
list(APPEND CMAKE_MESSAGE_INDENT "    ")

if (ODT_SYSTEM_SUNDIALS)
    set(SUNDIALS_DIR "/usr/local/lib")      # sundials default install location
    find_package(sundials REQUIRED COMPONENTS cvodes nvecserial sunmatrixdense)      # throws an error if required components not found
    if(sundials_FOUND)
        message(STATUS "Using system installation of Sundials")
    endif()
else()
    message(STATUS "Using Sundials installed by Cantera")
endif ()

list(POP_BACK CMAKE_MESSAGE_INDENT)