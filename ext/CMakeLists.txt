################### DEPENDENCIES #########################

include(FetchContent)
#find_package(PkgConfig REQUIRED)

set(FETCHCONTENT_SOURCE_DIR_RADLIB  "/home/vbstephens/codes/radlib"  CACHE STRING "local radlib development location")
set(FETCHCONTENT_SOURCE_DIR_SOOTLIB "/home/vbstephens/codes/sootlib" CACHE STRING "local sootlib development location")

# ---- FetchContent declarations ----
FetchContent_Declare(yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG 9a3624205e8774953ef18f57067b3426c1c5ada6)           # yaml-cpp-0.6.3
FetchContent_Declare(sundials
        GIT_REPOSITORY https://github.com/LLNL/sundials.git
        GIT_TAG 887af4374af2271db9310d31eaa9b5aeff49e829)           # sundials-5.3.0
FetchContent_Declare(radlib
        GIT_REPOSTIORY https://github.com/BYUignite/radlib.git
        GIT_TAG a77c8052d3b9373de0f1996800af2878d7d6ea9b)           # radlib-2.0.0
FetchContent_Declare(sootlib
        GIT_REPOSTIORY https://github.com/BYUignite/sootlib.git
        GIT_TAG cc137fb82107a422df84ecd826575aba2e585beb)           # sootlib-0.1.0





find_package(Boost REQUIRED COMPONENTS system)
find_package(Threads REQUIRED)



#--------------------- YAML ------------------------

option(ODT_SYSTEM_YAML ON)
message(STATUS "Looking for yaml-cpp")

if (ODT_SYSTEM_YAML)
    find_package(yaml-cpp)
#    pkg_check_modules(yaml-cpp IMPORTED_TARGET)
    if(yaml-cpp_FOUND)
        message(STATUS "Found---using system yaml-cpp")
    endif()
endif ()

if (!ODT_SYSTEM_YAML OR !yaml-cpp_FOUND)
    FetchContent_MakeAvailable(yaml-cpp)
    if(yaml-cpp_FOUND)
        message(STATUS "Found---using repository yaml-cpp")
    else()
        message(FATAL_ERROR "yaml-cpp not found")
    endif()
endif ()


#--------------------- SUNDIALS ------------------------
#message(STATUS "Looking for Sundials")
#set(SUNDIALS_DIR "/usr/local/lib")
#find_package(SUNDIALS REQUIRED)
#message(STATUS "Found Sundials")



#FetchContent_GetProperties(sundials)
#if(NOT sundials_POPULATED)
#    FetchContent_Populate(sundials)
#    set(BUILD_ARKODE             OFF     CACHE INTERNAL "")
#    set(BUILD_CVODE              ON      CACHE INTERNAL "")
#    set(BUILD_CVODES             OFF     CACHE INTERNAL "")
#    set(BUILD_IDA                OFF     CACHE INTERNAL "")
#    set(BUILD_IDAS               OFF     CACHE INTERNAL "")
#    set(BUILD_KINSOL             OFF     CACHE INTERNAL "")
#    set(BUILD_NVECTOR_SERIAL     ON      CACHE INTERNAL "")
#    set(BUILD_NVECTOR_MANYVECTOR OFF     CACHE INTERNAL "")
#    set(EXAMPLES_ENABLE          OFF     CACHE INTERNAL "")
#    set(ENABLE_LAPACK            ON      CACHE INTERNAL "")
#    set(CMAKE_C_FLAGS            "-fPIC" CACHE INTERNAL "")
#
#    # Bring the populated content into the build
#    add_subdirectory(${sundials_SOURCE_DIR} ${sundials_BINARY_DIR})
#endif()

#--------------------- CANTERA ------------------------
#message(STATUS "Looking for Cantera")
#find_package(Cantera REQUIRED)
#message(STATUS "Found Cantera")