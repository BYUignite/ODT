################### CMAKE PROJECT ODT #########################

cmake_minimum_required(VERSION 3.14)

project(ODT
        VERSION 2.0
        LANGUAGES CXX C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)      # Release or Debug
endif()

#set(CMAKE_VERBOSE_MAKEFILE TRUE)        # vbs debug

################### BUILD & RUN DIRECTORIES #########################

set(CMAKE_BINARY_DIR  "${CMAKE_SOURCE_DIR}/build")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")       # tells find_package() where to look for Find[pkg].cmake files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/run")

# ---- C++ compiler flags -------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++11 -pthread"    CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3"           CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3"            CACHE STRING "" FORCE)

#################### EXECUTABLE ##########################

add_executable(odt.x)
target_compile_features(odt.x PUBLIC cxx_std_11)

################### PACKAGES & LIBRARIES #########################

include(FetchContent)

find_package(PkgConfig)
find_package(Boost REQUIRED COMPONENTS system)
find_package(Threads REQUIRED)

target_include_directories(odt.x PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(odt.x PUBLIC Threads::Threads -pthread)

#--------------------- YAML ------------------------
message(STATUS "Looking for yaml-cpp")
FIND_PACKAGE(yaml-cpp REQUIRED)
if(yaml-cpp_FOUND)
    message(STATUS "Found yaml-cpp")
else()
    message(STATUS "yaml-cpp not found---fetching from repository")
    FetchContent_Declare(yaml-cpp
            GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
            GIT_TAG yaml-cpp-0.6.3)
    FetchContent_MakeAvailable(yaml-cpp)
endif()

target_include_directories(odt.x PRIVATE ${YAML_INCLUDE_DIR})
target_link_libraries(odt.x PUBLIC yaml-cpp -lyaml-cpp)

#--------------------- SUNDIALS ------------------------
message(STATUS "Looking for Sundials")
set(SUNDIALS_DIR "/usr/local/lib")
find_package(SUNDIALS REQUIRED)
message(STATUS "Found Sundials")

target_include_directories(odt.x PRIVATE ${SUNDIALS_INCLUDE_DIR})
target_link_libraries(odt.x PUBLIC SUNDIALS::cvode SUNDIALS::nvecserial SUNDIALS::sunmatrixdense)

#FetchContent_Declare(sundials
#        GIT_REPOSITORY https://github.com/LLNL/sundials.git
#        GIT_TAG v5.8.0)
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
message(STATUS "Looking for Cantera")
find_package(Cantera REQUIRED)
message(STATUS "Found Cantera")

include_directories("${Cantera_INCLUDE_DIR}/cantera")
link_directories("${Cantera_LIBRARY}")

target_include_directories(odt.x PRIVATE ${CANTERA_INCLUDE_DIR})
target_link_libraries(odt.x PUBLIC cantera sundials_cvodes sundials_sunlinsollapackband sundials_sunlinsollapackdense lapack blas -lyaml-cpp)

#################### TESTS ###############################

#if (BUILD_TESTS)
#    include(Ctest)
#    add_subdirectory(test)
#endif()

#################### EXAMPLES ###############################

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

############################### DOCUMENTATION ###########################

if (BUILD_DOCS)
    message(STATUS "Looking for Doxygen")
    find_package(Doxygen COMPONENTS dot)
    if(Doxygen_FOUND)
        message(STATUS "Found Doxygen")
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        add_custom_target( docs_doxygen ALL
                COMMAND rm -rf ../docs/doxygen/html
                COMMAND mkdir  ../docs/doxygen/html
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating ODT documentation with Doxygen"
                VERBATIM )
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
add_definitions(-DDOCANTERA)