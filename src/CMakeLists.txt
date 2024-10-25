# Set an environment variable CMAKE_PREFIX_PATH with paths for find_library etc to search
# or add: HINTS "/path/to/lib/dir" at the end of the find_library call

add_executable(odt.x "")

target_sources(odt.x PRIVATE main.cc
                             inputoutput.cc
                             domain.cc
                             processor.cc
                             streams.cc
                             meshManager.cc
                             param.cc
                             micromixer.cc
                             eddy.cc
                             solver.cc
                             cvodeDriver.cc)

add_subdirectory(cvode)
add_subdirectory(domainvariables)
add_subdirectory(radiation)
add_subdirectory(domaincases)
add_subdirectory(user_chemical_mechanisms)

####################

set_property(TARGET odt.x PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)  # adds LC_RPATH to executable

find_library(YAML   yaml-cpp)
find_library(CANTERA NAMES cantera_shared cantera)

target_include_directories(odt.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(odt.x PRIVATE "$ENV{CONDA_PREFIX}/include")
target_include_directories(odt.x PRIVATE "/usr/local/include")

target_link_libraries(odt.x ${CANTERA} ${YAML})

####################

install(TARGETS odt.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
