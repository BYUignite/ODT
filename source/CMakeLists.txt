#################### EXECUTABLE ##########################

add_executable(odt.x)
target_compile_features(odt.x PUBLIC cxx_std_11)

#################### Local source files

target_include_directories(odt.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt.x PRIVATE main.cc
                               inputoutput.cc  inputoutput.h
                               domain.cc       domain.h
                               processor.cc    processor.h
                               streams.cc      streams.h
                               meshManager.cc  meshManager.h
                               param.cc        param.h
                               micromixer.cc   micromixer.h
                               eddy.cc         eddy.h
                               solver.cc       solver.h
                               cvodeDriver.cc  cvodeDriver.h
                                               interp_linear.h
              )

#################### Subdirectories

add_subdirectory(domainvariables)
add_subdirectory(radiation)
add_subdirectory(domaincases)
add_subdirectory(user_chemical_mechanisms)


#target_include_directories(odt.x PRIVATE ${Boost_INCLUDE_DIR})
#target_link_libraries(odt.x PUBLIC Threads::Threads -pthread)
#
#target_include_directories(odt.x PRIVATE ${YAML_INCLUDE_DIR})
#target_link_libraries(odt.x PUBLIC yaml-cpp -lyaml-cpp)
#
#target_include_directories(odt.x PRIVATE ${SUNDIALS_INCLUDE_DIR})
#target_link_libraries(odt.x PUBLIC SUNDIALS::cvode SUNDIALS::nvecserial SUNDIALS::sunmatrixdense)
#
#include_directories("${Cantera_INCLUDE_DIR}/cantera")
#link_directories("${Cantera_LIBRARY}")
#target_include_directories(odt.x PRIVATE ${CANTERA_INCLUDE_DIR})
#target_link_libraries(odt.x PUBLIC cantera sundials_cvodes sundials_sunlinsollapackband sundials_sunlinsollapackdense lapack blas -lyaml-cpp)