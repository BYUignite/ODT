
add_executable(odt.x)

#################### Compile options

target_compile_features(odt.x PUBLIC cxx_std_11)

#----------- Link libraries and flags ---------
target_include_directories(odt.x PRIVATE ${YAML_INCLUDE_DIR}
        ${CANTERA_INCLUDE_DIR}
        ${SUNDIALS_INCLUDE_DIR}
        ${Boost_INCLUDE_DIR}
        ${LAPACK_INCLUDE_DIR})

set(ODT_LINK_FLAGS "-pthread -lcantera -lyaml-cpp -lblas -llapack -lsundials_cvodes -lsundials_ida -lsundials_nvecserial -lgfortran")

target_link_libraries(odt.x PUBLIC cantera sundials_nvecserial
        sundials_cvodes lapack yaml-cpp blas ${ODT_LINK_FLAGS})

#################### Local source files

target_include_directories(odt.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/inputoutput.cc             ${CMAKE_CURRENT_SOURCE_DIR}/inputoutput.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domain.cc                  ${CMAKE_CURRENT_SOURCE_DIR}/domain.h
        ${CMAKE_CURRENT_SOURCE_DIR}/processor.cc               ${CMAKE_CURRENT_SOURCE_DIR}/processor.h
        ${CMAKE_CURRENT_SOURCE_DIR}/streams.cc                 ${CMAKE_CURRENT_SOURCE_DIR}/streams.h
        ${CMAKE_CURRENT_SOURCE_DIR}/meshManager.cc             ${CMAKE_CURRENT_SOURCE_DIR}/meshManager.h
        ${CMAKE_CURRENT_SOURCE_DIR}/param.cc                   ${CMAKE_CURRENT_SOURCE_DIR}/param.h
        ${CMAKE_CURRENT_SOURCE_DIR}/micromixer.cc              ${CMAKE_CURRENT_SOURCE_DIR}/micromixer.h
        ${CMAKE_CURRENT_SOURCE_DIR}/eddy.cc                    ${CMAKE_CURRENT_SOURCE_DIR}/eddy.h
        ${CMAKE_CURRENT_SOURCE_DIR}/solver.cc                  ${CMAKE_CURRENT_SOURCE_DIR}/solver.h
        ${CMAKE_CURRENT_SOURCE_DIR}/cvodeDriver.cc             ${CMAKE_CURRENT_SOURCE_DIR}/cvodeDriver.h
                                                               ${CMAKE_CURRENT_SOURCE_DIR}/interp_linear.h

)

#################### Subdirectories

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/domainvariables)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/radiation)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/domaincases)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/user_chemical_mechanisms)

#################### Doxygen documentation

find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_custom_target(doxygen
        COMMAND rm -rf ../doc/doxygen/html
        COMMAND mkdir  ../doc/doxygen/html
        COMMAND cp    ../doc/doxygen_build_pages/montage_odt.png ../doc/doxygen/html
        COMMAND doxygen -s ../doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building doxygen documentation"
        VERBATIM
    )
else()
    message("\n********************************************************")
    message("Cannot build documentation because doxygen was not found")
    message("********************************************************\n")
endif()