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