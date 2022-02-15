#################### EXECUTABLE ##########################

# -------- executable targets --------
add_executable(odt-run)

# -------- target properties --------
target_compile_features(odt-run PUBLIC cxx_std_11)

# -------- preprocessor directives --------

if(ODT_PARALLEL)
    message(STATUS "MPI support is on")
    add_compile_definitions(DOMPI)
endif()

if(NOT ODT_CHEMISTRY STREQUAL "CANTERARR")
    add_compile_definitions(PROBLEMSPECIFICRR)
endif()

if(ODT_VERBOSE)
#    message(STATUS "ODT will suppress some runtime output")
    add_compile_definitions(VERBOSE)
endif()

add_compile_definitions(_USESTDVECTOR_)
add_compile_definitions(DOCANTERA)

# -------- src files --------
target_include_directories(odt-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(odt-run PRIVATE main.cc
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

# -------- src subdirectories --------
add_subdirectory(domainvariables)
add_subdirectory(radiation)
add_subdirectory(domaincases)
add_subdirectory(user_chemical_mechanisms)

# -------- dependencies --------
get_target_property(cantera_INCLUDE_DIR cantera-ext INTERFACE_INCLUDE_DIRECTORIES)  # tell CMake where to find cantera headers
list(APPEND sundials_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/cantera/ext)

target_include_directories(odt-run PRIVATE ${Boost_INCLUDE_DIR} ${yaml-cpp_INCLUDE_DIR} ${sundials_INCLUDE_DIR} ${cantera_INCLUDE_DIR})

target_link_libraries(odt-run PUBLIC Threads::Threads -pthread)
target_link_libraries(odt-run PUBLIC yaml-cpp)
target_link_libraries(odt-run PUBLIC radlib::radProps sootlib::sootModel)
target_link_libraries(odt-run PUBLIC cantera-ext -lyaml-cpp)

############# INSTALLATION ################

include(GNUInstallDirs)

#---- allow odt to be used in other targets ----
target_include_directories(odt-run
        INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
        )

install(TARGETS odt-run
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )