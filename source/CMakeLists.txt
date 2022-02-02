#################### EXECUTABLE ##########################

# -------- executable targets --------
add_executable(odt-run)

# -------- target properties --------
target_compile_features(odt-run PUBLIC cxx_std_11)

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
list(APPEND sundials_INCLUDE_DIR ${FETCHCONTENT_BASE_DIR}/sundials-build/include)   # tell CMake where to find sundials headers
list(APPEND sundials_INCLUDE_DIR ${FETCHCONTENT_BASE_DIR}/sundials-src/include)

target_include_directories(odt-run PRIVATE ${Boost_INCLUDE_DIR} ${yaml-cpp_INCLUDE_DIR} ${sundials_INCLUDE_DIR} ${cantera_INCLUDE_DIR})

target_link_libraries(odt-run PUBLIC Threads::Threads -pthread)
target_link_libraries(odt-run PUBLIC yaml-cpp)
target_link_libraries(odt-run PUBLIC radlib::radProps sootlib::sootModel)
target_link_libraries(odt-run PUBLIC sundials_cvodes sundials_sunlinsollapackband sundials_sunlinsollapackdense lapack blas)
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
        EXPORT  odt
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )