#################### EXECUTABLE ##########################

# -------- executable targets --------
add_executable(odt-runOneRlz runOneRlz.cc)
add_executable(odt-runManyRlz runManyRlz.cc)

# -------- target properties --------
target_compile_features(odt-runOneRlz PUBLIC cxx_std_11)
target_compile_features(odt-runManyRlz PUBLIC cxx_std_11)

add_dependencies(odt-runOneRlz odt-run)
add_dependencies(odt-runManyRlz odt-run)

############# INSTALLATION ################

include(GNUInstallDirs)

install(TARGETS odt-runOneRlz odt-runManyRlz
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )