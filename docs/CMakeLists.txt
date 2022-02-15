############################### DOCUMENTATION ###########################

message(STATUS "Looking for Doxygen")
find_package(Doxygen COMPONENTS dot)
if(Doxygen_FOUND)
#    message(STATUS "Found Doxygen")
    set(DOXYGEN_IN  ${PROJECT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${PROJECT_BINARY_DIR}/docs/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    add_custom_target(odt-docs ALL
            file(REMOVE_RECURSE ${PROJECT_SOURCE_DIR}/docs/doxygen)
            file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/docs/doxygen)
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            COMMENT "Generating ODT documentation with Doxygen"
            VERBATIM )
else()
    message(STATUS "Doxygen not found---documentation will not be built")
endif()