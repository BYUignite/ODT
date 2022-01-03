############################### DOCUMENTATION ###########################

message(STATUS "Looking for Doxygen")
find_package(Doxygen COMPONENTS dot)
if(Doxygen_FOUND)
    message(STATUS "Found Doxygen")
    set(DOXYGEN_IN  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    add_custom_target(odt-docs ALL
            COMMAND rm -rf ../docs/doxygen/html
            COMMAND mkdir  ../docs/doxygen/html
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating ODT documentation with Doxygen"
            VERBATIM )
else()
    message(STATUS "Doxygen not found---documentation will not be built")
endif()