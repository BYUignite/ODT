#################### Doxygen documentation

find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_custom_target(docs
        COMMENT "Building doxygen documentation"
        COMMAND rm -rf doxygen/html
        COMMAND mkdir  doxygen/html
        COMMAND doxygen -s Doxyfile
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        VERBATIM
    )
else()
    message("Cannot build documentation because doxygen was not found")
endif()
