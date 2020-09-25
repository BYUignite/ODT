find_package(Git)
if(GIT_FOUND)
    add_custom_target(yaml 
        COMMAND echo "rm -rf include lib yaml-cpp"                       > delme.sh
        COMMAND echo "git clone https://github.com/jbeder/yaml-cpp.git" >> delme.sh
        COMMAND echo "cd yaml-cpp"                                      >> delme.sh
        COMMAND echo "git checkout tags/yaml-cpp-0.6.3"                 >> delme.sh
        COMMAND echo "mkdir build"                                      >> delme.sh
        COMMAND echo "cd build"                                         >> delme.sh
        COMMAND echo "cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../.."  >> delme.sh
        COMMAND echo "make"                                             >> delme.sh
        COMMAND echo "make install"                                     >> delme.sh
        COMMAND bash ./delme.sh
        COMMAND rm -f delme.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Downloading and building yaml-cpp"
        VERBATIM
    )
else()
    message("Cannot get and build yaml because git was not found")
endif()

