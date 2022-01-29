include(FindPkgConfig)
pkg_check_modules(PC_yaml-cpp QUIET yaml-cpp)

find_path(yaml-cpp_INCLUDE_DIR
        NAMES yaml-cpp
        PATHS ${PC_yaml-cpp_INCLUDE_DIRS} /usr/include
        PATH_SUFFIXES yaml-cpp
        )
find_library(yaml-cpp_LIBRARY
        NAMES libyaml-cpp.a
        PATHS ${PC_yaml-cpp_LIBRARY_DIRS} /usr/lib
        )

set(yaml-cpp_VERSION ${PC_yaml-cpp_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(yaml-cpp
        FOUND_VAR yaml-cpp_FOUND
        REQUIRED_VARS
        yaml-cpp_LIBRARY
        yaml-cpp_INCLUDE_DIR
        VERSION_VAR yaml-cpp_VERSION
        )

if(yaml-cpp_FOUND)
    set(yaml-cpp_LIBRARIES ${yaml-cpp_LIBRARY})
    set(yaml-cpp_INCLUDE_DIRS ${yaml-cpp_INCLUDE_DIR})
    set(yaml-cpp_DEFINITIONS ${PC_yaml-cpp_CFLAGS_OTHER})
endif()

if(yaml-cpp_FOUND AND NOT TARGET yaml-cpp::yaml-cpp)
    add_library(yaml-cpp::yaml-cpp UNKNOWN IMPORTED)
    set_target_properties(yaml-cpp::yaml-cpp PROPERTIES
            IMPORTED_LOCATION "${yaml-cpp_LIBRARY}"
            INTERFACE_COMPILE_OPTIONS "${PC_yaml-cpp_CFLAGS_OTHER}"
            INTERFACE_INCLUDE_DIRECTORIES "${yaml-cpp_INCLUDE_DIR}"
            )
endif()

mark_as_advanced(
        yaml-cpp_INCLUDE_DIR
        yaml-cpp_LIBRARY
)

# compatibility variables
set(yaml-cpp_VERSION_STRING ${yaml-cpp_VERSION})