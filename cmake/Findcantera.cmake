include(FindPkgConfig)
pkg_check_modules(PC_cantera QUIET cantera)

find_path(cantera_INCLUDE_DIR
        NAMES cantera Cantera
        PATHS ${PC_cantera_INCLUDE_DIRS} ${ODT_CANTERA_INCLUDE_LOCATION}
        PATH_SUFFIXES cantera Cantera
        )
find_library(cantera_LIBRARY
        NAMES libcantera.a
        PATHS ${PC_cantera_LIBRARY_DIRS} ${ODT_CANTERA_LIBRARY_LOCATION}
        )

set(cantera_VERSION ${PC_cantera_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cantera
        FOUND_VAR cantera_FOUND
        REQUIRED_VARS
        cantera_LIBRARY
        cantera_INCLUDE_DIR
        VERSION_VAR cantera_VERSION
        )

if(cantera_FOUND)
    set(cantera_LIBRARIES ${cantera_LIBRARY})
    set(cantera_INCLUDE_DIRS ${cantera_INCLUDE_DIR})
    set(cantera_DEFINITIONS ${PC_cantera_CFLAGS_OTHER})
endif()

if(cantera_FOUND AND NOT TARGET cantera::cantera)
    add_library(cantera::cantera UNKNOWN IMPORTED)
    set_target_properties(cantera::cantera PROPERTIES
            IMPORTED_LOCATION "${cantera_LIBRARY}"
            INTERFACE_COMPILE_OPTIONS "${PC_cantera_CFLAGS_OTHER}"
            INTERFACE_INCLUDE_DIRECTORIES "${cantera_INCLUDE_DIR}"
            )
endif()

mark_as_advanced(
        cantera_INCLUDE_DIR
        cantera_LIBRARY
)

# compatibility variables
set(cantera_VERSION_STRING ${cantera_VERSION})