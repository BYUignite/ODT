include(FindPkgConfig)
pkg_check_modules(PC_Cantera QUIET Cantera)

find_path(Cantera_INCLUDE_DIR
        NAMES cantera Cantera
        PATHS ${PC_Cantera_INCLUDE_DIRS} ${ODT_CANTERA_INCLUDE_LOCATION}
        PATH_SUFFIXES Cantera
        )
find_library(Cantera_LIBRARY
        NAMES libcantera.a
        PATHS ${PC_Cantera_LIBRARY_DIRS} ${ODT_CANTERA_LIBRARY_LOCATION}
        )

set(Cantera_VERSION ${PC_Cantera_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cantera
        FOUND_VAR Cantera_FOUND
        REQUIRED_VARS
        Cantera_LIBRARY
        Cantera_INCLUDE_DIR
        VERSION_VAR Cantera_VERSION
        )

if(Cantera_FOUND)
    set(Cantera_LIBRARIES ${Cantera_LIBRARY})
    set(Cantera_INCLUDE_DIRS ${Cantera_INCLUDE_DIR})
    set(Cantera_DEFINITIONS ${PC_Cantera_CFLAGS_OTHER})
endif()

if(Cantera_FOUND AND NOT TARGET Cantera::Cantera)
    add_library(Cantera::Cantera UNKNOWN IMPORTED)
    set_target_properties(Cantera::Cantera PROPERTIES
            IMPORTED_LOCATION "${Cantera_LIBRARY}"
            INTERFACE_COMPILE_OPTIONS "${PC_Cantera_CFLAGS_OTHER}"
            INTERFACE_INCLUDE_DIRECTORIES "${Cantera_INCLUDE_DIR}"
            )
endif()

mark_as_advanced(
        Cantera_INCLUDE_DIR
        Cantera_LIBRARY
)

# compatibility variables
set(Cantera_VERSION_STRING ${Cantera_VERSION})