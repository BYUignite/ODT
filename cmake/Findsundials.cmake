include(FindPkgConfig)
pkg_check_modules(PC_sundials QUIET sundials)

find_path(sundials_INCLUDE_DIR
        NAMES cantera
        PATHS ${PC_sundials_INCLUDE_DIRS} /usr/include
        PATH_SUFFIXES sundials
        )
find_library(sundials_LIBRARY
        NAMES libcantera.a
        PATHS ${PC_sundials_LIBRARY_DIRS} /usr/lib
        )

set(sundials_VERSION ${PC_sundials_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(sundials
        FOUND_VAR sundials_FOUND
        REQUIRED_VARS
        sundials_LIBRARY
        sundials_INCLUDE_DIR
        VERSION_VAR sundials_VERSION
        )

if(sundials_FOUND)
    set(sundials_LIBRARIES ${sundials_LIBRARY})
    set(sundials_INCLUDE_DIRS ${sundials_INCLUDE_DIR})
    set(sundials_DEFINITIONS ${PC_sundials_CFLAGS_OTHER})
endif()

if(sundials_FOUND AND NOT TARGET sundials::sundials)
    add_library(sundials::sundials UNKNOWN IMPORTED)
    set_target_properties(sundials::sundials PROPERTIES
            IMPORTED_LOCATION "${sundials_LIBRARY}"
            INTERFACE_COMPILE_OPTIONS "${PC_sundials_CFLAGS_OTHER}"
            INTERFACE_INCLUDE_DIRECTORIES "${sundials_INCLUDE_DIR}"
            )
endif()

mark_as_advanced(
        sundials_INCLUDE_DIR
        sundials_LIBRARY
)

# compatibility variables
set(sundials_VERSION_STRING ${sundials_VERSION})