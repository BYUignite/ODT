target_include_directories(sec.x PRIVATE .)

target_sources(sec.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/sundials_smalldense.c
        ${CMAKE_CURRENT_SOURCE_DIR}/sundials_nvector.c
        ${CMAKE_CURRENT_SOURCE_DIR}/sundials_math.c
        ${CMAKE_CURRENT_SOURCE_DIR}/sundials_dense.c
        ${CMAKE_CURRENT_SOURCE_DIR}/nvector_serial.c
        ${CMAKE_CURRENT_SOURCE_DIR}/cvode_io.c
        ${CMAKE_CURRENT_SOURCE_DIR}/cvode_dense.c
        ${CMAKE_CURRENT_SOURCE_DIR}/cvode.c
)
