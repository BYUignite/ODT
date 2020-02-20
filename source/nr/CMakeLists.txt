target_include_directories(odt.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/interp_linear.cc  ${CMAKE_CURRENT_SOURCE_DIR}/interp_linear.h
)
