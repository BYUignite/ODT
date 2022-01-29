#################### TARGET SOURCES ##########################

target_include_directories(odt-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(ODT_CHEMISTRY STREQUAL "SIMPLEDLR")
    target_sources(odt-run PRIVATE simple_dlr.cc)
else()
    message(FATAL_ERROR "ODT_CHEMISTRY=${ODT_CHEMISTRY} not supported")
endif()