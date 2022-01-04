#################### TARGET SOURCES ##########################

target_include_directories(odt-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(CHEMISTRY STREQUAL "SIMPLEDLR")
    target_sources(odt-run PRIVATE simple_dlr.cc)
else()
    message(FATAL_ERROR "CHEMISTRY=${CHEMISTRY} not supported")
endif()