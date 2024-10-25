target_include_directories(odt.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(CHEMISTRY STREQUAL "SIMPLEDLR")
    target_sources(odt.x PRIVATE simple_dlr.cc)
else()
    message("CHEMISTRY=${CHEMISTRY} not supported")
endif()

