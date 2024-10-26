target_include_directories(odt.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(CHEMISTRY STREQUAL "SIMPLEDLR")
    target_sources(odt.x PRIVATE simple_dlr.cc)
elseif(CHEMISTRY STREQUAL "P4_1987")
    target_sources(odt.x PRIVATE p4_1987.cc)
elseif(CHEMISTRY STREQUAL "P4")
    target_sources(odt.x PRIVATE p4.cc)
else()
    message("CHEMISTRY=${CHEMISTRY} not supported")
endif()

