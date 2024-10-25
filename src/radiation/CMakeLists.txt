target_include_directories(odt.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt.x PRIVATE radiation.cc
                             radiationProperties.cc
                             rad_opthin.cc
                             rad_twoflux.cc
                             rad_fvdom.cc)
