#################### TARGET SOURCES ##########################

target_include_directories(odt-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt-run PRIVATE
        radiation.cc            radiation.h
#        radiationProperties.cc  radiationProperties.h
        rad_opthin.cc           rad_opthin.h
        rad_twoflux.cc          rad_twoflux.h
        rad_fvdom.cc            rad_fvdom.h
)