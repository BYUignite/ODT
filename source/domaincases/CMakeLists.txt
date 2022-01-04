#################### TARGET SOURCES ##########################

target_include_directories(odt-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt-run PRIVATE
        domaincase.cc                      domaincase.h
        domaincase_odt_channel.cc          domaincase_odt_channel.h
        domaincase_odt_channelScalar.cc    domaincase_odt_channelScalar.h
        domaincase_odt_isothermalWall.cc   domaincase_odt_isothermalWall.h
        domaincase_odt_RT.cc               domaincase_odt_RT.h
        domaincase_odt_coldPropaneJet.cc   domaincase_odt_coldPropaneJet.h
        domaincase_odt_coldJet.cc          domaincase_odt_coldJet.h
        domaincase_odt_jetFlame.cc         domaincase_odt_jetFlame.h
        domaincase_odt_MFjetFlame.cc       domaincase_odt_MFjetFlame.h
        domaincase_odt_jetMixlRxn.cc       domaincase_odt_jetMixlRxn.h
)