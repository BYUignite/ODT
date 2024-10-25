target_include_directories(odt.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt.x PRIVATE domaincase.cc
                             domaincase_odt_channel.cc
                             domaincase_odt_channelScalar.cc
                             domaincase_odt_isothermalWall.cc
                             domaincase_odt_RT.cc
                             domaincase_odt_coldPropaneJet.cc
                             domaincase_odt_coldJet.cc
                             domaincase_odt_jetFlame.cc
                             domaincase_odt_MFjetFlame.cc
                             domaincase_odt_jetMixlRxn.cc)
