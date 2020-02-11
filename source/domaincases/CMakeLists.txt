target_include_directories(sec.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(sec.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase.cc                      ${CMAKE_CURRENT_SOURCE_DIR}/domaincase.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_channel.cc          ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_channel.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_channelScalar.cc    ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_channelScalar.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_isothermalWall.cc   ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_isothermalWall.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_RT.cc               ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_RT.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_coldPropaneJet.cc   ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_coldPropaneJet.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_coldJet.cc          ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_coldJet.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_jetFlame.cc         ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_jetFlame.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_MFjetFlame.cc       ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_MFjetFlame.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_jetMixlRxn.cc       ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_odt_jetMixlRxn.h
)
