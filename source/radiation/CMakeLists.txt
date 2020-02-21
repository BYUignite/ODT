target_include_directories(sec.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(sec.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/radiation.cc           ${CMAKE_CURRENT_SOURCE_DIR}/radiation.h
        ${CMAKE_CURRENT_SOURCE_DIR}/radiationProperties.cc ${CMAKE_CURRENT_SOURCE_DIR}/radiationProperties.h
        ${CMAKE_CURRENT_SOURCE_DIR}/rad_opthin.cc          ${CMAKE_CURRENT_SOURCE_DIR}/rad_opthin.h
        ${CMAKE_CURRENT_SOURCE_DIR}/rad_twoflux.cc         ${CMAKE_CURRENT_SOURCE_DIR}/rad_twoflux.h
        ${CMAKE_CURRENT_SOURCE_DIR}/rad_fvdom.cc           ${CMAKE_CURRENT_SOURCE_DIR}/rad_fvdom.h
)
