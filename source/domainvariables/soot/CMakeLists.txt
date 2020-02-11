target_include_directories(sec.x PRIVATE .)

target_sources(sec.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/pdgen.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/wheeler.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot.cc          ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_MONO.cc     ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_MONO.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_QMOM.cc     ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_QMOM.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_LOGN.cc     ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_LOGN.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_MOMIC.cc    ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_MOMIC.h
        ${CMAKE_CURRENT_SOURCE_DIR}/eispack.cc          ${CMAKE_CURRENT_SOURCE_DIR}/eispack.h
)

#CQMOM.cc
#dv_soot_CQMOM.cc
#wheeler_ARCHES.cc
#vandermonde.cc
#dv_soot_CQMOM.h
