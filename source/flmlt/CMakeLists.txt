target_include_directories(sec.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(sec.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_flmlt.cc     ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_flmlt.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_flmltX.cc    ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_flmltX.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_flmlt.cc        ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_flmlt.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_chi_flmlt.cc         ${CMAKE_CURRENT_SOURCE_DIR}/dv_chi_flmlt.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_enth_flmlt.cc        ${CMAKE_CURRENT_SOURCE_DIR}/dv_enth_flmlt.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_enth_flmltX.cc       ${CMAKE_CURRENT_SOURCE_DIR}/dv_enth_flmltX.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_flmlt.cc        ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_flmlt.h
        ${CMAKE_CURRENT_SOURCE_DIR}/solver_flmlt.cc         ${CMAKE_CURRENT_SOURCE_DIR}/solver_flmlt.h
        ${CMAKE_CURRENT_SOURCE_DIR}/micromixer_flmlt.cc     ${CMAKE_CURRENT_SOURCE_DIR}/micromixer_flmlt.h
                                                            ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_flmlt_LOGN.h
                                                            ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_flmlt_MOMIC.h
                                                            ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_flmlt_MONO.h
                                                            ${CMAKE_CURRENT_SOURCE_DIR}/dv_soot_flmlt_QMOM.h
)

#dv_soot_flmlt_CQMOM.h
