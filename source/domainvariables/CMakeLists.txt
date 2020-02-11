target_include_directories(sec.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(sec.x
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/dv.cc                 ${CMAKE_CURRENT_SOURCE_DIR}/dv.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_pos.cc             ${CMAKE_CURRENT_SOURCE_DIR}/dv_pos.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_posf.cc            ${CMAKE_CURRENT_SOURCE_DIR}/dv_posf.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_rho.cc             ${CMAKE_CURRENT_SOURCE_DIR}/dv_rho.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_rho_const.cc       ${CMAKE_CURRENT_SOURCE_DIR}/dv_rho_const.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_rho_mf.cc          ${CMAKE_CURRENT_SOURCE_DIR}/dv_rho_mf.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_dvisc.cc           ${CMAKE_CURRENT_SOURCE_DIR}/dv_dvisc.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_dvisc_const.cc     ${CMAKE_CURRENT_SOURCE_DIR}/dv_dvisc_const.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_uvw.cc             ${CMAKE_CURRENT_SOURCE_DIR}/dv_uvw.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_sca.cc             ${CMAKE_CURRENT_SOURCE_DIR}/dv_sca.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_enth.cc            ${CMAKE_CURRENT_SOURCE_DIR}/dv_enth.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_temp.cc            ${CMAKE_CURRENT_SOURCE_DIR}/dv_temp.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_mixf.cc            ${CMAKE_CURRENT_SOURCE_DIR}/dv_mixf.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_chi.cc             ${CMAKE_CURRENT_SOURCE_DIR}/dv_chi.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_chi_dmf.cc         ${CMAKE_CURRENT_SOURCE_DIR}/dv_chi_dmf.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_hr.cc              ${CMAKE_CURRENT_SOURCE_DIR}/dv_hr.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas.cc            ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_noRxn.cc      ${CMAKE_CURRENT_SOURCE_DIR}/dv_ygas_noRxn.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_aDL.cc             ${CMAKE_CURRENT_SOURCE_DIR}/dv_aDL.h
)
