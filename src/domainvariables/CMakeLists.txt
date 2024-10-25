target_include_directories(odt.x PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt.x PRIVATE dv.cc
                             dv_pos.cc
                             dv_posf.cc
                             dv_rho.cc
                             dv_rho_const.cc
                             dv_rho_mf.cc
                             dv_dvisc.cc
                             dv_dvisc_const.cc
                             dv_uvw.cc
                             dv_sca.cc
                             dv_enth.cc
                             dv_temp.cc
                             dv_mixf.cc
                             dv_chi.cc
                             dv_chi_dmf.cc
                             dv_hr.cc
                             dv_ygas.cc
                             dv_ygas_noRxn.cc
                             dv_aDL.cc)
