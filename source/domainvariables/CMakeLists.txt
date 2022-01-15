#################### TARGET SOURCES ##########################

target_include_directories(odt-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt-run PRIVATE
        dv.cc                 dv.h
        dv_pos.cc             dv_pos.h
        dv_posf.cc            dv_posf.h
        dv_rho.cc             dv_rho.h
        dv_rho_const.cc       dv_rho_const.h
        dv_rho_mf.cc          dv_rho_mf.h
        dv_dvisc.cc           dv_dvisc.h
        dv_dvisc_const.cc     dv_dvisc_const.h
        dv_uvw.cc             dv_uvw.h
        dv_sca.cc             dv_sca.h
        dv_enth.cc            dv_enth.h
        dv_temp.cc            dv_temp.h
        dv_mixf.cc            dv_mixf.h
        dv_chi.cc             dv_chi.h
        dv_chi_dmf.cc         dv_chi_dmf.h
        dv_hr.cc              dv_hr.h
        dv_ygas.cc            dv_ygas.h
        dv_ygas_noRxn.cc      dv_ygas_noRxn.h
        dv_aDL.cc             dv_aDL.h
        dv_soot.cc            dv_soot.h
)