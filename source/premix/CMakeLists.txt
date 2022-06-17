target_include_directories(odt-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt-run
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/solver_premix.cc  ${CMAKE_CURRENT_SOURCE_DIR}/solver_premix.h
        ${CMAKE_CURRENT_SOURCE_DIR}/dv_temp_profile.cc  ${CMAKE_CURRENT_SOURCE_DIR}/dv_temp_profile.h
        ${CMAKE_CURRENT_SOURCE_DIR}/micromixer_premix.cc  ${CMAKE_CURRENT_SOURCE_DIR}/micromixer_premix.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_premixedFlameBurner.cc  ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_premixedFlameBurner.h
        ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_premixedFlameBurner_fixed_T.cc  ${CMAKE_CURRENT_SOURCE_DIR}/domaincase_premixedFlameBurner_fixed_T.h
)