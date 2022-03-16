#################### TARGET SOURCES ##########################

target_include_directories(odt-run PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(odt-run PRIVATE chemMech.h
                               canteraRR.h
                               chemNone.h
                               simple_dlr.h     simple_dlr.cc
                               onestep_c2h4.h   onestep_c2h4.cc
                               onestep_ch4.h    onestep_ch4.cc
                               fourstep_ch4.h   fourstep_ch4.cc
                               ch4red.h         ch4red.cc
                               c2h4red.h        c2h4red.cc)