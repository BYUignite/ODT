target_include_directories(odt.x PRIVATE .)

target_sources(odt.x PRIVATE sundials_smalldense.c
                             sundials_nvector.c
                             sundials_math.c
                             sundials_dense.c
                             nvector_serial.c
                             cvode_io.c
                             cvode_dense.c
                             cvode.c)
