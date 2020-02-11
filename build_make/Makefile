
ifndef $(MACH)
	MACH=MAC
	#MACH=OMEN
	#MACH=SNL
	#MACH=FSL
endif

include make.$(MACH)

####################   DIRECTORIES   ##########################################

SRCDIR       = ../source/
CVODE_SRCDIR = $(SRCDIR)cvode/
PROG         = ../run/sec.x

YAML_LIB = ../source/yaml/lib/
YAML_INC = ../source/yaml/include/

####################   COMPILATION   ##########################################

includes += -I$(SRCDIR) -I$(SRCDIR)flmlt -I$(SRCDIR)hips -I$(SRCDIR)domainvariables/
includes += -I$(SRCDIR)domaincases/
includes += -I$(SRCDIR)domainvariables/soot/
includes += -I$(SRCDIR)nr/
includes += -I$(SRCDIR)user_chemical_mechanisms/

####################   Cantera flags   ########################################

ext_libs += -lcantera
includes += -I$(CANT_INC)
libs     += -L$(CANT_LIB)

####################   YAML flags   ########################################

ext_libs += -lyaml-cpp
includes += -I$(YAML_INC) -I/usr/local/include
libs     += -L$(YAML_LIB)

####################   other flags   ##########################################
ifneq ($(CHEMISTRY), CANTERARR)
    D1        = -DPROBLEMSPECIFICRR
endif
ifeq ($(PARALLEL), YES)
    D2        = -DDOMPI
    includes += -I$(MPI_PATH)include/
endif
ifeq ($(SILENT), YES)
    D3        = -DSILENT
endif

DIRECTIVES = $(D1) $(D2) $(D3) -D_USESTDVECTOR_


COMPCC   = $(CC_COMPILER_NAME) $(includes) $(options) $(cppoptions) $(DIRECTIVES) -c
COMPC    = $(C_COMPILER_NAME)  $(includes) $(options) $(DIRECTIVES) -c

COMPDEP  = $(CC_COMPILER_NAME) $(includes) $(options) $(cppoptions) $(DIRECTIVES)  -MM -MT
COMPDEPC = $(C_COMPILER_NAME)  $(includes) $(options) $(DIRECTIVES)  -MM -MT

LINK     = $(CC_COMPILER_NAME) $(options) $(linkoptions) $(libs)

######################## SOURCE FILES ##########################

SRC   = $(SRCDIR)inputoutput.cc                         \
        $(SRCDIR)main.cc                                \
        $(SRCDIR)domain.cc                              \
        $(SRCDIR)processor.cc                           \
        $(SRCDIR)streams.cc                             \
        $(SRCDIR)meshManager.cc                         \
        $(SRCDIR)param.cc                               \
        $(SRCDIR)micromixer.cc                          \
        $(SRCDIR)eddy.cc                                \
        $(SRCDIR)solver.cc                              \
        $(SRCDIR)cvodeDriver.cc                         \
        $(SRCDIR)radiation/radiation.cc                 \
        $(SRCDIR)radiation/radiationProperties.cc       \
        $(SRCDIR)radiation/rad_opthin.cc                \
        $(SRCDIR)radiation/rad_twoflux.cc               \
        $(SRCDIR)radiation/rad_fvdom.cc                 \
        $(SRCDIR)domainvariables/dv.cc                    \
        $(SRCDIR)domainvariables/dv_pos.cc                \
        $(SRCDIR)domainvariables/dv_posf.cc               \
        $(SRCDIR)domainvariables/dv_rho.cc                \
        $(SRCDIR)domainvariables/dv_rho_const.cc          \
        $(SRCDIR)domainvariables/dv_rho_mf.cc             \
        $(SRCDIR)domainvariables/dv_dvisc.cc              \
        $(SRCDIR)domainvariables/dv_dvisc_const.cc        \
        $(SRCDIR)domainvariables/dv_uvw.cc                \
        $(SRCDIR)domainvariables/dv_sca.cc                \
        $(SRCDIR)domainvariables/dv_enth.cc               \
        $(SRCDIR)domainvariables/dv_temp.cc               \
        $(SRCDIR)domainvariables/dv_mixf.cc               \
        $(SRCDIR)domainvariables/dv_chi.cc                \
        $(SRCDIR)domainvariables/dv_chi_dmf.cc            \
        $(SRCDIR)domainvariables/dv_hr.cc                 \
        $(SRCDIR)domainvariables/dv_ygas.cc               \
        $(SRCDIR)domainvariables/dv_ygas_noRxn.cc         \
        $(SRCDIR)domainvariables/dv_aDL.cc                \
        $(SRCDIR)domainvariables/soot/dv_soot.cc          \
        $(SRCDIR)domainvariables/soot/dv_soot_MONO.cc     \
        $(SRCDIR)domainvariables/soot/dv_soot_QMOM.cc     \
        $(SRCDIR)domainvariables/soot/dv_soot_LOGN.cc     \
		$(SRCDIR)domainvariables/soot/dv_soot_MOMIC.cc    \
        $(SRCDIR)domainvariables/soot/pdgen.cc            \
        $(SRCDIR)domainvariables/soot/wheeler.cc          \
        $(SRCDIR)domainvariables/soot/eispack.cc          \
        $(SRCDIR)nr/interp_linear.cc                      \
        $(SRCDIR)domaincases/domaincase.cc                        \
        $(SRCDIR)domaincases/domaincase_odt_channel.cc            \
        $(SRCDIR)domaincases/domaincase_odt_channelScalar.cc      \
        $(SRCDIR)domaincases/domaincase_odt_isothermalWall.cc     \
        $(SRCDIR)domaincases/domaincase_odt_RT.cc                 \
        $(SRCDIR)domaincases/domaincase_odt_coldPropaneJet.cc     \
        $(SRCDIR)domaincases/domaincase_odt_coldJet.cc            \
        $(SRCDIR)domaincases/domaincase_odt_jetFlame.cc           \
        $(SRCDIR)domaincases/domaincase_odt_MFjetFlame.cc         \
        $(SRCDIR)domaincases/domaincase_odt_jetMixlRxn.cc         \
        $(SRCDIR)flmlt/domaincase_flmlt.cc         \
        $(SRCDIR)flmlt/dv_ygas_flmlt.cc         \
        $(SRCDIR)flmlt/dv_chi_flmlt.cc          \
        $(SRCDIR)flmlt/dv_enth_flmlt.cc         \
        $(SRCDIR)flmlt/dv_enth_flmltX.cc        \
        $(SRCDIR)flmlt/dv_soot_flmlt.cc         \
        $(SRCDIR)flmlt/solver_flmlt.cc          \
        $(SRCDIR)flmlt/micromixer_flmlt.cc      \
        $(SRCDIR)hips/dv_enth_hips.cc           \
        $(SRCDIR)hips/dv_mixf_hips.cc           \
        $(SRCDIR)hips/dv_ygas_hips.cc           \
        $(SRCDIR)hips/dv_ygas_cold_hips.cc      \
        $(SRCDIR)hips/micromixer_hips.cc        \
        $(SRCDIR)hips/solver_hips.cc            \
        $(SRCDIR)hips/domaincase_hips.cc           \
        $(SRCDIR)hips/domaincase_hips_comb.cc      \
        $(SRCDIR)hips/domaincase_hips_simpleRxn.cc	
#        $(SRCDIR)domainvariables/soot/dv_soot_CQMOM.cc    \
#        $(SRCDIR)domainvariables/soot/CQMOM.cc            \
#        $(SRCDIR)domainvariables/soot/vandermonde.cc      \


CVSRC = $(CVODE_SRCDIR)sundials_smalldense.c     \
        $(CVODE_SRCDIR)sundials_nvector.c        \
        $(CVODE_SRCDIR)sundials_math.c           \
        $(CVODE_SRCDIR)sundials_dense.c          \
        $(CVODE_SRCDIR)nvector_serial.c          \
        $(CVODE_SRCDIR)cvode_io.c                \
        $(CVODE_SRCDIR)cvode_dense.c             \
        $(CVODE_SRCDIR)cvode.c

#-------Add your rxn rate file here: make another ifeq block

ifeq ($(CHEMISTRY), ONESTEPC2H4)
    SRC += $(SRCDIR)user_chemical_mechanisms/onestep_c2h4.cc
endif
ifeq ($(CHEMISTRY), ONESTEPCH4)
    SRC += $(SRCDIR)user_chemical_mechanisms/onestep_ch4.cc
endif
ifeq ($(CHEMISTRY), C2H4RED)
    SRC += $(SRCDIR)user_chemical_mechanisms/c2h4red.cc
endif
ifeq ($(CHEMISTRY), CH4RED)
    SRC += $(SRCDIR)user_chemical_mechanisms/ch4red.cc
endif
ifeq ($(CHEMISTRY), HEPTRED)
    SRC += $(SRCDIR)user_chemical_mechanisms/nHept_LiuHewsonPitsch_red.cc
endif
ifeq ($(CHEMISTRY), SIMPLEDLR)
    SRC += $(SRCDIR)user_chemical_mechanisms/simple_dlr.cc
endif
ifeq ($(CHEMISTRY), ABFRED)
    SRC += $(SRCDIR)user_chemical_mechanisms/abf_red62.cc
endif
ifeq ($(CHEMISTRY), ARAMCORED)
    SRC += $(SRCDIR)user_chemical_mechanisms/aramco_red86_nonstiff.cc
endif
ifeq ($(CHEMISTRY), ARAMCORED_STIFF)
    SRC += $(SRCDIR)user_chemical_mechanisms/aramco_red86_stiff.cc
endif

######################## OBJECT FILES ##########################

OBJS   = $(SRC:.cc=.o)
CVOBJS = $(CVSRC:.c=.o)

######################## BUILD TARGETS ##########################

default : $(PROG)

all : depend default documentation

$(PROG) : $(OBJS) $(CVOBJS)
	$(LINK) $(OBJS) $(CVOBJS) -o $@ $(ext_libs)

%.o: %.cc
	$(COMPCC) $< -o $@
%.o: %.c
	$(COMPC)  $< -o $@

clean :
	rm -f $(OBJS) $(CVOBJS) ./*.mod ./.depend
	rm -f ../run/sec.x
	rm -f ../run/*.log ../run/*.xml
	rm -rf ../doc/doxygen/html/*
	rm -f ../source/user_chemical_mechanisms/*.o

documentation :
	rm -rf ../doc/doxygen/html
	mkdir ../doc/doxygen/html
	cp ../doc/doxygen_build_pages/*.png ../doc/doxygen/html/
	doxygen -s ../doc/Doxyfile




####################### OBJECT DEPENDENCIES #####################

DEPENDS   = $(SRC:.cc=.d)
CVDEPENDS = $(CVSRC:.c=.cvd)

depend :  $(DEPENDS) $(CVDEPENDS)
	cat   $(DEPENDS) $(CVDEPENDS) > .depend
	rm -f $(DEPENDS) $(CVDEPENDS)

%.d :
	$(COMPDEP) $*.o $*.cc > $*.d
%.cvd :
	$(COMPDEPC) $*.o $*.c > $*.cvd

#ifeq ($(wildcard .depend), .depend)
-include .depend
#endif

######################## BUILD YAML  ##########################

yaml :
	./build_yaml.sh

