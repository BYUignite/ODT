#generates plots for multiple cases together
#see bottom section to specify output plots

from __future__ import division
import os
import sys
import numpy as np
from data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from data_tools import commentHdr

#plot functions
from uclu0_xD import uclu0_xD
from u0ucl_xD import u0ucl_xD
from wD_xD import wD_xD
from urmsucl_xD import urmsucl_xD
from uucl_rxx0 import uucl_rxx0
from urmsucl_rxx0 import urmsucl_rxx0
from logTKEdiss_xD import logTKEdiss_xD
from logTKEdissrms_xD import logTKEdissrms_xD
from logTKEdiss_rxx0 import logTKEdiss_rxx0
from CZLES_decay import CZLES_decay, C_decay, ZLES_decay
from CZLES_growth import CZLES_growth, C_growth, ZLES_growth
from CZLES_decay_x0 import CZLES_decay_x0, C_decay_x0, ZLES_decay_x0
from CZLES_growth_x0 import CZLES_growth_x0, C_growth_x0, ZLES_growth_x0
from intTKEdiss import intTKEdiss
from intTKEdiss_xD import intTKEdiss_xD
from TKEdiss_rxx0 import TKEdiss_rxx0, TKEdiss_rxx0_half
from TKEdiss_xD import TKEdiss_xD
from peak_urms_plot import peak_urms_plot

#------------------------------------------------------------------------

DI=[]


try:
    caseN1 = sys.argv[1]
except: 
    raise ValueError("Include case name in the call to plot.py")
DI.append({'pdir':'../../data/'+caseN1+'/post/',     \
       'ddir':'../../data/'+caseN1+'/data/',     \
       'cdir':'../../data/'+caseN1+'/',          \
       'cn':caseN1})
if not os.path.exists(DI[0]['pdir']):
    os.mkdir(DI[0]['pdir'])

for i in range(2,len(sys.argv)):

    caseN = sys.argv[i]
    DI.append({'pdir':'../../data/'+caseN+'/post/',     \
       'ddir':'../../data/'+caseN+'/data/',     \
       'cdir':'../../data/'+caseN+'/',          \
       'cn':caseN})
    if not os.path.exists(DI[i-1]['pdir']):
        os.mkdir(DI[i-1]['pdir'])

#-------------------------------------------------------------------------
##uncomment plot functions for desired plots--best not to generate all plots at once for run speed. Check plot function files for correct/desired figure directory/name before running to avoid overwriting previously made plots with different case combinations.

##velocity plots##
#uclu0_xD(DI)
#u0ucl_xD(DI)
#wD_xD(DI)
#urmsucl_xD(DI) 
#uucl_rxx0(DI)
#urmsucl_rxx0(DI)
#peak_urms_plot(DI) #use peak_urms.py to extract points, work in progress

##TKE and logTKE diss plots##
#logTKEdiss_xD(DI)
#logTKEdissrms_xD(DI)
#logTKEdiss_rxx0(DI)
#intTKEdiss(DI)
#TKEdiss_rxx0(DI)
#TKEdiss_rxx0_half(DI)
#intTKEdiss_xD(DI)

##C and ZLES parameter scatter plots##
#CZLES_decay(DI)
#C_decay(DI)
#ZLES_decay(DI)
#CZLES_growth(DI)
#C_growth(DI)
#ZLES_growth(DI)
#CZLES_decay_x0(DI)
#C_decay_x0(DI)
#ZLES_decay_x0(DI)
#CZLES_growth_x0(DI)
#C_growth_x0(DI)
#ZLES_growth_x0(DI)


