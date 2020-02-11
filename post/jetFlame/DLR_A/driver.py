
from __future__ import division
import os
import sys
import numpy as np
from data_py       import data_py
from basic_stats   import basic_stats
from fwhm          import fwhm, plot_fwhm_mixf, plot_fwhm_uvel
from eddyInfo      import eddyMaps, eddyStats
from contourPlot   import contourPlot
from plot_RorZprof import *
from cl            import cl, plot_cl
from checkBadRuns  import *

#--------------------------------------------------------------------------------------------

try :
    caseN = sys.argv[1]
except :
    raise ValueError("Include the case name in the call to driver.py")

DI = {'pdir':'../../../data/'+caseN+'/post/',       \
      'ddir':'../../../data/'+caseN+'/data/',       \
      'cdir':'../../../data/'+caseN+'/',            \
      'cn':caseN}

if not os.path.exists(DI['pdir']):
    os.mkdir(DI['pdir'])

#--------------------------------------------------------------------------------------------

checkForIncompleteRlzs(DI, caseN)

#--------------------------------------------------------------------------------------------

if not os.path.exists(DI['cdir']+"data_py/") :
    data_py(DI)

#--------------------------------------------------------------------------------------------

#checkForBlowout(DI)

#--------------------------------------------------------------------------------------------

basic_stats(DI, nfbins=60, do_yt=True)

#--------------------------------------------------------------------------------------------

eddyStats(DI)
eddyMaps(DI)

#--------------------------------------------------------------------------------------------

#fwhm(DI, profName="mixf")
#fwhm(DI, profName="uvel")
#plot_fwhm_uvel(DI)
#plot_fwhm_mixf(DI)

cl(DI, profName="mixf")
cl(DI, profName="uvel")
cl(DI, profName="temp")
plot_cl(DI,'uvel', 'exp/exp_v_cl.dat', [0,60],  [1,2])
plot_cl(DI,'temp', 'exp/exp_cl.dat',   [0,2500],[3,4])
plot_cl(DI,'mixf', 'exp/exp_cl.dat',   [0,1.2], [1,2])

expCol   = 1
expFlist = ['exp/exp_r_T_05.dat', 'exp/exp_r_T_10.dat', 'exp/exp_r_T_20.dat', 'exp/exp_r_T_40.dat', 'exp/exp_r_T_60.dat', 'exp/exp_r_T_80.dat' ]
plot_RorZprof(DI, 'R', 'means', 'temp', np.array([5,10,20,40,60,80]), expFlist, expCol, [-20,20], [292,2000], [-20,-10,0,10,20], [500,1000,1500,2000], 'r/D', 'T (K)')
expCol   = 2
expFlist = ['exp/exp_r_T_05.dat', 'exp/exp_r_T_10.dat', 'exp/exp_r_T_20.dat', 'exp/exp_r_T_40.dat', 'exp/exp_r_T_60.dat', 'exp/exp_r_T_80.dat' ]
plot_RorZprof(DI, 'R', 'sig', 'temp', np.array([5,10,20,40,60,80]), expFlist, expCol, [-20,20], [0,500], [-20,-10,0,10,20], [100,200,300,400,500], 'r/D', 'T (K)')
expCol   = 1
expFlist = ['exp/exp_r_Z_05.dat', 'exp/exp_r_Z_10.dat', 'exp/exp_r_Z_20.dat', 'exp/exp_r_Z_40.dat', 'exp/exp_r_Z_60.dat', 'exp/exp_r_Z_80.dat' ]
plot_RorZprof(DI, 'R', 'means', 'mixf', np.array([5,10,20,40,60,80]), expFlist, expCol, [-20,20], [0,1], [-20,-10,0,10,20], [0.25,0.5,0.75,1.0], 'r/D', 'mixf')
expCol   = 2
expFlist = ['exp/exp_r_Z_05.dat', 'exp/exp_r_Z_10.dat', 'exp/exp_r_Z_20.dat', 'exp/exp_r_Z_40.dat', 'exp/exp_r_Z_60.dat', 'exp/exp_r_Z_80.dat' ]
plot_RorZprof(DI, 'R', 'sig', 'mixf', np.array([5,10,20,40,60,80]), expFlist, expCol, [-20,20], [0,0.25], [-20,-10,0,10,20], [0.05,0.10,0.20,0.25], 'r/D', 'mixf')

expCol = 2
expFlist = ['exp/DLRA2_05.Ycnd', 'exp/DLRA2_10.Ycnd', 'exp/DLRA2_20.Ycnd', 'exp/DLRA2_40.Ycnd', 'exp/DLRA2_60.Ycnd', 'exp/DLRA2_80.Ycnd']
plot_RorZprof(DI, 'Z', 'cmeans', 'temp', np.array([5,10,20,40,60,80]), expFlist, expCol, [0,1], [292,2500], [0,0.2,0.4,0.6,0.8,1], [500,1000,1500,2000,2500], 'mixf', 'T (K)')

#--------------------------------------------------------------------------------------------

figName = DI['pdir']+"contour_temp_"+DI['cn']
levels = np.linspace(290,2000,50); ticks = [500,1000,1500,2000]
contourPlot(DI, 'temp', levels, ticks, figName, 'mixf', False, 'means')

figName = DI['pdir']+"contour_temp_rms_"+DI['cn']
levels = np.linspace(0,800,50); ticks = [0,200,400,600,800]
contourPlot(DI, 'temp', levels, ticks, figName, 'mixf', True, 'sig')

figName = DI['pdir']+"contour_mixf_"+DI['cn']
levels = np.linspace(-0.01,0.3,50); ticks = [0,0.1,0.2,0.3,0.4]
contourPlot(DI, 'mixf', levels, ticks, figName, 'mixf', True, 'means')

figName = DI['pdir']+"contour_rho_"+DI['cn']
levels = np.linspace(0.0,1.2,50); ticks = [0,0.4,0.8,1.2]
contourPlot(DI, 'rho', levels, ticks, figName, 'mixf', True, 'means')

figName = DI['pdir']+"contour_uvel_"+DI['cn']
levels = np.linspace(0.0,30,50); ticks = [0,10,20,30]
contourPlot(DI, 'uvel', levels, ticks, figName, 'mixf', True, 'means')
