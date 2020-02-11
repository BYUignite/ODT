
from __future__ import division
import os
import sys
import numpy as np
from data_py       import data_py
from basic_stats   import basic_stats
from fwhm          import fwhm, plot_fwhm_mixf, plot_fwhm_uvel
from eddyInfo      import eddyMaps
from eddyInfo      import eddyStats
from contourPlot   import contourPlot
from cl            import cl, plot_cl
from plot_RorZprof import *
from fv            import *
from checkBadRuns  import *


# USER OPTIONS ############################################################################

# sanity checks
doIncompleteRlzCheck = True
doBlowoutCheck = True
extinctionTemp = 800    #K; sets temperature threshold for extinguished flame

# data and plots to generate
doBasicStats = True  
doEddyStats = True
doEddyMaps = True
doFWHM = False
doCenterline = True
doContourPlots = False

#############################################################################################

# caseName check
try :
    caseN = sys.argv[1]
except :
    raise ValueError("Include the case name in the call to driver.py")

# establish file paths
DI = {'pdir':'../../../data/'+caseN+'/post/',       \
      'ddir':'../../../data/'+caseN+'/data/',       \
      'cdir':'../../../data/'+caseN+'/',            \
      'cn':caseN}

if not os.path.exists(DI['pdir']):
    os.mkdir(DI['pdir'])

#--------------------------------------------------------------------------------------------

# check for incomplete/unusable data
if doIncompleteRlzCheck:
    checkForIncompleteRlzs(DI, caseN)

# copies data into .npy files
if not os.path.exists(DI['cdir']+"data_py/") :
    data_py(DI)

# check for extinguished cases
if doBlowoutCheck:
    checkForBlowout(DI, extinctionTemp)

# appends soot volume fraction data to .npy files if applicable
sootFlag = get_inputFileParameter(DI, ("params","Lsoot"))
if sootFlag:
    fv(DI)

#--------------------------------------------------------------------------------------------

# generates basic statistics
if doBasicStats:
    basic_stats(DI, nfbins=60, do_yt=True, soot=sootFlag)

# generates eddy statistics and figures
if doEddyStats:
    eddyStats(DI)

if doEddyMaps:
    eddyMaps(DI)

# generates FWHM data and figures
if doFWHM:
    fwhm(DI, profName="mixf")
    fwhm(DI, profName="uvel")
    plot_fwhm_uvel(DI)
    plot_fwhm_mixf(DI)

# generates centerline data and plots
if doCenterline:
    cl(DI, profName="temp")
    cl(DI, profName="fv")
    plot_cl(DI,'temp', 'exp/DCL_temp.Yave', [0,2500], [1,2])
    plot_cl(DI,'fv',   'exp/DCL_soot.Yave', [1e-4,1], [1,2])

#--------------------------------------------------------------------------------------------

expCol   = 2    # temperature
expFlist = ["exp/D055_T.txt", "exp/D070_T.txt", "exp/D085_T.txt", "exp/D100_T.txt", "exp/D115_T.txt", "exp/D125_T.txt", "exp/D130_T.txt", "exp/D135_T.txt", "exp/D145_T.txt", "exp/D160_T.txt", "exp/D175_T.txt", "exp/D190_T.txt", "exp/D205_T.txt", "exp/D220_T.txt", "exp/D235_T.txt"]
plot_RorZprof(DI, 'R', 'means', 'temp', np.array([55,70,85,100,115,125,130,135,145,160,175,190,205,220,235]), expFlist, expCol, [-20,20], [292,2200], [-20,-10,0,10,20], [500,1000,1500,2000], 'r/D', 'T (K)')

expCol   = 3    # standard deviation of temperature
expFlist = ["exp/D055_T.txt", "exp/D070_T.txt", "exp/D085_T.txt", "exp/D100_T.txt", "exp/D115_T.txt", "exp/D125_T.txt", "exp/D130_T.txt", "exp/D135_T.txt", "exp/D145_T.txt", "exp/D160_T.txt", "exp/D175_T.txt", "exp/D190_T.txt", "exp/D205_T.txt", "exp/D220_T.txt", "exp/D235_T.txt"]
plot_RorZprof(DI, 'R', 'sig', 'temp', np.array([55,70,85,100,115,125,130,135,145,160,175,190,205,220,235]), expFlist, expCol, [-20,20], [0,1000], [-20,-10,0,10,20], [0,300,600,900], 'r/D', 'T_std (K)')

if sootFlag:
    expCol   = 2    # soot weighted temperature
    expFlist = ["exp/D055_T.txt", "exp/D070_T.txt", "exp/D085_T.txt", "exp/D100_T.txt", "exp/D115_T.txt", "exp/D125_T.txt", "exp/D130_T.txt", "exp/D135_T.txt", "exp/D145_T.txt", "exp/D160_T.txt", "exp/D175_T.txt", "exp/D190_T.txt", "exp/D205_T.txt", "exp/D220_T.txt", "exp/D235_T.txt"]
    plot_RorZprof(DI, 'R', 'SWmeans', 'temp', np.array([55,70,85,100,115,125,130,135,145,160,175,190,205,220,235]), expFlist, expCol, [-20,20], [292,2200], [-20,-10,0,10,20], [500,1000,1500,2000], 'r/D', 'T (K)')

    expCol   = 2    # soot volume fraction
    expFlist = ["exp/D023_fv.txt", "exp/D039_fv.txt", "exp/D055_fv.txt", "exp/D070_fv.txt", "exp/D086_fv.txt", "exp/D102_fv.txt", "exp/D117_fv.txt", "exp/D133_fv.txt", "exp/D148_fv.txt", "exp/D164_fv.txt", "exp/D180_fv.txt", "exp/D195_fv.txt", "exp/D211_fv.txt", "exp/D227_fv.txt", "exp/D242_fv.txt"]
    plot_RorZprof(DI, 'R', 'means', 'fv', np.array([23,39,55,70,86,102,117,133,148,164,180,195,211,227,242]), expFlist, expCol, [-25,25], [1e-4,3], [-25,-15,0,15,25], [1e-4,1e-3,1e-2,1e-1,1], 'r/D', 'fv (ppm)')

    expCol   = 3    # soot volume fraction RMS
    expFlist = ["exp/D023_fv.txt", "exp/D039_fv.txt", "exp/D055_fv.txt", "exp/D070_fv.txt", "exp/D086_fv.txt", "exp/D102_fv.txt", "exp/D117_fv.txt", "exp/D133_fv.txt", "exp/D148_fv.txt", "exp/D164_fv.txt", "exp/D180_fv.txt", "exp/D195_fv.txt", "exp/D211_fv.txt", "exp/D227_fv.txt", "exp/D242_fv.txt"]
    plot_RorZprof(DI, 'R', 'sig', 'fv', np.array([23,39,55,70,86,102,117,133,148,164,180,195,211,227,242]), expFlist, expCol, [-25,25], [1e-4,3], [-25,-15,0,15,25], [1e-4,1e-3,1e-2,1e-1,1], 'r/D', 'fv_rms (ppm)')

#--------------------------------------------------------------------------------------------

if doContourPlots:

    figName = DI['pdir']+"contour_temp_"+DI['cn']
    levels = np.linspace(295,2500,50); ticks = [500,1000,1500,2000,2500]
    contourPlot(DI, 'temp', levels, ticks, figName, 'mixf', True, 'means')
    
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
    
