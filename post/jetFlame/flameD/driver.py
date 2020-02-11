
from __future__ import division
import os
import sys
import numpy as np
from data_py     import data_py
from basic_stats import basic_stats
from fwhm        import fwhm, plot_fwhm_mixf
from eddyInfo    import eddyMaps
from contourPlot import contourPlot
from plot_radial import *
from plot_conditional import *
from cl            import cl, plot_cl

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

if not os.path.exists(DI['cdir']+"data_py/") :
    data_py(DI)

#--------------------------------------------------------------------------------------------

basic_stats(DI, nfbins=60, do_yt=True)

#--------------------------------------------------------------------------------------------

eddyMaps(DI)

#--------------------------------------------------------------------------------------------

fwhm(DI, profName="mixf")
plot_fwhm_mixf(DI)


cl(DI, profName="mixf")
cl(DI, profName="temp")
plot_cl(DI,'mixf', 'exp/DCL.Yave', [0,1.2],  [1,2])
plot_cl(DI,'temp', 'exp/DCL.Yave', [0,2500], [3,4])


expCol   = 3
expFlist = ["exp/D075.Yave", "exp/D15.Yave", "exp/D30.Yave", "exp/D45.Yave", "exp/D60.Yave", "exp/D75.Yave"]
plot_RorZprof(DI, 'R', 'means', 'temp', np.array([7.5,15,30,45,60,75]), expFlist, expCol, [-20,20], [292,2000], [-20,-10,0,10,20], [500,1000,1500,2000], 'r/D', 'T (K)')
expCol   = 4
expFlist = ["exp/D075.Yave", "exp/D15.Yave", "exp/D30.Yave", "exp/D45.Yave", "exp/D60.Yave", "exp/D75.Yave"]
plot_RorZprof(DI, 'R', 'sig', 'temp', np.array([7.5,15,30,45,60,75]), expFlist, expCol, [-20,20], [0,800], [-20,-10,0,10,20], [0,200,400,600,800], 'r/D', 'T rms (K)')
expCol    = 1
expFlist = ["exp/D075.Yave", "exp/D15.Yave", "exp/D30.Yave", "exp/D45.Yave", "exp/D60.Yave", "exp/D75.Yave"]
plot_RorZprof(DI, 'R', 'means', 'mixf', np.array([7.5,15,30,45,60,75]), expFlist, expCol, [-20,20], [0,1], [-20,-10,0,10,20], [0,0.25,0.5,0.75,1], 'r/D', 'mixf')
expCol    = 2
expFlist = ["exp/D075.Yave", "exp/D15.Yave", "exp/D30.Yave", "exp/D45.Yave", "exp/D60.Yave", "exp/D75.Yave"]
plot_RorZprof(DI, 'R', 'sig', 'mixf', np.array([7.5,15,30,45,60,75]), expFlist, expCol, [-20,20], [0,1], [-20,-10,0,10,20], [0,0.25,0.5,0.75,1], 'r/D', 'mixf rms')



expCol   = 2
expFlist = ["exp/D075.Ycnd", "exp/D15.Ycnd", "exp/D30.Ycnd", "exp/D45.Ycnd", "exp/D60.Ycnd", "exp/D75.Ycnd"]
plot_RorZprof(DI, 'Z', 'cmeans', 'temp', np.array([7.5,15,30,45,60,75]), expFlist, expCol, [0,1], [292,2500], [0,0.2,0.4,0.6,0.8,1], [500,1000,1500,2000,2500], 'mixf', 'T (K)')
expCol   = 3
expFlist = ["exp/D075.Ycnd", "exp/D15.Ycnd", "exp/D30.Ycnd", "exp/D45.Ycnd", "exp/D60.Ycnd", "exp/D75.Ycnd"]
plot_RorZprof(DI, 'Z', 'csig', 'temp', np.array([7.5,15,30,45,60,75]), expFlist, expCol, [0,1], [0,800], [0,0.2,0.4,0.6,0.8,1], [0,200,400,600,800], 'mixf', 'T (K)')


#--------------------------------------------------------------------------------------------

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

