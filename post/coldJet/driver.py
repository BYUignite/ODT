
from __future__ import division
import os
import sys
import numpy as np
from data_py       import data_py
from basic_stats   import basic_stats
from pdfs          import get_pdfs
from uvel        import uvel, plot_uvel
from fwhm          import fwhm, plot_fwhm_uvel
from eddyInfo      import eddyMaps
from eddyInfo      import eddyStats
from grid_stats  import get_dxpdfs
#from cl            import cl, plot_cl
#from plot_RorZprof import *
from logTKEdiss import logTKEdiss
from TKEdiss import TKEdiss

#--------------------------------------------------------------------------------------------

try :
    caseN = sys.argv[1]
except :
    raise ValueError("Include the case name in the call to driver.py")

DI = {'pdir':'../../data/'+caseN+'/post/',       \
      'ddir':'../../data/'+caseN+'/data/',       \
      'cdir':'../../data/'+caseN+'/',            \
      'cn':caseN}

if not os.path.exists(DI['pdir']):
    os.mkdir(DI['pdir'])

#--------------------------------------------------------------------------------------------

if not os.path.exists(DI['cdir']+"data_py/") :
    data_py(DI)

#--------------------------------------------------------------------------------------------

basic_stats(DI, nfbins=60, do_yt=True, soot=False)

get_pdfs(DI, nbins=200)
#--------------------------------------------------------------------------------------------

eddyMaps(DI)

eddyStats(DI)

#--------------------------------------------------------------------------------------------

uvel(DI, profName="uvel")

plot_uvel(DI)

# check for differences here -- it's nice to have the width
fwhm(DI, profName="uvel")
plot_fwhm_uvel(DI)

#--------------------------------------------------------------------------------------------

logTKEdiss(DI, profName="logTKEdiss")
TKEdiss(DI, profName="TKEdiss")

#--------------------------------------------------------------------------------------------

get_dxpdfs(DI, nbins=60)

