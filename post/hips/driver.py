
from __future__ import division
import os
import sys

from data_py    import data_py
from hips_stats import hips_stats
from energy_spectra_rlz  import energy_spectra_rlz
from energy_spectra_time import energy_spectra_time

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

#hips_stats(DI, LSR=False)
#hips_stats(DI, LSR=True)

#energy_spectra_rlz(DI)
energy_spectra_time(DI)


#--------------------------------------------------------------------------------------------

