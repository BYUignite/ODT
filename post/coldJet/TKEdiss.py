#does centerline processing of TKEdiss

from __future__ import division
import numpy as np
from data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from data_tools import commentHdr

#---------------------------------------------------------------------------------------------

def TKEdiss(DI, profName="TKEdiss"): 

    mfile = DI['pdir']+"means_" + profName + ".dat"
    sfile = DI['pdir']+"sig_"   + profName + ".dat"

    data = np.loadtxt(mfile, comments=commentHdr)
    y    = data[:,1:]
    data = np.loadtxt(sfile, comments=commentHdr)
    ys   = data[:,1:]

    npts, ntimes = np.shape(data)
    ntimes = ntimes - 1

    cLine  = np.zeros(ntimes)
    scLine = np.zeros(ntimes)
    imid   = int(npts/2) + 1

    for it in range(ntimes) :
        print("TKEdiss: time %i of %i" %(it+1, ntimes))
        cLine[it] = y[imid,it]
        scLine[it] = ys[imid,it]

    times = get_inputFileParameter(DI, ("dumpTimes",))
    D     = get_inputFileParameter(DI, ("initParams","djeti"))
    times = np.array(times)

    if len(times) != ntimes :
        raise ValueError("TKEdiss.py: wrong number of times in data file and in input file.")

    data = np.vstack([times/D, cLine, scLine]).T

    head = "  x/D,               centerline,       rms_cL"
    fname = DI['pdir']+ profName + "_cl.dat"
    np.savetxt(fname, data, header=head, fmt="%15.8e ", comments=commentHdr)


