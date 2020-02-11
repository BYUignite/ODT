from __future__ import division
import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, get_domainBounds, compute_pdf
from data_tools import commentHdr

#--------------------------------------------------------------------------------------------

def get_dxpdfs(DI, nbins=60):

    #--------------------------------------------------------------------------------------------
    
    dataFiles = glob.glob(DI['cdir']+"data_py/data_*.npy")
    ntimes    = len(dataFiles)
    nrlz      = get_nRlz(DI)
    varNames  = get_dataHeaderVars(DI)
    times     = get_inputFileParameter(DI, ("dumpTimes",))         # times or ypositions if spatial
    times     = times[0:ntimes]    # limit times to number of dataFiles, not all listed dumpTimes
    nvar      = len(varNames)
    dxmin     = get_inputFileParameter(DI, ("params", "dxmin"))
    dxmax     = get_inputFileParameter(DI, ("params", "dxmax"))
    L         = get_inputFileParameter(DI, ("params", "domainLength"))
    x0, xL    = get_domainBounds(DI) 

    dxmin *= L
    dxmax *= L
    
    try :
        iposf = varNames.index("posf")
    except :
        raise ValueError("In basic_stats: no posf variable found")
    try :
        ipos = varNames.index("pos")
    except :
        raise ValueError("In basic_stats: no pos variable found")

    dxpdfs = np.zeros((nbins, ntimes))

    #--------------------------------------------------------------------------------------------
    
    for itime in range(ntimes) :
        
        dx = np.empty(0)

        fname = DI['cdir']+"data_py/data_py_" + "{0:0>5}".format(itime) + ".npy"
        data_all = np.load(fname)
        posf_all   = data_all[:,1]
        dx_all = posf_all[1:]-posf_all[0:-1]

        istarts = np.where(dx_all <0.0)[0] + 1
        istarts = np.insert(istarts, 0, 0.0)

        iends = np.where(dx_all < 0.0)[0]
        iends = np.append(iends, len(posf_all)-1)
        
        nrlz = len(istarts)  # some sims end early so compute nrlz for each time

        for irlz in range(nrlz) : 
    
            print("Processing time # %i of %i; for realization %i of %i" %(itime+1, ntimes, irlz+1, nrlz))
    
            #data = get_data_realization(DI, itime, irlz)
            i_s = istarts[irlz]
            i_e = iends[irlz]

            data = data_all[i_s:i_e+1, :]
            x  = data[:,ipos]
            xf = data[:,iposf]
            xf = np.append(xf,x[-1]+(xf[-1]-x[-1]))       
            dx = np.append(dx, (xf[1:] - xf[0:-1]))

        #--------------

        #dxpdfs[:,itime], bins = compute_pdf(dx, dxmin, dxmax, nbins)
        dxpdfs[:,itime], bins = compute_pdf(np.log10(dx), np.log10(dxmin), np.log10(dxmax), nbins)

        bdx = np.vstack([bins,dxpdfs.T]).T

    #--------------

    head = " log(dx)[m]     "
    for i,time in enumerate(times) :
        hi = "loc_" + str(i+2) + "_" + str(time) 
        hi = hi + (22-len(hi))*" "
        head = head + hi
    np.savetxt(DI['pdir']+'dxpdfs.dat', bdx, header=head, fmt="%15.8e ", comments=commentHdr)

