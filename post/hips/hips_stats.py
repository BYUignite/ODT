from __future__ import division
import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, extrap1d, get_domainBounds

#--------------------------------------------------------------------------------------------

def hips_stats(DI, LSR=False) :

    # LSR : True for simple reaction, False otherwise

    #--------------------------------------------------------------------------------------------

    dataFiles = glob.glob(DI['cdir']+"data_py/data_*.npy")
    ntimes    = len(dataFiles)
    varNames  = get_dataHeaderVars(DI)
    times     = get_inputFileParameter(DI, ("dumpTimes",))         # times or ypositions if spatial
    nvar      = len(varNames)
    nLevels   = get_inputFileParameter(DI, ("params", "nLevels"))    
    nx = int(2**(nLevels-1))
    nrlz      = get_nRlz(DI, nx)


    #--------------------------------------------------------------------------------------------

    if LSR:
        iR = varNames.index("y_R")
        iP = varNames.index("y_P")

    #--------------------------------------------------------------------------------------------

    X      = np.arange(nx) + 1

    means  = np.zeros((ntimes, nvar, nx))
    mean2  = np.zeros((ntimes, nvar, nx))
    mT     = np.zeros((ntimes,nvar))
    mT2    = np.zeros((ntimes,nvar))

    if LSR:
        meansS = np.zeros((ntimes,nx))
        mean2S = np.zeros((ntimes,nx))
        mTS = np.zeros(ntimes)
        mT2S = np.zeros(ntimes)

    for itime in range(ntimes) :
        for irlz in range(nrlz) :

            print("Processing time # %i of %i; for realization %i of %i" %(itime+1, ntimes, irlz, nrlz))

            data = get_data_realization(DI, itime, irlz, nx)

            if LSR:
                y_R = data[:,iR]
                y_P = data[:,iP]
                S   = (y_R+1E-12)/(y_R+y_P+1.0E-12)

                meansS[itime,:] = meansS[itime,:] + S
                mean2S[itime,:] = mean2S[itime,:] + S**2
                mTS[itime]  = mTS[itime] + np.sum(S)
                mT2S[itime] = mTS[itime] + np.sum(S**2)

            for ivar in range(nvar) :
                Y = data[:,ivar]
                Y2 = Y**2.0
                means[itime,ivar,:] = means[itime,ivar,:] + Y
                mean2[itime,ivar,:] = mean2[itime,ivar,:] + Y2
                mT[itime,ivar]  = mT[itime,ivar] + np.sum(Y)
                mT2[itime,ivar] = mT2[itime,ivar] + np.sum(Y2)


    means /= nrlz
    mean2 /= nrlz
    mT    /= nrlz*nx
    mT2   /= nrlz*nx

    if LSR:
        meansS /= nrlz
        mean2S /= nrlz
        mTS    /= nrlz*nx
        mT2S   /= nrlz*nx

    sig2 = mean2 - means*means
    sig  = np.sqrt(np.abs(sig2))
    sigT2 = mT2 - mT*mT
    sigT = np.sqrt(np.abs(sigT2))

    if LSR:
        sigS2 = mean2S - meansS*meansS
        sigS  = np.sqrt(np.abs(sigS2))
        sigTS2 = mT2S - mTS*mTS
        sigTS = np.sqrt(np.abs(sigTS2))

    #--------------------------------------------------------------------------------------------
    # write the mean and std data

    head = "iParcel         "
    for i,time in enumerate(times) :
        hi = str(i+2) + "_" + str(time)
        hi = hi + (17-len(hi))*" "
        head = head + hi

    for i in range(nvar) :
        var = means[:,i,:]
        var = np.reshape(var,(ntimes,nx))
        var = np.vstack([X,var]).T
        fname = DI['pdir'] + "means_" + varNames[i] + ".dat"
        np.savetxt(fname, var, header=head, fmt="%15.8e ")

    for i in range(nvar) :
        var = sig[:,i,:]
        var = np.reshape(var,(ntimes,nx))
        var = np.vstack([X,var]).T
        fname = DI['pdir'] + "sig_" + varNames[i] + ".dat"
        np.savetxt(fname, var, header=head)

    if LSR:
        var = np.vstack([X,meansS]).T
        fname = DI['pdir'] + "means_S.dat"
        np.savetxt(fname, var, header=head, fmt="%15.8e ")

        var = np.vstack([X,sigS]).T
        fname = DI['pdir'] + "sig_S.dat"
        np.savetxt(fname, var, header=head, fmt="%15.8e ")

    #--------------------------------------------------------------------------------------------
    # write the TOTAL mean and std data (that is, averaged over the domains)

    head = "times  "
    for i in range(len(varNames)):
        head += varNames[i] + "  "

    fname = DI['pdir'] + "TotMeans.dat"
    var = np.column_stack((times, mT))
    np.savetxt(fname, var, header=head, fmt="%15.8e ")

    fname = DI['pdir'] + "TotSig.dat"
    var = np.column_stack((times, sigT))
    np.savetxt(fname, var, header=head, fmt="%15.8e ")

    if LSR:
        fname = DI['pdir'] + "TotSMean.dat"
        var = np.column_stack((times, mTS))
        np.savetxt(fname, var, header="times  S", fmt="%15.8e ")

        fname = DI['pdir'] + "TotSSig.dat"
        var = np.column_stack((times, sigTS))
        np.savetxt(fname, var, header="times  S", fmt="%15.8e ")
    

