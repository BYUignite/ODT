from __future__ import division
import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, extrap1d, get_domainBounds
from data_tools import commentHdr
#--------------------------------------------------------------------------------------------

def basic_stats(DI, nfbins=60, nx=-1, favre=False, do_yt=False, soot=False) :

    #--------------------------------------------------------------------------------------------

    dataFiles = glob.glob(DI['cdir']+"data_py/data_*.npy")
    ntimes    = len(dataFiles)
    nrlz      = get_nRlz(DI)
    varNames  = get_dataHeaderVars(DI)
    times     = get_inputFileParameter(DI, ("dumpTimes",))         # times or ypositions if spatial
    times     = times[0:ntimes]    # limit times to number of dataFiles, not all listed dumpTimes
    nvar      = len(varNames)
    cCoord    = get_inputFileParameter(DI, ("params", "cCoord"))
    df        = 1.0/nfbins
    dxmin     = get_inputFileParameter(DI, ("params", "dxmin"))
    nx        = int(1.0/dxmin/10) if nx < 0 else nx;
    x0, xL    = get_domainBounds(DI)
    if nx%2 == 0 :
        nx=nx+1   # make it odd for centerline
    
    # Get the colflow velocity
    try: 
        uCoflow = get_inputFileParameter(DI, ("params", "U_a"))
    except :
        uCoflow = 0.0

    try :
        imixf = varNames.index("mixf")
    except :
        imixf = -1
    try :
        irho = varNames.index("rho")
    except :
        irho = -1
    try :
        iposf = varNames.index("posf")
    except :
        raise ValueError("In basic_stats: no posf variable found")
    try :
        ipos = varNames.index("pos")
    except :
        raise ValueError("In basic_stats: no pos variable found")
    try :
        iuvel = varNames.index("uvel")
    except :
        iuvel = -1
    try :
        iM1 = varNames.index("M1")
    except :
        iM1 = -1

    doConditional = True if imixf > -1 else False

    if irho == -1 and favre :
        raise ValueError("In basic_stats: favre is true, but there is no rho in data file")

    if iM1 == -1 and soot :
        raise ValueError("In basic_stats: soot is true, but there is no rho or M1 in data file")

    if do_yt and (irho == -1 or iuvel == -1):
        raise ValueError("In basic_stats: do_yt is true but missing rho or uvel in data file")

    #--------------------------------------------------------------------------------------------

    X        = np.linspace(x0,xL,nx)
    fbins    = np.linspace(df/2.0, 1.0-df/2.0, nfbins)

    means    = np.zeros((ntimes, nvar, nx))
    mean2    = np.zeros((ntimes, nvar, nx))
    SWmeans  = np.zeros((ntimes, nvar, nx))         # soot weighted means
    SWmean2  = np.zeros((ntimes, nvar, nx))         # soot weighted means
    rhoM     = np.zeros((ntimes, nvar, nx))         # needed for favre to normalize
    rhoS     = np.zeros((ntimes, nvar, nx))         # needed for soot weighted to normalize
    cmeans   = np.zeros((ntimes,nvar,nfbins))
    cmean2   = np.zeros((ntimes,nvar,nfbins))
    binNrm   = np.zeros((ntimes,nvar,nfbins))       # to normalize conditional means
    rhouu    = np.zeros(ntimes)
    rhou     = np.zeros(ntimes)

    for itime in range(ntimes) :
        for irlz in range(nrlz) :

            print("basic_stats: Processing time # %i of %i; for realization %i of %i" %(itime+1, ntimes, irlz+1, nrlz))

            try :
                data = get_data_realization(DI, itime, irlz)
                x  = data[:,ipos]
                xf = data[:,iposf]
                xf = np.append(xf,x[-1]+(xf[-1]-x[-1]))
                
                if doConditional :
                    # weigth samples by dx or annular area
                    dx = np.abs(np.abs(xf[1:])**cCoord - np.abs(xf[0:-1])**cCoord)
                    i = np.where(xf[1:] * xf[0:-1] < 0)[0]
                    dx[i] = np.abs(np.abs(xf[i+1])**cCoord + np.abs(xf[i])**cCoord)
                    # this is only for even bins of width df
                    ibin = (data[:,imixf] / df).astype(int)
                    ibin[np.where(ibin < 0)] = 0
                    ibin[np.where(ibin > nfbins-1)] = nfbins-1
                    wt = dx*data[:,irho] if favre else dx.copy()
                    
                if favre :
                    fextrap = extrap1d(x,data[:,irho])
                    rho = fextrap(X)
                    rhoM[itime,:,:] = rhoM[itime,:,:] + rho
                    
                if soot :
                    fextrap = extrap1d(x,data[:,irho])
                    rho = fextrap(X)
                    fextrap = extrap1d(x,data[:,iM1])
                    M1 = fextrap(X) * rho
                    rhoS[itime,:,:] = rhoS[itime,:,:] + M1
                    
                for ivar in range(nvar) :

                    y = data[:,ivar]
                    fextrap = extrap1d(x,y)
                    Y = fextrap(X)
                    Y2 = Y**2.0

                    means[itime,ivar,:] = means[itime,ivar,:] + (Y* rho if favre else Y)
                    mean2[itime,ivar,:] = mean2[itime,ivar,:] + (Y2*rho if favre else Y2)
                    
                    SWmeans[itime,ivar,:] = SWmeans[itime,ivar,:] + (Y* M1 if soot else Y)
                    SWmean2[itime,ivar,:] = SWmean2[itime,ivar,:] + (Y2*M1 if soot else Y2)

                    if doConditional :
                        cmeans[itime,ivar,ibin] = cmeans[itime,ivar,ibin] + y   * wt
                        cmean2[itime,ivar,ibin] = cmean2[itime,ivar,ibin] + y*y * wt
                        binNrm[itime,ivar,ibin] = binNrm[itime,ivar,ibin] + wt

                if do_yt :
                    fextrap = extrap1d(x,data[:,irho])
                    rho = fextrap(X)
                    fextrap = extrap1d(x,data[:,iuvel])
                    uvel = fextrap(X)
                    rhouu[itime] += np.sum(rho*uvel*uvel)
                    rhou[itime]  += np.sum(rho*uvel)

            except IndexError:
                print(" Missing realizations for itime = " + str(itime))

    if favre:
        means /=  rhoM
        mean2 /=  rhoM
    else:
        means /= nrlz
        mean2 /= nrlz

    if soot :
        SWmeans /=  rhoS
        SWmean2 /=  rhoS

    sig2 = mean2 - means*means
    sig  = np.sqrt(np.abs(sig2))

    SWsig2 = SWmean2 - SWmeans*SWmeans
    SWsig  = np.sqrt(np.abs(SWsig2))

    if doConditional :
        cmeans = cmeans / binNrm
        cmean2 = cmean2 / binNrm
        cmeans[np.where(binNrm==0)]=0

        csig2 = cmean2 - cmeans*cmeans
        csig  = np.sqrt(np.abs(csig2))
        csig[np.where(binNrm==0)]=0

    if do_yt:
        uavg     = rhouu/rhou
        uavg_mid = 0.5*(uavg[1:]+uavg[0:-1])
        Lspatial = get_inputFileParameter(DI, ("params", "Lspatial"))
        if Lspatial:
            print("setting ytu.dat for spatial case")
            ypos = np.array(times)
            dy   = ypos[1:]-ypos[0:-1]
            tpos = np.cumsum(dy/uavg_mid)
            tpos = np.insert(tpos,0,0.0)
        else:    # temporal
            print("setting ytu.dat for temporal case")
            tpos = np.array(times)
            dt   = tpos[1:]-tpos[0:-1]
            ypos = np.cumsum(dt*uavg_mid)
            ypos = np.insert(ypos,0,0.0)
        ytu = np.vstack([ypos,tpos,uavg]).T
        np.savetxt(DI['pdir']+'ytu.dat', ytu, header=" y(m) t(s) u(m/s)", fmt="%15.8e ", comments=commentHdr)
           
    #--------------------------------------------------------------------------------------------
    # write the mean and std data

    head = "x_(m)           "
    for i,time in enumerate(times) :
        hi = "time_" + str(i+2) + "_" + str(time) 
        hi = hi + (22-len(hi))*" "
        head = head + hi

    for i in range(nvar) :
        var = means[:,i,:]
        var = np.reshape(var,(ntimes,nx))
        var = np.vstack([X,var]).T
        fname = DI['pdir'] + "means_" + varNames[i] + ".dat"
        np.savetxt(fname, var, header=head, fmt="%15.8e ", comments=commentHdr)
    
    for i in range(nvar) :
        var = sig[:,i,:]
        var = np.reshape(var,(ntimes,nx))
        var = np.vstack([X,var]).T
        fname = DI['pdir'] + "sig_" + varNames[i] + ".dat"
        np.savetxt(fname, var, header=head, comments=commentHdr)

    #--------------------------------------------------------------------------------------------
    # write the conditional mean and std data

    if doConditional :

        head = "mixf_           "
        for i,time in enumerate(times) :
            hi = "time_" + str(i+2) + "_" + str(time) 
            hi = hi + (22-len(hi))*" "
            head = head + hi

        for i in range(nvar) :
            var = cmeans[:,i,:]
            var = np.reshape(var,(ntimes,nfbins))
            var = np.vstack([fbins,var]).T
            fname = DI['pdir'] + "cmeans_" + varNames[i] + ".dat"
            np.savetxt(fname, var, header=head, fmt="%15.8e ", comments=commentHdr)
        
        for i in range(nvar) :
            var = csig[:,i,:]
            var = np.reshape(var,(ntimes,nfbins))
            var = np.vstack([fbins,var]).T
            fname = DI['pdir'] + "csig_" + varNames[i] + ".dat"
            np.savetxt(fname, var, header=head, comments=commentHdr)
        
    #--------------------------------------------------------------------------------------------
    # write the soot weighted mean and std data

    if soot :

        head = "x_(m)           "
        for i,time in enumerate(times) :
            hi = "time_" + str(i+2) + "_" + str(time) 
            hi = hi + (22-len(hi))*" "
            head = head + hi

        for i in range(nvar) :
            var = SWmeans[:,i,:]
            var = np.reshape(var,(ntimes,nx))
            var = np.vstack([X,var]).T
            fname = DI['pdir'] + "SWmeans_" + varNames[i] + ".dat"
            np.savetxt(fname, var, header=head, fmt="%15.8e ", comments=commentHdr)

        for i in range(nvar) :
            var = SWsig[:,i,:]
            var = np.reshape(var,(ntimes,nx))
            var = np.vstack([X,var]).T
            fname = DI['pdir'] + "SWsig_" + varNames[i] + ".dat"
            np.savetxt(fname, var, header=head, fmt="%15.8e ", comments=commentHdr)
        
