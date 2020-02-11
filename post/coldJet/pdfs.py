from __future__ import division
import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, get_domainBounds, compute_pdf, compute_wpdf
from data_tools import commentHdr

#--------------------------------------------------------------------------------------------

def get_pdfs(DI, favre=False, nbins=60):

    #--------------------------------------------------------------------------------------------
    
    dataFiles = glob.glob(DI['cdir']+"data_py/data_*.npy")
    ntimes    = len(dataFiles)
    nrlz      = get_nRlz(DI)
    varNames  = get_dataHeaderVars(DI)
    times     = get_inputFileParameter(DI, ("dumpTimes",))         # times or ypositions if spatial
    times     = times[0:ntimes]    # limit times to number of dataFiles, not all listed dumpTimes
    nvar      = len(varNames)
    cCoord    = get_inputFileParameter(DI, ("params", "cCoord"))
    dxmin     = get_inputFileParameter(DI, ("params", "dxmin"))
    dxmax     = get_inputFileParameter(DI, ("params", "dxmax"))
    L         = get_inputFileParameter(DI, ("params", "domainLength"))
    umax      =  get_inputFileParameter(DI, ("initParams", "vel_max"))
    umin      =  get_inputFileParameter(DI, ("initParams", "vel_min"))
    vmin      = -0.05*np.abs(umax-umin)
    vmax      = 0.05*np.abs(umax-umin)
    wmin      = -0.05*np.abs(umax-umin)
    wmax      = 0.05*np.abs(umax-umin)
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
    try :
        irho = varNames.index("rho")
    except :
        irho = -1
    try :
        idvisc = varNames.index("dvisc")
    except :
        idvisc = -1
    try :
        imixf = varNames.index("mixf")
    except :
        imixf = -1
    try :
        iuvel = varNames.index("uvel")
    except :
        iuvel = -1
    try :
        ivvel = varNames.index("vvel")
    except :
        ivvel = -1
    try :
        iwvel = varNames.index("wvel")
    except :
        iwvel = -1

    if irho == -1 and favre :
        raise ValueError("In basic_stats: favre is true, but there is no rho in data file")

    dxpdfs = np.zeros((nbins, ntimes))

    P_uvel = np.zeros([nbins,ntimes])
    P_vvel = np.zeros([nbins,ntimes])
    P_wvel   = np.zeros([nbins,ntimes])
    P_diss = np.zeros([nbins,ntimes])
    P_logDiss = np.zeros([nbins,ntimes])
    P_diffU  =  np.zeros([nbins,ntimes])
    #P_logDiffUpos = np.zeros([nbins,ntimes])
    #P_logDiffUneg = np.zeros([nbins,ntimes])

    dumpTimesString = ''

    #--------------------------------------------------------------------------------------------
    
    for itime in range(ntimes) :
        
        dx = np.empty(0)

        dumpTimesString = dumpTimesString + '"P(t='+'{:.2e}'.format(times[itime])+' s)" '
        fname = DI['cdir']+"data_py/data_py_" + "{0:0>5}".format(itime) + ".npy"
        data = np.load(fname)

        print("Processing time # %i of %i" %(itime+1, ntimes))
    
        x  = data[:,ipos]
        xf = data[:,iposf]
        xf = np.append(xf,x[-1]+(xf[-1]-x[-1]))    
        rho = data[:,irho]
        uvel = data[:,iuvel]
        vvel = data[:,ivvel]
        wvel = data[:,iwvel]
        dvisc = data[:,idvisc]

        wt1 = np.abs(np.abs(xf[1:])**cCoord - np.abs(xf[0:-1])**cCoord)
        i = np.where(xf[1:] * xf[0:-1] < 0)[0]
        wt1[i] = np.abs(np.abs(xf[i+1])**cCoord + np.abs(xf[i])**cCoord) 
        # the crossing between realizations is computed correctly, so we don't need to mask end points
        #maskCrossing = (xf[1:] > xf[0:-1])
        wt = ( wt1*rho if favre else wt1 ) #  don't need this: (wt1*rho if favre else wt1)*maskCrossing
        
        #--------------
        j = np.where(  ( uvel - umin ) > 0.0001*(umax-umin) )
        P_uvel[:,itime] , uvel_bins = compute_wpdf(uvel[j], wt[j], umin, umax, nbins)
        P_vvel[:,itime] , vvel_bins = compute_wpdf(vvel[j], wt[j], vmin, vmax, nbins)
        P_wvel[:,itime] , wvel_bins = compute_wpdf(wvel[j], wt[j], wmin, wmax, nbins)

        # dissipation computation
        uvel  = np.append(uvel,uvel[-1])  # get same number of entries as dx
        dvisc = np.append(dvisc,dvisc[-1]) 
        dx = (xf[1:] - xf[0:-1])
        du = ( uvel[1:] - uvel[0:-1] )
        dvisc = 0.5 * ( dvisc[1:] + dvisc[0:-1] )
        diss = dvisc * ( du / dx )**2.0 # uses dx from cell centers and averaged dvisc
        #don't compute dissipation at crossing and don't count points with negligible dissipation
        j = np.where( np.logical_and( dx > 0.0 , diss > 10**(-6) ) )
        P_diss[:,itime] , diss_bins = compute_wpdf(diss[j], wt[j], 0, 1000, nbins)
        P_logDiss[:,itime] , logDiss_bins = compute_wpdf(np.log10(diss[j]), wt[j], -6, 8, nbins)
        #--------------

        #dxpdfs[:,itime], bins = compute_pdf(dx, dxmin, dxmax, nbins)
        dxpdfs[:,itime], bins = compute_pdf(np.log10(dx[j]), np.log10(dxmin), np.log10(dxmax), nbins)

        bdx = np.vstack([bins,dxpdfs.T]).T

    #--------------

#    head = " log(dx)[m]     "  + dumpTimesString
#    for i,time in enumerate(times) :
#        hi = "time_" + str(i+2) + "_" + str(time) 
#        hi = hi + (22-len(hi))*" "
#        head = head + hi

    var = np.vstack([uvel_bins,P_uvel.T]).T
    head = " u[m/s]   "  + dumpTimesString
    np.savetxt(DI['pdir']+'pdfs_uvel.dat', var, header=head, fmt="%15.8e ", comments=commentHdr)

    var = np.vstack([vvel_bins,P_vvel.T]).T
    head = " v[m/s]   "  + dumpTimesString
    np.savetxt(DI['pdir']+'pdfs_vvel.dat', var, header=head, fmt="%15.8e ", comments=commentHdr)

    var = np.vstack([wvel_bins,P_wvel.T]).T
    head = " w[m/s]   "  + dumpTimesString
    np.savetxt(DI['pdir']+'pdfs_wvel.dat', var, header=head, fmt="%15.8e ", comments=commentHdr)

    var = np.vstack([diss_bins,P_diss.T]).T
    head = " TKEdiss   "  + dumpTimesString
    np.savetxt(DI['pdir']+'pdfs_TKEdiss.dat', var, header=head, fmt="%15.8e ", comments=commentHdr)

    var = np.vstack([logDiss_bins,P_logDiss.T]).T
    head = " log10(TKEdiss)   "  + dumpTimesString
    np.savetxt(DI['pdir']+'pdfs_logTKEdiss.dat', var, header=head, fmt="%15.8e ", comments=commentHdr)

