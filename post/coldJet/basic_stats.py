
from __future__ import division
import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, extrap1d, get_domainBounds
from data_tools import commentHdr

#--------------------------------------------------------------------------------------------
#COLDJET
def basic_stats(DI, nfbins=60, nx=-1, favre=False, do_yt=False, filter=True) :

    favre = False
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
    L         = get_inputFileParameter(DI, ("params", "domainLength"))
    Dx        = L/nx
    coflow    = get_inputFileParameter(DI, ("initParams", "vel_min"))   
    
    if nx%2 == 0 :
        nx=nx+1   # make it odd for centerline

    try :
        imixf = varNames.index("mixf")
    except :
        imixf = -1
    try :
        irho = varNames.index("rho")
    except :
        irho = -1
    try :
        idvisc = varNames.index("dvisc")
    except :
        idvisc = -1
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

    doConditional = True if imixf > -1 else False

    if irho == -1 and favre :
        raise ValueError("In basic_stats: favre is true, but there is no rho in data file")

    if do_yt and (irho == -1 or iuvel == -1):
        raise ValueError("In basic_stats: do_yt is true but missing rho or uvel in data file")

    #--------------------------------------------------------------------------------------------

    X      = np.linspace(x0,xL,nx)
    fbins  = np.linspace(df/2.0, 1.0-df/2.0, nfbins)

    # Add the turbulence dissipation rate statistics and its logarithm, so we have nvar+2
    means  = np.zeros((ntimes, nvar+2, nx))
    mean2  = np.zeros((ntimes, nvar+2, nx))
    rhoM   = np.zeros((ntimes, nvar+2, nx))         # needed for favre to normalize
    cmeans = np.zeros((ntimes,nvar,nfbins))
    cmean2 = np.zeros((ntimes,nvar,nfbins))
    binNrm = np.zeros((ntimes,nvar,nfbins))       # to normalize conditional means
    rhouu  = np.zeros(ntimes)
    rhou   = np.zeros(ntimes)

    for itime in range(ntimes) :
        
        # instead of using the get_data_realization() which reads each data_py*.npy again for each realization,
        # we read and hold the data for the dump time.
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

            print("basic_stats: Processing time # %i of %i; for realization %i of %i" %(itime+1, ntimes, irlz+1, nrlz))

            #data = get_data_realization(DI, itime, irlz)
            i_s = istarts[irlz]
            i_e = iends[irlz]
            
            data = data_all[i_s:i_e+1, :]
            x  = data[:,ipos]
            xf = data[:,iposf]
            xf = np.append(xf,x[-1]+(xf[-1]-x[-1]))

            if doConditional :
                dx = np.abs(np.abs(xf[1:])**cCoord - np.abs(xf[0:-1])**cCoord)
                i = np.where(xf[1:] * xf[0:-1] < 0)[0]
                dx[i] = np.abs(np.abs(xf[i+1])**cCoord + np.abs(xf[i])**cCoord)
                ibin = (data[:,imixf] / df).astype(int)
                ibin[np.where(ibin < 0)] = 0
                ibin[np.where(ibin > nfbins-1)] = nfbins-1
                wt = dx*data[:,irho] if favre else dx.copy()

            if favre :
                fextrap = extrap1d(x,data[:,irho])
                rho = fextrap(X)
                rhoM[itime,:,:] = rhoM[itime,:,:] + rho

            for ivar in range(nvar) :

                y = data[:,ivar]
                fextrap = extrap1d(x,y)
                Y = fextrap(X)
                Y2 = Y**2
                
                means[itime,ivar,:] = means[itime,ivar,:] + (Y* rho if favre else Y)
                mean2[itime,ivar,:] = mean2[itime,ivar,:] + (Y2*rho if favre else Y2)
                                                
                if doConditional :
                    cmeans[itime,ivar,ibin] = cmeans[itime,ivar,ibin] + y     * wt
                    cmean2[itime,ivar,ibin] = cmean2[itime,ivar,ibin] + y*y   * wt
                    binNrm[itime,ivar,ibin] = binNrm[itime,ivar,ibin] + wt

            # here compute the turbulent dissipation rate and other derived variables
            #  need to recompute dx using the cell centers and not area weighting
            dx = (x[1:] - x[0:-1])
            dU = ( data[1:,iuvel] - data[0:-1,iuvel] ) # ( Uvel[1:] - Uvel[0:-1] )
            dVisc = 0.5 * ( data[1:,idvisc] + data[0:-1,idvisc] )
            #i = np.where(( x[1:] * x[0:-1] < 0 ) and ( x[1:] > x[0:-1] ) )[0]   # center origin point
            i = np.where( dx < 0)[0] #crossing realizations
            dx[i] = dx[i-1]
            #dU = np.append(dU,0.0)  # append a null value for crossing to next rlzn
            #dVisc = np.append(dVisc,0.0)
            diss = dVisc * ( dU / dx )**2.0 # uses dx from cell centers and averaged dvisc
            diss = np.append(diss,0.0)
            fextrap = extrap1d(x,diss)
            yDiss = fextrap(X)
            yDiss2 = yDiss**2.0
            
            means[itime,nvar,:] = means[itime,nvar,:] + (yDiss* rho if favre else yDiss)
            mean2[itime,nvar,:] = mean2[itime,nvar,:] + (yDiss2*rho if favre else yDiss2)
            logDiss = np.log(np.maximum(yDiss,1e-20))
            logDiss2 = logDiss**2.0
            means[itime,nvar+1,:] = means[itime,nvar+1,:] + (logDiss* rho if favre else logDiss)
            mean2[itime,nvar+1,:] = mean2[itime,nvar+1,:] + (logDiss2*rho if favre else logDiss2)
            
            if do_yt :
                fextrap = extrap1d(x,data[:,irho])
                rho = fextrap(X) 
                fextrap = extrap1d(x,data[:,iuvel])
                uvel = fextrap(X)
                rhouu[itime] += np.sum(rho*uvel*uvel)
                rhou[itime]  += np.sum(rho*uvel)

        if favre: 
            means[itime,:,:] /=  rhoM
            mean2[itime,:,:] /=  rhoM
            
        else:
            means[itime,:,:] /= nrlz
            mean2[itime,:,:] /= nrlz
            

    sig2 = mean2 - means*means
    sig  = np.sqrt(np.abs(sig2))

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
            ypos = np.array(times[0:ntimes])
            dy   = ypos[1:]-ypos[0:-1]
            tpos = np.cumsum(dy/uavg_mid)
            tpos = np.insert(tpos,0,0.0)
        else:    # temporal
            print("setting ytu.dat for temporal case")
            tpos = np.array(times[0:ntimes])
            dt   = tpos[1:]-tpos[0:-1]
            ypos = np.cumsum(dt*uavg_mid)
            ypos = np.insert(ypos,0,0.0)
        ytu = np.vstack([ypos,tpos,uavg]).T
        np.savetxt(DI['pdir']+'ytu.dat', ytu, header=" y(m) t(s) u(m/s)", fmt="%15.8e ", comments=commentHdr)
           
    #--------------------------------------------------------------------------------------------
    #get jet width data to use for filter--this uses width from the fit line for the coldJet_base_gDens120 case, but could use a fit line from experimental data as long as the downstream distances match
    if filter:
        wfile = '/home/abaumga/odt2.0/post/coldJet/base_width.dat'
        #wfile = DI['pdir']+"fwhm_cl_uvel.dat"  #format same as fwhm_cl_uvel.dat
        jetWidth = np.loadtxt(wfile, comments=commentHdr)
        D = get_inputFileParameter(DI, ("initParams", "djeti"))
        width = np.zeros(ntimes)
        for itime in range(ntimes):
            width[itime] = (jetWidth[itime,1])*D

    #--------------------------------------------------------------------------------------------
    # write the mean and std data

    head = "x_(m)           "
    for i,time in enumerate(times) :
        hi = "loc_" + str(i+2) + "_" + str(time) 
        hi = hi + (22-len(hi))*" "
        head = head + hi

    if filter:
        wndws = np.zeros(ntimes)                                           
        for itime in range(ntimes):                                        
            wndws[itime] = (0.1*width[itime])/Dx #window size is 10% of the unfiltered fwhm value at that dumptime  

    mav_means = means
    mav_sig = sig
    if filter:
        for i in range(nvar+2): #nvar+2 so TKEdiss/logTKEdiss filtered                       
            for j in range(ntimes):                                                          
                wndw = int(wndws[j])                                                         
                for k in range(nx):                                                          
                    if k>=int(wndw/2):                                                       
                        mav_means[j,i,k] = np.mean(means[j,i,k-int(wndw/2):k+int(wndw/2)])   
                        mav_sig[j,i,k] = np.mean(sig[j,i,k-int(wndw/2):k+int(wndw/2)])       
               

    for i in range(nvar) :
        var = mav_means[:,i,:]
        var = np.reshape(var,(ntimes,nx))
        var = np.vstack([X,var]).T
        fname = DI['pdir'] + "means_" + varNames[i] + ".dat"
        np.savetxt(fname, var, header=head, fmt="%15.8e ", comments=commentHdr)
    
    for i in range(nvar) :
        var = mav_sig[:,i,:]
        var = np.reshape(var,(ntimes,nx))
        var = np.vstack([X,var]).T
        fname = DI['pdir'] + "sig_" + varNames[i] + ".dat"
        np.savetxt(fname, var, header=head, comments=commentHdr)

    var = mav_means[:,nvar,:]
    var = np.reshape(var,(ntimes,nx))
    var = np.vstack([X,var]).T
    fname = DI['pdir'] + "means_TKEdiss.dat"
    np.savetxt(fname, var, header=head, fmt="%15.8e ", comments=commentHdr)

    var = mav_means[:,nvar+1,:]
    var = np.reshape(var,(ntimes,nx))
    var = np.vstack([X,var]).T
    fname = DI['pdir'] + "means_logTKEdiss.dat"
    np.savetxt(fname, var, header=head, fmt="%15.8e ", comments=commentHdr)

    var = mav_sig[:,nvar,:]
    var = np.reshape(var,(ntimes,nx))
    var = np.vstack([X,var]).T
    fname = DI['pdir'] + "sig_TKEdiss.dat"
    np.savetxt(fname, var, header=head, fmt="%15.8e ", comments=commentHdr)

    var = mav_sig[:,nvar+1,:]
    var = np.reshape(var,(ntimes,nx))
    var = np.vstack([X,var]).T
    fname = DI['pdir'] + "sig_logTKEdiss.dat"
    np.savetxt(fname, var, header=head, fmt="%15.8e ", comments=commentHdr)

    

    #--------------------------------------------------------------------------------------------
    # write the conditional mean and std data

    if doConditional :

        head = "x_(m)           "
        for i,time in enumerate(times) :
            hi = "loc_" + str(i+2) + "_" + str(time) 
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
            
       
