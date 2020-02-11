from __future__ import division
import numpy as np
from   data_tools import get_inputFileParameter, get_axialLocations
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt
from data_tools import commentHdr

#--------------------------------------------------------------------------------------------

def cl(DI, soot=False, profName="mixf") :

    mfile = DI['pdir']+"means_"   + profName + ".dat"
    sfile = DI['pdir']+"sig_"     + profName + ".dat"

    data = np.loadtxt(mfile, comments=commentHdr)
    y    = data[:,1:]
    data = np.loadtxt(sfile, comments=commentHdr)
    ys   = data[:,1:]

    if soot :
        wfile = DI['pdir']+"SWmeans_" + profName + ".dat"
        zfile = DI['pdir']+"SWsig_"   + profName + ".dat"
        data = np.loadtxt(wfile, comments=commentHdr)
        yw   = data[:,1:]
        data = np.loadtxt(zfile, comments=commentHdr)
        yz   = data[:,1:]

    npts, ntimes = np.shape(data)
    ntimes = ntimes - 1

    cLine  = np.zeros(ntimes)
    scLine = np.zeros(ntimes)
    if soot:
        wLine  = np.zeros(ntimes)
        swLine = np.zeros(ntimes)
    imid   = int(npts/2) + 1

    for it in range(ntimes) :
        print("cl: time %i of %i" %(it, ntimes))
        cLine[it]  = y[imid,it]
        scLine[it] = ys[imid,it]
        if soot: 
            wLine[it]  = yw[imid,it]
            swLine[it] = yz[imid,it]

    times = get_axialLocations(DI, forceGetDumpTimes=False)
    times     = times[0:ntimes]    # limit times to number of dataFiles, not all listed dumpTimes
    D     = get_inputFileParameter(DI, ("initParams","d_f"))

    if len(times) != ntimes :
        raise ValueError("cl.py: wrong number of times in data file and in input file.")

    if soot:
        data = np.vstack([times/D, cLine, scLine, wLine, swLine]).T
        
        head = "  y/D,               cl,       rms_cL,      SW_cl,     SW_rms_cL"
        fname = DI['pdir']+"cl_" + profName + ".dat"
        np.savetxt(fname, data, header=head, fmt="%15.8e ", comments=commentHdr)

    else :
        data = np.vstack([times/D, cLine, scLine]).T
        
        head = "  y/D,               cl,       rms_cL"
        fname = DI['pdir']+"cl_" + profName + ".dat"
        np.savetxt(fname, data, header=head, fmt="%15.8e ", comments=commentHdr)

#--------------------------------------------------------------------------------------------

def plot_cl(DI, varName, expFname, ylim, expCols, soot=False) :

    fname = DI['pdir']+"cl_"+varName+".dat"
    odt = np.loadtxt(fname, comments=commentHdr)
    exp_cl = np.loadtxt(expFname)

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, ax = plt.subplots()

    ax.plot(odt[:,0],odt[:,1],'k-',label='ODT cL')
    ax.plot(odt[:,0],odt[:,2],'b-',label='ODT cL rms')
    if soot:
        ax.plot(odt[:,0],odt[:,3],'r-',label='ODT SW cL')
        ax.plot(odt[:,0],odt[:,4],'g-',label='ODT SW cL rms')
    ax.plot(exp_cl[:,0],exp_cl[:,expCols[0]],'ko:',label='EXP cL')
    ax.plot(exp_cl[:,0],exp_cl[:,expCols[1]],'bo:',label='EXP cL rms')
    ax.set_ylabel(varName+' cL, rms', fontsize=22)
    ax.set_xlabel("y/D", fontsize=22)
    ax.legend(loc='upper right', frameon=False, fontsize=16)
    if (varName == "fv") :
        ax.set_yscale('log')
    ax.set_ylim(ylim)

    plt.savefig(DI['pdir']+"cl_"+varName+"_"+DI['cn'].replace(".","o"))


#--------------------------------------------------------------------------------------------

