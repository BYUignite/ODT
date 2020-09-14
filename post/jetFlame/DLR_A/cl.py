from __future__ import division
import numpy as np
from   data_tools import get_inputFileParameter, get_axialLocations
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------------------------

def cl(DI, profName="mixf") :

    mfile = DI['pdir']+"means_" + profName + ".dat"
    sfile = DI['pdir']+"sig_"   + profName + ".dat"

    data = np.loadtxt(mfile)
    y    = data[:,1:]
    data = np.loadtxt(sfile)
    ys   = data[:,1:]

    npts, ntimes = np.shape(data)
    ntimes = ntimes - 1

    cLine  = np.zeros(ntimes)
    scLine = np.zeros(ntimes)
    imid   = int(npts/2) + 1

    for it in range(ntimes) :
        print("cl: time %i of %i" %(it, ntimes))
        cLine[it] = y[imid,it]
        scLine[it] = ys[imid,it]

    times = get_axialLocations(DI, forceGetDumpTimes=False)
    D     = get_inputFileParameter(DI, ("initParams","d_f"))

    if len(times) != ntimes :
        raise ValueError("cl.py: wrong number of times in data file and in input file.")

    data = np.vstack([times/D, cLine, scLine]).T

    head = "  y/D,               cl,       rms_cL"
    fname = DI['pdir']+"cl_" + profName + ".dat"
    np.savetxt(fname, data, header=head, fmt="%15.8e ")

#--------------------------------------------------------------------------------------------

def plot_cl(DI, varName, expFname, ylim, expCols) :

    fname = DI['pdir']+"cl_"+varName+".dat"
    odt = np.loadtxt(fname)
    exp_cl = np.loadtxt(expFname)

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, ax = plt.subplots()

    ax.plot(odt[:,0],odt[:,1],'k-',label='ODT')
    ax.plot(odt[:,0],odt[:,2],'k-')
    ax.plot(exp_cl[:,0],exp_cl[:,expCols[0]],'ko:',label='EXP')
    ax.plot(exp_cl[:,0],exp_cl[:,expCols[1]],'ko:')
#    ax.plot(odt[:,0],odt[:,1],'k-',label='ODT cL')
#    ax.plot(odt[:,0],odt[:,2],'b-',label='ODT cL rms')
#    ax.plot(exp_cl[:,0],exp_cl[:,expCols[0]],'ko:',label='EXP cL')
#    ax.plot(exp_cl[:,0],exp_cl[:,expCols[1]],'bo:',label='EXP cL rms')
#    ax.set_ylabel(varName+' cL, rms', fontsize=22)
    ax.set_xlabel("y/D", fontsize=22)
    ax.legend(loc='upper right', frameon=False, fontsize=16)
    ax.set_ylim(ylim)

    if varName == "uvel" : 
        ax.set_ylabel("v (m/s)", fontsize=22)
        ax.text(25, 35,  'CL',  fontsize=22, rotation=0)
        ax.text(25, 8,   'RMS', fontsize=22, rotation=0)
    if varName == "temp" : 
        ax.set_ylabel("T (K)", fontsize=22)
        ax.text(25, 1600,  'CL',  fontsize=22, rotation=0)
        ax.text(25, 300,   'RMS', fontsize=22, rotation=0)
    if varName == "mixf" : 
        ax.set_ylabel("$\\xi$", fontsize=22)
        ax.text(25, 0.65,  'CL',  fontsize=22, rotation=0)
        ax.text(25, 0.15,   'RMS', fontsize=22, rotation=0)
    
    plt.savefig(DI['pdir']+"cl_"+varName+"_"+DI['cn'].replace(".","o"))


#--------------------------------------------------------------------------------------------

