
from __future__ import division
import numpy as np
from   data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt
from data_tools import commentHdr

#--------------------------------------------------------------------------------------------

def uvel(DI, profName="uvel") :

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
        print("uvel: time %i of %i" %(it+1, ntimes))
        cLine[it] = y[imid,it]
        scLine[it] = ys[imid,it]

    times = get_inputFileParameter(DI, ("dumpTimes",))
    D     = get_inputFileParameter(DI, ("initParams","djeti"))
    times = np.array(times)

    if len(times) != ntimes :
        raise ValueError("uvel.py: wrong number of times in data file and in input file.")

    data = np.vstack([times/D, cLine, scLine]).T

    head = "  x/D,               centerline,       rms_cL"
    fname = DI['pdir']+ profName + "_cl.dat"
    np.savetxt(fname, data, header=head, fmt="%15.8e ", comments=commentHdr)

#--------------------------------------------------------------------------------------------

def plot_uvel(DI, profName='uvel') :

    fname = DI['pdir']+"uvel_cl.dat"
    odt = np.loadtxt(fname, comments=commentHdr)
    U0 = get_inputFileParameter(DI, ("initParams","vel_max"))
    fit = 1.0/5.8*(odt[:,0]-4.0)
    Exp = np.loadtxt("exp_u_cl.dat")
    ua = get_inputFileParameter(DI, ("initParams","vel_min"))

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()

    #axL.plot(odt[:,0],(U0-ua)/(odt[:,1]-ua),'ko-')
    axL.plot(odt[:,0],(U0-ua)/(odt[:,1]-ua),'k-')
    axL.plot(Exp[:,0],Exp[:,1],'^')
    axL.plot(odt[:,0],fit,'k--')
    axL.set_ylabel(r'$v_0/v_{cL}$', fontsize=22)
    axL.set_xlabel("y/D", fontsize=22)
    axL.legend(("ODT", "EXP", "Fit"), loc='upper right', frameon=False, fontsize=16)
    axL.set_ylim([0,25])

    plt.savefig(DI['pdir']+"uvel_cl_" + DI['cn'].replace(".","o"))

    #-----------------------------------------------

    D = get_inputFileParameter(DI, ("initParams","djeti"))
    Exp = np.loadtxt("exp_u_rad.dat")
    mfile = DI['pdir'] + "/means_" + profName + ".dat"
    data = np.loadtxt(mfile, comments=commentHdr)
    times = get_inputFileParameter(DI, ("dumpTimes",))
    ua = 0.0
    #ua = get_inputFileParameter(DI, ("initParams","vel_min"))

    npts = len(data[:,0])
    ntimes = len(times)
    rnorm = np.empty( (npts, ntimes) )
    U = data[:,1:]
    icl  = int(npts/2)+1
    Ucl = U[icl,:]
    for i in range(ntimes):
        rnorm[:,i] = data[:,0]/(times[i] - 4.0*D)

    plt.cla()
    #i = 8;  axL.plot(rnorm[:,i], (U[:,i]-ua)/(Ucl[i]-ua), 'k:')
    #i = 12; axL.plot(rnorm[:,i], (U[:,i]-ua)/(Ucl[i]-ua), 'b-.')
    #i = 16; axL.plot(rnorm[:,i], (U[:,i]-ua)/(Ucl[i]-ua), 'g--')
    #i = 20; axL.plot(rnorm[:,i], (U[:,i]-ua)/(Ucl[i]-ua), 'r-')
    i = 12; axL.plot(rnorm[:,i], (U[:,i]-ua)/(Ucl[i]-ua), 'k:')
    i = 16; axL.plot(rnorm[:,i], (U[:,i]-ua)/(Ucl[i]-ua), 'b--')
    i = 20; axL.plot(rnorm[:,i], (U[:,i]-ua)/(Ucl[i]-ua), 'r-')
    axL.plot(Exp[:,0], Exp[:,1], 'o')
    axL.set_ylim([0,1.2])
    axL.set_xlim([0,0.25])
    axL.set_xlabel(r"$r/(y-y_0)$", fontsize=22)
    axL.set_ylabel(r"$v/v_{cL}$", fontsize=22)
    #axL.legend(("y=30D", "y=50D", "y=70D", "y=90D", 'EXP'), frameon=False, fontsize=16)
    axL.legend(("y=50D", "y=70D", "y=90D", 'EXP'), frameon=False, fontsize=16)

    plt.savefig(DI['pdir']+"uvel_r_" + DI['cn'].replace(".","o"))

    #-----------------------------------------------

    mfile = DI['pdir']+"sig_uvel.dat"
    sig_uvel = np.loadtxt(mfile, comments=commentHdr)
    Exp = np.loadtxt("exp_uu_rad.dat")

    plt.cla()
    #i = 8;  axL.plot(rnorm[:,i], sig_uvel[:,i]/Ucl[i], 'k:')
    #i = 12; axL.plot(rnorm[:,i], sig_uvel[:,i]/Ucl[i], 'b-.')
    #i = 16; axL.plot(rnorm[:,i], sig_uvel[:,i]/Ucl[i], 'g--')
    #i = 20; axL.plot(rnorm[:,i], sig_uvel[:,i]/Ucl[i], 'r-')
    i = 12; axL.plot(rnorm[:,i], sig_uvel[:,i]/Ucl[i], 'k:')
    i = 16; axL.plot(rnorm[:,i], sig_uvel[:,i]/Ucl[i], 'b--')
    i = 20; axL.plot(rnorm[:,i], sig_uvel[:,i]/Ucl[i], 'r-')
    axL.plot(Exp[:,0], Exp[:,1]**0.5, 'o', ms=4)
    axL.set_ylim([0,0.6])
    axL.set_xlim([0,0.25])
    axL.set_xlabel(r"$r/(y-y_0)$", fontsize=22)
    axL.set_ylabel(r"$v_{rms}/v_{cL}$", fontsize=22)
    #axL.legend(("y=30D", "y=50D", "y=70D", "y=90D", 'EXP'), frameon=False, fontsize=16)
    axL.legend(("y=50D", "y=70D", "y=90D", 'EXP'), frameon=False, fontsize=16)

    plt.savefig(DI['pdir']+"uu_r_" + DI['cn'].replace(".","o"))

