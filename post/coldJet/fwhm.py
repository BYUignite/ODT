from __future__ import division
import numpy as np
from   data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt
from data_tools import commentHdr

#--------------------------------------------------------------------------------------------

def fwhm(DI, profName="mixf") :

    mfile = DI['pdir']+"means_" + profName + ".dat"
    sfile = DI['pdir']+"sig_"   + profName + ".dat"

    data = np.loadtxt(mfile, comments=commentHdr)
    x    = data[:,0]
    y    = data[:,1:]
    data = np.loadtxt(sfile, comments=commentHdr)
    ys   = data[:,1:]

    npts, ntimes = np.shape(data)
    ntimes = ntimes - 1

    width  = np.zeros(ntimes)
    cLine  = np.zeros(ntimes)
    scLine = np.zeros(ntimes)
    imid   = int(npts/2) + 1

    for it in range(ntimes) :
        print("fwhm: time %i of %i" %(it, ntimes))
        width[it] = local_fwhm(x,y[:,it])
        cLine[it] = y[imid,it]
        scLine[it] = ys[imid,it]

    times = get_inputFileParameter(DI, ("dumpTimes",))
    #D     = get_inputFileParameter(DI, ("initParams","d_f"))
    D     = get_inputFileParameter(DI, ("initParams","djeti"))
    times = np.array(times)[0:ntimes]

    if len(times) != ntimes :
        raise ValueError("fwhm.py: wrong number of times in data file and in input file.")

    data = np.vstack([times/D, width/D, cLine, scLine]).T

    head = "  x/D,               fwhm/D,         centerline,       rms_cL"
    fname = DI['pdir']+"fwhm_cl_" + profName + ".dat"
    np.savetxt(fname, data, header=head, fmt="%15.8e ", comments=commentHdr)

#--------------------------------------------------------------------------------------------

def local_fwhm(x,y) :

    ymin = 0.5*(y[0]+y[-1])

    m = np.max(y)
    im = np.where(y==m)[0]
    im = int((im[0]+im[-1])/2)
    hm = ymin + 0.5*(m-ymin)

    i = np.where(y>=hm)[0]

    ii = i[0]
    x1 = x[ii-1]
    x2 = x[ii]
    y1 = y[ii-1]
    y2 = y[ii]
    xlo = x1 + (x2-x1)/(y2-y1)*(hm-y1)

    ii = i[-1]
    x1 = x[ii+1]
    x2 = x[ii]
    y1 = y[ii+1]
    y2 = y[ii]
    xhi = x1 + (x2-x1)/(y2-y1)*(hm-y1)

    return xhi-xlo

#--------------------------------------------------------------------------------------------

def plot_fwhm_mixf(DI) :

    #------------ do mixf

    fname = DI['pdir']+"fwhm_cl_mixf.dat"
    odt = np.loadtxt(fname, comments=commentHdr)
    #exp_fwhm = np.loadtxt("SANDC3H8.JET/ray/raystat/proc/fwhm_mixf.dat")
    #exp_cl   = np.loadtxt("SANDC3H8.JET/ray/raystat/paxray.txt", comments="CC")

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()

    axL.plot(odt[:,0],odt[:,2],'k-',label='ODT cL')
    axL.plot(odt[:,0],odt[:,3],'b-',label='ODT cL rms')
    #axL.plot(exp_cl[:,0],exp_cl[:,3],'ko:',label='EXP cL')
    #axL.plot(exp_cl[:,0],exp_cl[:,5],'bs:',label='EXP cL rms')
    axL.set_ylabel('cL, rms', fontsize=22)
    axL.set_xlabel("x/D", fontsize=22)
    axL.legend(loc='upper right', frameon=False, fontsize=16)
    axL.set_ylim([0,1.2])

    axR = axL.twinx()

    axR.plot(odt[:,0],odt[:,1],'r-',label='ODT fwhm')
    #axR.plot(exp_fwhm[:,0],exp_fwhm[:,1],'r^:',label='EXP fwhm')
    axR.set_ylabel("FWHM/D", fontsize=22)
    axR.legend(loc='upper center', frameon=False, fontsize=16)
    axR.set_ylim([0,20])

    plt.savefig(DI['pdir']+"fwhm_mixf_" + DI['cn'].replace(".","o"))


#--------------------------------------------------------------------------------------------

def plot_fwhm_uvel(DI) :

    #------------ do uvel

    fname = DI['pdir']+"fwhm_cl_uvel.dat"
    odt = np.loadtxt(fname, comments=commentHdr)
    #exp_fwhm_air = np.loadtxt("SANDC3H8.JET/vel/velstat/proc/fwhm_air_u.dat")
    #exp_fwhm_jet = np.loadtxt("SANDC3H8.JET/vel/velstat/proc/fwhm_jet_u.dat")
    #exp_cl_air   = np.loadtxt("SANDC3H8.JET/vel/velstat/paxv.air.txt", comments="CC")   # 0 2 4
    #exp_cl_jet   = np.loadtxt("SANDC3H8.JET/vel/velstat/paxv.jet.txt", comments="CC")   # 0 2 4

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()

    axL.plot(odt[:,0],odt[:,2],'k-')
    axL.plot(odt[:,0],odt[:,3],'b-')
    #axL.plot(exp_cl_air[:,0],exp_cl_air[:,2],'ko:')
    #axL.plot(exp_cl_jet[:,0],exp_cl_jet[:,2],'ko:')
    #axL.plot(exp_cl_air[:,0],exp_cl_air[:,4],'bs:')
    #axL.plot(exp_cl_jet[:,0],exp_cl_jet[:,4],'bs:')
    axL.set_ylabel('cL, rms', fontsize=22)
    axL.set_xlabel("x/D", fontsize=22)
    axL.legend(("ODT cL", "ODT cL rms", "EXP cL", "EXP cL rms"), loc='upper right', frameon=False, fontsize=16)
    axL.set_ylim([0,80])

    axR = axL.twinx()

    axR.plot(odt[:,0],odt[:,1],'r-')
    #axR.plot(exp_fwhm_air[:,0],exp_fwhm_air[:,1],'r^:')
    #axR.plot(exp_fwhm_jet[:,0],exp_fwhm_jet[:,1],'r^:')
    axR.set_ylabel("FWHM/D", fontsize=22)
    #axR.legend(("ODT fwhm"),loc='upper center', frameon=False, fontsize=16)
    axR.legend(("ODT fwhm", "EXP fwhm"),loc='upper center', frameon=False, fontsize=16)
    axR.set_ylim([0,20])

    plt.savefig(DI['pdir']+"fwhm_uvel_" + DI['cn'].replace(".","o"))



