from __future__ import division
import numpy as np
from   data_tools import get_inputFileParameter, get_axialLocations
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------------------------

def fwhm(DI, profName="mixf") :

    mfile = DI['pdir']+"means_" + profName + ".dat"
    sfile = DI['pdir']+"sig_"   + profName + ".dat"

    data = np.loadtxt(mfile)
    x    = data[:,0]
    y    = data[:,1:]
    data = np.loadtxt(sfile)

    npts, ntimes = np.shape(data)
    ntimes = ntimes - 1

    width  = np.zeros(ntimes)

    for it in range(ntimes) :
        print("fwhm: time %i of %i" %(it, ntimes))
        width[it] = local_fwhm(x,y[:,it])

    times = get_axialLocations(DI, forceGetDumpTimes=False)
    D     = get_inputFileParameter(DI, ("initParams","d_f"))

    if len(times) != ntimes :
        raise ValueError("fwhm.py: wrong number of times in data file and in input file.")

    data = np.vstack([times/D, width/D]).T

    head = "  y/D,               fwhm/D"
    fname = DI['pdir']+"fwhm_" + profName + ".dat"
    np.savetxt(fname, data, header=head, fmt="%15.8e ")

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

    fname = DI['pdir']+"fwhm_mixf.dat"
    odt = np.loadtxt(fname)
    exp_fwhm = np.loadtxt("exp/fwhm_Z.dat")

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, ax = plt.subplots()

    ax.plot(odt[:,0],odt[:,1],'ro-',label='ODT')
    ax.plot(exp_fwhm[:,0],exp_fwhm[:,1],'r^:',label='EXP')
    ax.set_ylabel("FWHM/D", fontsize=22)
    ax.set_xlabel("y/D", fontsize=22)
    ax.legend(loc='upper center', frameon=False, fontsize=16)
    ax.set_ylim([0,20])

    plt.savefig(DI['pdir']+"fwhm_mixf_" + DI['cn'].replace(".","o"))

#--------------------------------------------------------------------------------------------

def plot_fwhm_uvel(DI) :

    #------------ do uvel

    fname = DI['pdir']+"fwhm_uvel.dat"
    odt = np.loadtxt(fname)
    exp_fwhm = np.loadtxt("exp/fwhm_u.dat")

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, ax = plt.subplots()

    ax.plot(odt[:,0],odt[:,1],'ro-')
    ax.plot(exp_fwhm[:,0],exp_fwhm[:,1],'r^:')
    ax.set_ylabel("FWHM/D", fontsize=22)
    ax.set_xlabel("y/D", fontsize=22)
    ax.legend(("ODT", "EXP"),loc='upper center', frameon=False, fontsize=16)
    ax.set_ylim([0,20])

    plt.savefig(DI['pdir']+"fwhm_uvel_" + DI['cn'].replace(".","o"))



