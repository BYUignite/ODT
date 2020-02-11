#plots normalized rms velocity fluctuations vs. r/(x-x0) and allows user to select x-axis ranges on that plot to search for peak values
#outputs file with downstream distance, left peak urms/ucl value, and right peak urms/ucl value
#run for one case at a time

from __future__ import division
import os
import sys
import numpy as np
from data_tools import get_inputFileParameter, commentHdr
import matplotlib
import matplotlib.pyplot as plt
#-----------------------------------------------------------------------------
#click functions
def onclick(event):
        
    if event.dblclick:
        global pt
        pt = event.xdata
        return pt
#--------------------------------------------------------------------------------
#make plot
def urmsucl_rxx0(DI,instructions,dist):

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})
    fig, axL = plt.subplots()

    mfile = DI['pdir']+"means_uvel.dat"
    data = np.loadtxt(mfile, comments=commentHdr)
    times = get_inputFileParameter(DI, ("dumpTimes",))
    ua = get_inputFileParameter(DI, ("initParams","vel_min"))
    D = get_inputFileParameter(DI, ("initParams","djeti"))
    sfile = DI['pdir']+"sig_uvel.dat"
    sig_uvel = np.loadtxt(sfile, comments=commentHdr)

    npts = len(data[:,0])
    ntimes = len(times)
    rnorm = np.empty((npts,ntimes))
    U = data[:,1:]
    icl = int(npts/2)+1
    Ucl = U[icl,:]
    rnorm[:,dist] = data[:,0]/(times[dist] - 4.0*D)

    axL.plot(rnorm[:,dist],(sig_uvel[:,dist]/Ucl[dist]))
    axL.set_ylabel(r"$v_{rms}/v_{cL}$", fontsize=22)
    axL.set_xlabel(r"$r/(x-x_0)$", fontsize=22)
    axL.set_xlim([-0.5,0.5])
    axL.set_title(instructions)

    data = np.vstack([rnorm[:,dist],(sig_uvel[:,dist]/Ucl[dist])]).T
    head = "  rnorm                urms/ucL"
    np.savetxt(DI['pdir']+"urmsucl_rxx0.dat", data, header=head, comments=commentHdr)

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    return pt

#-----------------------------------------------------------------------------------
#find peaks
def peaks(DI,pts,dist):
    
    file = DI['pdir']+"urmsucl_rxx0.dat"
    data = np.loadtxt(file)
    times = get_inputFileParameter(DI, ("dumpTimes",))
    rnorm = data[:,0]
    urmsucl = data[:,1]
    idx1 = (np.abs(rnorm-pts[0]).argmin())
    idx2 = (np.abs(rnorm-pts[1]).argmin())
    print idx1,idx2
    leftPeak = max(urmsucl[:idx1])
    rightPeak = max(urmsucl[idx2:])
    print times[dist], leftPeak, rightPeak
    return ([times[dist], leftPeak, rightPeak])

#----------------------------------------------------------------------------------
#get peaks vs distance data
def all_peaks(DI):
    
    dist = np.array([3,6,9,12,15]) #indicate indices of downstream distance to find peaks for
    allPeaks = np.zeros((len(dist),3))

    for i in range(0,len(dist)):
        pt1 = urmsucl_rxx0(DI, "double click to right of left peak, then exit plot",dist[i])
        pt2 = urmsucl_rxx0(DI, "double click to left of right peak, then exit plot",dist[i])
        pts = np.array([pt1,pt2])
        peak = peaks(DI,pts,dist[i])
        allPeaks[i,0] = peak[0]
        allPeaks[i,1] = peak[1]
        allPeaks[i,2] = peak[2]

    head = "  x(m)             left peak urms          right peak urms"
    np.savetxt(DI['pdir']+"urmspeaks.dat", allPeaks, header=head, comments=commentHdr)
#-----------------------------------------------------------------------------------------
  
try :
    caseN = sys.argv[1]
except :
    raise ValueError("Include the case name in the call to peak_urms.py")

DI = {'pdir':'../../data/'+caseN+'/post/',       \
      'ddir':'../../data/'+caseN+'/data/',       \
      'cdir':'../../data/'+caseN+'/',            \
      'cn':caseN}

all_peaks(DI)

