#plots TKEdiss vs. r/(x-x0) and allows user to select x-axis ranges on that plot to search for peak values
#outputs file with downstream distance, left peak TKEdiss value, and right peak TKEdiss value
#run for one case at a time

from __future__ import division
import os
import sys
import numpy as np
from data_tools import get_inputFileParameter, commentHdr
import matplotlib
import matplotlib.pyplot as plt

#--------------------------------------------------------------------------------#click functions
def onclick(event):
        
    if event.dblclick:
        global pt
        pt = event.xdata
        return pt
#-------------------------------------------------------------------------------

#make plot
def TKEdiss_rxx0(DI,instructions,dist):

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})
    fig, axL = plt.subplots()

    mfile = DI['pdir']+"means_TKEdiss.dat"
    data = np.loadtxt(mfile, comments=commentHdr)
    times = get_inputFileParameter(DI, ("dumpTimes",))
    D = get_inputFileParameter(DI, ("initParams","djeti"))
    ufile = DI['pdir']+"uvel_cl.dat"
    ucl = np.loadtxt(ufile, comments=commentHdr)
    u = ucl[:,1]

    npts = len(data[:,0])
    ntimes = len(times)
    rnorm = np.empty((npts,ntimes))
    TKE = data[:,1:]
    
    rnorm[:,dist] = data[:,0]/(times[dist] - 4.0*D)

    axL.plot(rnorm[:,dist],TKE[:,dist]/((u[dist]**3)/times[dist]))
    axL.set_ylabel(r"$TKEdiss \bullet \frac{ x}{(u_{cL})^3}$", fontsize=22)
    axL.set_xlabel(r"$r/(x-x_0)$", fontsize=22)
    axL.set_xlim([-0.5,0.5])
    axL.set_title(instructions)

    data = np.vstack([rnorm[:,dist],(sig_uvel[:,dist]/Ucl[dist])]).T
    head = "  rnorm                TKEdiss x/ucL3"
    np.savetxt(DI['pdir']+"TKE_rxx0.dat", data, header=head, comments=commentHdr)

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    return pt

#-----------------------------------------------------------------------------------
#find peaks
def peaks(DI,pts,dist):
    
    file = DI['pdir']+"TKE_rxx0.dat"
    data = np.loadtxt(file)
    times = get_inputFileParameter(DI, ("dumpTimes",))
    rnorm = data[:,0]
    TKEdiss = data[:,1]
    idx1 = (np.abs(rnorm-pts[0]).argmin())
    idx2 = (np.abs(rnorm-pts[1]).argmin())
    print idx1,idx2
    leftPeak = max(TKEdiss[:idx1])
    rightPeak = max(TKEdiss[idx2:])
    print times[dist], leftPeak, rightPeak
    return ([times[dist], leftPeak, rightPeak])

#----------------------------------------------------------------------------------
#get peaks vs distance data
def all_peaks(DI):
    
    dist = np.array([3,6,9,12,15]) #indicate indices of downstream distance to find peaks for
    allPeaks = np.zeros((len(dist),3))

    for i in range(0,len(dist)):
        pt1 = TKE_rxx0(DI, "double click to right of left peak, then exit plot",dist[i])
        pt2 = TKE_rxx0(DI, "double click to left of right peak, then exit plot",dist[i])
        pts = np.array([pt1,pt2])
        peak = peaks(DI,pts,dist[i])
        allPeaks[i,0] = peak[0]
        allPeaks[i,1] = peak[1]
        allPeaks[i,2] = peak[2]

    head = "  x(m)             left peak TKEdiss         right peak TKEdiss"
    np.savetxt(DI['pdir']+"TKEpeaks.dat", allPeaks, header=head, comments=commentHdr)
#-----------------------------------------------------------------------------------------
  
try :
    caseN = sys.argv[1]
except :
    raise ValueError("Include the case name in the call to peak_TKEdiss.py")

DI = {'pdir':'../../data/'+caseN+'/post/',       \
      'ddir':'../../data/'+caseN+'/data/',       \
      'cdir':'../../data/'+caseN+'/',            \
      'cn':caseN}

all_peaks(DI)
