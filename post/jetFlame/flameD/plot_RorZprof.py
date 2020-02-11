
from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt
from data_tools import get_axialLocations, extrap1d, get_inputFileParameter

#----------------------------------------------------------

def plot_RorZprof(DI, RorZ, fileType, varName, yD, expFlist, expCol, xlim, ylim, xtics, ytics, xlabel, ylabel):


    D     = get_inputFileParameter(DI, ("initParams","d_f"))
    y     = yD*D

    if RorZ=='R':
        odt  = np.loadtxt(DI['pdir']+fileType+"_"+varName+".dat")
        scaleX = D
        figName = DI['pdir']+"radial_"+varName+"_" + DI['cn'].replace(".","o")
    else:
        odt  = np.loadtxt(DI['pdir']+fileType+"_"+varName+".dat")
        scaleX = 1
        figName = DI['pdir']+"conditional_"+varName+"_" + DI['cn'].replace(".","o")
    yodt = get_axialLocations(DI, forceGetDumpTimes=False)
    npts = len(odt[:,0])
    var  = np.zeros((npts, len(yD)))

    for i in range(len(yD)):
        for k in range(npts):
            fextrap = extrap1d(yodt, odt[k,1:])
            var[k,i] = fextrap(y[i])

    plt.figure()
    plt.rc("font", size=14)

    #for i,ifile in enumerate(expFlist):
    for i in range(len(yD)):
        plt.subplot(3,2,i+1)
        plt.plot(odt[:,0]/scaleX,var[:,i],'k-')
        if expFlist != []:
            exp = np.loadtxt(expFlist[i])
            plt.plot(exp[:,0],exp[:,expCol],'r-')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xticks([])
        plt.yticks([])
        plt.title("y/D = %i" %yD[i])

    plt.subplot(3,2,1)
    if expFlist == []:
        plt.legend(("ODT"), frameon=False)
    else:
        plt.legend(("ODT","EXP"), frameon=False)
    plt.yticks(ytics)
    plt.ylabel(ylabel)

    plt.subplot(3,2,3)
    plt.yticks(ytics)
    plt.ylabel(ylabel)

    plt.subplot(3,2,5)
    plt.yticks(ytics)
    plt.ylabel(ylabel)
    plt.xticks(xtics)
    plt.xlabel(xlabel)

    plt.subplot(3,2,6)
    plt.xticks(xtics)
    plt.xlabel(xlabel)

    #plt.show()

    plt.savefig(figName)

