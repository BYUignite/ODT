
from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt
from data_tools import get_axialLocations, extrap1d, get_inputFileParameter
from data_tools import commentHdr

#----------------------------------------------------------

def plot_RorZprof(DI, RorZ, fileType, varName, yD, expFlist, expCol, xlim, ylim, xtics, ytics, xlabel, ylabel):


    D     = get_inputFileParameter(DI, ("initParams","d_f"))
    y     = yD*D

    if RorZ=='R':
        odt  = np.loadtxt(DI['pdir']+fileType+"_"+varName+".dat", comments=commentHdr)
        scaleX = D
        figName = DI['pdir']+"radial_"+varName+"_"+fileType+"_" + DI['cn'].replace(".","o")
    else:
        odt  = np.loadtxt(DI['pdir']+fileType+"_"+varName+".dat", comments=commentHdr)
        scaleX = 1
        figName = DI['pdir']+"conditional_"+varName+"_"+fileType+"_"+ DI['cn'].replace(".","o")
    # we might not have all of the dump locations from the dump list
    yodt = get_axialLocations(DI, forceGetDumpTimes=False)[:len(odt[1,1:])]
    npts = len(odt[:,0])
    var  = np.zeros((npts, len(yD)))

    for i in range(len(yD)):
        for k in range(npts):
            fextrap = extrap1d(yodt, odt[k,1:])
            var[k,i] = fextrap(y[i])

    plt.figure(figsize=(12,16))
    plt.rc("font", size=14)

    #for i,ifile in enumerate(expFlist):
    for i in range(len(yD)):
        ax = plt.subplot(5,3,i+1)
        plt.plot(odt[:,0]/scaleX,var[:,i],'k-')
        if expFlist != []:
            exp = np.loadtxt(expFlist[i])
            #plt.plot(exp[:,0],exp[:,expCol],'r-')   # r (mm)
            plt.plot(exp[:,1],exp[:,expCol],'r-')    # r/D
        if (varName=="fv") :
            ax.set_yscale('log')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xticks([])
        plt.yticks([])
        plt.title("y/D = %i" %yD[i])

    plt.subplot(5,3,1)
    if expFlist == []:
        plt.legend(("ODT"), frameon=False)
    else:
        plt.legend(("ODT","EXP"), frameon=False)
    plt.yticks(ytics)
    plt.ylabel(ylabel)

    plt.subplot(5,3,4)
    plt.yticks(ytics)
    plt.ylabel(ylabel)

    plt.subplot(5,3,7)
    plt.yticks(ytics)
    plt.ylabel(ylabel)

    plt.subplot(5,3,10)
    plt.yticks(ytics)
    plt.ylabel(ylabel)

    plt.subplot(5,3,13)
    plt.yticks(ytics)
    plt.ylabel(ylabel)
    plt.xticks(xtics)
    plt.xlabel(xlabel)

    plt.subplot(5,3,14)
    plt.xticks(xtics)
    plt.xlabel(xlabel)

    plt.subplot(5,3,15)
    plt.xticks(xtics)
    plt.xlabel(xlabel)

    #plt.show()

    plt.savefig(figName)

    print("plot_RorZprof: " + figName)

