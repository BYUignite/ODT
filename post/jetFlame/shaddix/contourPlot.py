
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from   data_tools import get_inputFileParameter, get_fstoic
from data_tools import commentHdr

#----------------------------------------------------------

def contourPlot(DI, varName, levels, ticks, figName, ZvarName='mixf', addZstContour=False, means_rms='means'):
    """ Contour plot with input based on previously computed ?????_varName.dat
    varName    -- variable for plot using either means_varName.dat or rms_varName.dat depending on the value of 
    means_rms  -- using either means_varName.dat or rms_varName.dat depending on this arg (????? above)
    """

    print("Making contour plot for figure: ", figName)

    fname = DI['pdir']+means_rms+"_"+varName+".dat"
    
    data  = np.loadtxt(fname, comments=commentHdr)
    x     = data[:,0]
    #x     = x-np.average(x)
    var   = data[:,1:]

    y = get_inputFileParameter(DI, ("dumpTimes",))
    # we might not have all of the dump locations from the dump list 
    # so truncate at length of data
    y = np.array(y)[:len(data[1,1:])]
    D     = get_inputFileParameter(DI, ("initParams","d_f"))
    x = x/D
    y = y/D

    X,Y = np.meshgrid(x,y)

    #------------------ plot

    aspectRatio = (np.max(y) - np.min(y))/(np.max(x)-np.min(x))
    plt.figure(figsize = (4,2.2*aspectRatio))
    plt.rc("font", size=16)

    if levels != [] :
        C = plt.contourf(X,Y,var.T,levels=levels,extend='max', cmap=cm.viridis)
    else :
        C = plt.contourf(X,Y,var.T,50, cmap=cm.viridis)

    #C.cmap.set_over('k')     # color for high out of range values
    #C.cmap.set_under('k')    # color for low out of range values
    ax = plt.gca()
    ax.axis((np.min(x),np.max(x),np.min(y),np.max(y)))
    ax.set_xbound(lower=np.min(x), upper=np.max(x))
    ax.set_ybound(lower=np.min(y), upper=np.max(y))
    ax.set_aspect('equal')
    plt.xlabel('position (m)',fontsize=16)
    plt.ylabel('axial position (m)',fontsize=16)
    plt.colorbar(ticks = ticks)
    #plt.xticks((-0.2,0,0.2))
    #plt.yticks(())
    #ax.grid(True, color='w', linewidth=1, linestyle = ':')

    if addZstContour:
        zstoic = get_fstoic(DI)
        Z = np.loadtxt(DI['pdir']+"means_" + ZvarName + '.dat', comments=commentHdr)
        Z = Z[:,1:]
        plt.contour(X,Y,Z.T,levels=np.array([zstoic]), colors='white', linewidth=0.1)

    #plt.show()
    plt.savefig(figName)
