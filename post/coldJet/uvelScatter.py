#plots scatter plot of u-velocity vs. position at specified dmp time for a single case
#or, plots a group of 4 scatter plots of u-velocity vs. position at 4 specified dmp times for a single case (easier to compare different dmp times, but each subplot lower quality)
#plot file directory/name at end of function definition

from __future__ import division
import os
import sys
import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, extrap1d, get_domainBounds
from data_tools import commentHdr
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
#--------------------------------------------------------------------------------------------

def uvelScatter(DI): #uses one case, NOT a list of cases

    dataFiles = glob.glob(DI['cdir']+"data_py/data_*.npy")
    ntimes    = len(dataFiles)
    nrlz      = get_nRlz(DI)
    varNames  = get_dataHeaderVars(DI)
    times     = get_inputFileParameter(DI, ("dumpTimes",))         # times or ypositions if spatial
    times     = times[0:ntimes]    # limit times to number of dataFiles, not all listed dumpTimes
    nvar      = len(varNames)
 
    
    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

#comment or uncomment the following sections for either 4 plots or 1 plot (do not have both uncommented)

####FOR 4 SIDE-BY-SIDE PLOTS FOR 4 DUMPTIMES####

#    fig,axL=plt.subplots(2,2)
#    dmps=['00005','00010','00015','00020']
#    for dmp in range(len(dmps)):

#        fname = DI['cdir']+'data_py/data_py_'+dmps[dmp]+'.npy'
#        data_all = np.load(fname)
#        pos_all = data_all[:,0]
#        uvel_all = data_all[:,4]
#        iplt = np.random.choice(pos_all.shape[0],10000,replace=False)
        
#        for i in iplt:
#            if dmp==0:
#                axL[0,0].scatter(pos_all[i],uvel_all[i],s=0.5,c='b')
#                axL[0,0].set_title('00005')
#            if dmp==1:
#                axL[0,1].scatter(pos_all[i],uvel_all[i],s=0.5,c='b')
#                axL[0,1].set_title('00010')
#            if dmp==2:
#                axL[1,0].scatter(pos_all[i],uvel_all[i],s=0.5,c='b')
#                axL[1,0].set_title('00015')
#            if dmp==3:
#                axL[1,1].scatter(pos_all[i],uvel_all[i],s=0.5,c='b')
#                axL[1,1].set_title('00020')
#    for ax in axL.flat:
#        ax.set(xlabel='pos', ylabel='uvel')
#    for ax in axL.flat:
#        ax.label_outer()
#    plt.savefig(DI['pdir']+'uvelScatter'.replace(",","o"))

####FOR 1 PLOT FOR 1 DUPMTIME####

    fig,axL=plt.subplots()

    fname = DI['cdir']+'data_py/data_py_00010.npy' #change dmp file name here
    data_all = np.load(fname)
    pos_all = data_all[:,0]
    uvel_all = data_all[:,4]
    iplt = np.random.choice(pos_all.shape[0],10000,replace=False)

    for i in iplt:
        axL.scatter(pos_all[i],uvel_all[i],s=0.5,c='b')
    axL.set_ylabel('uvel')
    axL.set_xlabel('pos')
    axL.set_title(DI['cn'][8:])
    plt.savefig('../../data/plots_coldJet/uvelScatter_00010'+DI['cn'][8:].replace(",","o")) #change dmp file name here

#-------------------------------------------------------------------
    
try :
    caseN = sys.argv[1]
except :
    raise ValueError("Include the case name in the call to uvelScatter.py")

DI = {'pdir':'../../data/'+caseN+'/post/',       \
      'ddir':'../../data/'+caseN+'/data/',       \
      'cdir':'../../data/'+caseN+'/',            \
      'cn':caseN}

if not os.path.exists(DI['pdir']):
    os.mkdir(DI['pdir'])


uvelScatter(DI)
