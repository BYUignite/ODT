#plot of centerline logTKEdiss fluctuations vs. x/D
#plot file directory/name at end of function definition

from __future__ import division
import numpy as np
from data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from data_tools import commentHdr

#---------------------------------------------------------------------------------------------

def logTKEdissrms_xD(DI): #uses list of cases

    cases=[]
    plotname="logTKEdissrms_xD"

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()
    color = {'coldJet_base':'k', 'coldJet_base_gDens60':'m', 'coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b--', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r--'}
    
    for i in range(0,len(DI)):

        cases.append(DI[i]['cn'][8:])
        plotname+="__"+DI[i]['cn'][8:]

        fname = DI[i]['pdir']+"logTKEdiss_cl.dat"

        data = np.loadtxt(fname, comments=commentHdr)
       
        axL.plot(data[:,0],data[:,2],color[DI[i]['cn']])

    axL.set_ylabel(r"$logTKEdiss_{rms,cL}$", fontsize=22)
    axL.set_xlabel("$x/D$", fontsize=22)
    axL.legend(cases, loc='best', frameon=False, fontsize=8)

    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig("../../data/plots_coldJet/"+"logTKEdissrms_xD__ALL".replace(".","o"))
