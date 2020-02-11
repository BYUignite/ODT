#plots centerline TKEdiss vs. x/D
#plot file directory/name at end of function definition

from __future__ import division
import numpy as np
from   data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt
from data_tools import commentHdr

#-------------------------------------------------------------------------------------

def TKEdiss_xD(DI, profName="TKEdiss"): #uses list of cases
    
    cases=[]
    plotname="TKEdiss_xD"

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()
    color = {'coldJet_base':'k','coldJet_base_gDens60':'m', 'coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b--', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r--'}

    for i in range(0,len(DI)):

        cases.append(DI[i]['cn'][8:])
        plotname+="__"+DI[i]['cn'][8:]
        fname = DI[i]['pdir']+"TKEdiss_cl.dat"
        odt = np.loadtxt(fname, comments=commentHdr)
        axL.plot(odt[:,0],odt[:,1], color[DI[i]['cn']])

    #axL.set_ylim([0,1.2])
    #axL.set_xlim([-0.2,0.2])
    axL.set_xlabel(r"$x/D$", fontsize=22)
    axL.set_ylabel(r"$TKEdiss_{cL}$", fontsize=22)
    axL.set_title("coldJet",fontsize=22)
    axL.legend(cases,loc="best", frameon=False, fontsize=8)

    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'TKEdiss_xD__ALL'.replace(".","o"))
