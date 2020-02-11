#plots U0/UcL vs. x/D, and finds x-intercept
#plot file directory/name at end of function definition

from __future__ import division
import numpy as np
from data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from data_tools import commentHdr

#---------------------------------------------------------------------------------------------

def u0ucl_xD(DI): #use list of cases

    cases=[]
    plotname="u0ucl_xD"
    x0values=[]
    
    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()
    color = {'coldJet_base':'k', 'coldJet_base_gDens60':'m','coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b--', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r--'}
  

    for i in range(0,len(DI)):
        
        cases.append(DI[i]['cn'][8:])
        plotname+="__"+DI[i]['cn'][8:]
        fname = DI[i]['pdir']+"uvel_cl.dat"
        odt = np.loadtxt(fname, comments=commentHdr)
        U0 = get_inputFileParameter(DI[i], ("initParams", "vel_max"))
        ua = get_inputFileParameter(DI[i], ("initParams", "vel_min"))
        axL.plot(odt[:,0],(U0-ua)/(odt[:,1]-ua),color[DI[i]['cn']])
        m,b = np.polyfit(odt[4:,0],(U0-ua)/(odt[4:,1]-ua),1) #fit line starts at x/d=20
        fit = m*odt[:,0]+b
        axL.plot(odt[:,0],fit,color[DI[i]['cn']][0]+':')
        D = get_inputFileParameter(DI[i], ("initParams", "djeti"))
        x0 = (-b/m)*D
        cases.append('fit, $x_0$ = '+str(x0)[:7]+' m')

    axL.set_ylabel(r"$v_0/v_{cL}$", fontsize=22)
    axL.set_xlabel("$x/D$", fontsize=22)
    axL.legend(cases, loc='best', frameon=False, fontsize=8)
    axL.set_title("coldJet", fontsize=22)
    axL.set_ylim([0,25])
       
    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'u0ucl_xD__ALL'.replace(".","o"))
