#plots FWHM/nozzle diameter vs. x/D
#plot file directory/name at end of function definition

from __future__ import division
import numpy as np
from data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from data_tools import commentHdr

#-----------------------------------------------------------------------------------------------------------

def wD_xD(DI): #use list of cases

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    cases=[]
    plotname="wD_xD"
    
    fig, axL = plt.subplots()
    color = {'coldJet_base':'k','coldJet_base_gDens60':'m', 'coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b--', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r--'}
    
    for i in range(0,len(DI)):
        
        cases.append(DI[i]['cn'][8:])
        plotname+="__"+DI[i]['cn'][8:]
        fname = DI[i]['pdir']+"fwhm_cl_uvel.dat"
        odt = np.loadtxt(fname, comments=commentHdr)
        axL.plot(odt[:,0], odt[:,1], color[DI[i]['cn']])
        m,b = np.polyfit(odt[4:,0],odt[4:,1],1) #start fit line at x/d=20
        fit = m*odt[:,0]+b
        axL.plot(odt[:,0],fit,color[DI[i]['cn']][0]+':')
        D = get_inputFileParameter(DI[i], ("initParams","djeti"))
        x0 = (-b/m)*D
        cases.append("fit, $x_0$ = "+str(x0)[:7]+' m')

    axL.set_ylabel(r"$w/D$", fontsize=22)
    axL.set_xlabel("$x/D$", fontsize=22)
    axL.set_title("coldJet", fontsize=22)
    axL.legend(cases, loc="best", frameon=False, fontsize=8)
    axL.set_ylim([0,25])

    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'wD_xD__ALL'.replace(".","o"))

    #make filter, use on only one case--coldJet_base_gDens120?
#    for i in range(0, len(fit)):
#        if fit[i] < 0:
#            if fit[i+1] > 0:
#                fit[i] = fit[i+1]
#            elif fit[i+2] > 0:
#                fit[i] = fit[i+2]
#            elif fit[i+3] > 0:
#                fit[i] = fit[i+3]

#    fitLine = np.vstack([odt[:,0], fit]).T
#    fname = "/home/abaumga/odt2.0/post/coldJet/base_width.dat"
#    head = "x/D               width             "
#    np.savetxt(fname, fitLine, header=head, fmt="%15.8e ", comments=commentHdr)
