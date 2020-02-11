#plots U/UcL vs. r/(x-x0)
#plot file directory/name at end of function definition

from __future__ import division
import numpy as np
from   data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt
from data_tools import commentHdr

#-------------------------------------------------------------------------------------

def uucl_rxx0(DI, profName="uvel"): #use list of cases

    cases=[]
    #plotname="uucl_rxx0"
    plotname="uucl_rxx0_xD" #for plotting multiple distances
    
    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()
    color = {'coldJet_base':'k', 'coldJet_base_gDens60':'m','coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b--', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r--'}
    #color['coldJet_LplanarTau_gDens60']='k' #use if plotting LPlanarTau cases at multiple distances
    #color['coldJet__LPlanarTau_C_10_ZLES_7_gDens60']='k'     #use if plotting LPlanarTau cases at multiple distances

    for i in range(0,len(DI)):

        D = get_inputFileParameter(DI[i], ("initParams","djeti"))
        mfile = DI[i]['pdir'] + "/means_" + profName + ".dat"
        data = np.loadtxt(mfile, comments=commentHdr)
        times = get_inputFileParameter(DI[i], ("dumpTimes",))
        ua = get_inputFileParameter(DI[i], ("initParams","vel_min"))

        npts = len(data[:,0])
        ntimes = len(times)
        rnorm = np.empty( (npts, ntimes) )
        U = data[:,1:]
        icl  = int(npts/2)+1
        Ucl = U[icl,:]
        
        for j in range(ntimes):
            rnorm[:,j] = data[:,0]/(times[j] - 4.0*D)    

        #if plotting several cases together, best not to plot multiple downstream distances
        #j = 2; axL.plot(rnorm[:,j], (U[:,j]-ua)/(Ucl[j]-ua), color[DI[i]['cn']]+':') #x/D=10
        #j = 4; axL.plot(rnorm[:,j], (U[:,j]-ua)/(Ucl[j]-ua), color[DI[i]['cn']]+'--') #x/D=20
        j = 10; axL.plot(rnorm[:,j], (U[:,j]-ua)/(Ucl[j]-ua), color[DI[i]['cn']]) #x/D=50
        
        #cases.append(DI[i]['cn'][8:]+', x=10D')
        #cases.append(DI[i]['cn'][8:]+', x=20D')
        cases.append(DI[i]['cn'][8:]+', x=50D')
        plotname+="__"+DI[i]['cn'][8:]

    axL.set_ylim([0,1.2])
    axL.set_xlim([-0.25,0.25])
    axL.set_xlabel(r"$r/(x-x_0)$", fontsize=22)
    axL.set_ylabel(r"$v/v_{cL}$", fontsize=22)
    axL.set_title("coldJet",fontsize=22)
    axL.legend(cases,loc="best", frameon=False, fontsize=8)

    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'uucl_rxx0__ALL'.replace(".","o"))
