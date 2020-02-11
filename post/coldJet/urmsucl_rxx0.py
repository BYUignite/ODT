#plots rms velocity fluctuations normalized by centerline velocity vs. r/(x-x0)
#plot file directory/name at end of function definition

from __future__ import division
import numpy as np
from data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from data_tools import commentHdr


#---------------------------------------------------------------------------------------------

def urmsucl_rxx0(DI): #use list of cases

    cases=[]
    #plotname="urmsucl_rxx0"
    plotname="urmsucl_rxx0_xD"

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()
    color = {'coldJet_base':'k','coldJet_base_gDens60':'m', 'coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b--', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r--'}
    #color['coldJet_LplanarTau_gDens60']='k' #use if plotting LPlanarTau cases at multiple distances
    #color['coldJet__LPlanarTau_C_10_ZLES_7_gDens60']='k' #use if plotting LPlanarTau cases at multiple distances

    for i in range(0,len(DI)):

        #cases.append(DI[i]['cn'][8:]+' ,$x$ = 10$D$')
        #cases.append(DI[i]['cn'][8:]+' ,$x$ = 20$D$')
        cases.append(DI[i]['cn'][8:]+' ,$x$ = 50$D$')
        plotname += '__'+DI[i]['cn'][8:]
        mfile = DI[i]['pdir']+"means_uvel.dat"
        data = np.loadtxt(mfile, comments=commentHdr)
        times = get_inputFileParameter(DI[i], ("dumpTimes",))
        ua = get_inputFileParameter(DI[i], ("initParams","vel_min"))
        D = get_inputFileParameter(DI[i], ("initParams", "djeti"))
        sfile = DI[i]['pdir']+"sig_uvel.dat"
        sig_uvel = np.loadtxt(sfile, comments=commentHdr)

        npts = len(data[:,0])
        ntimes = len(times)
        rnorm = np.empty((npts,ntimes))
        
        U = data[:,1:]
        icl = int(npts/2)+1
        Ucl = U[icl,:]

        for j in range(ntimes):
            rnorm[:,j] = data[:,0]/(times[j] - 4.0*D)
        
        #if plotting several cases together, best not to plot multiple downstream distance
        #j=2; axL.plot(rnorm[:,j],(sig_uvel[:,j]/Ucl[j]), color[DI[i]['cn']]+':') #x/D=10
        #j=4; axL.plot(rnorm[:,j],(sig_uvel[:,j]/Ucl[j]), color[DI[i]['cn']]+'--') #x/D=20
        #j=4; axL.plot(rnorm[:,j],(sig_uvel[:,j]/Ucl[j]), color[DI[i]['cn']]) #x/D=20
        j=10; axL.plot(rnorm[:,j],(sig_uvel[:,j]/Ucl[j]), color[DI[i]['cn']]) #x/D=50

    axL.set_ylabel(r"$v_{rms}/v_{cL}$", fontsize=22)
    axL.set_xlabel(r"$r/(x-x_0)$", fontsize=22)
    axL.legend(cases, loc='best', frameon=False, fontsize=8)
    axL.set_title("coldJet", fontsize=22)
    #axL.set_ylim([-0.1,0.5])
    axL.set_xlim([-0.5,0.5])
           
    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'urmsucl_rxx0__ALL'.replace(".","o"))
    #plt.savefig('../../data/plots_coldJet/'+'urmsucl_rxx0__20_ALL'.replace(".","o"))
