#plots TKEdiss vs. r/(x-x0); TKEdiss_rxx0_half cuts the vertical range in half to zoom in on cases with lower peaks
#plot file directories/names at end of function definitions

from __future__ import division
import numpy as np
from   data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt
from data_tools import commentHdr

#-------------------------------------------------------------------------------------

def TKEdiss_rxx0(DI, profName="TKEdiss"): #uses list of cases
    
    cases=[]
    plotname="TKEdiss_rxx0"
    #plotname="TKEdiss_rxx0_xD"

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()
    color = {'coldJet_base':'k','coldJet_base_gDens60':'m', 'coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b--', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r--'}
    #color['coldJet_LplanarTau_gDens60']='k' #use if plotting LPlanarTau cases at multiple distances
    #color['coldJet__LPlanarTau_C_10_ZLES_7_gDens60']='k' #use if plotting LPlanarTau cases at multiple distances

    for i in range(0,len(DI)):

        D = get_inputFileParameter(DI[i],("initParams","djeti"))
        mfile = DI[i]['pdir']+"/means_"+profName+".dat"
        data = np.loadtxt(mfile, comments=commentHdr)
        times = get_inputFileParameter(DI[i],("dumpTimes",))
        ucl = DI[i]['pdir']+"/uvel_cl.dat"
        ucl = np.loadtxt(ucl, comments=commentHdr)
        u = ucl[:,1]
        
        npts = len(data[:,0])
        ntimes = len(times)
        rnorm = np.empty((npts,ntimes))
        TKE = data[:,1:]
        for j in range(ntimes):
            rnorm[:,j] = data[:,0]/(times[j]-4.0*D)

        #if plotting several cases at once, best not to plot multiple downstream distances
        #j = 2; axL.plot(rnorm[:,j], TKE[:,j]/((u[j]**3)/times[j]), color[DI[i]['cn']]+':') #x/D=10
        #j = 4; axL.plot(rnorm[:,j], TKE[:,j]/((u[j]**3)/times[j]), color[DI[i]['cn']]+'--') #x/D=20
        j = 10; axL.plot(rnorm[:,j], TKE[:,j]/((u[j]**3)/times[j]), color[DI[i]['cn']]) #x/D=50
        
        #cases.append(DI[i]['cn'][8:]+', x=10D')
        #cases.append(DI[i]['cn'][8:]+', x=20D')
        cases.append(DI[i]['cn'][8:]+', x=50D')
        plotname+="__"+DI[i]['cn'][8:]

    #axL.set_ylim([0,4000])
    axL.set_xlim([-0.4,0.4])
    axL.set_xlabel(r"$r/(x-x_0)$", fontsize=22)
    axL.set_ylabel(r"$TKEdiss \bullet \frac{ x}{(u_{cL})^3}$", fontsize=22)
    axL.set_title("coldJet",fontsize=22)
    axL.legend(cases,loc="best", frameon=False, fontsize=8)

    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'TKEdiss_rxx0__ALL'.replace(".","o"))

#-------------------------------------------------------------------------
   
def TKEdiss_rxx0_half(DI, profName="TKEdiss"): #uses list of cases
    
    cases=[]
    plotname="TKEdiss_rxx0_half"
    #plotname="TKEdiss_rxx0_half_xD"

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()
    color = {'coldJet_base':'k','coldJet_base_gDens60':'m', 'coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b--', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r--'}
    #color['coldJet_LplanarTau_gDens60']='k' #use if plotting LPlanarTau cases at multiple distances
    #color['coldJet__LPlanarTau_C_10_ZLES_7_gDens60']='k' #use if plotting LPlanarTau cases at multiple distances

    for i in range(0,len(DI)):

        D = get_inputFileParameter(DI[i],("initParams","djeti"))
        mfile = DI[i]['pdir']+"/means_"+profName+".dat"
        data = np.loadtxt(mfile, comments=commentHdr)
        times = get_inputFileParameter(DI[i],("dumpTimes",))
        ucl = DI[i]['pdir']+"/uvel_cl.dat"
        ucl = np.loadtxt(ucl, comments=commentHdr)
        u = ucl[:,1]
        
        npts = len(data[:,0])
        ntimes = len(times)
        rnorm = np.empty((npts,ntimes))
        TKE = data[:,1:]
        for j in range(ntimes):
            rnorm[:,j] = data[:,0]/(times[j]-4.0*D)

        #if plotting several cases at once, best not to plot multiple downstream distances
        #j = 2; axL.plot(rnorm[:,j], TKE[:,j]/((u[j]**3)/times[j]), color[DI[i]['cn']]+':') #x/D=10
        #j = 4; axL.plot(rnorm[:,j], TKE[:,j]/((u[j]**3)/times[j]), color[DI[i]['cn']]+'--') #x/D=20
        j = 10; axL.plot(rnorm[:,j], TKE[:,j]/((u[j]**3)/times[j]), color[DI[i]['cn']]) #x/D=50
        
        #cases.append(DI[i]['cn'][8:]+', x=10D')
        #cases.append(DI[i]['cn'][8:]+', x=20D')
        cases.append(DI[i]['cn'][8:]+', x=50D')
        plotname+="__"+DI[i]['cn'][8:]

    axL.set_ylim([0,0.18])
    axL.set_xlim([-0.3,0.3])
    axL.set_xlabel(r"$r/(x-x_0)$", fontsize=22)
    axL.set_ylabel(r"$TKEdiss \bullet \frac{ x}{(u_{cL})^3}$", fontsize=22)
    axL.set_title("coldJet",fontsize=22)
    axL.legend(cases,loc="best", frameon=False, fontsize=8)

    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'TKEdiss_rxx0_half__ALL'.replace(".","o"))
