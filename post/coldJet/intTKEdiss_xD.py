#plot of integral of TKEdiss times r with respect to r, vs. x/D
#plot file directory/name at end of function definition

from __future__ import division
import numpy as np
from   data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt
from data_tools import commentHdr
from scipy import integrate

#-------------------------------------------------------------------------------------

def intTKEdiss_xD(DI, profName="TKEdiss"): #uses list of cases
    
    cases = []
    plotname = "intTKEdiss_xD"
    
    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()
    color = {'coldJet_base':'k','coldJet_base_gDens60':'m', 'coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b--', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r--'}

    for i in range(0,len(DI)):

        D = get_inputFileParameter(DI[i],("initParams","djeti"))
        L = get_inputFileParameter(DI[i], ("params","domainLength"))
        mfile = DI[i]['pdir']+"/means_"+profName+".dat"
        data = np.loadtxt(mfile, comments=commentHdr)
        times = get_inputFileParameter(DI[i],("dumpTimes",))

        npts = len(data[:,0])
        ntimes = len(times)
        icl = int(npts/2)+1
        r = data[:,0]
        times = np.array(times)
        xD = times/D

        TKErdr = np.empty((npts-1,ntimes))

        TKEr = np.empty((npts,ntimes))
        TKE = data[:,1:]

        for ipt in range(npts):
            for itime in range(ntimes):
                TKEr[ipt,itime] = TKE[ipt,itime]*r[ipt]
        
        #integrated from centerline out to maximum positive radial distance
        for itime in range(ntimes):
            TKErdr[icl:,itime] = integrate.cumtrapz(TKEr[icl:,itime],r[icl:])
           
        
        axL.plot(xD, TKErdr[-1,:], color[DI[i]['cn']])
       
        cases.append(DI[i]['cn'][8:])
        plotname+="__"+DI[i]['cn'][8:]

    axL.set_xlabel(r"$x/D$", fontsize=22)
    axL.set_ylabel(r"$\int \epsilon rdr$", fontsize=22)
    axL.set_title("coldJet", fontsize=22)
    axL.legend(cases,loc='best', frameon=False, fontsize=8)
    #axL.set_xlim([-.3,.3])
    #axL.set_ylim([-0.05,0.4])

    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'intTKEdiss_xD__ALL'.replace(".","o"))
