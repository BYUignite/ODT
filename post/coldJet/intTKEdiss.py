#plots the integral of TKEdiss times r with respect to r, vs. r/(x-x0)
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

def intTKEdiss(DI, profName="TKEdiss"): #uses list of cases
    
    cases = []
    plotname = "intTKEdiss"
    
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
        rnorm = np.empty((npts,ntimes))
        intTKEdiss = np.empty((npts-1,ntimes))

        TKEr = np.empty((npts,ntimes))
        TKE = data[:,1:]

        for ipt in range(npts):
            for itime in range(ntimes):
                TKEr[ipt,itime] = TKE[ipt,itime]*r[ipt]
        
        for itime in range(ntimes):
            intTKEdiss[icl:,itime] = integrate.cumtrapz(TKEr[icl:,itime],r[icl:])
            rnorm[:,itime] = data[:,0]/(times[itime]-4.0*D)
        
        j = 10; axL.plot(rnorm[icl+1:,j], intTKEdiss[icl:,j], color[DI[i]['cn']]) #x=50D
        #plot fit line to check that it starts after asymptote reached to ensure correct y-intercept value
        m,b=np.polyfit(rnorm[icl+8001:,j], intTKEdiss[icl+8000:,j], 1) #change index from icl+8001/icl+8000 depending on start of asymptote
        #fit = rnorm[icl+1:,j]*m+b 
        #axL.plot(rnorm[icl+1:,j],fit,color[DI[i]['cn']][0]+':')

        #cases.append(DI[i]['cn'][8:]+', x=50D') #use if plotting fit
        #cases.append('$y_{int} = $'+str(b)[:7]) #use if plotting fit
        cases.append(DI[i]['cn'][8:]+', x=50D, $y_{int}$ ='+str(b)[:7]) #use if not plotting fit
        plotname+="__"+DI[i]['cn'][8:]

    axL.set_xlabel(r"$r/(x-x_0)$", fontsize=22)
    axL.set_ylabel(r"$\int \epsilon rdr$", fontsize=22)
    axL.set_title("coldJet", fontsize=22)
    axL.legend(cases,loc='best', frameon=False, fontsize=8)
    #axL.set_xlim([-.3,.3]) 
    #axL.set_ylim([-0.05,0.4])

    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'intTKEdiss__ALL'.replace(".","o"))
