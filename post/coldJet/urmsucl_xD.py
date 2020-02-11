#plots rms velocity fluctuations normalized by centerline velocity vs. x/D
#plot file directory/name at end of function definition

from __future__ import division
import numpy as np
from data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from data_tools import commentHdr

#---------------------------------------------------------------------------------------------

def urmsucl_xD(DI): #use list of cases

    cases=[]
    plotname="urmsucl_xD"

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()
    color = {'coldJet_base':'k', 'coldJet_base_gDens60':'m','coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b--', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r--'}
    
    for i in range(0,len(DI)):

        cases.append(DI[i]['cn'][8:])
        plotname+="__"+DI[i]['cn'][8:]


        mfile = DI[i]['pdir']+"means_uvel.dat"
        sfile = DI[i]['pdir']+"sig_uvel.dat"

        data = np.loadtxt(mfile, comments=commentHdr)
        y    = data[:,1:]
        data = np.loadtxt(sfile, comments=commentHdr)
        ys   = data[:,1:]

        npts, ntimes = np.shape(data)
        ntimes = ntimes - 1

        cLine  = np.zeros(ntimes)
        scLine = np.zeros(ntimes)
        imid   = int(npts/2) + 1

        for it in range(ntimes) :
            cLine[it] = y[imid,it]
            scLine[it] = ys[imid,it]

        times = get_inputFileParameter(DI[i], ("dumpTimes",))
        D     = get_inputFileParameter(DI[i], ("initParams","djeti"))
        times = np.array(times)

        data = np.vstack([times/D, cLine, scLine]).T
  
        axL.plot(data[:,0], data[:,2]/data[:,1], color[DI[i]['cn']]) 

    axL.set_ylabel(r"$v_{rms,cL}/v_{cL}$", fontsize=22)
    axL.set_xlabel("$x/D$", fontsize=22)
    axL.legend(cases, loc='best', frameon=False, fontsize=8)
    axL.set_title("coldJet", fontsize=22)
    #axL.set_ylim([-3,3])
    #axL.set_xlim([-0.1,0.1])

  
    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig("../../data/plots_coldJet/"+"urmsucl_xD__ALL".replace(".","o"))
