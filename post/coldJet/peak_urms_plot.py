#plot of peak values of normalized centerline rms velocity fluctuations vs. x/D
#plot file directory/name at end of function definition
#need to run peak_urms.py first for all cases to be plotted here

from __future__ import division
import numpy as np
from data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from data_tools import commentHdr

#------------------------------------------------------------------------
#plot peaks vs distance
def peak_urms_plot(DI): #uses list of cases
    
    cases=[]
        
    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    color = {'coldJet_base':'k', 'coldJet_base_gDens60':'m','coldJet_base_gDens120':'b', 'coldJet_C_10_gDens120':'g', 'coldJet_C_10_ZLES_7_gDens120':'r', 'coldJet_LplanarTau_gDens60':'b', 'coldJet_ZLES_7_gDens120':'c', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'r'}
    marker = {'coldJet_base':'o', 'coldJet_base_gDens60':'o','coldJet_base_gDens120':'o', 'coldJet_C_10_gDens120':'o', 'coldJet_C_10_ZLES_7_gDens120':'o', 'coldJet_LplanarTau_gDens60':'v', 'coldJet_ZLES_7_gDens120':'o', 'coldJet__LPlanarTau_C_10_ZLES_7_gDens60':'v'}

    for i in range(0,len(DI)):
        
        cases.append(DI[i]['cn'][8:])
        allPeaks = np.loadtxt(DI[i]['pdir']+"urmspeaks.dat",comments=commentHdr)

        plt.scatter(allPeaks[:,0], allPeaks[:,1], color=color[DI[i]['cn']], marker=marker[DI[i]['cn']])
    
    plt.ylabel("peak $u_{rms}/u_{cL}$", fontsize=22)
    plt.xlabel("x/D", fontsize=22)
    plt.title('coldJet')
    plt.legend(cases, loc='lower left', fontsize=8, frameon=False)

    plt.savefig('../../data/plots_coldJet/'+"peak_urms_ALL".replace(".","o"))

