import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, extrap1d, get_domainBounds

import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

#--------------------------------------------------------------------------------------------

def energy_spectra_time(DI):

    dataFiles = glob.glob(DI['cdir']+"data_py/data_*.npy")
    ntimes    = len(dataFiles)
    varNames  = get_dataHeaderVars(DI)
    times     = get_inputFileParameter(DI, ("dumpTimes",))         # times or ypositions if spatial
    nvar      = len(varNames)
    nLevels   = get_inputFileParameter(DI, ("params", "nLevels"))
    if get_inputFileParameter(DI, ("params", "LScHips")):
        scalarSc = get_inputFileParameter(DI, ("scalarSc",))
        maxSc = np.max(scalarSc)
        if maxSc > 1.0:
            nLevels += int(np.ceil(np.log2(maxSc)/2.0))
    nx        = int(2**(nLevels-1))
    nrlz      = get_nRlz(DI, nx)
    L0        = get_inputFileParameter(DI, ("params", "domainLength"))
    A         = get_inputFileParameter(DI, ("params", "Afac"))


    
    #--------------------------------------------------------------------------------------------
    
    def get_scalar_variance_at_all_levels(scalarName):
    
        sVar  = np.zeros(nLevels)              # scalar variance at each level
        iScal = varNames.index(scalarName)     # scalar index
        irlz  = 0
    
        itimes = range(50,ntimes)
        for itime in itimes:
    
            data = get_data_realization(DI, itime, irlz, nx)
    
            for iLevel in range(nLevels):
    
                nSubtrees          = 2**iLevel
                nParcelsPerSubtree = 2**(nLevels-1-iLevel)    # int(nx/nSubtrees)
    
                subtreeVar = np.zeros(nSubtrees)              # variance in each subtree
    
                for iTree in range(nSubtrees):
                    istart = iTree*nParcelsPerSubtree
                    iend   = istart + nParcelsPerSubtree
                    subtreeVar[iTree] = np.var(data[istart:iend, iScal])
    
                sVar[iLevel] += np.mean(subtreeVar)
        
        return sVar / len(itimes)

    #--------------------------------------------------------------------------------------------

    E = np.zeros((nvar, nLevels))
    L_levels = L0*A**np.arange(nLevels)
    kwave = 2*np.pi/L_levels            # wavenumbers at levels

    for ivar in range(nvar):
        print(f'energy spectra: {varNames[ivar]}')
        scalVar = get_scalar_variance_at_all_levels(varNames[ivar])
        E[ivar, 1:] = (scalVar[:-1] - scalVar[1:])/kwave[1:]

    #--------------------------------------------------------------------------------------------

    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True})

    head = "kwave            E"
    for ivar in range(nvar):
        data = np.column_stack((kwave, E[ivar,:]))
        fname = DI['pdir'] + "Espec_" + varNames[ivar] + ".dat"
        np.savetxt(fname, data, header=head, fmt="%15.8e ")

        #-------------- plot

        fig, ax = plt.subplots()
        plt.cla()

        X = kwave[1:-1]/kwave[1]
        Y = E[ivar,1:-1]/E[ivar,1]
        
        x1 = 1                      # -5/3 line
        x2 = 2
        y1 = np.log10(Y[3]) + 1
        y2 = y1 - 5/3*(x2-x1)
        x53 = np.array([x1,x2])
        y53 = np.array([y1,y2])

        ax.loglog(X,Y, 'ko-')
        ax.plot(10**(x53), 10**(y53), '--', color='gray')

        ax.set_xlim([0.6, 10**(np.ceil(np.log10(np.max(X))))])
        ax.set_ylim([10**(np.floor(np.log10(np.min(Y)))), 2])

        ax.minorticks_off()
        ax.tick_params(top=True, right=True, which='both', direction='in')
        ax.xaxis.set_major_locator(LogLocator(base=10, numticks=15))
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=15))
        #locmin = LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=1000)
        #ax.xaxis.set_minor_locator(locmin)
        #ax.yaxis.set_minor_locator(locmin)
        
        ax.set_xlabel(r"$\kappa/\kappa_0$", fontsize=22)
        ax.set_ylabel(r"$E/E_0$",      fontsize=22)
        
        pname = DI['pdir'] + "Espec_" + varNames[ivar]
        plt.savefig(pname)
