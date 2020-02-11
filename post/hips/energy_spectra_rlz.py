import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, extrap1d, get_domainBounds

#--------------------------------------------------------------------------------------------

def energy_spectra_rlz(DI):

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
    
    def get_scalar_variance_at_all_levels(scalarName, itime):
    
        sVar  = np.zeros(nLevels)              # scalar variance at each level
        iScal = varNames.index(scalarName)     # scalar index
    
        for irlz in range(nrlz):
    
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
        
        return sVar / nrlz

    #--------------------------------------------------------------------------------------------

    E = np.zeros((nvar, nLevels, ntimes))
    L_levels = L0*A**np.arange(nLevels)
    kwave = 2*np.pi/L_levels            # wavenumbers at levels

    for itime in range(ntimes):
        for ivar in range(nvar):
            print(f'energy spectra: var, time = {ivar}, {itime}')
            scalVar = get_scalar_variance_at_all_levels(varNames[ivar], itime)
            E[ivar, 1:, itime] = (scalVar[:-1] - scalVar[1:])/kwave[1:]

    #--------------------------------------------------------------------------------------------

    head = "kwave            E"
    for i,time in enumerate(times) :
        hi = str(i+2) + "_" + str(time)
        hi = hi + (17-len(hi))*" "
        head = head + hi

    for ivar in range(nvar):
        data = np.column_stack((kwave, E[ivar, :, :]))
        fname = DI['pdir'] + "Espec_" + varNames[ivar] + ".dat"
        np.savetxt(fname, data, header=head, fmt="%15.8e ")

