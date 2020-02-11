from __future__ import division
import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, extrap1d, get_domainBounds

#--------------------------------------------------------------------------------------------

def fv(DI) :
    """
    Adds a column to data_py_*.npy files for fv (ppmv) and adds 'fv' to header.dat.
    After it's added, fv post data is calculated automatically with the other variables.
    """
    #--------------------------------------------------------------------------------------------

    varNames  = get_dataHeaderVars(DI)

    fname = DI['cdir']+"data_py/header.dat"
    with open(fname,'r') as f :
        line = f.readline().strip()

    try :
        varNames.index("fv")
    except :
        with open(fname,'w') as f :
            f.writelines(line + "    " + str(len(varNames)+1) + "_fv")

    #--------------------------------------------------------------------------------------------

        dataFiles = glob.glob(DI['cdir']+"data_py/data_*.npy")
        ntimes    = len(dataFiles)
        varNames  = get_dataHeaderVars(DI)

        M1_col    = varNames.index('M1')
        rho_col   = varNames.index('rho')
        rho_soot  = get_inputFileParameter(DI, ("sootParams","rho_s"))

        for itime in range(ntimes) :

            print("fv: appending fv to data_py_" + "{0:0>5}".format(itime) + ".npy")
            fname    = DI['cdir']+"data_py/data_py_" + "{0:0>5}".format(itime) + ".npy"
            data     = np.load(fname)
            f_v      = np.empty([len(data),1])

            for row in range(len(f_v)) :
                f_v[row] = data[row,M1_col] * data[row,rho_col] / rho_soot * 1e6    # m3-soot / m3-total, ppmv

            data_all = np.hstack([data,f_v])

            np.save(fname,data_all)

