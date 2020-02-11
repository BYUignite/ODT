from __future__ import division
import os
import sys
import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, extrap1d, get_domainBounds
import subprocess

#--------------------------------------------------------------------------------------------

def checkForIncompleteRlzs(DI, caseN) :

    #--------------------------------------------------------------------------------------------
    
    nrlz      = len(glob.glob(DI['ddir']+"data_*")) 
    times     = get_inputFileParameter(DI, ("dumpTimes",))         # times or ypositions if spatial
    
    # Incomplete runs ---------------------------------------------------------------------------
    
    print("Checking for incomplete realizations.")
    
    incompleteRlzs = []
    
    for irlz in range(nrlz) :

        irlz_str = format(irlz, '05') 
    
        files = glob.glob(DI['ddir']+"data_"+irlz_str+"/dmp_*.dat")
        if (len(files) < len(times)) :
            print("INCOMPLETE: rlz " + str(irlz))
            incompleteRlzs += [irlz]
    
    if len(incompleteRlzs) != 0 :
    
        copyOverOldFlag = input(str(len(incompleteRlzs))+" realization(s) have incomplete data. Copy over them from adjacent realizations? (y/n) ")
    
        if copyOverOldFlag == "y" :
    
            copying_args    = ['']*3
            copying_args[0] = "./copyOverIncompleteRlz.sh"
            copying_args[1] = caseN
    
            for i in incompleteRlzs :
                copying_args[2] = str(i)
                subprocess.call(copying_args)
    
    else :
    
        print("No incomplete realizations.")
    
#--------------------------------------------------------------------------------------------

def checkForBlowout(DI) :
    
    dataFiles = glob.glob(DI['cdir']+"data_py/data_*.npy")
    ntimes    = len(dataFiles)
    nrlz      = get_nRlz(DI)
    varNames  = get_dataHeaderVars(DI)
    times     = get_inputFileParameter(DI, ("dumpTimes",))         # times or ypositions if spatial

    # Blowout -----------------------------------------------------------------------------------
    
    extinction_temp = 1800.0
    
    itemp = varNames.index("temp")
    blowout_num = 0
    
    print("Checking for extinguished/blowout cases.")
    
    for irlz in range(nrlz) :
        data = get_data_realization(DI,ntimes-1,irlz)
        T = data[:,itemp]
    
        Tmax = np.amax(T)
        if Tmax < extinction_temp :
            print("BLOWOUT: rlz " + str(irlz))
            blowout_num += 1
    
    if blowout_num >= 1 :
        warning_flag = input("WARNING: Blowout in %i realizations. Process at your own risk. Quit? (y/n) " %(blowout_num))
        
        if warning_flag == "y" :
            sys.exit()
    
    else :
        print("No blowout cases.")
    
    
