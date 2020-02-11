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
    dumpTimes = get_inputFileParameter(DI, ("dumpTimes",))         # times or ypositions if spatial
    tEnd      = get_inputFileParameter(DI, ("params","tEnd"))      # specified end time
    
    # Incomplete runs ---------------------------------------------------------------------------
    
    print("Checking for incomplete realizations.")
    
    incompleteRlzs = []
    times = []

    # removes dump times bigger than tEnd to avoid discrepancies in the input file information
    for t in dumpTimes :
        if t < tEnd :
            times.append(t)
    
    # identifies rlzs that don't have enough dump files
    for irlz in range(nrlz) :

        irlz_str = format(irlz, '05') 
    
        files = glob.glob(DI['ddir']+"data_"+irlz_str+"/dmp_*.dat")
        if (len(files) < len(times)) :
            print("INCOMPLETE: rlz " + str(irlz))
            incompleteRlzs += [irlz]
    
    # if rlzs are incomplete, copy over them with adjacent rlzs according to user input
    if len(incompleteRlzs) > 0 :
    
        copyOverOldFlag = input(str(len(incompleteRlzs))+" realization(s) with incomplete data. Further post-processing scripts may not function properly. Copy over incomplete realizations with adjacent realizations? (y/n/q) ")
    
        if copyOverOldFlag == "q" :
            sys.exit()

        elif copyOverOldFlag == "y" :
    
            copying_args    = ['']*3
            copying_args[0] = "./copyOverIncompleteRlz.sh"
            copying_args[1] = caseN
    
            for i in incompleteRlzs :
                copying_args[2] = str(i)
                subprocess.call(copying_args)
        
        elif copyOverOldFlag == "n" :
            print("Continuing with incomplete realizations.")

    else :
        print("No incomplete realizations.")
    
#--------------------------------------------------------------------------------------------

def checkForBlowout(DI,extinction_temp) :
    
    dataFiles = glob.glob(DI['cdir']+"data_py/data_*.npy")
    ntimes    = len(dataFiles)
    nrlz      = get_nRlz(DI)
    varNames  = get_dataHeaderVars(DI)
    times     = get_inputFileParameter(DI, ("dumpTimes",))         # times or ypositions if spatial

    # Blowout -----------------------------------------------------------------------------------
    
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
        warning_flag = input("WARNING: Blowout in %i realizations. Process at your own risk. Continue? (y/n) " %(blowout_num))
        
        if warning_flag == "n" :
            sys.exit()
    
    else :
        print("No blowout cases.")
    
    
