
from __future__ import division
import glob
import os
import numpy as np

def data_py(DI) :
    """ Make directory for python numpy formatted data caseN/data_py.
    From data directories, get header by scanning contents of line 
    containing 'pos' and save this in data_py/header.dat.  Save each dump
    time #### from all realization into file caseN/data_py/data_py_####.npy
    """

    os.mkdir(DI['cdir']+'data_py')

    #-------- get the header

    fname = DI['ddir'] + "data_00000/dmp_00000.dat"
    ifile  = open(fname,'r')
    lines = ifile.readlines()
    ifile.close()
    for header in lines :
        if "pos" in header or "mixf" in header :
            break
    ofile  = open(DI['cdir'] + "/data_py/header.dat", 'w')
    ofile.write(header)
    ofile.close()

    #---------

    nrlz   = len(glob.glob(DI['ddir']+"data_*"))
    ntimes = len(glob.glob(DI['ddir']+"data_00000/dmp_*.dat"))

    for it in range(ntimes) :
        for irlz in range(nrlz) :
            print("data_py: it, irlz = ", it, irlz)

            fname = DI['ddir']+"data_"+ "{0:0>5}".format(irlz) +"/dmp_"+"{0:0>5}".format(it)+".dat"
            try:
                A = np.loadtxt(fname)
    
                if irlz==0 :
                    data = A
                else :
                    data = np.vstack([data,A])
            except IOError:
                print("Missing file " + fname)

        if np.ndim(A) == 1:
            data = np.reshape(data, np.size(data))
            data = data[:,np.newaxis]

        fname = DI['cdir']+"data_py/data_py_" + "{0:0>5}".format(it) + ".npy"
        np.save(fname,data)

