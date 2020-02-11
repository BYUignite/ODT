#!/usr/bin/env python

import sys
import glob
import os
import numpy as np
sys.path.append('../tools')
from compute_wpdf import *
from compute_pdf import *
from getInputFileParameter_func import *

Ldomain = float(getInputFileParameter_func('../../input/odtParam.inp', 'domainLength'))
nxbins = 40
dfbin  = 0.1
fbin   = np.array([0.0319,0.06374, 0.2, 0.4])
nfbins = len(fbin)
fbinlo = fbin - dfbin/2.0
fbinhi = fbin + dfbin/2.0

caseName = sys.argv[1];

nrlz   = len(glob.glob("../../data/"+caseName+"/data_*"));
ntimes = len(glob.glob("../../data/"+caseName+"/data_0/dmp_*.dat"));

for it in range(ntimes) :

    print "it = ", it

    fname = "../../data_py_" + caseName + "/data_py_" + "{0:0>4}".format(it+1) + ".npy"
    data = np.load(fname)

    posf   = data[:,1]
    rho  = data[:,2]
    mixf = data[:,5]
    temp = data[:,7]
    Ys   = data[:,36]

    fv   = rho*Ys

    posf = np.append(posf,Ldomain)

    dx = posf[1:]-posf[0:-1]
    dx = dx.clip(0)

    Ptemp = np.zeros([nxbins,nfbins])
    Pfv   = np.zeros([nxbins,nfbins])

    for ifbin in range(nfbins) :

        i = np.where( np.logical_and( mixf >= fbinlo[ifbin], mixf <= fbinhi[ifbin] ) )

        w = dx[i]

        X = np.array(temp[i])
        P,xbins_temp = compute_wpdf(X,w,300,2500,nxbins)
        #P,xbins_temp = compute_pdf(X,300,2500,nxbins)
        Ptemp[:,ifbin] = P

        X = np.array(fv[i])
        P = np.zeros(nxbins)
        if len(X) > 0 :
            P,xbins_fv = compute_wpdf(X,w,0,0.02,nxbins)
            #P,xbins_fv = compute_pdf(X,0,0.02,nxbins)
        Pfv[:,ifbin] = P

    var = np.column_stack((xbins_temp,Ptemp))
    fname = "cpdfs/P_temp_"+str(it+1)+".dat"
    header = "# T(K) P(f=0.0319+-0.05), P(f=0.06374+-0.05) P(f=0.2+-0.05) P(f=0.4+-0.05)"
    np.savetxt(fname,var)

    var = np.column_stack((xbins_fv,Pfv))
    fname = "cpdfs/P_fv_"+str(it+1)+".dat"
    header = "# fv P(f=0.0319+-0.05), P(f=0.06374+-0.05) P(f=0.2+-0.05) P(f=0.4+-0.05)"
    np.savetxt(fname,var)






