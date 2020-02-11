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
dTbin  = 200
Tbin   = np.array([800,1000,1200,1400,1600,1800,2000])
nTbins = len(Tbin)
Tbinlo = Tbin - dTbin/2.0
Tbinhi = Tbin + dTbin/2.0

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

    Pfv   = np.zeros([nxbins,nTbins])

    for iTbin in range(nTbins) :

        i = np.where( np.logical_and( temp >= Tbinlo[iTbin], temp <= Tbinhi[iTbin] ) )

        w = dx[i]

        X = np.array(fv[i])
        P = np.zeros(nxbins)
        if len(X) > 0 :
            P,xbins_fv = compute_wpdf(X,w,0,0.02,nxbins)
            #P,xbins_fv = compute_pdf(X,0,0.02,nxbins)
        Pfv[:,iTbin] = P

    var = np.column_stack((xbins_fv,Pfv))
    fname = "cpdfs_T/P_fv_"+"{0:0>3}".format(it+1)+".dat"
    header = "fv "
    for i in range(len(Tbin)):
        header = header + "P(T=" + str(Tbin[i]) + "+-" + str(dTbin/2) + "), "
    np.savetxt(fname,var,header=header)






