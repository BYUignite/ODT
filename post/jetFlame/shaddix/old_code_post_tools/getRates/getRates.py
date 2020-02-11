#!/usr/bin/env python

import cantera as ct
import glob
import os
import numpy as np
import sys
sys.path.append('../tools')
from compute_wpdf import *
from compute_pdf import *
from getInputFileParameter_func import *
from scipy.interpolate import *

class reactionRates :

   ###########################################################################################

   def __init__(self, cname) :

       self.gas      = ct.Solution('../../input/c2h4red.xml')
       self.caseName = cname
       self.Ldomain  = float(getInputFileParameter_func('../../input/odtParam.inp', 'domainLength'))
       self.nrlz     = len(glob.glob("../../data/"+self.caseName+"/data_*"))
       self.ntimes   = len(glob.glob("../../data/"+self.caseName+"/data_0/dmp_*.dat"))
       self.pressure = float(getInputFileParameter_func('../../input/streams.inp', 'pres'))
       self.rho_soot = 1850.

       dxmin         = float(getInputFileParameter_func('../../input/odtParam.inp', 'dxmin'))
       self.ngu      = np.floor(1.0/dxmin)     # dxmin = dx/Ldomain
       self.pos_u    = np.linspace(0,self.Ldomain,self.ngu)

   ###########################################################################################

   def computeReactionRates(self) :

       """
       """

       rs_m_all     = np.zeros([self.ngu, self.nmom, self.ntimes])
       rs_rms_all   = np.zeros([self.ngu, self.nmom, self.ntimes])

       for it in range(self.ntimes) :

           print "it = ", it

           fname = "../../data_py_" + self.caseName + "/data_py_" + "{0:0>4}".format(it+1) + ".npy"
           data = np.load(fname)

           posf   = data[:,1]
           posf = np.append(posf,self.Ldomain)

           dx = posf[1:]-posf[0:-1]
           i = np.where( dx < 0.0 )[0]
           dx[i] = self.Ldomain - posf[i]

           self.irlz_start = np.where( posf == 0.0 )[0]
           self.irlz_end   = self.irlz_start[1:]-1
           self.irlz_end   = np.append(self.irlz_end, len(dx-1))

           #-----------

           Yg  = data[:,14:33]
           T   = data[:,7]
           rho = data[:,2]
           Ys  = data[:,36]
           pos = data[:,1]

           #-----------

           rs = np.zeros(len(dx))

           for irlz in range(len(self.irlz_start)) :

               irs = self.irlz_start[irlz]
               ire = self.irlz_end[irlz] + 1
               rs[irs:ire], qp[irs:ire] = self.getRadiativeFluxes_minus_plus(kgs[irs:ire], T[irs:ire], posf[irs:ire])

           #----------- compute mean and rms quantities

           rs_m, rs_m2     = self.computeMeanAndMean2OverRealizations(rs,   pos)

           rs_rms   = (rs_m2    - rs_m**2.0)**0.5

           #----------- save results

           rs_m_all[:,it]   = rs_m

           rs_rms_all[:,it]   = rs_rms

        #----------------------------------------------------------------------
        #------------ output results

       times = np.loadtxt("../../input/dumpTimes.inp", skiprows=1)
       header = " x(m), time="
       for i in range(len(times)) :
           header = header + str(times[i]) + ", "

       var = np.column_stack((self.pos_u, rs_m_all))
       fname = "rs_m.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((self.pos_u, rs_rms_all))
       fname = "rs_rms.dat"
       np.savetxt(fname, var, header=header)


   ###########################################################################################

   def computeMeanAndMean2OverRealizations(self, var, pos) :

       var_m  = np.zeros(self.ngu)
       var_m2 = np.zeros(self.ngu)

       for irlz in range(self.nrlz) :
           irs      = self.irlz_start[irlz]
           ire      = self.irlz_end[irlz] + 1
           f_linear = interp1d(pos[irs:ire], var[irs:ire], bounds_error=False, fill_value=0.0)
           var_u    = f_linear(self.pos_u)
           var_m    = var_m  + var_u
           var_m2   = var_m2 + var_u * var_u

       var_m  = var_m  / self.nrlz
       var_m2 = var_m2 / self.nrlz
       return (var_m, var_m2)

   ###########################################################################################


###########################################################################################

rt = reactionRates(sys.arsv[1])
rt.computeReactionRates()














