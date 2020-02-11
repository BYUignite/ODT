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

class radiationProperties :

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

   def computeRadiationProperties(self) :

       """
       """

       qm_m_all   = np.zeros([self.ngu, self.ntimes])
       qp_m_all   = np.zeros([self.ngu, self.ntimes])
       qnet_m_all = np.zeros([self.ngu, self.ntimes])
       kg_m_all   = np.zeros([self.ngu, self.ntimes])
       ks_m_all   = np.zeros([self.ngu, self.ntimes])

       qm_rms_all   = np.zeros([self.ngu, self.ntimes])
       qp_rms_all   = np.zeros([self.ngu, self.ntimes])
       qnet_rms_all = np.zeros([self.ngu, self.ntimes])
       kg_rms_all   = np.zeros([self.ngu, self.ntimes])
       ks_rms_all   = np.zeros([self.ngu, self.ntimes])

       Lopt_m_all   = np.zeros(self.ntimes)
       Lopt_rms_all = np.zeros(self.ntimes)

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

           #----------- get gas and soot absorption coefficients

           kg  = np.zeros(len(dx))
           ks  = np.zeros(len(dx))
           kgs = np.zeros(len(dx))

           Yg  = data[:,14:33]
           T   = data[:,7]
           rho = data[:,2]
           Ys  = data[:,36]
           pos = data[:,1]

           for i in range(len(dx)) :
               fv = rho[i]*Ys[i]/self.rho_soot
               kg[i], ks[i] = self.getAbsorptionCoefficient_gas_soot(Yg[i,:], T[i], fv)
               kgs[i] = kg[i] + ks[i]

           #----------- get radiative heat flux fields

           self.TBClo  = float(getInputFileParameter_func('../../input/bc.inp', 'freeStreamTemp'))
           self.TBChi  = float(getInputFileParameter_func('../../input/bc.inp', 'freeStreamTemp'))

           qp = np.zeros(len(dx))
           qm = np.zeros(len(dx))

           for irlz in range(len(self.irlz_start)) :

               irs = self.irlz_start[irlz]
               ire = self.irlz_end[irlz] + 1
               qm[irs:ire], qp[irs:ire] = self.getRadiativeFluxes_minus_plus(kgs[irs:ire], T[irs:ire], posf[irs:ire])

           qnet = qp - qm

           #----------- get optical thickness (nondimensional)

           Lopt = np.zeros(self.nrlz)

           for irlz in range(self.nrlz) :
               irs = self.irlz_start[irlz]
               ire = self.irlz_end[irlz] + 1

               Lopt[irlz] = np.sum(dx[irs:ire]*kgs[irs:ire])

           #----------- compute mean and rms quantities

           qm_m, qm_m2     = self.computeMeanAndMean2OverRealizations(qm,   pos)
           qp_m, qp_m2     = self.computeMeanAndMean2OverRealizations(qp,   pos)
           qnet_m, qnet_m2 = self.computeMeanAndMean2OverRealizations(qnet, pos)
           Lopt_m          = np.mean(Lopt)
           kg_m, kg_m2     = self.computeMeanAndMean2OverRealizations(kg,   pos)
           ks_m, ks_m2     = self.computeMeanAndMean2OverRealizations(ks,   pos)

           qm_rms   = (qm_m2    - qm_m**2.0)**0.5
           qp_rms   = (qp_m2    - qp_m**2.0)**0.5
           qnet_rms = (qnet_m2  - qnet_m**2.0)**0.5
           Lopt_rms = np.std(Lopt)
           kg_rms   = (kg_m2    - kg_m**2.0)**0.5
           ks_rms   = (ks_m2    - ks_m**2.0)**0.5

           #----------- save results

           qm_m_all[:,it]   = qm_m
           qp_m_all[:,it]   = qp_m
           qnet_m_all[:,it] = qnet_m
           Lopt_m_all[it]   = Lopt_m
           kg_m_all[:,it]   = kg_m
           ks_m_all[:,it]   = ks_m

           qm_rms_all[:,it]   = qm_rms
           qp_rms_all[:,it]   = qp_rms
           qnet_rms_all[:,it] = qnet_rms
           Lopt_rms_all[it]   = Lopt_m
           kg_rms_all[:,it]   = kg_rms
           ks_rms_all[:,it]   = ks_rms

        #----------------------------------------------------------------------
        #------------ output results

       times = np.loadtxt("../../input/dumpTimes.inp", skiprows=1)
       header = " x(m), time="
       for i in range(len(times)) :
           header = header + str(times[i]) + ", "

       var = np.column_stack((self.pos_u, qm_m_all))
       fname = "qm_m.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((self.pos_u, qp_m_all))
       fname = "qp_m.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((self.pos_u, qnet_m_all))
       fname = "qnet_m.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((self.pos_u, kg_m_all))
       fname = "kg_m.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((self.pos_u, ks_m_all))
       fname = "ks_m.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((self.pos_u, qm_rms_all))
       fname = "qm_rms.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((self.pos_u, qp_rms_all))
       fname = "qp_rms.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((self.pos_u, qnet_rms_all))
       fname = "qnet_rms.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((self.pos_u, kg_rms_all))
       fname = "kg_rms.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((self.pos_u, ks_rms_all))
       fname = "ks_rms.dat"
       np.savetxt(fname, var, header=header)

       var = np.column_stack((times[0:self.ntimes], Lopt_m_all, Lopt_rms_all))
       header = "time Lopt_mean Lopt_rms (nondim)"
       fname = "Lopt.dat"
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

   def getRadiativeFluxes_minus_plus(self, kgs, T, pf) :

       """
       Comput q- and q+ for two-flux radiation model for a single realization.
       kgs is the gas + soot absorption coefficient at ODT grid points for the realization
       T   is the gas temperature
       pf  is the posf of the odt cells (but without the last face, which is appended below
       """

       sigmaSB = 5.670E-8             # W/m2*K4
       npts    = len(kgs+2)

       kabs = kgs.copy()
       kabs = np.insert(kabs,0,kabs[0])
       kabs = np.append(kabs,kabs[-1])

       posf = pf.copy()
       posf = np.append(posf, self.Ldomain)
       dx = np.zeros(npts-1)
       dx[0] = 0.5*(posf[1]-posf[0])
       for i in range(1,npts-2) :
           dx[i]  = 0.5*(posf[i+1]-posf[i-1])
       dx[npts-2] = 0.5*(posf[-1]-posf[-2])

       emmTerm = np.zeros(npts)
       for i in range(npts-2) :
           emmTerm[i+1] = 2.0*kabs[i+1] * sigmaSB*T[i]**4.0
       emmTerm[0]  = emmTerm[1]
       emmTerm[-1] = emmTerm[-2]

       qpBClo = sigmaSB * self.TBClo**4.0
       qmBChi = sigmaSB * self.TBChi**4.0

       qp = np.zeros(npts)
       qm = np.zeros(npts)

       qp[0] = qpBClo
       for i in range(1,npts) :
           qp[i] = ( qp[i-1] + dx[i-1]*emmTerm[i] ) / ( 1. + 2.*kabs[i]*dx[i-1] )

       qm[-1] = qmBChi
       for i in range(npts-2, -1, -1) :
           qm[i] = ( qm[i+1] + dx[i]*emmTerm[i] ) / ( 1. + 2.*kabs[i]*dx[i] )

       return (qm, qp)



   ###########################################################################################

   def getAbsorptionCoefficient_gas_soot(self, Ysp, T, fvSoot=-1) :

       """
       Comput the gas and soot absorption coeffients for a single point
       Ysp is the species array
       T   is the temperature
       fvSoot is the volume fraction (raw, not ppm (without the *1E6 factor))
       returns kg and ks (units are 1/m)
       """

       nRadSp     = 4                    # CH4, CO2, H2O, CO
       radCoefs   = np.zeros([4,6])
       sootFactor = 1863.0

       radCoefs[0,0] =  1.017015E+1    # ch4; kp=2.798 at 1150K
       radCoefs[0,1] = -7.947312E-03
       radCoefs[0,2] =  4.342446E-7
       radCoefs[0,3] =  1.048611E-9
       radCoefs[0,4] = -2.287861E-13
       radCoefs[0,5] =  0.000000E+0
       radCoefs[1,0] =  3.24442E+1     # co2; kp=29.197 at 925K
       radCoefs[1,1] =  7.537513E-02
       radCoefs[1,2] = -1.535140E-04
       radCoefs[1,3] =  9.48794E-8
       radCoefs[1,4] = -2.509259E-11
       radCoefs[1,5] =  2.447995E-15
       radCoefs[2,0] =  6.86948E+1     # h2o; kp=4.474  at 1119K
       radCoefs[2,1] = -1.52349E-01
       radCoefs[2,2] =  1.417848E-04
       radCoefs[2,3] = -6.620996E-8
       radCoefs[2,4] =  1.52415E-11
       radCoefs[2,5] = -1.373456E-15
       radCoefs[3,0] =  1.56536E+0     # co; kp=2.501 at 1007 K
       radCoefs[3,1] =  1.483914E-02
       radCoefs[3,2] = -2.656035E-05
       radCoefs[3,3] =  1.68798E-8
       radCoefs[3,4] = -4.674473E-12
       radCoefs[3,5] =  4.767887E-16

       fmissing = False
       iRadIndx = np.zeros(nRadSp)

       isp = self.gas.species_index("CH4")
       isp if (isp > 0) else self.gas.species_index("ch4")
       iRadIndx[0] = isp
       if(isp < 0) : fmissing = True

       isp = self.gas.species_index("CO2")
       isp if (isp > 0) else self.gas.species_index("co2")
       iRadIndx[1] = isp
       if(isp < 0) : fmissing = True

       isp = self.gas.species_index("H2O")
       isp if (isp > 0) else self.gas.species_index("h2o")
       iRadIndx[2] = isp
       if(isp < 0) : fmissing = True

       isp = self.gas.species_index("CO")
       isp if (isp > 0) else self.gas.species_index("co")
       iRadIndx[3] = isp
       if(isp < 0) : fmissing = True

       if(fmissing) :
           print "\nWarning one or more radiating species missing from mechanism\n"

       self.gas.TPY = T, self.pressure, Ysp
       xMole   = self.gas.X

       Kgas = 0.0

       for k in range(nRadSp) :

           if(iRadIndx[k] < 0) :         # this radiation species k is not in the mechanism
               continue
           Kabs = radCoefs[k][5]
           for j in range(4,-1,-1) :
               Kabs = Kabs * T + radCoefs[k][j]

           Kgas += xMole[iRadIndx[k]]*self.pressure/101325.0*Kabs

       if(fvSoot > 0.0) :
           Ksoot = sootFactor * fvSoot * T
       else :
           Ksoot = 0.0

       return (Kgas, Ksoot)

###########################################################################################

rp = radiationProperties(sys.argv[1])
rp.computeRadiationProperties()














