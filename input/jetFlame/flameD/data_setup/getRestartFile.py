
# coding: utf-8

# ## Setup restart file for Flame D initialization

# In[4]:

import cantera as ct
import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate   import interp1d


# ### Velocity profile

# In[12]:

ngrid = 2000
dj = 0.0072
dp = 0.0182
Ld = dj*100
rj = dj/2
rp = dp/2
print("rp = ", rp)
print("Ap = ", np.pi*(rp**2-rj**2))

Fac = 1.5
rp = np.sqrt(rj**2 + Fac*(rp**2-rj**2))
dp = 2*rp
print("rp = ", rp)
print("Ap = ", np.pi*(rp**2-rj**2))

bc = np.loadtxt('exp_U_bc.dat')
rbc = bc[:,0] * dj
rbc = np.append(rbc,Ld/2)
Ubc = bc[:,1]
Ubc = np.append(Ubc,0.9)
rbc = np.hstack([-np.flipud(rbc),rbc])
Ubc = np.hstack([ np.flipud(Ubc),Ubc])

dx = Ld/ngrid
r   = np.linspace(-Ld/2+dx/2,Ld/2-dx/2,ngrid)
rf  = np.linspace(-Ld,Ld-dx,ngrid)
finterp = interp1d(rbc,Ubc)
U   = finterp(r)
V = np.zeros(ngrid)
W = np.zeros(ngrid)

#plt.plot(r,U,'-')
#plt.xlim([-0.05,0.05])


# ### Composition profile

# In[6]:

gas = ct.Solution("ch4red.xml")

xj  = "CH4:0.25, O2:0.1575, N2:0.5925"
xp  = "N2:0.7342, O2:0.054, O:7.47E-4, H2:1.29E-4, H:2.48E-5, H2O:0.0942, CO:4.07E-3, CO2:0.1098, OH:0.0028"
xair = "N2:0.79, O2:0.21"

Tj = 294
Tp = 1880
Tair = 291

P  = float(int(0.993*101325))

WO = gas.atomic_weight(0)
WH = gas.atomic_weight(1)
WC = gas.atomic_weight(2)
WN = gas.atomic_weight(3)

gas.TPX = Tj,P,xj
hj = gas.enthalpy_mass
rhoj = gas.density
muj  = gas.viscosity
xj = gas.X
yj = gas.Y
yHj = gas.elemental_mass_fraction(1)
yCj = gas.elemental_mass_fraction(2)

gas.TPX = Tair,P,xair
hair = gas.enthalpy_mass
rhoair = gas.density
muair  = gas.viscosity
xair = gas.X
yair = gas.Y
yHair = gas.elemental_mass_fraction(1)
yCair = gas.elemental_mass_fraction(2)

gas.TPX = Tp,P,xp
hp = gas.enthalpy_mass
rhop = gas.density
mup  = gas.viscosity
xp = gas.X
yp = gas.Y
yHp = gas.elemental_mass_fraction(1)
yCp = gas.elemental_mass_fraction(2)
mixfp = (0.5*(yHp-yHair)/WH + 2*(yCp-yCair)/WC)/         (0.5*(yHj-yHair)/WH + 2*(yCj-yCair)/WC)



x = np.zeros((ngrid,gas.n_species))
y = np.zeros((ngrid,gas.n_species))
T = np.zeros(ngrid)
h = np.zeros(ngrid)
rho = np.zeros(ngrid)
mu  = np.zeros(ngrid)
mixf = np.zeros(ngrid)

for i in range(ngrid):
    if np.abs(r[i]) <= rj:
        x[i,:] = xj
        y[i,:] = yj
        T[i]   = Tj
        h[i]   = hj
        rho[i] = rhoj
        mu[i]  = muj
        mixf[i] = 1
    elif np.abs(r[i]) > rj and np.abs(r[i]) <= rp:
        x[i,:] = xp
        y[i,:] = yp
        T[i]   = Tp
        h[i]   = hp
        rho[i] = rhop
        mu[i]  = mup
        mixf[i] = mixfp
    else:
        x[i,:] = xair
        y[i,:] = yair
        T[i]   = Tair
        h[i]   = hair
        rho[i] = rhoair
        mu[i]  = muair
        mixf[i] = 0
print(mixfp)


# In[13]:

#plt.plot(r,mixf,'.-')
#plt.xlim([-0.02,0.02])
#plt.ylim([0,1.1])


# ### Write Restart File

# In[14]:

chi = np.zeros(ngrid)
hr = np.zeros(ngrid)
data = np.vstack([r, rf, rho, mu, U, V, W, T, mixf, chi, hr, y.T, h]).T

header = "# time = 0\n# Grid points = " + str(ngrid) + "\n# Domain Size = " + str(Ld) + "\n"
header += "# Pressure (Pa) = " + str(int(P)) + "\n"
header += "# 1_pos        2_posf          3_rho          4_dvisc        5_uvel         6_vvel         7_wvel         8_temp         9_mixf         10_chi         11_hr          "
for i in range(12, 12+gas.n_species):
    header += str(i) + "_" + gas.species_name(i-12) + "         "
header += "     " + str(12+gas.n_species) + "_h"

np.savetxt("restart.dat", data, header=header, comments="", fmt="%14.7E")


# ### Dump times

# In[15]:

xoverd = np.arange(5,85,5)
dmb = np.array([1,2,3,7.5])
xoverd = np.hstack([xoverd,dmb])
xoverd = np.sort(xoverd)
x = xoverd * dj
for a in x:
    print(a, "  #", a/dj)


# In[ ]:




# In[ ]:



