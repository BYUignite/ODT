import numpy as np
from scipy.interpolate import interp1d

x = np.array([5,10,20,40,60,80])
nx = len(x)
flist = ('exp_r_Z_05.dat', 'exp_r_Z_10.dat', 'exp_r_Z_20.dat', 'exp_r_Z_40.dat', 'exp_r_Z_60.dat', 'exp_r_Z_80.dat')
print(flist)
fwhm_y = np.zeros(nx)
for i,file in enumerate(flist):
    data = np.loadtxt(file)
    r = data[:,0]
    y = data[:,1]
    n = len(y)

    ymid = 0.5*(np.max(y) + np.min(y))
    f_interp = interp1d(y[int(n/2):],r[int(n/2):])
    rmid = f_interp(ymid)
    print(rmid, ymid)
    fwhm_y[i] = 2*rmid

    data = np.vstack((x,fwhm_y)).T
    np.savetxt("fwhm_Z.dat",data)

#-----------------------------------------------------

x = np.array([0,5,10,20,40,60,80])
nx = len(x)
flist = ('exp_r_v_0.dat', 'exp_r_v_05.dat', 'exp_r_v_10.dat', 'exp_r_v_20.dat', 'exp_r_v_40.dat', 'exp_r_v_60.dat', 'exp_r_v_80.dat')
print(flist)
fwhm_y = np.zeros(nx)
for i,file in enumerate(flist):
    data = np.loadtxt(file)
    r = data[:,0]
    y = data[:,1]
    n = len(y)

    ymid = 0.5*(np.max(y) + np.min(y))
    f_interp = interp1d(y[int(n/2):],r[int(n/2):])
    rmid = f_interp(ymid)
    print(rmid, ymid)
    fwhm_y[i] = 2*rmid/0.008

    data = np.vstack((x,fwhm_y)).T
    np.savetxt("fwhm_u.dat",data)

