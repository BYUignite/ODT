import numpy as np
from scipy.interpolate import interp1d

x = np.array([1,2,3,7.5,15,30,45,60,75])
nx = len(x)
flist = ["D01.Yave", "D02.Yave", "D03.Yave", "D075.Yave", "D15.Yave", "D30.Yave", "D45.Yave", "D60.Yave", "D75.Yave"]
print(flist)
fwhm_y = np.zeros(nx)
for i,file in enumerate(flist):
    print(file)
    data = np.loadtxt(file)
    r = data[:,0]
    y = data[:,1]
    n = len(y)
    ymax = np.max(y)
    ymin = np.min(y)
    ymid = 0.5*(ymax - ymin)
    imax = np.where(y==ymax)[0][0]

    f_interp = interp1d(y[imax:],r[imax:])
    rmid = f_interp(ymid)
    fwhm_y[i] = 2*rmid

data = np.vstack((x,fwhm_y)).T
np.savetxt("fwhm_Z.dat",data, header=" x/d, fwhm(m)")


