from __future__ import division
import numpy as np
from scipy.interpolate import interp1d
import yaml

#commentHdr = "VARIABLES = "
commentHdr = "#"

##########################################################################

def get_data_realization(DI, itime, iRlz) :
    """
    Get the data for a single realization for a single time
    DI: input: used to find the right directory
    itime: input: which time we want
    iRlz: input: which time we want
    """

    fname = DI['cdir']+"data_py/data_py_" + "{0:0>5}".format(itime) + ".npy"
    data = np.load(fname)

    posf   = data[:,1]
    dx = posf[1:]-posf[0:-1]

    istarts = np.where(dx <0.0)[0] + 1
    istarts = np.insert(istarts, 0, 0.0)

    iends = np.where(dx < 0.0)[0]
    iends = np.append(iends, len(posf)-1)

    i_s = istarts[iRlz]
    i_e = iends[iRlz]

    return data[i_s:i_e+1, :]

##########################################################################

def get_nRlz(DI) :

    fname = DI['cdir']+"data_py/data_py_00000.npy"

    data = np.load(fname)

    posf   = data[:,1]
    dx = posf[1:]-posf[0:-1]

    istarts = np.where(dx <0.0)[0] + 1
    istarts = np.insert(istarts, 0, 0.0)

    return np.size(istarts)

##########################################################################

def get_dataHeaderVars(DI) :

    fname = DI['cdir']+"data_py/header.dat"
    ifile = open(fname,'r')
    varNames = ifile.readline().split()[1:]
    ifile.close()

    for i,name in enumerate(varNames) :
        varNames[i] = varNames[i][name.index("_")+1:]

    return varNames

##########################################################################

def get_axialLocations(DI, forceGetDumpTimes=False) :

    if forceGetDumpTimes:
        return np.array(get_inputFileParameter(DI, ("dumpTimes",)))

    Lspatial = get_inputFileParameter(DI, ("params", "Lspatial"))
    if Lspatial:
        yt = get_inputFileParameter(DI, ("dumpTimes",))
        return np.array(yt)
    else:
        yt = np.loadtxt(DI['pdir']+"ytu.dat", comments=commentHdr)
        return yt[:,0]

##########################################################################

def get_inputFileParameter(DI, paramNest) :

    #dirname = os.getcwd().split("/")[-1] + "/"
    with open(DI['cdir']+"input/input.yaml") as yfile :
        yfile = yaml.load(yfile, Loader=yaml.FullLoader)
    if len(paramNest) == 1 :
        return yfile[paramNest[0]]
    elif len(paramNest) == 2 :
        return yfile[paramNest[0]][paramNest[1]]
    elif len(paramNest) == 3 :
        return yfile[paramNest[0]][paramNest[1]][paramNest[2]]
    else :
        raise ValueError("paramNest is too big in get_inputFileParameter")

##########################################################################

def get_fstoic(DI) :

    fname = DI['cdir']+"runtime/runtime_00000"
    with open(fname) as ifile :
        lines = ifile.readlines()
        for line in lines :
            if "mixfStoic" in line :
                return float(line.split()[4])

##########################################################################

def get_domainBounds(DI) :

    data = get_data_realization(DI, 0, 0)
    x0 = data[0,1]
    Ld = get_inputFileParameter(DI, ("params", "domainLength"))
    xL = x0 + Ld

    return x0, xL


##########################################################################


def extrap1d(x,y) :

    xs = x[0]-1000.0*(x[-1]-x[0])
    xe = x[-1]+1000.0*(x[-1]-x[0])
    ys = y[0]+(xs-x[0])*(y[1]-y[0])/(x[1]-x[0])
    ye = y[-1]+(xe-x[-1])*(y[-1]-y[-2])/(x[-1]-x[-2])

    x = np.insert(x,0,xs)
    x = np.append(x,xe)
    y = np.insert(y,0,ys)
    y = np.append(y,ye)

    return interp1d(x,y)

# def extrap1d(interpolator) :
#     '''
#     Found this online. Its clever, but SLOW
#     '''
#
#     x = interpolator.x
#     y = interpolator.y
#
#     def interp_extrap_1pt(xpt):
#         if xpt < x[0]:
#             return y[0]+(xpt-x[0])*(y[1]-y[0])/(x[1]-x[0])
#         elif xpt > x[-1]:
#             return y[-1]+(xpt-x[-1])*(y[-1]-y[-2])/(x[-1]-x[-2])
#         else:
#             return interpolator(xpt)
#
#     def my_interp_extrap(X):
#         return np.array(map(interp_extrap_1pt, np.array(X)))
#
#     return my_interp_extrap


##########################################################################

def compute_pdf(x, xmin, xmax, nbins=100) :

    dxbin = (xmax-xmin)/nbins;
    xbins = np.linspace(xmin+dxbin/2, xmax-dxbin/2, num=nbins)

    ibin = np.floor((x-xmin)/dxbin)
    np.clip(ibin,0,nbins-1)

    P = np.zeros(nbins)

    for i in range(nbins) :
        P[i] = np.size(np.where(ibin==i))

    if len(x) > 0:
        P = P/len(x)/dxbin;
    else:
        P = 0.0;

    return(P,xbins)

##########################################################################

def compute_wpdf(x, w, xmin, xmax, nbins=100) :

    dxbin = (xmax-xmin)/nbins;
    xbins = np.linspace(xmin+dxbin/2, xmax-dxbin/2, num=nbins)

    ibin = np.floor((x-xmin)/dxbin)
    np.clip(ibin,0,nbins-1)

    P = np.zeros(nbins)

    for i in range(nbins) :
        ii = np.where(ibin==i)
        if ( np.transpose(ii).size > 0 ):
            P[i] = np.sum(w[ii])
        else:
            P[i] = 0.0

    if (np.sum(P) > 0):
        P = P/np.sum(P)/dxbin;
    else:
        P = np.zeros(nbins)

    return(P,xbins)


