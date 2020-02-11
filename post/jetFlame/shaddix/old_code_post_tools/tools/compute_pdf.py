import numpy as np

def compute_pdf(x, xmin, xmax, nbins=100) :

    dxbin = (xmax-xmin)/nbins;
    xbins = np.linspace(xmin+dxbin/2, xmax-dxbin/2, num=nbins)

    ibin = np.floor((x-xmin)/dxbin)
    np.clip(ibin,0,nbins-1)

    P = np.zeros(nbins)

    for i in range(nbins) :
        P[i] = np.size(np.where(ibin==i))

    P = P/len(x)/dxbin;

    return(P,xbins)

