
from __future__ import division
import os
import glob
import numpy as np
import subprocess as sp
from data_tools import get_nRlz, get_domainBounds, compute_pdf, commentHdr

##############################################################################################

def eddyMaps(DI) :

    nToDo = 1 
    nRlz  = get_nRlz(DI)
    x0, xL = get_domainBounds(DI)

    if not os.path.exists(DI['pdir']+"eddyMaps") :
        os.mkdir(DI['pdir']+"eddyMaps")

    gpfile = open(DI['pdir']+"eddyMaps/plotEddyMaps.gnu",'w')
    gpfile.write("set bar 0\n")
    gpfile.write("set xlabel('Position (m)')\n")
    gpfile.write("set ylabel('Time or space (s or m)')\n")
    gpfile.write("set xrange [%g:%g]\n" %(x0,xL))

    flist = glob.glob(DI['cdir']+"runtime/runtime*")

    for ii,filename in enumerate(flist[0:nToDo]) :

        print("eddyMaps: writing for realization %i of %i (with %i total possible)" %(ii+1,nToDo,nRlz))

        gpfile.write("plot 'eddy_map_"+"{0:0>5}".format(ii)+".dat' us 1:2:3 with xerrorbars ti 'R_%i'; pause -1\n" %(ii))

        lines = []
        ##     the below doesn't work for Hewson's version
        theGrep = sp.getoutput("grep -h '^ *[1-9]' " + filename).split('\n')
        ## Instead the below seems to be a new way to do this
        #p_out = sp.Popen(["grep", "-h", "^ *[1-9]" , filename], stdout=sp.PIPE)
        #theGrep = p_out.communicate()[0].split('\n')[0:-1]  # the [0] arg is output, [1] is stderr

        for i in range(len(theGrep)) :
            lines.append(theGrep[i])

        ### check for zero length in 'lines' in case there were no eddies
        if len(lines) > 1:
            x = np.empty((len(lines),len(lines[0].split()) ))
            for i,line in enumerate(lines) :
                x[i,:] = np.fromstring(line, sep=" ")

            x = x[:,(6,1,5)]         # edPos, time, edSize
            x[:,2] = x[:,2] / 2.0    # fix "error bar" for gnuplot plotting

            fname = DI['pdir']+"eddyMaps/eddy_map_"+"{0:0>5}".format(ii)+".dat"
            np.savetxt(fname,x,fmt="%12.5e", comments=commentHdr)

    gpfile.close()

##############################################################################################

def eddyStats(DI) :

#    grp = sp.check_output(["grep -h '^ *[1-9]' "+ DI['cdir']+"runtime/runtime_*"],shell=True).split(str="'\n'")[:-1]
    grp = sp.getoutput("grep -h '^ *[1-9]' " + DI['cdir']+"runtime/runtime_*").split('\n')

    arr = np.empty( ( len(grp), len(grp[0].split()) ) )

    for i in range(len(grp)) :
        arr[i,:] = np.fromstring(grp[i], sep=" ")

    L = arr[:,5]
    y0 = arr[:,6]

    neddies = np.size(L) / len(glob.glob(DI['cdir']+"runtime/runtime_*"))
    Lavg = np.mean(L)
    Lrms = np.std(L)
    y0avg = np.mean(y0)
    y0rms = np.std(y0)

    nbins = 100
    P_y0, y0_bins = compute_pdf(y0, np.min(y0), np.max(y0), nbins)
    P_logL,  logL_bins = compute_pdf(np.log10(L), np.min(np.log10(L)), np.max(np.log10(L)), nbins)
    #P_L,  L_bins  = compute_pdf(L,  np.min(L),  np.max(L),   nbins)

    #fname = "../../data/"+caseN + "/post/eddyStats.dat"
    fname = DI['pdir']+"eddyStats.dat"
    with open(fname, 'w') as ofile :
        ofile.write("# neddies = %i\n" %(neddies))
        ofile.write("# Lavg    = %f\n" %(Lavg))
        ofile.write("# Lrms    = %f\n" %(Lrms))
        ofile.write("# y0avg   = %f\n" %(y0avg))
        ofile.write("# y0rms   = %f\n" %(y0rms))
        ofile.write("#\n")
        ofile.write("# y0, PDF_y0, log10L, PDF_log10L\n")
        for i in range(nbins) :
            ofile.write("%12.5e  %12.5e  %12.5e  %12.5e \n" %(y0_bins[i],P_y0[i],logL_bins[i],P_logL[i]))


