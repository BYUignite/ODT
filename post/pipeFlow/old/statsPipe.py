
#from __future__ import division
import numpy as np
import glob as gb
import yaml
import sys
import matplotlib
matplotlib.use('PDF')       # or Agg (for png), SVG, PS
import matplotlib.pyplot as plt
import os


#--------------------------------------------------------------------------------------------

try :
    caseN = sys.argv[1]
except :
    raise ValueError("Include the case name in the call")

if not os.path.exists("../../data/"+caseN+"/post") :
    os.mkdir("../../data/"+caseN+"/post")

#-----------------------------------

def extrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    y = np.interp(x, xp, yp)
    y[x<xp[ 0]] = yp[ 0] + (x[x<xp[ 0]]-xp[ 0])*(yp[ 0]-yp[ 1])/(xp[ 0]-xp[ 1])
    y[x>xp[-1]] = yp[-1] + (x[x>xp[-1]]-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2])
    return y

#-----------------------------------

with open("../../data/"+caseN+"/input/odt_input.yaml") as ifile :
    y = yaml.load(ifile, Loader=yaml.FullLoader)
kvisc = y["params"]["kvisc0"]
delta = y["params"]["domainLength"] * 0.5
Retau = 1.0/kvisc

#flist = gb.glob('../../data/'+caseN+'/data/data_00000/odt_[1-9]????_*.dat')
flist = gb.glob('../../data/'+caseN+'/data/data_00000/dmp_*.dat')
print flist

nunif = 4000       # should be 1/dxmin

nfiles = len(flist)
xu = np.linspace(-delta,delta,nunif)
um  = np.zeros(nunif)
u2m = np.zeros(nunif)
vm  = np.zeros(nunif)
v2m = np.zeros(nunif)
wm  = np.zeros(nunif)
w2m = np.zeros(nunif)


for ifile in flist :
    print ifile
    data = np.loadtxt(ifile);
    x = data[:,0]
    u = data[:,2]
    v = data[:,3]
    w = data[:,4]

    uu = extrap(xu,x,u)
    vv = extrap(xu,x,v)
    ww = extrap(xu,x,w)

    um += uu
    vm += vv
    wm += ww
    u2m += uu*uu
    v2m += vv*vv
    w2m += ww*ww

um /= nfiles
vm /= nfiles
wm /= nfiles
um = 0.5*(um[:nunif/2] + np.flipud(um[nunif/2:]))
vm = 0.5*(vm[:nunif/2] + np.flipud(vm[nunif/2:]))
wm = 0.5*(wm[:nunif/2] + np.flipud(wm[nunif/2:]))

u2m /= nfiles
v2m /= nfiles
w2m /= nfiles
u2m = 0.5*(u2m[:nunif/2] + np.flipud(u2m[nunif/2:]))
v2m = 0.5*(v2m[:nunif/2] + np.flipud(v2m[nunif/2:]))
w2m = 0.5*(w2m[:nunif/2] + np.flipud(w2m[nunif/2:]))

urms = np.sqrt(u2m - um*um)
vrms = np.sqrt(v2m - vm*vm)
wrms = np.sqrt(w2m - wm*wm)

xu += delta
xu = xu[:nunif/2]

dudx = (um[1]-um[0])/(xu[1]-xu[0])
utau = np.sqrt(kvisc * np.abs(dudx))
RetauOdt = utau * delta / kvisc

xu *= utau    # OPTIONAL: utau should be 1
xu /= kvisc

data = np.vstack([xu,um,vm,wm,urms,vrms,wrms]).T
fname = "../../data/"+caseN+"/post/ODTstat.dat"
np.savetxt(fname, data, header="x, um, vm, wm, urms, vrms, urms")

print "nominal Retau: ", Retau
print "nominal  utau: ", 1
print "actual  Retau: ", RetauOdt
print "actual   utau: ", utau

#--------------------------------------------------------------------------------------------
# make the plot 1
print "MAKING DNS PLOT FOR VARIOUS RE"

dns_1 = np.loadtxt("DNS_raw/180_Re_1.dat", comments='%')
dns_2 = np.loadtxt("DNS_raw/360_Re_1.dat", comments='%')
dns_3 = np.loadtxt("DNS_raw/550_Re_1.dat", comments='%')
dns_4 = np.loadtxt("DNS_raw/1000_Re_1.dat", comments='%')

matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

fig, ax = plt.subplots()

ax.semilogx(dns_1[:,1],dns_1[:,2], 'r:',  label=r'$\mathit{Re}_{\tau}=180$')
ax.semilogx(dns_2[:,1],dns_2[:,2], 'g-.', label=r'$\mathit{Re}_{\tau}=360$')
ax.semilogx(dns_3[:,1],dns_3[:,2], 'b--', label=r'$\mathit{Re}_{\tau}=550$')
ax.semilogx(dns_4[:,1],dns_4[:,2], 'k-',  label=r'$\mathit{Re}_{\tau}=1000$')

ax.set_title(r'DNS')
ax.set_xlabel(r'$(R-r)^+$') #, fontsize=22)
ax.set_ylabel(r'$U_z^+$') #, fontsize=22)
ax.legend(loc='upper left', frameon=True, fontsize=16)
ax.set_ylim([0.,30.])
ax.set_xlim([0.1,1000.])

plt.savefig("../../data/"+caseN+"/post/UmeanDns")

#--------------------------------------------------------------------------------------------
# make the plot 2

print("MAKING UMEAN PLOT FOR DNS RE=1000: make sure that is the Re of your case.")
dns_4 = np.loadtxt("DNS_raw/1000_Re_1.dat", comments='%')
odt   = np.loadtxt("../../data/"+caseN+"/post/ODTstat.dat")

matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

fig, ax = plt.subplots()

ax.semilogx(  odt[:,0],  odt[:,1], 'k-',  label=r'ODT')
ax.semilogx(dns_4[:,1],dns_4[:,2], 'k--', label=r'DNS')

ax.set_title(r'$\mathit{Re}_\tau=1000$')
ax.set_xlabel(r'$(R-r)^+$') #, fontsize=22)
ax.set_ylabel(r'$U_z^+$') #, fontsize=22)
ax.legend(loc='upper left', frameon=True, fontsize=16)
ax.set_ylim([0.,30.])
ax.set_xlim([0.1,1000.])

plt.savefig("../../data/"+caseN+"/post/Umean_Re1000")

#--------------------------------------------------------------------------------------------
# make the plot 3

print("MAKING UPROF PLOT FOR DNS RE=1000: make sure that is the Re of your case.")
dns_4 = np.loadtxt("DNS_raw/1000_Re_1.dat", comments='%')
odt   = np.loadtxt("../../data/"+caseN+"/post/ODTstat.dat")

matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

fig, ax = plt.subplots()

Xodt   = np.append(  odt[:,0]*kvisc/utau-delta, delta-np.flipud(  odt[:,0])*kvisc/utau)
Xdns_4 = np.append(dns_4[:,1]*kvisc/utau-delta, delta-np.flipud(dns_4[:,1])*kvisc/utau)

Uodt   = np.append(  odt[:,1], np.flipud(  odt[:,1]))*utau
Udns_4 = np.append(dns_4[:,2], np.flipud(dns_4[:,2]))*utau

ax.plot(Xodt  , Uodt  , 'k-',  label=r'ODT')
ax.plot(Xdns_4, Udns_4, 'k--', label=r'DNS')

##ax.plot(  odt[:,0]*kvisc/utau,   odt[:,1]*utau, 'k-',  label=r'ODT')
##ax.plot(dns_4[:,1]*kvisc/utau, dns_4[:,2]*utau, 'k--', label=r'DNS')

ax.set_title(r'pipe flow')
ax.set_xlabel(r'$r$') #, fontsize=22)
ax.set_ylabel(r'$U_z$') #, fontsize=22)
ax.legend(loc='upper left', frameon=True, fontsize=16)
ax.set_ylim([0., 30.])
ax.set_xlim([-delta, delta])

plt.savefig("../../data/"+caseN+"/post/Uprof_Re1000")

#--------------------------------------------------------------------------------------------
# make the plot 4

print("MAKING RMS PLOT FOR DNS RE=1000: make sure that is the Re of your case.")
dns = np.loadtxt("DNS_raw/DNS_1000_reversed.dat", comments='%')
odt = np.loadtxt("../../data/"+caseN+"/post/ODTstat.dat")

matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

fig, ax = plt.subplots()

ax.plot(odt[:,0], odt[:,4], 'k-' ,label='')
ax.plot(odt[:,0], odt[:,5], 'b--',label='')
ax.plot(odt[:,0], odt[:,6], 'r:' ,label='')
ax.plot(dns[:,1], dns[:,6], 'k-' ,label=r'$u_{z,\mathrm{rms}}$')
ax.plot(dns[:,1], dns[:,4], 'b--',label=r'$u_{r,\mathrm{rms}}$')
ax.plot(dns[:,1], dns[:,5], 'r:' ,label=r'$u_{\theta,\mathrm{rms}}$')

color = (0.6,0.6,0.6)
ax.plot((0,0),(0,3.5), c=color, label='')
#ax.lines[6].set_color(color)
ax.arrow(-450,1.7,-300,0, head_width=.1, head_length=50, fc=color, ec=color)
ax.arrow( 450,1.7, 300,0, head_width=.1, head_length=50, fc=color, ec=color)
ax.text(-700,1.9,r'DNS', fontsize=20, color=color)
ax.text(500, 1.9,r'ODT', fontsize=20, color=color)

ax.set_xlabel(r'$(R-r)^+$') #, fontsize=22)
ax.set_ylabel(r'$u_{RMS}^+$') #, fontsize=22)
ax.legend(loc='upper right', frameon=True, fontsize=16)
ax.set_ylim([0.,3.5])
ax.set_xlim([-1000.,1000.])

plt.savefig("../../data/"+caseN+"/post/uvwRms_Re1000")

