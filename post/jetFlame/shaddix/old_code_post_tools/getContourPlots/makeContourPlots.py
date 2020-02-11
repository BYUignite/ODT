#!/usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys


caseN   = sys.argv[1]
ivar    = 8
ivarZ   = 6                # mixf
figName = 'test.png'

#----------------------------------------------------------

fname = '../' + caseN + '/' + '/means_' + str(ivar) + '.dat'

data  = np.loadtxt(fname)
x     = data[:,0]
x     = x-np.average(x)
var   = data[:,1:]

Z = np.loadtxt('../' + caseN + '/means_' + str(ivarZ) + '.dat')
Z = Z[:,1:]

times = np.loadtxt('../../input/dumpTimes.inp', skiprows=1)

uxt   = np.loadtxt('../' + caseN + '/uxt.dat')
uxt   = (uxt[1:] + uxt[0:-1])/2
dt    = times[1:] - times[0:-1]
ypos  = np.cumsum(dt*uxt)
ypos  = np.insert(ypos, 0, 0.0)

X,Y = np.meshgrid(x,ypos)

#------------------ plot

zstoic = 0.06374

aspectRatio = (np.max(ypos) - np.min(ypos))/(np.max(x)-np.min(x))
fig = plt.figure(figsize = (4,2.2*aspectRatio))
plt.rc("font", size=16)

C  = plt.contourf(X,Y,var.T,50)
ax = plt.gca()
ax.axis((np.min(x),np.max(x),np.min(ypos),np.max(ypos)))
ax.set_xbound(lower=np.min(x), upper=np.max(x))
ax.set_ybound(lower=np.min(ypos), upper=np.max(ypos))
ax.set_aspect('equal')
plt.xlabel('position (m)',fontsize=16)
plt.ylabel('axial position (m)',fontsize=16)
plt.colorbar()
plt.xticks((-0.2,0,0.2))
#plt.yticks(())
#ax.grid(True, color='w', linewidth=1, linestyle = ':')

plt.contour(X,Y,Z.T,levels=np.array([zstoic]), colors='black', linewidth=0.5)

#plt.show()
plt.savefig(figName)

