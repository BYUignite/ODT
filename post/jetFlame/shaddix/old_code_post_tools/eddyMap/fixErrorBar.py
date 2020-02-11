#!/usr/bin/python

import glob
import os
import sys

os.chdir("eddySpaceTime")
for file in glob.glob("*.dat"):

    print "fixing error bars for file ", file
    ifile = open(file, 'r')
    lines = ifile.readlines()
    ifile.close()
    ofile = open(file,'w')
    for line in lines :
        x, y, z = line.split()
        z = float(z) / 2.0
        ofile.write('%12.5e %12.5e %12.5e\n' % (float(x), float(y), z))

    ofile.close()

