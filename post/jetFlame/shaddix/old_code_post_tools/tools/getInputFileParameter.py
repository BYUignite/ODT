#!/usr/bin/python

import sys, math, string

try:
    ifile  = open(sys.argv[1], 'r')
    param  = sys.argv[2]

except:
    print 'useage: %s ifile' % (sys.argv[0])

param = "(" + param
lines = ifile.readlines()

value = None
for line in lines :
    if line.find(param) >= 1 :
        value = string.split(line)[0]
if value is None :
    exit(1)

print value
