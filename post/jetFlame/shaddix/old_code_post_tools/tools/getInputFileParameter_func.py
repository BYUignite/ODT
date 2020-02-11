
import sys, math, string

def getInputFileParameter_func(fname, param) :

    try:
        ifile  = open(fname, 'r')

    except:
        print 'error opening file ', fname

    param = "(" + param
    lines = ifile.readlines()

    value = None
    for line in lines :
        if line.find(param) >= 1 :
            value = string.split(line)[0]
    if value is None :
        print 'error finding param: ', param
        return -1

    return value
