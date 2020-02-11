# Call as: 
#      python changeInputParam.py caseName paramName paramValue paramName paramValue ...
# For example:
#      python changeInputParam.py case1 C_param 10
#      ...This will edit the input file ../data/case1/odt_input.yaml,
#      ...which is what the odt code reads

import sys

#-----------------------------

with open('../data/'+sys.argv[1]+'/input/odt_input.yaml', 'r') as ifile:
    lines = ifile.readlines()

for i,line in enumerate(lines):
    for iarg in range(0,len(sys.argv),2):
        if sys.argv[iarg] in line :
            print " In changeParam ",len(sys.argv), iarg , sys.argv[iarg] 
            words = line.split(None,2)
            print " Looking for ",sys.argv[iarg] #,sys.argv[iarg+1] 
            print " Found ",words
            words[1] = sys.argv[iarg+1]
            line = "         ".join(words) 
            line = "    " + line  
            lines[i] = line

with open('../data/'+sys.argv[1]+'/input/odt_input.yaml', 'w') as ofile:
    ofile.writelines(lines)

