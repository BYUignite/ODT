# Call as: 
#      python changeInputParam.py caseName paramName paramValue paramName paramValue ...
# For example:
#      python changeInputParam.py case1 C_param 10
#      ...This will edit the input file ../data/case1/odt_input.yaml,
#      ...which is what the odt code reads

import sys

#-----------------------------

if __name__ == __main__:
    changeInputParam(sys.argv)

#-----------------------------

def changeInputParam(args):

    with open('../data/'+argv[1]+'/input/odt_input.yaml', 'r') as ifile:
        lines = ifile.readlines()
    
    for i,line in enumerate(lines):
        for iarg in range(0,len(argv),2):
            if argv[iarg] in line :
                print " In changeParam ",len(argv), iarg , argv[iarg] 
                words = line.split(None,2)
                print " Looking for ",argv[iarg] #,argv[iarg+1] 
                print " Found ",words
                words[1] = argv[iarg+1]
                line = "         ".join(words) 
                line = "    " + line  
                lines[i] = line
    
    with open('../data/'+argv[1]+'/input/odt_input.yaml', 'w') as ofile:
        ofile.writelines(lines)

