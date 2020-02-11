#need data/dmp files to have original timestamps (if copying files from another directory, use rsync, not cp)

import os 
import numpy as np
import sys
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

#---------------------------------------------------------------------
#returns number of files
def filecount(dir_name):
    # return the number of files in directory dir_name
    try:
        return len([f for f in os.listdir(dir_name) if os.path.isfile(os.path.join(dir_name, f))])
    except Exception, e:
        return None
#-----------------------------------------------------------------------
#called on a case and list of realizations (allRlz, or specify)
#lists of avg. time between specific dmp file #s, calculates avg. time to a specified dmp file
def dump_times(DI): #call on a case followed by 'allRlz', or a list of specific realizations

    nrlz = len(DI)
    n = np.zeros((nrlz-1,1)) #number of dmp files in each realization
    dmps = [] #will be list  of lists of dmp file names for each realization
    
    for rlz in range(1, nrlz):
        n[rlz-1,0] = filecount(DI[rlz]['ddir']) - 4 #3 subtract off number of files that aren't dmp.dat (odt_*, *gnu...)
        dmps.append([])
        rn = int(n[rlz-1,0])
        for i in range(0,rn):
            dmps[rlz-1].append(str(i))
            if len(dmps[rlz-1][i])==1:
                dmps[rlz-1][i]='0'+dmps[rlz-1][i]
    nmin = int(np.amin(n))
    nmax = int(np.amax(n))
    dmp_times = np.zeros((nrlz, nmax))
    dmp_dtimes = np.zeros((nrlz, nmax))
    total_dtimes = np.zeros((nrlz,1))

    for rlz in range(1,nrlz):
        ndmps = len(dmps[rlz-1])
        for dmp in range(0,ndmps):
            try:
                dmp_times[rlz-1,dmp] = float(os.stat(DI[rlz]['ddir']+'dmp_000' + dmps[rlz-1][dmp] + '.dat').st_mtime) 
            except:
                dmp_times[rlz-1,dmp] = dmp_times[rlz-1,dmp-1]

            dmp_dtimes[rlz-1,dmp] = dmp_times[rlz-1,dmp]-dmp_times[rlz-1,dmp-1]

        ###select definition of total_dtimes based on desired output
        #total_dtimes[rlz-1,0] = dmp_times[rlz-1,ndmps-1]-dmp_times[rlz-1,0] #total time to end dmp file
        #total_dtimes[rlz-1,0] = dmp_times[rlz-1,nmin-2]-dmp_times[rlz-1,0] #total time up to min number of dmp files
        total_dtimes[rlz-1,0] = dmp_times[rlz-1,5]-dmp_times[rlz-1,0] #total time up to specified  dmp file
    
    #get average of times between each  dmp file, excluding realizations without those dmp times
    avg_dmp_dtimes = np.zeros((nmax-1,1))

    for i in range(0,nmax-1):
        nonzero_dmps=[]
        for rlz in range(1,nrlz+1):
            if dmp_dtimes[rlz-1,i+1]!=0:
                nonzero_dmps.append(dmp_dtimes[rlz-1,i+1])
        avg_dmp_dtimes[i,0] = np.mean(nonzero_dmps)

    #get average of total time to either the end dmp file, the minimum dmp file reached by any realization, or up to a specified dmp file #
    avg_total_time = np.mean(total_dtimes[:,0])

    print avg_dmp_dtimes
    #print "average total time: ",avg_total_time
    #print "average time to dmp",nmin-1,":", avg_total_time
    print "average time to dmp",6,":", avg_total_time

#-----------------------------------------------------------------------
#called on a case and list of realizations (allRlz, or specify)
#calculates mean and variance of times between dmp*.dat files (averaged for all dmp times), generates pdf of times between dmp files
def time_stats(DI): #DI has dictionary of directories for case followed by dictionaries for rlz data directories
    
    nrlz = len(DI)
    mean_dts = np.zeros((nrlz,1))
    var_dts = np.zeros((nrlz,1))
    ndmps = np.zeros((nrlz,1))
    dmps = []

    for rlz in range(1,nrlz):
        ndmps[rlz-1,0] = filecount(DI[rlz]['ddir'])
        n = filecount(DI[rlz]['ddir'])
        dmps.append([])
        for i in range(0,n):
            dmps[rlz-1].append(str(i))
            if len(dmps[rlz-1][i])==1:
                dmps[rlz-1][i]='0'+dmps[rlz-1][i]
    
    maxdmp = np.amax(ndmps[:,0])
    times = np.zeros((nrlz,int(maxdmp),2)) #for last entry, index 0 is timestamp, index 1 is time difference
    
    for rlz in range(1,nrlz):

        for i in range(0, int(ndmps[rlz-1,0])):
            try:
                times[rlz-1,i,0] = float(os.stat(DI[rlz]['ddir']+'dmp_000' + dmps[rlz-1][i] + '.dat').st_mtime) 
            except:
                times[rlz-1,i,0] = times[rlz-1,i-1,0]
            
        for i in range(1, int(ndmps[rlz-1,0])):
            times[rlz-1,i,1] = times[rlz-1,i,0] - times[rlz-1,i-1,0] #use this for pdfs

        mean_dts[rlz,0] = np.mean(times[rlz-1,:,1])
        var_dts[rlz,0] = 0
        for i in range(0, int(ndmps[rlz-1,0])):
            var_dts[rlz,0] += ((mean_dts[rlz,0] - times[rlz-1,i,1])**2)/int(ndmps[rlz-1,0])
    
    
    mean_time = np.mean(mean_dts[:,0])
    var_time = np.mean(var_dts[:,0])

    print "mean total: ",mean_time, " seconds"
    print "variance total: ",var_time, " seconds"

    ###make pdfs------------------------------------------------------------------------------
    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})
    plt.hist(times[:,:,1], bins=20, range=(100,10000), normed=True) #start range sometime after 0--for realizations that don't reach the max dmp time, the time differences for dmp times not reached are set to 0, but don't want to include these
    plt.xlabel("$t$", fontsize=22)
    plt.savefig("../../data/"+caseN+"/post/time_pdf".replace(".","o"))
    
#-----------------------------------------------------------------------
#get case and list of all realizations to be used in statistics

DI=[]


try:
    caseN = sys.argv[1]
except: 
    raise ValueError("Include data directories.")
DI.append({'pdir':'../../data/'+caseN+'/post/',     \
       'ddir':'../../data/'+caseN+'/data/',     \
       'cdir':'../../data/'+caseN+'/',          \
       'cn':caseN})
if not os.path.exists(DI[0]['pdir']):
    os.mkdir(DI[0]['pdir'])

if sys.argv[2] == 'allRlz':
    data_rlz = []
    for i in range(0,10):
        data_rlz.append('data_0000'+str(i))
    for i in range(10,100):
        data_rlz.append('data_000'+str(i))
    for i in range(100,1000):
        data_rlz.append('data_00'+str(i))
    for rlz in data_rlz:
        if os.path.exists('../../data/'+caseN+'/data/'+rlz+'/'):
            DI.append({'ddir':'../../data/'+caseN+'/data/'+rlz+'/'})

else: #if typing out data_00* directories individually
    for i in range(2,len(sys.argv)):
        data_rlz = sys.argv[i]
        DI.append({'ddir':'../../data/'+caseN+'/data/'+data_rlz+'/'})
#----------------------------------------------------------------------

time_stats(DI)
dump_times(DI)
