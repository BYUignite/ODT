#scatter plots of the FWHM growth rate (slope) vs. C/ZLES, C, and ZLES 
#plot file directory/name at end of each function definition

from __future__ import division
import numpy as np
from data_tools import get_inputFileParameter
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from data_tools import commentHdr

#---------------------------------------------------------------------------------------------

def CZLES_growth(DI): #uses list of cases

    cases=[]
    plotname="CZLES_growth"
    
    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()

    growth = np.zeros((len(DI), 1))
    CZLES = np.zeros((len(DI),1))

    for i in range(0, len(DI)):
        
        cases.append(DI[i]['cn'][8:])
        plotname+='__'+DI[i]['cn'][8:]
        fname = DI[i]['pdir']+"fwhm_cl_uvel.dat"
        odt = np.loadtxt(fname, comments=commentHdr)
        m,b = np.polyfit(odt[4:,0],odt[4:,1],1) #line of fit starting at x/d=20
        D = get_inputFileParameter(DI[i], ("initParams", "djeti"))
        x0 = (-b/m)*D
        C = get_inputFileParameter(DI[i], ("params", "C_param"))
        Z_LES = get_inputFileParameter(DI[i], ("params","Z_LES"))
        
        growth[i,0] = m
        CZLES[i,0] = C/Z_LES

    
    plt.scatter(CZLES[:,0], growth[:,0])
    plt.ylabel("FWHM growth")
    plt.xlabel("$C/Z_{ZLES}$")
    plt.title("coldJet")
    
    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'CZLES_growth__ALL'.replace(".","o"))

#------------------------------------------------------------------------------------------------

def C_growth(DI): #uses list of cases

    cases=[]
    plotname="C_growth"
    
    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()

    growth = np.zeros((len(DI), 1))
    C_params = np.zeros((len(DI),1))

    for i in range(0, len(DI)):
        
        cases.append(DI[i]['cn'][8:])
        plotname+='__'+DI[i]['cn'][8:]
        fname = DI[i]['pdir']+"fwhm_cl_uvel.dat"
        odt = np.loadtxt(fname, comments=commentHdr)
        m,b = np.polyfit(odt[4:,0],odt[4:,1],1) #line of fit starting at x/d=20
        D = get_inputFileParameter(DI[i], ("initParams", "djeti"))
        x0 = (-b/m)*D
        C = get_inputFileParameter(DI[i], ("params", "C_param"))
               
        growth[i,0] = m
        C_params[i,0] = C

    
    plt.scatter(C_params[:,0], growth[:,0])
    plt.ylabel("FWHM growth")
    plt.xlabel("$C$")
    plt.title("coldJet")
    
    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'C_growth__ALL'.replace(".","o"))

#---------------------------------------------------------------------------------------------

def ZLES_growth(DI): #uses list of cases

    cases=[]
    plotname="ZLES_growth"
    
    matplotlib.rcParams.update({'font.size':20, 'figure.autolayout': True}) #, 'font.weight':'bold'})

    fig, axL = plt.subplots()

    growth = np.zeros((len(DI), 1))
    ZLES_params = np.zeros((len(DI),1))

    for i in range(0, len(DI)):
        
        cases.append(DI[i]['cn'][8:])
        plotname+='__'+DI[i]['cn'][8:]
        fname = DI[i]['pdir']+"fwhm_cl_uvel.dat"
        odt = np.loadtxt(fname, comments=commentHdr)
        m,b = np.polyfit(odt[4:,0],odt[4:,1],1) #line of fit starting at x/d=20
        D = get_inputFileParameter(DI[i], ("initParams", "djeti"))
        x0 = (-b/m)*D
        Z_LES = get_inputFileParameter(DI[i], ("params","Z_LES"))
        
        growth[i,0] = m
        ZLES_params[i,0] = Z_LES

    
    plt.scatter(ZLES_params[:,0], growth[:,0])
    plt.ylabel("FWHM growth")
    plt.xlabel("$Z_{LES}$")
    plt.title("coldJet")
    
    #plt.savefig('../../data/plots_coldJet/'+plotname.replace(".","o"))
    plt.savefig('../../data/plots_coldJet/'+'ZLES_growth__ALL'.replace(".","o"))
