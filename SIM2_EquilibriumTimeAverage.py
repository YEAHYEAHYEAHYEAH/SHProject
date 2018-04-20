#========== Script to find equilibrium time from txt files =============
#_______________________________________________________________________


#import modules
import math as m
import numpy as np
import matplotlib.pyplot as pl
from numpy import random as rd
import re
import operator
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import glob
import scipy.odr as odr


class AnalysedFiles(object):
    def __init__(self,dirPath,equib=False,equibVals = None,equibSigs=None):
        self.equibTimes = []
        self.Lks        = []
        self.equib      = equib
        self.dirPath    = dirPath
        if self.equib == False:    
            self.equibVals = []
            self.equibSigs = []
        else:                      
            self.equibVals = equibVals
            self.equibSigs = equibSigs
        self.getData(dirPath)
            
    def getData(self,dirPath):
        for fname in sorted(glob.glob(dirPath+'/*.txt')):
            title = fname[fname.find("/")+1:-4]
            dataFile = open(fname,'r')
            
            print dataFile
            rigidity    =  float(title[-8:-4])
            print rigidity
            exit()
            plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)

            initialLink =  int(title[5:-4])
            frameNo = plotArray[:,0]
            twists = plotArray[:,1]
            writhes = plotArray[:,2]        
            total = plotArray[:,1]+plotArray[:,2]
            
            equibTime = self.equilibriumRange(twists,initialLink)
            self.equibTimes.append(equibTime)

            self.Lks.append(initialLink)
            
    def equilibriumRange(self,xData,initialLink):
        if self.equib==False: 
            equibVal = np.mean(xData[20000:])
            self.equibVals.append(equibVal)
            
            equibSig = np.std(xData[20000:])
            self.equibSigs.append(equibSig)

        else: 
            equibVal = self.equibVals[initialLink-1]
            equibSig = self.equibSigs[initialLink-1]
        endRange = 0
        
        for i in xrange(len(xData)):
            if abs(equibVal - xData[i])<1.0*equibSig:
                endRange = i
                break
        if False:
            pl.plot(xData[:endRange])
            pl.title(self.dirPath + "\ Initial Twist: "+str(initialLink))
            pl.xlabel("Frame Number (50 Natural Time Units)")
            pl.ylabel("Twist Value")
            pl.savefig(self.dirPath+"/Equib{}InLk".format(initialLink), dpi=300, 
                       facecolor='w', edgecolor='w', orientation='landscape', 
                       papertype=None, format=None, transparent=False, 
                       bbox_inches='tight', pad_inches=0.01, frameon=None)
            pl.close()
        return endRange

longSim = AnalysedFiles("LongMain",equib=False)

equibSim1 = AnalysedFiles("Equib1",equib=True,equibVals=longSim.equibVals,equibSigs=longSim.equibSigs)
equibSim2 = AnalysedFiles("Equib2",equib=True,equibVals=longSim.equibVals,equibSigs=longSim.equibSigs)
equibSim3 = AnalysedFiles("Equib3",equib=True,equibVals=longSim.equibVals,equibSigs=longSim.equibSigs)


print longSim.Lks

exit()
averageTimes = []
error = []



for i in xrange(5):
    simTimes = np.array((longSim.equibTimes[indices[i]], equibSim1.equibTimes[i], equibSim2.equibTimes[i], equibSim3.equibTimes[i]))
    averageTimes.append(np.mean(simTimes))
    stdev = np.std(simTimes)
    error.append(stdev/(4.0**0.5))

pl.errorbar(equibSim1.Lks,averageTimes,yerr=error,color ="r", marker = '.',ls="None", label="Averaged Values")

pl.legend(loc="best")

pl.title("Equilibrium Time as a Function of Initial Linking Number")
pl.xlabel("Initial Linking Number")
pl.ylabel("Equilibrium Time")
pl.savefig("InitialLkEqTime", dpi=300, 
           facecolor='w', edgecolor='w', orientation='landscape', 
           papertype=None, format=None, transparent=False, 
           bbox_inches='tight', pad_inches=0.01, frameon=None)
pl.close()


