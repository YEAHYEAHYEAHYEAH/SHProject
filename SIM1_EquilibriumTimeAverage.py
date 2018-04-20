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
from scipy.odr import *
import matplotlib

matplotlib.rcParams.update({'font.size': 22})

class AnalysedFiles(object):
    def __init__(self,dirPath,equib=False,equibVals = None,equibSigs=None):
        self.equibTimes = []
        self.Lks        = []
        self.equib      = equib
        print "\n{}\n".format(dirPath)
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
            
            plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)

            initialLink =  int(title[5:-4])
            frameNo = plotArray[:,0]
            twists = plotArray[:,1]
            writhes = plotArray[:,2]        
            total = plotArray[:,1]+plotArray[:,2]
            
            equibTime = self.equilibriumRange(twists,initialLink)
            self.equibTimes.append(equibTime)

            self.Lks.append(initialLink)
            
            dataFile.close()
            
    def equilibriumRange(self,xData,link):
        if self.equib==False: 
            equibVal = [link,np.mean(xData[20000:])]
            self.equibVals.append(equibVal)
            
            equibSig = [link,np.std(xData[20000:])]
            self.equibSigs.append(equibSig)

        else:
            for i in self.equibVals:
                if i[0] == link: equibVal = i
            for i in self.equibSigs:
                if i[0] == link: equibSig = i
            
        endRange = 0
        
        for i in xrange(len(xData)):
            
            if abs(equibVal[1] - xData[i])<equibSig[1]:
                endRange = i
                print "Lk: {}, EqVal: {}, endFrame: {}".format(link,equibVal[1],endRange)
                break
        if False:
            
            pl.plot(xData[:endRange])
            pl.title(self.dirPath + "\Initial Twist: "+str(link))
            pl.xlabel("Frame Number (50 Natural Time Units)")
            pl.ylabel("Twist Value")
            pl.savefig(self.dirPath+"/Equib{}InLk".format(link), dpi=300, 
                       facecolor='w', edgecolor='w', orientation='landscape', 
                       papertype=None, format=None, transparent=False, 
                       bbox_inches='tight', pad_inches=0.01, frameon=None)
            pl.close()
        return endRange

longSim = AnalysedFiles("LongMain",equib=False)
equibSim1 = AnalysedFiles("EquilibriumTime1",equib=True,equibVals=longSim.equibVals,equibSigs=longSim.equibSigs)
equibSim2 = AnalysedFiles("EquilibriumTime2",equib=True,equibVals=longSim.equibVals,equibSigs=longSim.equibSigs)
equibSim3 = AnalysedFiles("EquilibriumTime3",equib=True,equibVals=longSim.equibVals,equibSigs=longSim.equibSigs)

averageTimes = []
error = []
indices = [0,5,11,19]


for i in xrange(4):
    simTimes = np.array((longSim.equibTimes[indices[i]], equibSim1.equibTimes[i], equibSim2.equibTimes[i], equibSim3.equibTimes[i]))
    averageTimes.append(np.mean(simTimes))
    stdev = np.std(simTimes)
    error.append(stdev/(4.0**0.5))


goodPointsY = []
goodPointsX = []
for i in xrange(4):
    goodPointsX.append(longSim.Lks[indices[i]])
    goodPointsY.append(longSim.equibTimes[indices[i]])

def str8Line(A,x): return A[0]*x + A[1]
def exponent(c,averages): return c[0]*np.exp(-(c[1]*averages))

xFit = np.linspace(0,25, num=1000) #Prepare to plot fit
dataIn = mydata = RealData(equibSim1.Lks, averageTimes, sy=error)
model = Model(exponent)
varODR = ODR(dataIn, model, beta0=[2000.0,2.0])
modOut = varODR.run()
yFit = exponent(modOut.beta, xFit)
chisq = modOut.res_var

#yFit += exponent([-1700.0,1.0],xFit)

print chisq

pl.plot(xFit[120:,], yFit[120:,], "k--",label = ("Exp ($\chi^{2}$="+"{:.2f})".format(chisq)),linewidth=2)


pl.plot(longSim.Lks,longSim.equibTimes, color="0.7", marker = ".",ls="None")
pl.plot(goodPointsX,goodPointsY, color="k", marker = ".",ls="None")
pl.plot(equibSim1.Lks,equibSim1.equibTimes, color="k", marker = ".",ls="None")
pl.plot(equibSim2.Lks,equibSim2.equibTimes, color="k", marker = ".",ls="None")
pl.plot(equibSim3.Lks,equibSim3.equibTimes, color="k", marker = ".",ls="None")


pl.errorbar(equibSim1.Lks,averageTimes,yerr=error,color ="g", marker = '.',ls="None", label="Average Time")

pl.legend(loc="best")

#pl.title("Equilibrium Time as a Function of Initial Linking Number")
pl.xlabel("Initial Linking Number")
pl.ylabel("Equilibrium Time")
pl.savefig("InitialLkEqTime", dpi=300, 
           facecolor='w', edgecolor='w', orientation='landscape', 
           papertype=None, format=None, transparent=False, 
           bbox_inches='tight', pad_inches=0.01, frameon=None)
pl.close()
