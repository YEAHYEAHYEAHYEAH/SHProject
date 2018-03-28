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
    def __init__(self,dirPath,equib=False,equibVals = None):
        self.equibTimes = []
        self.Lks        = []
        self.equib      = equib
        print "\n{}\n".format(dirPath)
        self.dirPath    = dirPath
        if self.equib == False:    self.equibVals = []
        else:                      self.equibVals = equibVals
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
            
    def equilibriumRange(self,xData,pos,link):
        if self.equib==False: 
            equibVal = np.mean(xData[20000:])
            self.equibVals.append([link,equibVal])

        else:
            for i in equibVals:
                if i[0] == link: equibVal = i
        endRange = 0
        
        for i in xrange(len(xData)):
            if abs(equibVal[1] - xData[i])<0.4:
                endRange = i
                print "{},{}".format(link,endRange)
                break
        if True:
            
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
equibSim1 = AnalysedFiles("EquilibriumTime1",equib=True,equibVals=longSim.equibVals)



#equibSim2 = AnalysedFiles("EquilibriumTime2",equib=True,equibVals=longSim.equibVals)
#equibSim3 = AnalysedFiles("EquilibriumTime3",equib=True,equibVals=longSim.equibVals)
#equibSim4 = AnalysedFiles("EquilibriumTime4",equib=True,equibVals=longSim.equibVals)


