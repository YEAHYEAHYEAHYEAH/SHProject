#======= Script to find equilibrium values w/ error from txt files =====
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
from supercoiling_cython import calculate_Writhe
import glob
import scipy.odr as odr


def str8Line(A,x): return A[0]*x + A[1]

def autoCorrelation(xData,maxRange):                                    
    """
    A function to calculate the auto correlation time of a data set.
    :xData: the data for which the time is found.
    :maxRange: the range over which an exponential decay is fitted(=100).
    """
    maxTp = maxRange                                                    # the maximum distance between values in a list, eg List[0] and List[100]
    xMean = np.average(xData)                                           # finds the mean of all values
    C = np.zeros(maxTp)                                                 # numpy arrays to hold Correlation
    tPr = np.zeros(maxTp)                                               # different Correlation times, t', No of elements between values
    for i in xrange(maxTp):
        xLen = xData.shape[0]
        runTot = 0
        for j in xrange(xLen-i):
            runTot+=xData[j]*xData[j+i]                                 # checks how correlated each bit is
        C[i]=(float(runTot)/float(xLen-i))
        tPr[i]=(i)
    C = C - xMean**2.
    
    lnC = np.log(C)
    myodr = odr.ODR(odr.Data(tPr,lnC), odr.Model(str8Line), beta0=[-0.00558, -1.0])
    result = myodr.run()
    fitLine = result.beta[0]*tPr + result.beta[1]
    
    #pl.plot(tPr,lnC)
    #pl.plot(tPr,fitLine)
    #pl.show()
    
    correlationTime = -1.0/result.beta[0]
    return correlationTime

def computeError(vals,corrT):
    corrT = int(corrT)+1
    uncorrVals = []
    resampVals = []
    for i in xrange(len(vals)):
        if (i%corrT==0):
            uncorrVals.append(vals[i])
    error = np.std(uncorrVals)/(float(len(uncorrVals))**0.5)
    return error


def readDataFiles():
    
    resFile = open("InitialLkWritheData.txt",'w')
    resFile.write("500 beads\nTotal Lk,Wr,StdDev,Error\n\n")
    
    for fname in glob.glob('LongMain/*.txt'):
        title = fname[fname.find("/")+1:-4]
        dataFile = open(fname,'r')
        plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=5003)
        
        frameNo = plotArray[:,0]
        twists = plotArray[:,1]
        writhes = plotArray[:,2]        
        total = plotArray[:,1]+plotArray[:,2]
        
        corrTime = autoCorrelation(writhes,100)
        
        totalLk  = int(np.average(total))+1
        avTw     = np.average(writhes)
        stddevTw = np.std(writhes)
        errTw    = computeError(writhes,corrTime)
        resFile.write("{},{},{},{}\n".format(totalLk,avTw,stddevTw,errTw))
    resFile.close()

readDataFiles()
