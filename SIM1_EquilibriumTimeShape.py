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

f, axarr = pl.subplots(2, sharex=True)

class AnalysedFiles(object):
    def __init__(self,fName):
        fName+=".txt"
        plotArray = np.genfromtxt(fName,delimiter=',',skip_header=4)

        initialLink =  20
        self.frameNo = plotArray[:,0][0:4999]/100.0
        self.twists = plotArray[:,1][0:4999]
        self.writhes = plotArray[:,2][0:4999]
        self.total = plotArray[:,1]+plotArray[:,2]


longSim = AnalysedFiles("twist15DataL")
equibSim1 = AnalysedFiles("twist15Data")
equibSim2 = AnalysedFiles("twist15Data1")
equibSim3 = AnalysedFiles("twist15Data2")

averageTwists = []
error = []


for i in xrange(4999):
    avTw = np.array((longSim.twists[i], equibSim1.twists[i], equibSim2.twists[i], equibSim3.twists[i]))
    averageTwists.append(np.mean(avTw))
    stdev = np.std(avTw)
    error.append(stdev/(4.0**0.5))


def str8Line(A,x): return A[0]*x + A[1]
def exponent(c,averages): return c[0]*np.exp(-(c[1]*averages))

def exponent2(c,data): return c[0]*np.exp(-(c[1]*(data)**c[3]))+c[2]

def exponent1(c,data):  return (c[0]*np.exp(-(c[1]*data))+c[2])+\
                               (c[3]*np.exp(-c[4]*data))

def powerlaw(c,data):   return c[1] + np.power(data,c[0])
def invlogistic(c,data): return c[0]/(c[0]+np.exp(-c[1]*data))
def log(c,data): return c[0]*(np.log(c[1]*data)/np.log(c[1]))

xFit = np.linspace(0,50, num=150000) #Prepare to plot fit
dataIn = mydata = RealData(equibSim1.frameNo,averageTwists, sy=error)
model = Model(exponent1)
varODR = ODR(dataIn, model, beta0=[1.0,1.0,3.0,14.0,1.0])
modOut = varODR.run()
yFit = exponent1(modOut.beta, xFit)

print varODR.beta0
print modOut.beta
chisq = modOut.res_var


averageTwists=np.array(averageTwists)

print chisq


axarr[0].errorbar(equibSim1.frameNo,averageTwists,yerr=error,color ="r", marker = '.',ls="None")

axarr[0].plot(xFit, yFit, "k--",label = ("Exponential ($\chi^{2}$="+"{:.2f})".format(chisq)),linewidth=3)
axarr[0].legend(loc="best")

axarr[0].set_ylabel("Twist Value")

longSim = AnalysedFiles("twist20DataL")
equibSim1 = AnalysedFiles("twist20Data")
equibSim2 = AnalysedFiles("twist20Data1")
equibSim3 = AnalysedFiles("twist20Data2")

averageTwists = []
error = []


for i in xrange(4999):
    avTw = np.array((longSim.twists[i], equibSim1.twists[i], equibSim2.twists[i], equibSim3.twists[i]))
    averageTwists.append(np.mean(avTw))
    stdev = np.std(avTw)
    error.append(stdev/(4.0**0.5))


xFit = np.linspace(0,50, num=150000) #Prepare to plot fit
dataIn = mydata = RealData(equibSim1.frameNo,averageTwists, sy=error)
model = Model(exponent1)
varODR = ODR(dataIn, model, beta0=[1.0,1.0,3.0,17.0,1.0,10.0])
modOut = varODR.run()
yFit = exponent1(modOut.beta, xFit)

print varODR.beta0
print modOut.beta
chisq = modOut.res_var


averageTwists=np.array(averageTwists)

print chisq


axarr[1].errorbar(equibSim1.frameNo,averageTwists,yerr=error,color ="r", marker = '.',ls="None", label="Averaged Values")

axarr[1].plot(xFit, yFit, "k--",label = ("Exponential ($\chi^{2}$="+"{:.2f})".format(chisq)),linewidth=3)



axarr[1].legend(loc="best")

#pl.title("Equilibrium Time as a Function of Initial Linking Number")
pl.xlabel("Frame Number (100s)")
axarr[1].set_ylabel("Twist Value")
pl.savefig("InitialLkEqTime", dpi=300, 
           facecolor='w', edgecolor='w', orientation='landscape', 
           papertype=None, format=None, transparent=False, 
           bbox_inches='tight', pad_inches=0.01, frameon=None)
pl.close()
