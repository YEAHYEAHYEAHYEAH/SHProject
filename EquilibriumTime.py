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
from supercoiling_cython import calculate_Writhe
import glob
import scipy.odr as odr


def equilibriumRange(xData):
    equibVal = np.mean(xData[20000:])
    
    counter = 0
    endRange = 0
    for i in xrange(len(xData)):
        if abs(equibVal - xData[i])<0.04:
            endRange = i
            break
    #pl.plot(xData[:endRange])
    #pl.show()
    return endRange

def findEquilibriumTime():
    
    Lks = []
    equibTimes = []
    
    for fname in glob.glob('data/*'):
        title = fname[fname.find("/")+1:-4]
        dataFile = open(fname,'r')
        plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)

        initialLk =  int(title[5:-4])
        frameNo = plotArray[:,0]
        twists = plotArray[:,1]
        writhes = plotArray[:,2]        
        total = plotArray[:,1]+plotArray[:,2]
        
        equibTime = equilibriumRange(twists)
        equibTimes.append(equibTime)
        Lks.append(initialLk)
    pl.plot(Lks,equibTimes,  color="k", label="Twists" ,marker=".", linestyle=" ")
    pl.xlabel("Initial Twist")
    pl.ylabel("Equilibrium Time")
    pl.title("Average Equilibrium Times for Different Initial Twists")

    pl.savefig("EquilibriumTime", dpi=300, facecolor='w', edgecolor='w',
               orientation='landscape', papertype=None, format=None,
               transparent=False, bbox_inches='tight', pad_inches=0.01,
               frameon=None)
    pl.close()
        
        
findEquilibriumTime()
