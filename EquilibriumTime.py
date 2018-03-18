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
    equibVal = np.mean(xData[10000:])
    endRange = 0
    for i in xrange(len(xData)):
        if abs(equibVal - xData[i])<0.05:
            endRange = i
            break
    pl.plot(xData[:endRange])
    pl.show()


def findEquilibriumTime():
    
    resFile = open("EquilibriumTime.txt",'w')
    resFile.write("500 beads\nTotal Lk,Equilibrium Time,Error\n\n")
    
    for fname in glob.glob('data/*'):
        title = fname[fname.find("/")+1:-4]
        dataFile = open(fname,'r')
        plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)
        
        frameNo = plotArray[:,0]
        twists = plotArray[:,1]
        writhes = plotArray[:,2]        
        total = plotArray[:,1]+plotArray[:,2]
        
        equibTime = equilibriumRange(twists)
        
        #pl.plot(frameNo,twists,  color="r", label="Twists" ,marker="o", linestyle=" ")
        #pl.plot(frameNo,writhes, color="b", label="Writhes",marker="o", linestyle=" ")
        #pl.plot(frameNo,total,   color="k", label="Total"  ,marker="o", linestyle=" ")
        #pl.xlabel("Frame Number")
        #pl.ylabel("Total Linking Number")
        #pl.title(title)
        #pl.legend(loc="best")
        #pl.savefig(title, dpi=300, facecolor='w', edgecolor='w',
        #        orientation='landscape', papertype=None, format=None,
        #        transparent=False, bbox_inches='tight', pad_inches=0.01,
        #        frameon=None)
        #pl.close()
