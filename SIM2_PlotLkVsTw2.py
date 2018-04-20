#============== Script to plot Analysed Data ===========================
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


matplotlib.rcParams.update({'font.size': 18})



dataFile = open("AnalysedData.txt",'r')
plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)
    
total = plotArray[:,0]
rigidity = plotArray[:,1]
writhe = plotArray[:,2]        
stdev = plotArray[:,3]
error = plotArray[:,4]

initialLks = [10,20,30]


f, axarr = pl.subplots(2, sharex=True)


#axarr[0].errorbar(rigidity[0:5],  writhe[0:5]/total[0:5],  yerr=error[0:5]/total[0:5],   color="#cc0099", marker=".", linestyle=" ",label="10-Lk")
axarr[0].errorbar(rigidity[5:10], writhe[5:10]/total[5:10], yerr=error[5:10]/total[5:10],  color="g", marker=".", linestyle=" ",label="20-Lk")
axarr[0].errorbar(rigidity[10:15],writhe[10:15]/total[10:15],yerr=error[10:15]/total[10:15], color="k", marker=".", linestyle=" ",label="30-Lk")

dataFile = open("EquilibriumWr4.txt",'r')
plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)
    
total = plotArray[:,0]
rigidity = plotArray[:,1]
writhe = plotArray[:,2]        
stdev = plotArray[:,3]
error = plotArray[:,4]

initialLks = [10,20,30]

#axarr[1].errorbar(rigidity[0:5],  writhe[0:5]/total[0:5],  yerr=error[0:5]/total[0:5],   color='#cc0099', marker=".", linestyle=" ")
axarr[1].errorbar(rigidity[5:10], writhe[5:10]/total[5:10], yerr=error[5:10]/total[5:10],  color='g', marker=".", linestyle=" ")
axarr[1].errorbar(rigidity[10:15],writhe[10:15]/total[10:15],yerr=error[10:15]/total[10:15], color="k", marker=".", linestyle=" ")


def str8Line(A,x): return A[0]*x + A[1]
def exponent(c,averages): return c[0]*np.exp(-(c[1]*averages))

xFit = np.linspace(5,55, num=1000) #Prepare to plot fit
dataIn = mydata = RealData(rigidity[5:10], writhe[5:10]/total[5:10], sy=error[5:10]/total[5:10])
model = Model(str8Line)
varODR = ODR(dataIn, model, beta0=[2000.0,2.0])
modOut = varODR.run()
yFit = str8Line(modOut.beta, xFit)
chisq = modOut.res_var

#yFit += exponent([-1700.0,1.0],xFit)

print chisq

axarr[1].plot(xFit, yFit, "k--",label = ("Linear ($\chi^{2}$="+"{:.2f})".format(chisq)),linewidth=2,color='#333300')





axarr[0].legend(loc='best')
axarr[1].legend(loc='best')

pl.xlim(5,55)

pl.xlabel("Rigidity Coefficient")
axarr[1].set_ylabel("Wr Fraction")
axarr[0].set_ylabel("Tw Fraction")

pl.savefig("analysedData", dpi=300, facecolor='w', edgecolor='w',
        orientation='landscape', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.01,
        frameon=None)
pl.close()
