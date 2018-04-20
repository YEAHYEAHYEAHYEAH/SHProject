#========== Script to plot histogram from txt files ====================
#_______________________________________________________________________


#import modules
import math as m
import numpy as np
import matplotlib.pyplot as pl
from numpy import random as rd
import matplotlib.mlab as mlab
import re
import operator
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import glob
import scipy.odr as odr
from scipy.stats import norm
from scipy.odr import *
import matplotlib

matplotlib.rcParams.update({'font.size': 22})

lks = [5,10,15,20]


def str8Line(A,x): return A[0]*x + A[1]
def exponent(c,data): return c[0]*np.exp(-(c[1]*(data)**c[3]))+c[2]
def power(c,data): return c[0]*data**c[1]
def sin(c,data):  return c[0]*np.sin(c[1]*data + c[2])

def fitALine(xD,yD,yE):
    xFit = np.linspace(xD[0],xD[-1], num=len(xD)) #Prepare to plot fit
    
    dataIn = RealData(xD,yD,sy=yE)
    model = Model(exponent)
    varODR = ODR(dataIn, model, beta0=[ yD[0],0.01526902,yD[-1],0.33685753])
    modOut = varODR.run()
    yFit = exponent(modOut.beta, xFit)
    chisq = modOut.res_var

    print chisq
    print modOut.beta
    print modOut.sd_beta
    print "\n\n"
    return [xFit,yFit,chisq]
    
def fitASin(xD,yD,yE):
    xFit = np.linspace(xD[0],xD[-1], num=len(xD)) #Prepare to plot fit
    
    dataIn = RealData(xD,yD,sy=yE)
    model = Model(sin)
    varODR = ODR(dataIn, model, beta0=[1.0,1.0,45.0])
    modOut = varODR.run()
    yFit = sin(modOut.beta, xFit)
    chisq = modOut.res_var
    
    print chisq
    print modOut.beta
    print "\n\n"
    return [xFit,yFit,chisq]

for i in lks:
    
    frame = np.arange(10000)
    
    
    twists  = np.zeros((5, 10000))
    twistAv = np.zeros(10000)
    twistEr = np.zeros(10000)
    
    writhes  = np.zeros((5,10000))
    writheAv = np.zeros(10000)
    writheEr = np.zeros(10000)
    
    for j in xrange(1,6):
        dataFile = open("twist{}Data{}.txt".format(i,j),'r')
        plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)
        twists[j-1] += plotArray[:,1]
        #pl.plot(plotArray[:,1][0:400],color="0.9",linestyle="  ",marker=".")
        writhes[j-1] += plotArray[:,2]
    
    twistAv = np.average(twists,0)
    twistEr = np.std(twists,0)/(5.0**0.5)
    
    writheAv = np.average(writhes,0)
    writheEr = np.std(writhes,0)/(5.0**0.5)
    
    twistAv,twistEr,frameTw = twistAv[0:400],twistEr[0:400],frame[0:400]
    writheAv,writheEr,frameWr = writheAv[0:4000],writheEr[0:4000],frame[0:4000]


    pl.errorbar(frameTw,twistAv,yerr=twistEr, color="b", marker=".", linestyle=" ")
    print "TWIST \n"
    
    fitDataTw = fitALine(frameTw,twistAv,twistEr)
    pl.plot(fitDataTw[0],fitDataTw[1], "--",color="k",label = ("Exponential ($\chi^{2}$="+"{:.2f})".format(fitDataTw[2])),linewidth=3)
    pl.legend(loc="best")
    pl.ylim(-0.5,3.5)
    pl.savefig("TwAv{}".format(i), dpi=300, 
               facecolor='w', edgecolor='w', orientation='landscape', 
               papertype=None, format=None, transparent=False, 
               bbox_inches='tight', pad_inches=0.01, frameon=None)
    pl.close()    
        
    
    pl.errorbar(frameWr/100.0,writheAv,yerr=writheEr, color="r", marker=".", linestyle=" ")
    print "WRITHE \n"
    
    fitDataWr = fitALine(frameWr,writheAv,writheEr)
    pl.ylim(0,20)
    pl.plot(fitDataWr[0]/100.0,fitDataWr[1], "k--",label = ("Exponential ($\chi^{2}$="+"{:.2f})".format(fitDataWr[2])),linewidth=3)
    pl.legend(loc="best")
    pl.savefig("WrAv{}".format(i), dpi=300, 
               facecolor='w', edgecolor='w', orientation='landscape', 
               papertype=None, format=None, transparent=False, 
               bbox_inches='tight', pad_inches=0.01, frameon=None)
    pl.close()
