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
from scipy.odr import *
import matplotlib


equibVals = np.zeros((6,4))
equibSigs = np.zeros((6,4))
matplotlib.rcParams.update({'font.size': 20})
#threshold = [0.4,1.05,2.0]

counter = 0
#percent = (float(input("%: "))/100.0)



def equilibriumRange(xData,rig,lk):
    if equibVals[rig][lk]==0.0: 
        equibVals[rig][lk] = np.mean(xData[10000:])
        equibSigs[rig][lk]  = (np.std(xData[10000:]))
    equibVal = equibVals[rig][lk]
    equibSig = equibSigs[rig][lk]
    endRange = 0
    for i in xrange(len(xData)):
        
        if abs(equibVal - xData[i])<3.0*equibSig:
            endRange = i
            break
    #pl.plot(xData[:endRange])
    #pl.show()

    return endRange

def findEquilibriumTime(sim,num):
    
    rigidities = []
    equibTimes = []
    
    
    
    for fname in sorted(glob.glob('{}/{}/*'.format(sim,num))):
        title = fname[fname.find("/")+1:-4]
        #print title
        
        initialLink =  int(title[8:-11])
        lk = int(title[8:9])
        
        rigidity    =  float(title[-8:-4])
        rig         = int(title[-8:-7])
        
        dataFile = open(fname,'r')
        plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)


        frameNo = plotArray[:,0]
        twists = plotArray[:,1]
        writhes = plotArray[:,2]        
        total = plotArray[:,1]+plotArray[:,2]
        
        equibTime = equilibriumRange(writhes,rig,lk)
        equibTimes.append(equibTime)
        rigidities.append(rigidity)
    
    return np.array([rigidities,equibTimes])
    #pl.plot(Lks,equibTimes,   label="Initial Lk: {}".format(num) ,marker=".", linestyle=" ")


nms = ['LongMain','Equib1','Equib2','Equib3']
        
long10 = findEquilibriumTime(nms[0],10)
equib1_10 = findEquilibriumTime(nms[1],10)
equib2_10 = findEquilibriumTime(nms[2],10)
equib3_10 = findEquilibriumTime(nms[3],10)

counter+=1

long20 = findEquilibriumTime(nms[0],20)
equib1_20 = findEquilibriumTime(nms[1],20)
equib2_20 = findEquilibriumTime(nms[2],20)
equib3_20 = findEquilibriumTime(nms[3],20)

counter+=1

long30 = findEquilibriumTime(nms[0],30)
equib2_30 = findEquilibriumTime(nms[2],30)
equib1_30 = findEquilibriumTime(nms[1],30)
equib3_30 = findEquilibriumTime(nms[3],30)



average10 = (long10+equib1_10+equib2_10+equib3_10 )/ 4.0
average20 = (long20+equib1_20+equib2_20+equib3_20 )/ 4.0
average30 = (long30+equib1_30+equib2_30+equib3_30 )/ 4.0

error10 = []
error20 = []
error30 = []



for i in xrange(5):
    val10 = np.array((long10[1][i],equib1_10[1][i],equib2_10[1][i],equib3_10[1][i]))
    stdev10 = np.std(val10)
    error10.append(stdev10/(4.0**0.5))
    
    val20 = np.array((long20[1][i],equib1_20[1][i],equib2_20[1][i],equib3_20[1][i]))
    stdev20 = np.std(val20)
    error20.append(stdev20/(4.0**0.5))
    
    val30 = np.array((long30[1][i],equib1_30[1][i],equib2_30[1][i],equib3_30[1][i]))
    stdev30 = np.std(val30)
    error30.append(stdev30/(4.0**0.5))






def str8Line(A,x): return A[0]*x + A[1]
def exponent(c,averages): return c[0]*np.exp(-(c[1]*averages))

xFit10 = np.linspace(5,55, num=1000) #Prepare to plot fit
dataIn10 = RealData(average10[0],average10[1],sy=error10)
model10 = Model(str8Line)
varODR10 = ODR(dataIn10, model10, beta0=[10.0,2.0])
modOut10 = varODR10.run()
yFit10 = str8Line(modOut10.beta, xFit10)
chisq10 = modOut10.res_var

xFit20 = np.linspace(5,55, num=1000) #Prepare to plot fit
dataIn20 = RealData(average20[0],average20[1],sy=error20)
model20 = Model(str8Line)
varODR20 = ODR(dataIn20, model20, beta0=[10.0,2.0])
modOut20 = varODR20.run()
yFit20 = str8Line(modOut20.beta, xFit20)
chisq20 = modOut20.res_var

xFit30 = np.linspace(5,55, num=1000) #Prepare to plot fit
dataIn30 = RealData(average30[0],average30[1],sy=error30)
model30 = Model(str8Line)
varODR30 = ODR(dataIn30, model30, beta0=[10.0,2.0])
modOut30 = varODR30.run()
yFit30 = str8Line(modOut30.beta, xFit30)
chisq30 = modOut30.res_var

pl.plot(xFit10,yFit10, color="#cc0099",linestyle="--",linewidth=2)
pl.plot(xFit20,yFit20, color="g",linestyle="--",linewidth=2)
pl.plot(xFit30,yFit30, color="k",linestyle="--",linewidth=2)




pl.errorbar(average10[0],average10[1],yerr=error10,color ="#cc0099", marker = '.',ls="None", label=("10Lk ($\chi^{2}$="+"{:.2f})".format(chisq10)))
pl.errorbar(average20[0],average20[1],yerr=error20,color ="g", marker = '.',ls="None", label=("20Lk ($\chi^{2}$="+"{:.2f})".format(chisq20)))
pl.errorbar(average30[0],average30[1],yerr=error30,color ="k", marker = '.',ls="None", label=("30Lk ($\chi^{2}$="+"{:.2f})".format(chisq30)))

pl.xlabel("Rigidity Coefficient")
pl.ylabel("Equilibrium Time")
#pl.title("Average Equilibrium Times for Different Rigidty Coefficients")
pl.xlim(0,60)
pl.ylim(0,2250)



pl.legend(loc='best')

pl.savefig("EquilibriumTime All", dpi=300, facecolor='w', edgecolor='w',
           orientation='landscape', papertype=None, format=None,
           transparent=False, bbox_inches='tight', pad_inches=0.01,
           frameon=None)
pl.close()

print sum([chisq10,chisq20,chisq30])
