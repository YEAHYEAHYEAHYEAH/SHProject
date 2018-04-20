#============== Script to plot Analysed Data ===========================
#_______________________________________________________________________


#import modules
import math as m
import numpy as np
import matplotlib.pyplot as pl
import matplotlib
import glob
from scipy.odr import *

matplotlib.rcParams.update({'font.size': 26})


dataFile = open("InitialLkWritheData.txt",'r')
plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)
    
total1 = plotArray[:,0]
twists1 = plotArray[:,1]
stddev1 = plotArray[:,2]        
error1 = plotArray[:,3]

def str8Line(A,x): return A[0]*x + A[1]
def exponent(c,averages): return c[0]*np.exp(-(c[1]*averages))

xFit1 = np.linspace(1,24, num=1000) #Prepare to plot fit
dataIn = mydata = RealData(total1, twists1, sy=error1)
model = Model(str8Line)
varODR = ODR(dataIn, model, beta0=[2000.0,2.0])
modOut = varODR.run()
yFit1 = str8Line(modOut.beta, xFit1)
chisq = modOut.res_var

print modOut.beta
print modOut.sd_beta

dataFile = open("InitialLkTwistData.txt",'r')
plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)
    
total2 = plotArray[:,0]
twists2 = plotArray[:,1]
stddev2 = plotArray[:,2]        
error2 = plotArray[:,3]

def str8Line(A,x): return A[0]*x + A[1]
def exponent(c,averages): return c[0]*np.exp(-(c[1]*averages))

xFit2 = np.linspace(1,24, num=1000) #Prepare to plot fit
dataIn = mydata = RealData(total2, twists2, sy=error2)
model = Model(str8Line)
varODR = ODR(dataIn, model, beta0=[2000.0,2.0])
modOut = varODR.run()
yFit2 = str8Line(modOut.beta, xFit2)
chisq = modOut.res_var

print modOut.beta
print modOut.sd_beta


pl.errorbar(total1, twists1/twists2,yerr=(((error1/twists1)**2.0 + (error2/twists2)**2.0)**0.5),linestyle=" ",marker=".",color="g")
pl.plot(xFit1, yFit1/yFit2,"k--",linewidth=3)
pl.xlim(1,25)
#pl.ylim(0,1)
pl.ylabel("Ratio of Writhe to Twist")
pl.xlabel("Initial Linking Number")



#pl.plot(xFit, yFit, "k--",label = ("Linear ($\chi^{2}$="+"{:.2f})".format(chisq)),linewidth=1)

#pl.legend(loc="best")


#pl.errorbar(total,twists,yerr=error, color="b", marker=".", linestyle=" ")
#pl.xlabel("Initial Linking Number")
#pl.ylabel("Equilibrium Writhe")
#pl.title("Equilibrium Twist Value with Varying Initial Lk")
pl.savefig("analysedData", dpi=300, facecolor='w', edgecolor='w',
        orientation='landscape', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.01,
        frameon=None)
pl.close()



