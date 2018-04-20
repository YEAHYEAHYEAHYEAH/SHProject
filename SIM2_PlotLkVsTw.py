#============== Script to plot Analysed Data ===========================
#_______________________________________________________________________


#import modules
import math as m
import numpy as np
import matplotlib.pyplot as pl
import glob



dataFile = open("EquilibriumTw2.txt",'r')
plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)
    
total = plotArray[:,0]
rigidity = plotArray[:,1]
writhe = plotArray[:,2]        
stdev = plotArray[:,3]
error = plotArray[:,4]

initialLks = [10,20,30]



pl.errorbar(rigidity[0:5],  writhe[0:5]/total[0:5],  yerr=error[0:5]/total[0:5],   color="r", marker=".", linestyle=" ",label="Initial Lk: 10")
pl.errorbar(rigidity[5:10], writhe[5:10]/total[5:10], yerr=error[5:10]/total[5:10],  color="b", marker=".", linestyle=" ",label="Initial Lk: 20")
pl.errorbar(rigidity[10:15],writhe[10:15]/total[10:15],yerr=error[10:15]/total[10:15], color="k", marker=".", linestyle=" ",label="Initial Lk: 30")

pl.legend(loc='best')

pl.xlim(5,55)

pl.xlabel("Rigidity Coefficient")
pl.ylabel("Equilibrium Twist Fraction")
pl.title("Twists with Varying Initial Lk and Rigidity Coefficient")
pl.savefig("analysedData", dpi=300, facecolor='w', edgecolor='w',
        orientation='landscape', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.01,
        frameon=None)
pl.close()

