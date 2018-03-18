#============== Script to plot Analysed Data ===========================
#_______________________________________________________________________


#import modules
import math as m
import numpy as np
import matplotlib.pyplot as pl
import glob


def plotLkVsTw():
    dataFile = open("AnalysedData.txt",'r')
    plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)
        
    total = plotArray[:,0]
    twists = plotArray[:,1]
    stddev = plotArray[:,2]        
    error = plotArray[:,3]

    pl.errorbar(total,twists,yerr=error, color="k", marker=".", linestyle=" ")
    pl.xlabel("Total Linking Number")
    pl.ylabel("Twists")
    pl.title("Analysed Data")
    pl.savefig("analysedData", dpi=300, facecolor='w', edgecolor='w',
            orientation='landscape', papertype=None, format=None,
            transparent=False, bbox_inches='tight', pad_inches=0.01,
            frameon=None)
    pl.close()

plotLkVsTw()
