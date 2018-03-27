#============== Script to plot raw dumps ===============================
#_______________________________________________________________________


#import modules
import math as m
import numpy as np
import matplotlib.pyplot as pl
import glob


for fname in glob.glob('LongMain/*.txt'):
    title = fname[fname.find("/")+1:-10]
    print title
    dataFile = open(fname,'r')
    plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)
        
    frameNo = plotArray[:,0]
    twists = plotArray[:,1]
    writhes = plotArray[:,2]
    total = twists+writhes
    
    pl.plot(frameNo,twists,  color="r", label="Twists" ,marker="o", linestyle=" ")
    pl.plot(frameNo,writhes, color="b", label="Writhes",marker="o", linestyle=" ")
    pl.plot(frameNo,total,   color="k", label="Total"  ,marker="o", linestyle=" ")
    pl.xlabel("Frame Number")
    pl.ylabel("Total Linking Number")
    pl.title(title)
    pl.legend(loc="best")
    pl.savefig(title, dpi=300, facecolor='w', edgecolor='w',
            orientation='landscape', papertype=None, format=None,
            transparent=False, bbox_inches='tight', pad_inches=0.01,
            frameon=None)
    pl.close()
