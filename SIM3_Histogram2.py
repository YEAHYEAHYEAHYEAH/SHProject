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
import matplotlib

f, axarr = pl.subplots(2)


def str8Line(A,x): return A[0]*x + A[1]

def autoCorrelation(xData,maxRange,plot,col):                                    
    """
    A function to calculate the auto correlation time of a data set.
    :xData: the data for which the time is found.
    :maxRange: the range over which an exponential decay is fitted(=100).
    """
    maxTp = maxRange                                                    # the maximum distance between values in a list, eg List[0] and List[100]
    xMean = np.average(xData)                                           # finds the mean of all values
    C = np.zeros(maxTp)                                                 # numpy arrays to hold Correlation
    tPr = np.zeros(maxTp)                                               # different Correlation times, t', No of elements between values
    for i in xrange(maxTp):
        xLen = xData.shape[0]
        runTot = 0
        for j in xrange(xLen-i):
            runTot+=xData[j]*xData[j+i]                                 # checks how correlated each bit is
        C[i]=(float(runTot)/float(xLen-i))
        tPr[i]=(i)
    C = C - xMean**2.
    
    lnC = np.log(np.absolute(C))
    
    myodr = odr.ODR(odr.Data(tPr,lnC), odr.Model(str8Line), beta0=[-0.00558, -1.0])
    print myodr.beta0
    
    result = myodr.run()
    print result.beta
    print result.sd_beta
    
    print result.sd_beta/(result.beta**2.0)
    
    chisq = result.res_var
    print chisq
    
    if chisq<0.001: chisq = 1.09
    
    fitLine = result.beta[0]*tPr + result.beta[1]
    correlationTime = -1.0/result.beta[0]
        
    axarr[plot].plot(tPr,lnC,color=col,label="t*={:.1f}".format(correlationTime))
    axarr[plot].plot(tPr,fitLine,"k--",label=("$\chi^{2}$="+"{:.2f}".format(chisq)),linewidth=3)
    axarr[plot].legend(loc="best")

    return correlationTime

axarr[1].set_xlabel("Time (Frames)")
axarr[1].set_ylabel("                          log(Correlation Value)")

matplotlib.rcParams.update({'font.size': 22})

lks = [5,10,15,20]

writhes = np.zeros(1)
frames = np.zeros(1)
twists = np.zeros(1)

for i in lks:
    dataFile = open("twist{}Data.txt".format(i),'r')

    plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=20003)

    frame = np.append(frames,(plotArray[:,0]))
    twists= np.append(twists,(plotArray[:,1]))
    writhes=np.append(writhes,(plotArray[:,2]))

twCorr = autoCorrelation(twists,200,0,'r')
wrCorr = autoCorrelation(writhes,500,1,'b')

pl.savefig("Correlation Times", dpi=300, 
           facecolor='w', edgecolor='w', orientation='landscape', 
           papertype=None, format=None, transparent=False, 
           bbox_inches='tight', pad_inches=0.01, frameon=None)
pl.close()


exit()

n, twBins, patches = pl.hist(twists, 40, facecolor='r', alpha=0.5, normed=True,label="Twists",edgecolor='none')
n, wrBins, patches = pl.hist(writhes, 40, facecolor='b', alpha=0.5, normed=True,label="Writhes",edgecolor='none')

twMu, twSigma = norm.fit(twists)
wrMu, wrSigma = norm.fit(writhes)

wrFit = mlab.normpdf(wrBins, wrMu, wrSigma)
wrLine = pl.plot(wrBins, wrFit, 'k--', linewidth=2)

twFit = mlab.normpdf(twBins, twMu, twSigma)
twLine = pl.plot(twBins, twFit, 'k--', linewidth=2)



#pl.title("Histogram of Writhes and Twists in Long Simulation")
pl.xlabel("Twist and Writhes Values")
pl.ylim(0,1.0)
pl.xlim(-4,4)
pl.text(1.2, .17, r'$\mu={:.2f}${}$\sigma={:.2f}$'.format(wrMu,',\n',wrSigma))
pl.text(0.35, 0.55, r'$\mu={:.2f},\ \sigma={:.2f}$'.format(twMu,twSigma))
pl.legend(loc="best")
pl.ylabel("Normalised Probability")
pl.savefig("Histogram2", dpi=300, 
           facecolor='w', edgecolor='w', orientation='landscape', 
           papertype=None, format=None, transparent=False, 
           bbox_inches='tight', pad_inches=0.01, frameon=None)
pl.close()
