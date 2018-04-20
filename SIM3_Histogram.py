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

matplotlib.rcParams.update({'font.size': 22})

dataFile = open("twist20Data.txt",'r')

plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=21003)

end = 1

frameNo = plotArray[:,0][:-end]
twists = plotArray[:,1][:-end]
writhes = plotArray[:,2][:-end]
total = twists+writhes

n, twBins, patches = pl.hist(twists, 40, facecolor='r', alpha=0.5, normed=True,label="Twists")
n, wrBins, patches = pl.hist(writhes, 40, facecolor='b', alpha=0.5, normed=True,label="Writhes")

twMu, twSigma = norm.fit(twists)
wrMu, wrSigma = norm.fit(writhes)

wrFit = mlab.normpdf(wrBins, wrMu, wrSigma)
wrLine = pl.plot(wrBins, wrFit, 'k--', linewidth=2)

twFit = mlab.normpdf(twBins, twMu, twSigma)
twLine = pl.plot(twBins, twFit, 'k--', linewidth=2)

pl.legend(loc="best")

#pl.title("Histogram of Writhes and Twists in Long Simulation")
#pl.xlabel("Twist and Writhes Values")
pl.ylim(0,0.7)
pl.xlim(-4,4)
pl.text(1.2, .17, r'$\mu={:.2f}${}$\sigma={:.2f}$'.format(wrMu,',\n',wrSigma))
pl.text(0.35, 0.55, r'$\mu={:.2f},\ \sigma={:.2f}$'.format(twMu,twSigma))
#pl.ylabel("Normalised Probability")
pl.savefig("Discontinuity Histogram - legend Initial", dpi=300, 
           facecolor='w', edgecolor='w', orientation='landscape', 
           papertype=None, format=None, transparent=False, 
           bbox_inches='tight', pad_inches=0.01, frameon=None)
pl.close()
