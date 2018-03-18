#=============== Script to analyse LAMMPS dump files ===================
#_______________________________________________________________________


#import modules
import math as m
import numpy as np
import matplotlib.pyplot as pl
import re
import operator
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from supercoiling_cython import calculate_Writhe
import glob
import scipy.odr as odr


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def progress(count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush() #


def makeUnit(vec):
    norm = np.linalg.norm(vec)
    #norm = np.dot(vec,vec)**0.5
    return vec/norm
    
def rotationMatrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = makeUnit(np.asarray(axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

class Frame(object):
    """
    A class to hold frame data, a sorted list of atom objects which have
    been scaled and unwrapped.
    """
    def __init__(self,inputS,noFr=0,skipFrames = False):
        if skipFrames==True: self.skipFr(inputS,noFr)
        else: self.analyseFrame(inputS)
    
    def skipFr(self,inputStream,noF):
        noLines = noF*509
        for i in xrange(noLines):
            next(inputStream)
            if i==noLines-1: print i/509.

    def analyseFrame(self,inputStream):
        frameStart = r'c_quat[4] '
        frameEnd   = r'ITEM: TIMESTEP'
        inFrame = False
        frame = []
        atoms = []
        lnNum = 0
        bounds = np.empty([3],dtype=object)
        
        for _ in xrange(3):
            next(inputStream)

        for line in inputStream:
            lnNum += 1
            if line=="ITEM: BOX BOUNDS pp pp pp\n": lnNum = 2
            #if lnNum==1: self.noAtoms = int(line.strip(" \n"))
            if lnNum==3 or lnNum==4 or lnNum==5: 
                bounds[lnNum-3] = [int(s) for s in line.strip(" \n").split(" ")]
            
            if frameEnd in line:
                inFrame = False
                break
    
            if inFrame == True:
                lineFloat = [float(s) for s in line.strip(" \n").split(" ")]
                frame.append(lineFloat)
                newAtom = Atom(lineFloat,bounds)
                atoms.append(newAtom)

            if frameStart in line:
                inFrame = True

        atoms.sort(key=operator.attrgetter('_atmID'))

        self.atoms = atoms
        self.xyzPos = np.array(frame)[:,2:5]
        del frame[:]
        self.noAtms = len(atoms)        
        #self.calcTwists()
        #self.calcWrithes()

    def plot(self):
        """
        Plots the location of all the beads in a 3D box.
        """
        print self.totalTwists
        fig = pl.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax._axis3don = False
        ax.scatter(self.xyzPos[:,0],self.xyzPos[:,1],self.xyzPos[:,2])
        ax.text(1,0,0,'Twists: {0:.2f}'.format(self.totalTwists)+
                      '\nWrithes: {0:.2f}'.format(self.totalWrithes),color="k")
        pl.show()

    def radiusOfGyration(self):
        """ 
        Calculate the radius of gytation:
        Rg^2 = (1/N) sum ( r_k - r_mean)^2
        """

        # get mean position of all atoms
        meanRadius = np.zeros(3,dtype=np.float64)
        for i in range(self.noAtms):
            meanRadius[0] += self.atoms[i].getPos()[0]
            meanRadius[1] += self.atoms[i].getPos()[1]
            meanRadius[2] += self.atoms[i].getPos()[2]
        meanRadius = meanRadius/self.noAtms

        # get radius of gyration squared
        sqRadiusGy = 0.0
        for i in range(self.noAtms):
            sqRadiusGy += np.sum(np.square(self.atoms[i].getPos()-meanRadius))
        sqRadiusGy = sqRadiusGy/self.noAtms

        return np.sqrt(sqRadiusGy)


    def calcTwists(self):                                               # start by creating a set of vectors which point from one atom to the next
        
        angA     = np.empty(self.noAtms,dtype=float)
        binomVec = np.empty(self.noAtms,dtype=object)

        for i in xrange(self.noAtms):                                   # the -1 starts by connecting the last atom to the first, then goes as normal (total-1)
            self.atoms[i-1].toNext = makeUnit(self.atoms[i].getPos()-self.atoms[i-1].getPos())
        
        for i in xrange(self.noAtms):
            self.atoms[i].generateMVec()
        
        for i in xrange(self.noAtms):
            binomVec[i] = makeUnit(np.cross(self.atoms[i-1].toNext,self.atoms[i].toNext))
            
        for i in xrange(self.noAtms):
            angA[i] = np.arccos(np.dot(self.atoms[i-1].toNext,self.atoms[i].toNext))
        
        for i in xrange(self.noAtms):
            self.atoms[i].mShift = np.dot(rotationMatrix(binomVec[i],angA[i]),self.atoms[i-1].mVec)
        
        Tw = 0.0
        for i in xrange(self.noAtms):
            Tw += self.atoms[i].calcSinTheta()
        Tw=Tw/(2.0*m.pi)
        
        self.totalTwists = Tw
        return self.totalTwists
    
    def calcWrithes(self):
        self.totalWrithes = calculate_Writhe(self.atoms)
        return self.totalWrithes



class Atom(object):
    """
    An class to represent atom objects which takes line from dumpfile and
    assigns to different properties. 
    """
    def __init__(self,atomLine,bConds):
        self._atmID   = int(atomLine[0])
        self._atmType = int(atomLine[1])
        self._atmXYZ  = np.array(atomLine[2:5])
        self._atmIm   = np.array(atomLine[5:8],dtype=int)
        for i in range(3):
            bcDiff = bConds[i][1]-bConds[i][0]
            self._atmXYZ[i] = (self._atmXYZ[i]-0.5)*bcDiff # scale
            self._atmXYZ[i] = self._atmXYZ[i] + self._atmIm[i]*bcDiff #unwrap
        self._atmQTN  = np.array(atomLine[8:12])
        self.qtnAxisF()
        #self.qtnAxisU()
        #self.qtnAxisV()
        self.toNext = np.zeros(3)
        self.mVec      = np.zeros(3)
        self.mShift    = np.zeros(3)
        
        
    def getID(self):     return self._atmID
    def getType(self):   return self._atmType
    def getPos(self):    return self._atmXYZ
    def getQuat(self):   return self._atmQTN
    def getImFlgs(self): return self._atmIm
    def getF(self):      return self._fAxis
    #def getU(self):      return self._uAxis
    #def getV(self):      return self._vAxis

    def magSep(self,B):
        """ 
        Find magnitudinal separation of this atom and other atom B 
        """
        sepSq = 0
        for i in range(3):
            sepSq += (self.getPos()[i]-B.getPos()[i])**2.
        return sepSq**(0.5)
    
    def vecSep(self,B):
        """
        Find vector separation of this atom and other atom B
        """
        return self.getPos() - B.getPos()


    def dotProd(self,B):
        """
        Find dot product of atom positions
        """
        dotProd = 0
        for i in range(3):
            dotProd += (self.getPos()[i]*B.getPos()[i])
        return dotProd

        
    def qtnAxisU(self):
        """ 
        Return the unit vector corresponding to the z-axis of the bead 
        """
        u = np.zeros(3)
        q = self.getQuat()
        u[0] = 2.0*q[1]*q[3] + 2.0*q[0]*q[2]
        u[1] = 2.0*q[2]*q[3] - 2.0*q[0]*q[1]
        u[2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]
        self._uAxis = u
        return u
        
    def qtnAxisF(self):        
        """ 
        Return the unit vector corresponding to the z-axis of the bead 
        """
        f = np.zeros(3)
        q = self.getQuat()
        f[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3]
        f[1] = 2.0*q[1]*q[2] + 2.0*q[0]*q[3]
        f[2] = 2.0*q[1]*q[3] - 2.0*q[0]*q[2]
        self._fAxis = f
        return f

    def qtnAxisV(self):        
        """ 
        Return the unit vector corresponding to the V-axis of the bead 
        """
        v = np.zeros(3)
        q = self.getQuat()
        v[0] = 2.0*q[1]*q[2] - 2.0*q[0]*q[3]
        v[1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3]
        v[2] = 2.0*q[2]*q[3] + 2.0*q[0]*q[1]
        self._vAxis = v
        return v

    def generateMVec(self):
        self.mVec = makeUnit(np.cross(self.toNext,self.getF()))

    def calcSinTheta(self):
        return np.dot(self.toNext, np.cross(self.mShift,self.mVec))



def plotFigs():
    for fname in glob.glob('dna2/*'):
        title = fname[fname.find("/")+1:-4]
        lammpsDump = open(fname, 'r')

        frameNo = []
        twists = []
        writhes = []
        total = []
        
        for i in range(10000):
            frame = Frame(lammpsDump)
            progress(i, 10000, status='{}'.format(title))
            if (i%100==0):

                frameNo.append(i)
                
                Tw = frame.calcTwists()
                Wr = frame.calcWrithes()
                
                twists.append(Tw)
                writhes.append(Wr)
                total.append(Tw+Wr)

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
        print "Complete["
        
def writeFiles():
    for fname in glob.glob('dna2/*'):
        title = fname[fname.find("/")+1:-4]
        lammpsDump = open(fname, 'r')
        
        writeFile = title+"Data.txt"
        dataFile = open(writeFile,'w')
        dataFile.write("500 beads: initial {}\nFrame,Tw,Wr\n\n".format(title))
        frame = Frame(lammpsDump,noFr=10000,skipFrames=True)
        for i in range(40000):
            frame = Frame(lammpsDump)
            progress(i, 40000, status='{}'.format(title))
            dataFile.write("{},{},{}\n".format(i,frame.calcTwists(),frame.calcWrithes()))

def str8Line(A,x): return A[0]*x + A[1]


def autoCorrelation(xData):
    maxTp = 400

    xMean = np.average(xData)
    C = np.zeros(maxTp)
    tPr = np.zeros(maxTp)
    for i in xrange(maxTp):
        xLen = xData.shape[0]
        runTot = 0
        for j in xrange(xLen-i):
            runTot+=xData[j]*xData[j+i]
        C[i]=(float(runTot)/float(xLen-i))
        tPr[i]=(i)
    C = C - xMean**2.
    
    lnC = np.log(C)
    myodr = odr.ODR(odr.Data(tPr,lnC), odr.Model(str8Line), beta0=[-0.00558, -1.0])
    result = myodr.run()
    fitLine = result.beta[0]*tPr + result.beta[1]
    
    pl.plot(tPr,lnC)
    pl.plot(tPr,fitLine)
    pl.show()
    
    correlationTime = -1.0/result.beta[0]
    
    print correlationTime
    return correlationTime

def readDataFiles():
    for fname in glob.glob('data/*'):
        title = fname[fname.find("/")+1:-4]
        dataFile = open(fname,'r')
        plotArray = np.genfromtxt(dataFile,delimiter=',',skip_header=3)
        
        frameNo = plotArray[:,0]
        twists = plotArray[:,1]
        writhes = plotArray[:,2]        
        total = plotArray[:,1]+plotArray[:,2]
        
        corrTime = autoCorrelation(writhes)
        
        
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
        

readDataFiles()
