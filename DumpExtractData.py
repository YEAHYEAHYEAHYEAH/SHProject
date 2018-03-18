#=============== Script to extract features from LAMMPS dump files =====
#_______________________________________________________________________


#import modules
import glob
import LAMMPSInterface as lammps


def writeFiles():
    
    for fname in glob.glob('dna2/*'):
        title = fname[fname.find("/")+1:-4]
        lammpsDump = open(fname, 'r')
        
        print title
        
        writeFile = title+"Data.txt"
        dataFile = open(writeFile,'w')
        dataFile.write("500 beads: initial {}\nFrame,Tw,Wr\n\n".format(title))
        frame = lammps.Frame(lammpsDump)
        for i in range(40000):
            frame = Frame(lammpsDump)
            dataFile.write("{},{},{}\n".format(i,frame.calcTwists(),frame.calcWrithes()))
        
        dataFile.close()

writeFiles()
