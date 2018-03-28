#=============== Script to extract features from LAMMPS dump files =====
#_______________________________________________________________________


#import modules
import glob
import LAMMPSInterface as lammps

numbers = [5,10,15,20]

for x in numbers:
    startText = " LAMMPS data file\n\n 500 atoms\n 500 ellipsoids\n 500 bonds\n 1000 angles\n\n 1 atom types\n 1 bond types\n 2 angle types\n\n -100 100 xlo xhi\n -100 100 ylo yhi\n -100 100 zlo zhi\n\n\n Masses\n\n1 1\n\n Atoms\n\n"

    inputFile = open("500beads{}twistEq.in".format(x), 'w')

    defaults = open("500beads20twist.in",'r')

    for i in xrange(1025): next(defaults)


    dumpFile = open("twist{}.DNA".format(x))

    #frame = lammps.Frame(dumpFile)
    frame = lammps.Frame(dumpFile,noFr=4999,skipFrames = True)
    inputFile.write(startText)

    for i in xrange(frame.noAtms):
        atom = frame.atoms[i]
        
        x = atom.getPos()[0]
        y = atom.getPos()[1]
        z = atom.getPos()[2]
        
        q1 = atom.getQuat()[0]
        q2 = atom.getQuat()[1]
        q3 = atom.getQuat()[2]
        q4 = atom.getQuat()[3]
        
        inputFile.write(" {} 1 {} {} {} 1 1 1.90986\n".format(atom.getID(),x,y,z))

    inputFile.write("\n Ellipsoids\n\n")

    for i in xrange(frame.noAtms):
        atom = frame.atoms[i]
        
        x = atom.getPos()[0]
        y = atom.getPos()[1]
        z = atom.getPos()[2]
        
        q1 = atom.getQuat()[0]
        q2 = atom.getQuat()[1]
        q3 = atom.getQuat()[2]
        q4 = atom.getQuat()[3]
        
        inputFile.write(" {} 1 1 1 {} {} {} {}\n".format(atom.getID(),q1,q2,q3,q4))


    for line in defaults: inputFile.write(line)
