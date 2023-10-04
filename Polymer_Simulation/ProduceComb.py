#Comb polymer Data
#_|_|_|_|_|_


import numpy as np
import math  as m
import matplotlib.pyplot as plt
import random

file1 = open ("Comb.txt","w")

f=3  #the number of functional point
nt=5  #number of teeth of comb
N=nt+nt+2  #number of monomers
M=N #number of ions
L=120  #box size

file1.write("LAMMPS Description\n\n")

file1.write(str(int(N+M))+" atoms\n")
file1.write(str(int(N-1))+" bonds\n")
file1.write("0 angles\n")
file1.write("0 dihedrals\n")
file1.write("0 impropers\n\n")

file1.write("2 atom types\n")
file1.write("1 bond types\n")
file1.write("0 angle types\n")
file1.write("0 dihedral types\n")
file1.write("0 improper types\n\n")

file1.write(str(-L/2)+" "+str(L/2)+" xlo xhi\n")
file1.write(str(-L/2)+" "+str(L/2)+" ylo yhi\n")
file1.write(str(-L/2)+" "+str(L/2)+" zlo zhi\n\n")

file1.write("Masses\n\n")

file1.write("1 1.0\n")
file1.write("2 1.0\n")



file1.write("\nAtoms\n\n")


### Coordinates of Polymerchains
for a in range(1,int(N+1)):
    x=round(random.uniform(-L/2,L/2),3)
    y=round(random.uniform(-L/2,L/2),3)
    z=round(random.uniform(-L/2,L/2),3)
    file1.write(str((int(a)))+" "+str(1)+" 1 "+"14.97 "+str(x)+" "+str(y)+" "+str(z)+"\n")
    
### Coordinates of Countions
for m in range(1,M+1):
    x=round(random.uniform(-L/2,L/2),3)
    y=round(random.uniform(-L/2,L/2),3)
    z=round(random.uniform(-L/2,L/2),3)
    file1.write(str((int(N+m)))+" "+str(int(m+1))+" 2 "+"-14.97 "+str(x)+" "+str(y)+" "+str(z)+"\n")
 
### Bonds ###
file1.write("\nBonds\n\n")
     
    
file1.write(str(1)+" 1 "+str(int(1))+" "+str(int(2))+"\n") #bondnumer bondkind atoms (head)


#body part
for b in range(2,int(2*nt+1)):
    print(b)
    if b%2 == 1:
      file1.write(str(int(b))+" 1 "+str(int(b-1))+" "+str(int(b+1))+"\n") #bondnumer bondkind atoms  
    if b%2 == 0:
      file1.write(str(int(b))+" 1 "+str(int(b))+" "+str(int(b+1))+"\n") #bondnumer bondkind atoms
    
    
file1.write(str(int(2*nt+1))+" 1 "+str(int(N-2))+" "+str(int(N))+"\n") #bondnumer bondkind atoms (tail)

file1.close()

#print(file1.read())
