#Complete Dendrimer Data

import numpy as np
import math  as m
import matplotlib.pyplot as plt
import random

file1 = open ("Dendrimer.txt","w")

f=3  #the number of functional point
n=1  #number of generation - 1
N=((f-1)**(n+1)-1)/(f-2)*(f)+1  #number of monomers
M=10 #number of ions
L=120  #boxxize

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
     
    

for branch in range(f):
    file1.write(str(int(1+branch*(N-1)/f))+" 1 "+str(int(1+branch*(N-1)/f))+" "+str(int(N))+"\n") #bondnumer bondkind atoms
    for ge in range(n):
        nf = int((f-1)**(ge))
        ns = int(f-1)
        for father in range(1,nf+1):
          for son in range(1,ns+1):
            indexf=int(father+((f-1)**(ge)-1)/(f-2)+branch*(N-1)/f)
            indexs=int(son+(f-1)*(father-1)+((f-1)**(ge+1)-1)/(f-2)+branch*(N-1)/f)
            print(indexf,indexs)
            file1.write(str(int(indexs))+" 1 "+str(indexf)+" "+str(indexs)+"\n") #bondnumer bondkind atoms

    
file1.close()

#print(file1.read())


