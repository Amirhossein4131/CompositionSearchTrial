#-------------------------------------- Importing Needed Libraries -------------------------------------#

from lammps import lammps
lmp = lammps()

from mpi4py import MPI                            # To run parallel

import numpy as np
import pandas as pd 
import statistics

#-------------------------------------- Define What Compositions are Concerned -------------------------#

elements = ["Fe", "Ni", "Cr", "Co", "Mn"]
print ("Available Elements are: ", elements, "You Can Choose Only Three")

element1 = input ("Choose 1st Element: ")
print("1st Element is: " + element1)

element2 = input ("Choose 2nd Element: ")
print("2nd Element is: " + element2)

element3 = input ("Choose 3rd Element: ")
print("3rd Element is: " + element3)

range_element1=np.arange(10,100,10) 
range_element2=np.arange(10,100,10)
range_element3=np.arange(10,100,10)

compositions=[]

for i in range (9):
	for j in range (9):
		for k in range (9):
			if (range_element1[i] + range_element2[j] + range_element3[k]) == 100:
				compositions.append(element1+str(range_element1[i])+element2+str(range_element2[j])+element3+str(range_element3[k]))

compositions.append (element1+"33"+element2+"33"+element3+"33")	 #Equatomic Composition			
print (compositions)

#------------------------------------------- Define Masses ------------------------------------#

if element1 == "Fe":
	mass1 = 55.845
elif element1 == "Ni":
	mass1 = 58.6934
elif element1 == "Cr":
	mass1 = 51.9961
elif element1 == "Co":
	mass1 = 58.933195
elif element1 == "Mn":
	mass1 = 54.938044	

if element2 == "Fe":
	mass2 = 55.845
elif element2 == "Ni":
	mass2 = 58.6934
elif element2 == "Cr":
	mass2 = 51.9961
elif element2 == "Co":
	mass2 = 58.933195
elif element2 == "Mn":
	mass2 = 54.938044	

if element3 == "Fe":
	mass3 = 55.845
elif element3 == "Ni":
	mass3 = 58.6934
elif element3 == "Cr":
	mass3 = 51.9961
elif element3 == "Co":
	mass3 = 58.933195
elif element3 == "Mn":
	mass3 = 54.938044	

#-----------------------------Create Initial Super Cell with One Atom Type------------------------------#

# Simulation Box Boundries  
xlo = float (input ("Choose xlo: "))
xhi = float (input ("Choose xhi: "))
ylo = float (input ("Choose ylo: "))
yhi = float (input ("Choose yhi: "))
zlo = float (input ("Choose zlo: "))
zhi = float (input ("Choose zhi: "))

datainitcmds = ["units metal",
                "atom_style atomic",
                "boundary p p p",
                "lattice fcc 3.5",
                "region mybox block %d %d %d %d %d %d" %(xlo,xhi,ylo,yhi,zlo,zhi),
                "create_box 1 mybox",
                "create_atoms 1 region mybox",
                "mass 1 2000000",
                "run 0",
                "variable N equal count(all)",
                "write_data initial.data"
               ]

lmp.commands_list(datainitcmds)

Number_of_Atoms = lmp.get_natoms()
lmp.command ("clear")
   
# Delete "Mass" and Change "1 atom types" to "3 atom types" in initial.data datafile (To prevent errors)

fin = open("initial.data", "rt")

data = fin.read()
data = data.replace('1 atom types', '3 atom types')
data = data.replace('1 2000000', '1 %s  #%s \n2 %s  #%s \n3 %s  #%s ' %(mass1, element1, mass2, element2, mass3, element3) )

fin.close()
fin = open("initiall.data", "wt")
fin.write(data)
fin.close()

print ("mass1: ", mass1, "mass2: ", mass2, "mass3: ", mass3)  
print ("Number of Atoms are: ", Number_of_Atoms)

#------------------------------------------ Making Data Files of Differenet Compositions ---------------------------------#

def compositioncreation(elementtwo,elementthree,elementonepercentage,elementtwopercentage,elementthreepercentage):
	cmd = [
	"units metal",
	"atom_style atomic",
	"boundary p p p",
	"read_data initiall.data",
	"group kind1 type 1",
	"set group kind1 type/subset 2 %d  12345" % elementtwo,
	"group kind2 type 2",
	"group remain1 subtract all  kind2", 
	"set group remain1 type/subset 3 %d 74" % elementthree,
	"group kind3 type 3",
	"write_data %s%d%s%d%s%d.data" %(element1,elementonepercentage,element2,elementtwopercentage,element3,elementthreepercentage),
	"clear"
	]
	return cmd

for i in range (9):
	for j in range (9):
		for k in range (9):
			if (range_element1[i] + range_element2 [j] + range_element3[k]) == 100:
				elem2 = int(range_element2[j]*(Number_of_Atoms/100))
				elem3 = int(range_element3[k]*(Number_of_Atoms/100))
				lmp.commands_list(compositioncreation(elem2,elem3,range_element1[i],range_element2[j],range_element3[k]))

lmp.commands_list(compositioncreation(int (Number_of_Atoms/3),int(Number_of_Atoms/3),33,33,33))	#Equatomic Composition		


#------------------------------------ Runung the Simulation for 37 Different Compositions ------------------------------# 
BasicProperties = []

LattParms = ["a_1=", "a_2=", "a_3="]
NumAtoms  = ["Number of %s Atoms="%element1, "Number of %s Atoms="%element2, "Number of %s Atoms="%element3]	

steps = int(input ("Please Enter Number of Steps: "))
maxiter = int(input ("Please Enter Number Steps for Minimization: "))
swap = int (Number_of_Atoms / 4)
swapiter = int(input ("Please Enter Number Steps for Swap Command: "))
f = int ((steps*0.09))

for i in range (37):
	lmp.command (" units metal")
	lmp.command ("read_data %s.data" %compositions[i])
	lmp.command ("variable t_eq equal %d" %steps)
	lmp.command ("variable maxiter equal %d" %maxiter)
	lmp.command ("variable swaps equal %d" %swap)
	lmp.command ("variable swapiter equal %d" %swapiter)
	lmp.command ("pair_style meam/c")
	lmp.command ("pair_coeff * * library_CoNiCrFeMn.meam Co Ni Cr Fe Mn parameters.meam %s %s %s" %(element1, element2, element3))
	lines = open('equilibrium.in','r').readlines()
	for line in lines: lmp.command(line)
	data = np.loadtxt(fname="strainvalues.txt")[f :]
	stdd = np.std (data)
	mean = statistics.mean(data)
	STD = open ("STD.txt","a+")
	Mean = open ("Mean.txt", "a+")
	STD.write("%f \n" %stdd)
	Mean.write("%f \n" %mean)
	prop = open('properties.txt','r').read().split(" ")
	a1 = float(prop[1]) / (xhi-xlo)
	a2 = float(prop[2]) / (yhi-ylo)
	a3 = float(prop[3]) / (zhi-zlo)
	Numele1Atoms = prop[4]
	Numele2Atoms = prop[5]
	Numele3Atoms = prop[6]
	Basic = open ("Basic_properties.txt", "a+")
	Basic.write(str(compositions[i])+"\n"+str(LattParms[0])+str(a1)+"  "+str(LattParms[1])+str(a2)+"  "+str(LattParms[2])+str(a3)+"\n"+str(NumAtoms[0])+str(Numele1Atoms)+"  "+str(NumAtoms[1])+str(Numele2Atoms)+"  "+str(NumAtoms[2])+str(Numele3Atoms)+"STD="+str(stdd)+"   Mean="+str(mean)+"\n\n")
	lmp.command ("clear")
	print ("%dth Simulation is Finished"%i)
#----------------------------------------- Plots -------------------------------------------------#
STD.close()
Mean.close()
Basic.close()

from PostProcessor import*
STDvsComp(STD, comp=compositions)
MeanvsComp(Mean,compositions)

#------------------------------------------ End --------------------------------------------------#   
print ("DONE")
