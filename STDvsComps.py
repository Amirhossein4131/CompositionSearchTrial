#-------------------------------------- Importing Needed Libraries -------------------------------------#

from lammps import lammps
lmp = lammps()

from mpi4py import MPI                            # To run parallel

import numpy as np
import pandas as pd 
import statistics
import os
import json

#---------------------------------------- Open json File that Contains Inputs --------------------------#

with open('input.json') as json_file:
	Json = json.load(json_file)
	element1 = Json['element1']
	element2 = Json['element2']
	element3 = Json['element3']
	xlo = Json['xlo']
	xhi = Json['xhi']
	ylo = Json['ylo']
	yhi = Json['yhi']
	zlo = Json['zlo']
	zhi = Json['zhi']
	steps = Json['steps']
	maxiter = Json['maxiter']
	swapiter = Json['swapiter']

# Check input data
print ("\n\n\n\n\n","element1=",element1,"element2=",element2,"element3=",element3,
        "xlo=",xlo,"xhi=",xhi,"ylo=",ylo,"yhi=",yhi,"zlo=",zlo,"zhi=",zhi
        ,"steps=",steps,"maxiter=",maxiter,"swapiter=",swapiter,"\n\n\n\n\n\n")	
#-------------------------------------- Define What Compositions are Concerned -------------------------#

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

STD = []
Mean= []
LattParms = ["a_1=", "a_2=", "a_3="]
NumAtoms  = ["Number of %s Atoms="%element1, "Number of %s Atoms="%element2, "Number of %s Atoms="%element3]	
swap = int (Number_of_Atoms / 4)
f = int ((steps*0.09))
b = 1

#-------------------------- Initial lammps commands

def initial_lmp_commands():
	cmd = [
	"units metal",
	"read_data %s.data" %compositions[i],
	"variable t_eq equal %d" %steps,
	"variable maxiter equal %d" %maxiter,
	"variable swaps equal %d" %swap,
	"variable swapiter equal %d" %swapiter,
	"pair_style meam/c",
	"pair_coeff * * library_CoNiCrFeMn.meam Co Ni Cr Fe Mn parameters.meam %s %s %s" %(element1, element2, element3),
	]
	return cmd

#------------------------- Write basic data information foe EACH composition explicitly

dir_folder="post-data_NumberofAtoms_LatticeConstants_STD_Mean"
def write_to_file ():
	Basic=open ("%s_BasicProperties.txt"%compositions[i], "a")
	Basic.write (
			"Simulation "+str(b)+" is Finished"+"\n\n"+str(compositions[i])+"\n\n"
			+str(LattParms[0])+str(a1)+"\n"+str(LattParms[1])+str(a2)+"\n"
			+str(LattParms[2])+str(a3)+"\n\n"+str(NumAtoms[0])+str(Numele1Atoms)+
			"\n"+str(NumAtoms[1])+str(Numele2Atoms)+"\n"+str(NumAtoms[2])+
			str(Numele3Atoms)+"\n"+"STD="+str(stdd)+"\n"+"Mean="+str(mean)
		    )
	Basic.close()
	try:
		os.mkdir(dir_folder)
	except:
		print('dir is known')
		os.system('mv *.txt '+dir_folder)
	return None

#--------------------------- A function for flushing data on the terminal screen

def Print_to_Terminal():
	EsentialProperties = open ("essential.data", "a")
	EsentialProperties.write(
	"\n\n\n\n\n\n"+"Simulation "+str(b)+"/%d is Finished"%len(compositions)+"\n\n"+str(compositions[i])+"\n\n"+str(LattParms[0])
	+str(a1)+"\n"+str(LattParms[1])+str(a2)+"\n"+str(LattParms[2])+str(a3)+"\n\n"+str(NumAtoms[0])+str(Numele1Atoms)
	+"\n"+str(NumAtoms[1])+str(Numele2Atoms)+"\n"+str(NumAtoms[2])+str(Numele3Atoms)+"\n"+"STD="+str(stdd)+
	"\n"+"Mean="+str(mean)+"\n\n\n\n\n\n"
	                        )
	EsentialProperties.close()
	fileObject = open ("essential.data", "r")
	filecontent = fileObject.read()
	print (filecontent)
	fileObject.flush()
	fileObject.close()
	return None

#--------------------------- Final "for" loop to do the simulation
for i in range (len(compositions)):
	lmp.commands_list(initial_lmp_commands())
	lines = open('equilibrium.in','r').readlines()
	for line in lines: lmp.command(line)
	data = np.loadtxt(fname="strainvalues.data")[f :]
	stdd = np.std (data)
	mean = statistics.mean(data)
	# save numbers for plots
	STD.append (stdd)
	Mean.append (mean)
	# Calculate basic data to be written in file for each composition
	prop = open('properties.data','r').read().split(" ")
	a1 = float(prop[1]) / (xhi-xlo)
	a2 = float(prop[2]) / (yhi-ylo)
	a3 = float(prop[3]) / (zhi-zlo)
	Numele1Atoms = prop[4]
	Numele2Atoms = prop[5]
	Numele3Atoms = prop[6]
	# Write data to file
	write_to_file()
	lmp.command ("clear")
	# Print Data on the Terminal Window
	Print_to_Terminal()
	b += 1
#----------------------------------------- Plots -------------------------------------------------#
np.savetxt("STD.txt", STD)
np.savetxt("Mean.txt", Mean)

from PostProcessor import*
STDvsComp(STD, comp=compositions)
MeanvsComp(Mean,compositions)

#------------------------------------------ End --------------------------------------------------#   

print ("DONE")
