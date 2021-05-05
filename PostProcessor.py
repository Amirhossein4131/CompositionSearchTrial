import numpy as np
import matplotlib.pyplot as plt
import json

with open('input.json') as json_file:
	Json = json.load(json_file)
	element1 = Json['element1']
	element2 = Json['element2']
	element3 = Json['element3']

range_element1=np.arange(10,100,10) 
range_element2=np.arange(10,100,10)
range_element3=np.arange(10,100,10)

compositions=[]

for i in range (9):
	for j in range (9):
		for k in range (9):
			if (range_element1[i] + range_element2[j] + range_element3[k]) == 100:
				compositions.append(element1+str(range_element1[i])+element2+str(range_element2[j])+element3+str(range_element3[k]))

compositions.append ("equimolar")	 #Equatomic Composition			

print (compositions)
STD = np.loadtxt(fname="STD.txt")
print (STD)
Mean = np.loadtxt(fname="Mean.txt")
print (Mean)

def STDvsComp(std, comp):
	plt.xlabel("%s%s%s Compositions" %(element1,element2,element3))
	plt.ylabel("Standard Deviation")
	plt.title ("Standard Deviation vs Different Compositions")
	plt.xticks(rotation = 90)
	plt.scatter(comp, std)
	plt.show()
	return None

def MeanvsComp(mean, comp):
	plt.xlabel("%s%s%s Compositions"%(element1,element2,element3))
	plt.ylabel("Mean of Strain")
	plt.title ("Mean of Strains vs Different Compositions")
	plt.xticks(rotation = 45)
	plt.scatter(comp, mean)
	plt.show()
	return None

STDvsComp(STD, comp=compositions)
MeanvsComp(Mean,compositions)
