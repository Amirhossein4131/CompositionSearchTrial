import numpy as np
import matplotlib.pyplot as plt

elements = ["Fe", "Ni", "Cr", "Co", "Mn"]

print ("Elements are: ", elements)

element1 = input ("Choose 1st Element: ")
print("1st Element is: " + element1)

element2 = input ("Choose 1st Element: ")
print("2nd Element is: " + element2)

element3 = input ("Choose 1st Element: ")
print("3rd Element is: " + element3)

range_element1=np.arange(10,100,10) 
range_element2=np.arange(10,100,10)
range_element3=np.arange(10,100,10)

compositions=[]

for i in range (9):
	for j in range (9):
		for k in range (9):
			if (range_element1[i] + range_element2[j] + range_element3[k]) == 100:
				compositions.append(element1+str(range_element1[i])+element2+str(range_element2[j])         +element3+str(range_element3[k]))

compositions.append ("equimolar")	 #Equatomic Composition			

print (compositions)
STD = np.loadtxt(fname="STD_%s%s%s.txt" %(element1, element2, element3))
print (STD)
Mean = np.loadtxt(fname="Mean_%s%s%s.txt"%(element1, element2, element3))
print (Mean)

def STDvsComp(std, comp):
	plt.xlabel("%s%s%s Compositions" %(element1,element2,element3))
	plt.ylabel("Standard Deviation")
	plt.title ("Standard Deviation vs Different Compositions")
	plt.xticks(rotation = 90)
	plt.scatter(comp, std)
	plt.show()
	
def MeanvsComp(mean, comp):
	plt.xlabel("%s%s%s Compositions"%(element1,element2,element3))
	plt.ylabel("Mean of Strain")
	plt.title ("Mean of Strains vs Different Compositions")
	plt.xticks(rotation = 90)
	plt.scatter(comp, mean)
	plt.show()

STDvsComp(std=STD, comp=compositions)
MeanvsComp(Mean,compositions)

