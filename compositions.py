import numpy as np


elements=["Ni","Fe","Cr"]
range_Ni=np.arange(10,100,5) 
range_Fe=np.arange(10,100,5)
range_Cr=np.arange(10,100,5)

compositions=[]

for i in range (9):
	for j in range (9):
		for k in range (9):
			if (range_Ni[i] + range_Fe [j] + range_Cr[k]) == 100:
				compositions.append("Ni"+str(range_Ni[i])+"Fe"+str(range_Fe[j])+"Cr"+str(range_Cr[k]))


#compositions.append("Ni"+str(range_Ni[0])+"Fe"+str(range_Fe[0])+"Cr"+str(range_Cr[0])) 
#print (range_Ni, range_Fe, range_Cr)          
print (compositions)
