# CompositionSearchTrial

STDvsCompos.py:

You are abale to make a high entropy alloy (HEA) with 3 elements using this python script. First, it asks you the size of simulation box.
After that, you are able to chose 3 elements out of (Mn, Cr, Fe, Co, Ni). The code will create 37 different compositions of your desired 
HEA. Chosing the number of time steps for the simulation and number of steps you prefer for minimization process, it will equilibrate your
systems and calculates absolute standard deviation and mean deviation of strains using last 10 percent of the calculted strain values.
Finally, calculated properties will be saved in STD.txt, Mean.txt and BasicProprties.txt text files. 

equilibrium.in:
STDvsCompos.py script uses this LAMMPS code to run the simulations.

PostProcessor:
This script uses STD.txt and Mean.txt to produce two scatter plots with x axsis to be the different compositions of the investigating HEA
and y axsis to be absolute standard deviation and absolute mean deviation of strains.
