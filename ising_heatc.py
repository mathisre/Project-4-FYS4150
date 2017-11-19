from sys import argv
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["savefig.directory"] = "/home/mathisre/Dropbox/Uni/Images/Proj4"

file = open(argv[4], 'r')
param = (file.readline().split())
nspins = int(param[1])
n = int(param[0]) - 1 #Sometimes need to remove the -1

with file as filename:
    lines = [line.split() for line in filename]

heatC = np.zeros(n)
temperature = np.zeros(n)
susceptibility = np.zeros(n)
meanE = np.zeros(n)
meanM = np.zeros(n)
varE = np.zeros(n)
varM = np.zeros(n)
for k in range(len(lines)): #read data into array
    temperature[k] = float(lines[k][0])
    heatC[k] = float(lines[k][1])
    susceptibility[k] = float(lines[k][2])
    meanE[k] = float(lines[k][3])
    meanM[k] = float(lines[k][4])
    varE[k] = float(lines[k][5])
    varM[k] = float(lines[k][6])
file.close()



plt.plot(temperature, heatC, label='Heat capacity')
plt.xlabel('Temperature in k / J')
plt.ylabel('Heat capacity in units of k / J^2 [J/K]')
plt.title('Heat capacity for 100x100  lattice')
plt.grid(True)
plt.axis([temperature[0], 2.45, 0.75, 2.5])
plt.show()


plt.plot(temperature, susceptibility, label='Susceptibility')
plt.xlabel('Temperature in k/J')
plt.ylabel('Susceptibility in units of 1/J [/K]')
plt.title('Susceptibility for 100x100 lattice')
plt.grid(True)
plt.axis([temperature[0], 2.45, 0, 0.025])
plt.show()

plt.plot(temperature, meanE, label='Mean energy')
plt.xlabel('Temperature in k / J')
plt.ylabel('Mean energy in units of J')
plt.title('Mean energy for 100x100 lattice')
plt.grid(True)
plt.axis([temperature[0], 2.45, -1.7, -1.0])
plt.show()

plt.plot(temperature, meanM, label='Mean magnetisation')
plt.xlabel('Temperature in k / J')
plt.ylabel('Mean absolute magnetic moment')
plt.title('Mean absolute magnetic moment for 100x100 lattice')
plt.grid(True)
plt.axis([temperature[0], 2.45, 0, 1])
plt.show()



