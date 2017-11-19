from sys import argv
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["savefig.directory"] = "/home/mathisre/Dropbox/Uni/Images/Proj4"

file = open(argv[1], 'r')
n = (file.readline().split())
nspins = int(n[1])
n = int(n[0])
with file as filename:
    lines = [line.split() for line in filename]

probE = np.zeros(n)
sum = 0
for k in range(len(lines)): #read data into array
    probE[k] = float(lines[k][0])
file.close()
print np.sum(probE)

max_E = (n-1) / (2*(nspins*nspins))
MCc = np.linspace(-max_E,max_E, n)


plt.plot(MCc, probE, 'r--', label='Measured mean energy')
plt.xlabel('Energy in units of J')
plt.ylabel('Probability of state having energy E')
plt.title('Probability of spin state having energy E for T = 2.4 for a 20x20 lattice')
plt.grid(True)
#plt.legend()
plt.axis([-2,-1.95, 0, 1])

plt.show()


