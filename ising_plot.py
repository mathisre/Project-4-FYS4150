from sys import argv
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["savefig.directory"] = "/home/mathisre/Dropbox/Uni/Images/Proj4"

ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)


file = open(argv[1], 'r')
n = int(float(file.readline())) -1 


with file as filename:
    lines = [line.split() for line in filename]

expE = np.zeros(n+1)
meanM = np.zeros(n+1)
expE2 = np.zeros(n+1)
accepted = np.zeros(n+1)
E = np.zeros(n+1)
for k in range(len(lines)): #read data into array
    expE[k] = float(lines[k][0])
    meanM[k] = float(lines[k][1])
    accepted[k] = float(lines[k][2])
    E[k] = float(lines[k][3])

file.close()

MCc = np.linspace(0,n, n+1)



plt.plot(MCc, expE, label='Measured mean energy')
#plt.plot(MCc, -7.983928344*np.ones(n+1)/4, label= 'Theoretical energy expectation value')
plt.xlabel('Monte Carlo cycles')
plt.ylabel('Energy in units of J')
plt.title('Mean energy for T = 2.4 for initial ordered 20x20 lattice')
#plt.legend(loc = 4)
plt.grid(True)
plt.axis([0,10**6, -1.26, -1.22])

plt.show()


ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)

plt.plot(MCc, meanM, label='Measured mean magnetic moment')
#plt.plot(MCc, 3.994642931*np.ones(n+1) / 4, label= 'Theoretical mean magnetic moment')
plt.xlabel('Monte Carlo cycles')
plt.ylabel('Mean magnetic moment')
plt.title('Mean magnetic moment for T = 2.4 for initial ordered 20x20 lattice')
plt.grid(True)
#plt.legend(loc = 'best')
plt.axis([0,10**6, 0.44, 0.5])

plt.show()


ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)

plt.plot(MCc, accepted)
plt.xlabel('Monte Carlo cycles')
plt.ylabel('Accepted configurations')
plt.title('Accepted configurations for T = 2.4 for 20x20 lattice')
plt.grid(True)
plt.show()

plt.plot(MCc, E)
plt.xlabel('Monte Carlo cycles')
plt.ylabel('Accepted configurations')
plt.title('Accepted configurations for T = 2.4 in initial ordered 20x20 lattice')
plt.grid(True)
#plt.show()

