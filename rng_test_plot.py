from sys import argv
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["savefig.directory"] = "/home/mathisre/Dropbox/Uni/Images/Proj4"

file = open(argv[1], 'r')
n = int(file.readline())+1

with file as filename:
    lines = [line.split() for line in filename]

x = np.zeros(n)
xx = np.zeros(n)
ind = np.zeros(n)
autocorr = np.zeros(n)

mean = 0
meanxx = 0
for k in range(len(lines)): #read data into array
    ind[k] = float(lines[k][0])
    autocorr[k] = float(lines[k][1])

    
    #x[k] = float(lines[k][0])
    #xx[k] = float(lines[k][1])
    #mean = mean + x[k]
    #meanxx = meanxx + x[k]*x[k]
file.close()


#mean = mean / n
#meanxx = meanxx / n
#print mean
#print np.sqrt(meanxx - mean*mean)
#plt.hist(x, facecolor = 'green')#, alpha = 0.75)

plt.plot(ind, autocorr)
plt.xlabel('i\'th number in RNG sequence')
plt.ylabel('Auto correlation function')
plt.title('Auto correlation function of Mersenne twister random numbers')
plt.axis([0, n, -0.02, 0.02])

#plt.xlabel('Random number values')
#plt.ylabel('Random number count')
#plt.title('Histogram of random numbers in [0,1]')

plt.grid(True)
plt.show()






