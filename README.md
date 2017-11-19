# Project-4-FYS4150

RNG testing includes main.cpp and rng_testing.txt. main.cpp reads a file from the command line and then it calculates the autocorrelation function of 10000 random numbers created by the Mersenne twister algorithm. rng_testing.txt can be plotted using rng_test_plot.py. 

Ising model contains main.cpp which is the program that does the calculations for the ising model. It runs using mpic++. It writes text file from command line. Two example text files are added. One that includes heat capacity etc. results from 100x100 lattice (can be plotted in ising_heatc.py) and one that contains probability of state for 2x2 system (can be plotted using ising_prob_plot.py). I would have added one with mean values as function of MC cycles but it was too large for github (35 MB). That one can be plotted using ising_plot.py. It can be generated from main.cpp by commenting away some things and uncommenting an ofile command. It is vaguely described in the .cpp file.
