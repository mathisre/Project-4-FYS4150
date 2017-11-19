#include <iostream>
#include "mpi.h"
#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <time.h>
#include <random>
#include <omp.h>
using namespace arma;
using namespace std;
ofstream ofile;

int periodic_bc(int x, int N, int plusorminus){
    return (x + N + plusorminus) % N;
}

int main(int argc, char* argv[])
{
    string filename;
    filename = argv[0];
    ofile.open(filename);
    int numprocs;
    int myrank;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Initialize uniform RNG un [0,1]
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Set up the uniform distribution for x \in [[0, 1]
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    double energy = 0; int magnetic_moment = 0;
    double local_mean_values[4];
    double total_mean_values[4];
    for (int i = 0; i < 4; i++){
        local_mean_values[i] = 0;
        total_mean_values[i] = 0;
    }

    double total_energy;
    double total_magnetic_moment;
    double heat_capacity;
    double susceptibility;
    double variance_E;
    double variance_M;
    int accepted;
    int cut_off = 300000;
    double EnergyDiff[5];


    // Set temperature
    int Nspins = 100;
    double MC_cycles = pow(10,6); double T_init = 2.25; double T_final = 2.3; double T_step = 0.005;
    double Tpoints = double((T_final - T_init) / T_step) +1;
    double norm = 1.0/((MC_cycles - cut_off)*numprocs);


    if (myrank == 0){
        cout << MC_cycles  << " Monte Carlo cycles for " << Nspins << "x" << Nspins << " spins: " << endl;
        cout << "Initial temperature = " << T_init << ", final temperature = " << T_final << ", temperature step =  " << T_step << ", " <<Tpoints << " temperature points" <<  endl;
    }

    // Create spin matrix with max (negative) energy
    // Used for counting probability of state
    mat configuration_max_energy = ones<mat>(Nspins, Nspins);
    int max_energy = 0;
    for (int i = 0; i<Nspins; i++){
        for (int k = 0; k<Nspins; k++){ // Calculate the maximum energy of the system
            max_energy += -configuration_max_energy(i,k)*(configuration_max_energy(periodic_bc(i,Nspins,-1) , k)  +
                                           configuration_max_energy(i, periodic_bc(k,Nspins,-1) ) );
        }}
    max_energy = -max_energy;
    vec N = zeros<vec>(2*max_energy +1 );

    ofile << Tpoints << "   " << Nspins <<  endl; //If doing heat capacity
    //ofile << 2*max_energy +1 << "   " << Nspins <<endl; // If doing probability of state
    //ofile << MC_cycles << endl; // If doing mean values with time
    double Timestart = MPI_Wtime();

    for (double T = T_init; T <= T_final; T+=T_step){

        // Precalculate boltzmann factor of energy difference to save computation time
        for (int dE = -8; dE <= 8; dE+=4) EnergyDiff[dE / 4 + 2 ] = exp(-dE/(T));

        mat configuration= ones<mat>(Nspins, Nspins);
        for (int i = 0; i<Nspins; i++){
            for (int k = 0; k<Nspins; k++){ //Create random spin lattice
                //configuration(i,k) = int (RandomNumberGenerator(gen)*2) * 2 -1;
            }
        }

        accepted = 0;
        total_energy = 0;
        total_magnetic_moment = 0;
        energy = 0;
        magnetic_moment = 0;
        for (int i = 0; i<Nspins; i++){
            for (int k = 0; k<Nspins; k++){ // Calculate initial spin energy
                energy += -configuration(i,k)*(configuration(periodic_bc(i,Nspins,-1) , k)  +
                                               configuration(i, periodic_bc(k,Nspins,-1) ) );
                magnetic_moment += configuration(i,k);
            }}

        // Zeroing out averages values for averages
        for (int i = 0; i < 4; i++){
            local_mean_values[i] = 0;
            total_mean_values[i] = 0;
        }

        // Monte Carlo cycles
        for (int cycle = 0; cycle < MC_cycles; cycle ++ ){
            for (int x = 0; x < Nspins; x++){
                for (int y = 0; y < Nspins; y++){
                    int sx = int (RandomNumberGenerator(gen)*Nspins);
                    int sy = int (RandomNumberGenerator(gen)*Nspins);
                    int deltaE = 2*configuration(sx, sy)*(configuration(sx, periodic_bc(sy, Nspins, 1)) + configuration(sx, periodic_bc(sy, Nspins, -1))
                                                        + configuration(periodic_bc(sx, Nspins, 1), sy) + configuration(periodic_bc(sx, Nspins, -1), sy));

                    //Metropolis algorithm for accepting new configuration
                    if (RandomNumberGenerator(gen) <=  EnergyDiff[deltaE / 4 + 2]){
                        accepted += 1;
                        configuration(sx, sy) = configuration(sx, sy) * -1;
                        magnetic_moment += 2 * configuration(sx, sy);
                        energy += deltaE;
                    }}}
            total_energy += energy;
            total_magnetic_moment += abs(magnetic_moment);

            if (cycle >= cut_off) {
                local_mean_values[0] += energy;
                local_mean_values[1] += energy*energy;
                local_mean_values[2] += fabs(magnetic_moment);
                local_mean_values[3] += magnetic_moment*magnetic_moment;
                N(energy + max_energy) += 1;
            }
            // Use this ofile to get mean energy, mean magnetic moment and accepted spins as function of monte carlo cycles
            //ofile << setprecision(8) <<total_energy /(Nspins*Nspins*(cycle+1)) <<  "    " << total_magnetic_moment /(Nspins*Nspins*(cycle+1)) << "    " << accepted <<  endl;
        } // End Monte Carlo sampling

        for (int k = 0; k < 4; k++) MPI_Reduce(&local_mean_values[k], &total_mean_values[k], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        for (int k = 0; k < 4; k++) total_mean_values[k] *= norm;

        variance_E = (total_mean_values[1] - total_mean_values[0]*total_mean_values[0])/(Nspins*Nspins);
        variance_M = (total_mean_values[3] - total_mean_values[2]*total_mean_values[2])/(Nspins*Nspins);

        heat_capacity = variance_E/(T*T);
        susceptibility = variance_M / (Nspins*Nspins*T);

        if (myrank == 0){            
            // Use this ofile to get temperature, heat capacity, susceptibility, mean energy
            // mean magnetic moment, energy variance and magnetic moment variance
            ofile << setprecision(8) <<    T << "  "  << heat_capacity <<  "  " << susceptibility << "  " <<  total_mean_values[0]/(Nspins*Nspins)  << "  " <<  total_mean_values[2]/(Nspins*Nspins)  << "  " <<  variance_E << "  " <<  variance_M << endl;
            cout << "Temperature = " << T << ", heat capacity = " << heat_capacity <<" after " << MPI_Wtime() - Timestart << "s" << endl;
        }
    } // End temperature loop
    double Timestop = MPI_Wtime();
    double Totaltime = Timestop - Timestart;

    MPI_Finalize();
    if (myrank == 0){
        //Write probability of state to file:
        N /= (MC_cycles - cut_off - 1 );
        //for (int i = 0; i< 2*max_energy +1; i++) ofile << N(i) << endl;


        //Closing statements
        cout << "Time  = " << Totaltime << "s on "<<numprocs << " processors " << endl;
        cout << "Energy variance = " << variance_E /(Nspins*Nspins) << "  Magnetic moment variance = " << variance_M/(Nspins*Nspins) << endl;
        cout << "Final energy = " << energy/ (Nspins*Nspins) << endl;
        cout << "Final magnetic moment = " << total_mean_values[2]/ (Nspins*Nspins) << endl;
        cout << "Mean energy = " <<  total_mean_values[0]/ (Nspins*Nspins) << endl;
    }

    return 0;
}
