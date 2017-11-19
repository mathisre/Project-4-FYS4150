#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <time.h>
#include <random>
using namespace arma;
using namespace std;
ofstream ofile;

int main(int argc, char* argv[])
{
    string filename;
    filename = argv[1];
    ofile.open(filename);

    int n = 10000;
    ofile << n << endl;

    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Set up the uniform distribution for x \in [[0, 1]
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    double *X_array;
    X_array = new double[n];
    double X_tot = 0;
    double X_tot_sqr = 0;
    double x;

    for (int i = 0; i < n; i++){
        X_array[i] = RandomNumberGenerator(gen);
        X_tot += X_array[i];
        X_tot_sqr += X_array[i]*X_array[i];
        //ofile << setprecision(20) << x << "   " << x*x << endl;
    }    

    double Mean = X_tot / n;
    double Mean_sqr = X_tot_sqr / n;
    double Variance = sqrt(Mean_sqr - Mean*Mean);


    // Now we compute the autocorrelation function
    double *autocor;  autocor = new double[n];
    for (int j = 0; j < n; j++){
      double sum = 0.0;
      for (int k = 0; k < (n-j); k++){
          sum  += (X_array[k]-Mean)*(X_array[k+j]-Mean);
      }
      autocor[j] = sum/Variance/((double) n );
      //ofile << setiosflags(ios::showpoint | ios::uppercase);
      ofile << setw(15) << setprecision(8) << j;
      ofile << setw(15) << setprecision(8) << autocor[j] << endl;
    }
    ofile.close();  // close output file
    return 0;

}
