#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <Eigen/Core>
#include "SpectralClustering.h"
#include <iostream>
#include <fstream>
#include "common.h"

using namespace std;

double calculate_similarity(double a1, double a2, double b1, double b2) {
    return exp(-1 * (pow(a1 - b1, 2) + pow(a2 - b2, 2)) / 10);
}

int main( int argc, char **argv )
{
    char *filename = read_string( argc, argv, "-f", NULL );
    std::vector<double> points;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream linestream(line);
        std::string x;
        std::string y;

        linestream >> x >> y;
        points.push_back(std::stod(x)); 
        points.push_back(std::stod(y));

    }

    int nrows = points.size()/2;
    int ncols = nrows;

    int n = nrows * ncols;
   
    //double* points = points_vec.data();
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    double simulation_time = read_timer( );
  
    int count = nrows / n_proc;
    int remainder = nrows % n_proc;
    int start, stop;

    if (rank < remainder) {
        // The first 'remainder' ranks get 'count + 1' tasks each
        start = rank * (count + 1);
        stop = start + count;
    } else {
        // The remaining 'size - remainder' ranks get 'count' task each
        start = rank * count + remainder;
        stop = start + (count - 1);
    }
     
    long local_size = (stop - start + 1) * ncols;
    double *local_sim_mat = (double*) malloc( local_size * sizeof(double) );
    int row = 0;
    int col = 0;
    for(int i=2*start; i<=2*stop; i+=2){
        for(int j=0; j<nrows*2; j+=2){
            long cur_pos = row*ncols + col;
            local_sim_mat[cur_pos] = calculate_similarity(points[i], points[i+1], points[j], points[j+1]);
            col++;
        }
        col = 0;
        row++;
    }

    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
      printf( "n = %d, simulation time = %g seconds \n", nrows, simulation_time);
    } 
    MPI_Finalize( ); 
    return 0;
}
