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

inline int min( int a, int b ) { return a < b ? a : b; }
using namespace std;

double calculate_similarity(double a1, double a2, double b1, double b2) {
    return exp(-1 * (pow(a1 - b1, 2) + pow(a2 - b2, 2)) / 1000);
}

int main( int argc, char **argv )
{
    std::vector<double> points_vec;
    std::ifstream file("points.txt");
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream linestream(line);
        std::string x;
        std::string y;

        linestream >> x >> y;
        points_vec.push_back(std::stod(x)); 
        points_vec.push_back(std::stod(y));

    }

    int nrows = points_vec.size()/2;
    int ncols = nrows;

    int n = nrows * ncols;
   
    double* points = points_vec.data();
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
     
    int nrows_per_proc = (nrows + n_proc - 1) / n_proc;

    int *partition_offsets = (int*) malloc((n_proc + 1) * sizeof(int) );
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );

    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( 2 * i * nrows_per_proc , 2 * nrows );
    
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];

    double *local_points = (double*) malloc( nrows_per_proc * 2 *  sizeof(double) );    
    //
    //  allocate storage for local partition
    //
    int nstart = partition_offsets[rank];
    int nlocal = partition_sizes[rank];

//    MPI_Scatterv(points, partition_sizes, partition_offsets, MPI_DOUBLE,\
        local_points, nlocal, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    double *local_sim_mat = (double*) malloc( nrows_per_proc * ncols *  sizeof(double) );
    int row = 0;

    for(int i=nstart; i<nstart+nlocal; i+=2){
        for(int j=0; j<nrows*2; j+=2){
            local_sim_mat[row * ncols + j]= calculate_similarity(points[i], points[i+1], points[j], points[j+1]);
        }
        row++;
    }

    free( local_sim_mat );
    free( partition_offsets);
    free( partition_sizes);
    
    MPI_Finalize( );
    
    return 0;
}
