#include <mpi.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include "SpectralClustering.h"
#include <iostream>
#include <fstream>
#include "jacobi_parallel.h"

int main( int argc, char **argv )
{
    char *filename = read_string( argc, argv, "-f", NULL );
    std::vector<std::pair<double, double> > points;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream linestream(line);
        std::string x;
        std::string y;

        linestream >> x >> y;
        points.push_back(std::make_pair(std::stod(x), std::stod(y)));
    }

    int nrows = points.size();
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
    
    int nrows_local = (stop - start + 1);
    Matrix local_sim_mat = EmpMatrix(nrows_local, ncols);
    Matrix v = EmpMatrix(nrows_local, ncols);

    int row = 0;
    for(int i=start; i<=stop; i++){
        for(int j=0; j<ncols; j++){
            local_sim_mat[row][j] = calculate_similarity(\
                points[i], points[j], 2.236);
            if (i == j){
                v[row][j] = 1;
            }
            else{
                v[row][j] = 0;

            }
        }
        row++;
    }

    std::vector<int> offsets;

    offsets.resize(n_proc);

    MPI_Allgather(&start,1,MPI_INT,&offsets[0],1,MPI_INT,MPI_COMM_WORLD);

    std::vector<int> nrows_vec;

    for(int i=0; i< offsets.size()-1; i++){
        nrows_vec.push_back(offsets[i+1]-offsets[i]);
    }

    nrows_vec.push_back(ncols - offsets.back());

    ParallelJacobi(local_sim_mat, v, offsets, nrows_vec, ncols, 1e-5, rank);
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
      printf( "n = %d, simulation time = %g seconds \n", nrows, simulation_time);
    }


    double* local_eigs = new double[nrows_vec[rank]];
    //double* eigs = new double[nrows];
    std::vector<double> eigs;

    eigs.resize(nrows);

    row = 0;
    for(int i=start; i<=stop; i++){
        for(int j=0; j<ncols; j++){
            if(i == j){
                local_eigs[row] = local_sim_mat[row][j];
            }
        }
        row++;
    }

    MPI_Allgatherv(local_eigs, nrows_vec[rank],MPI_DOUBLE,&eigs[0],&nrows_vec[0],&offsets[0],MPI_DOUBLE, MPI_COMM_WORLD);

    if(rank==0){

        for (auto i: sort_indexes(eigs)) {
          std::cout << eigs[i] << std::endl;
          }
    }

    MPI_Finalize( ); 
    return 0;
}
