#include <mpi.h>
#include <math.h>
#include <stdlib.h>

inline int min( int a, int b ) { return a < b ? a : b; }

using namespace std;

int main( int argc, char **argv )
{    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    int nrows;
    int ncols;

    nrows = ncols = 10;
    int n = nrows * ncols;
   
    double *sim_mat = (double*) malloc( n * sizeof(double) );

    int nrows_per_proc = (nrows + n_proc - 1) / n_proc;

    double *local_sim_mat = (double*) malloc( nrows_per_proc * ncols *  sizeof(double) );

    int *partition_offsets = (int*) malloc((n_proc + 1) * sizeof(int) );
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );

    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * nrows_per_proc * ncols, n );
    
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    
    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];

    MPI_Scatterv(sim_mat, partition_sizes, partition_offsets, MPI_DOUBLE, local_sim_mat, nlocal, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    free( sim_mat);
    free( local_sim_mat );
    free( partition_offsets);
    free( partition_sizes);
    
    MPI_Finalize( );
    
    return 0;
}
