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

typedef double** DMatrix;

typedef struct dbl_twoindex_struct {
    double val;
    int    rank;
    int    i;
    int    j;
} dbl_twoindex;


void maxloc_dbl_twoindex(void *in, void *inout, int *len, MPI_Datatype *type){
    /* ignore type, just trust that it's our dbl_twoindex type */
    dbl_twoindex *invals    = static_cast<dbl_twoindex *>(in);
    dbl_twoindex *inoutvals = static_cast<dbl_twoindex *>(inout);

    for (int i=0; i<*len; i++) {
        if (invals[i].val > inoutvals[i].val) {
            inoutvals[i].val  = invals[i].val;
            inoutvals[i].rank = invals[i].rank;
            inoutvals[i].i = invals[i].i;
            inoutvals[i].j = invals[i].j;
        }

        else if (invals[i].val == inoutvals[i].val){
            if (invals[i].rank > inoutvals[i].rank){

                inoutvals[i].val  = invals[i].val;
                inoutvals[i].rank = invals[i].rank;
                inoutvals[i].i = invals[i].i;
                inoutvals[i].j = invals[i].j;

              }

        }
    }

    return;
}

struct Ind2Dmax {
  Ind2Dmax(int _i=0, int _j=0, double _max=0) {
    i = _i;
    j = _j;
    max = _max;

  }
  int i;
  int j;
  double max;
};

double calculate_similarity(double a1, double a2, double b1, double b2) {
    return exp(-1 * (pow(a1 - b1, 2) + pow(a2 - b2, 2)) / 10);
}

double NormMatrix(DMatrix m, const int nrows, const int ncols, const int row_offset) {
  double norm = 0;
  for (int j = 0; j < nrows; ++j) {
    for (int k = 0; k < ncols; ++k) {
      if (j + row_offset != k) {
        norm += m[j][k] * m[j][k];
      }
    }
  }
  return std::sqrt(norm);
}

void find_abs_max(DMatrix m, dbl_twoindex &max_ind, const int nrows, const int ncols, const int row_offset) {
  max_ind.val = m[1][0];
  max_ind.i = 1;
  max_ind.j = 0;
  //for (int j = 0; j < ncols; ++j) {
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      if (j != i+row_offset && abs(m[i][j]) > abs(max_ind.val)) {
        max_ind.val = m[i][j];
        max_ind.j = i;
        max_ind.i = j;
    }
  }
}
}

void ParallelJacobiRotate(DMatrix m, dbl_twoindex max, std::vector<int> &offsets, std::vector<int> &nrows_vec, int rank) {

  double* mj = new double[nrows_vec[rank]];
  double* mk = new double[nrows_vec[rank]];

  std::vector<int>::iterator high;
  high = std::upper_bound (offsets.begin(), offsets.end(), max.j);
  int rank_j = high - offsets.begin() - 1;

  //std::cout<<"maxj"<<max.j<<std::endl;
  //std::cout<<"ooo  "<<rank_j<<" "<< max.j<<std::endl;
  int local_j = max.j - offsets[rank_j];

  int rank_k = max.rank;
  int local_k = max.i;

  int kc = offsets[max.rank] + max.i;
  int jc = max.j;

  double mjj;
  double mjk;
  if(rank == rank_j)
  {
  mjj = m[local_j][jc];
  mjk = m[local_j][kc];
  }

  MPI_Bcast(&mjj, 1, MPI_DOUBLE, rank_j, MPI_COMM_WORLD);
  MPI_Bcast(&mjk, 1, MPI_DOUBLE, rank_j, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  double mkk;
  if(rank == rank_k)
  {
  mkk = m[local_k][kc];
  }

  MPI_Bcast(&mkk, 1, MPI_DOUBLE, rank_k, MPI_COMM_WORLD);
  MPI_Scatterv(m[local_j],&nrows_vec[0],&offsets[0],MPI_DOUBLE, mj, nrows_vec[rank], MPI_DOUBLE, rank_j, MPI_COMM_WORLD);
  MPI_Scatterv(m[local_k],&nrows_vec[0],&offsets[0],MPI_DOUBLE, mk, nrows_vec[rank], MPI_DOUBLE, rank_k, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  float c, s;
  if (mjj == mkk) {
    c = cos(M_PI / 4);
    s = sin(M_PI / 4);
  }
  else {
    float tau = (mjj - mkk) / (2 * mjk);
    float t = ((tau > 0) ? 1 : -1) / (abs(tau) + sqrt(1 + tau * tau));
    c = 1 / sqrt(1 + t * t);
    s = c * t;
  }

  if(rank = rank_k){  
      mj[local_k] = (c * c - s * s) * mjk + s * c * (mkk - mjj);
      mk[local_k] = s * s * mjj - 2 * s * c * mjk + c * c * mkk;
   }

  if(rank = rank_j)
  {
      mk[local_j] = (c * c - s * s) * mjk + s * c * (mkk - mjj);
      mj[local_j] = c * c * mjj + 2 * s * c * mjk + s * s * mkk;
  }

  double tmp_jl;
  int row_start = offsets[rank];
  int row_end = row_start + nrows_vec[rank];
  int row = 0;
  for (int l = row_start; l < row_end; ++l) {
    if (l != jc && l != kc) {
      tmp_jl = mj[row];
      mj[row] = c * tmp_jl + s * mk[row];
      mk[row] = s * tmp_jl - c * mk[row];
      m[row][jc] = mj[row];
      m[row][kc] = mk[row];
    }
    row++;
  }

  MPI_Gatherv(mj, nrows_vec[rank],MPI_DOUBLE,m[local_j],&nrows_vec[0],&offsets[0],MPI_DOUBLE, rank_j, MPI_COMM_WORLD);
  MPI_Gatherv(mk, nrows_vec[rank],MPI_DOUBLE,m[local_k],&nrows_vec[0],&offsets[0],MPI_DOUBLE, rank_k, MPI_COMM_WORLD);
}

void ParallelJacobi(DMatrix mat, std::vector<int> &offsets, const int ncols, const float eps, const int rank) {
    
    std::vector<int> nrows_vec;

    for(int i=0; i< offsets.size()-1; i++){
        nrows_vec.push_back(offsets[i+1]-offsets[i]);
    }

    nrows_vec.push_back(ncols - offsets.back());

    int nrows = nrows_vec[rank];
    dbl_twoindex local_max, global_max;
    local_max.rank = rank;
    find_abs_max(mat, local_max, nrows, ncols, offsets[rank]);

    /* create our new data type */
    MPI_Datatype mpi_dbl_twoindex;
    MPI_Datatype types[4] = { MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT };
    MPI_Aint disps[4] = { offsetof(dbl_twoindex, val),
                     offsetof(dbl_twoindex, rank),
                     offsetof(dbl_twoindex, i),
                     offsetof(dbl_twoindex, j),  };
    int lens[4] = {1,1,1,1};
    MPI_Type_create_struct(4, lens, disps, types, &mpi_dbl_twoindex);
    MPI_Type_commit(&mpi_dbl_twoindex);

   /* create our operator */
    MPI_Op mpi_maxloc_dbl_twoindex;
    MPI_Op_create(maxloc_dbl_twoindex, 1, &mpi_maxloc_dbl_twoindex);
 
    MPI_Allreduce(&local_max, &global_max, 1, mpi_dbl_twoindex, mpi_maxloc_dbl_twoindex, MPI_COMM_WORLD);
//   std::cout<<"rank "<< rank << " max " << local_max.val <<" grank "<< global_max.rank << " gmax " << global_max.val << " i " << global_max.i << " j " << global_max.j << std::endl;
    double local_norm = NormMatrix(mat, nrows, ncols, offsets[rank]);
    double global_norm;
    MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double tol = eps * global_norm;

//    ParallelJacobiRotate(mat, global_max, offsets, nrows_vec, rank);
//    std::cout<<"rank "<< rank<< " lnorm " << local_norm << " norm " << norm << std::endl;
//    if(rank==global_max.rank)
//        std::cout<<"rank " << rank << " max " << mat[global_max.i][global_max.j] <<  std::endl;

    int num_iter = 0;   
    while (global_norm > tol) {
            std::cout << "rank "<< rank << "Norm: " << global_norm << "/" << tol << std::endl;
            std::cout << "rank "<< rank << "Iter: " << num_iter << std::endl;
            std::cout<<"rank "<< rank << " max " << local_max.val <<" grank "<< global_max.rank << " gmax " << global_max.val << " i " << global_max.i << " j " << global_max.j << std::endl;
        if (rank == 0) {
            if(num_iter==1000000){
            MPI_Abort(MPI_COMM_WORLD, MPI_SUCCESS);
            }/*
            std::cout << "rank "<< rank << "Norm: " << global_norm << "/" << tol << std::endl;
            std::cout << "rank "<< rank << "Iter: " << num_iter << std::endl;
            std::cout<<"rank "<< rank << " max " << local_max.val <<" grank "<< global_max.rank << " gmax " << global_max.val << " i " << global_max.i << " j " << global_max.j << std::endl;
           */
            //num_iter++;
        }
        num_iter++;
        ParallelJacobiRotate(mat, global_max, offsets, nrows_vec, rank);
        MPI_Barrier(MPI_COMM_WORLD);
        local_norm = NormMatrix(mat, nrows, ncols, offsets[rank]);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        find_abs_max(mat, local_max, nrows, ncols, offsets[rank]);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&local_max, &global_max, 1, mpi_dbl_twoindex, mpi_maxloc_dbl_twoindex, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);  
  }


}


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
    DMatrix local_sim_mat = new double*[nrows_local];
    for (int i = 0; i < nrows_local; i++) {
        local_sim_mat[i] = new double[ncols];
      }

    int row = 0;
    for(int i=start; i<=stop; i++){
        for(int j=0; j<ncols; j++){
            local_sim_mat[row][j] = calculate_similarity(\
                points[i].first, points[i].second, points[j].first, points[j].second);
        }
        row++;
    }

    std::vector<int> offsets;

    offsets.resize(n_proc);

    MPI_Allgather(&start,1,MPI_INT,&offsets[0],1,MPI_INT,MPI_COMM_WORLD);

    //double local_norm = NormMatrix(local_sim_mat, nrows_local, ncols);
    //double norm;
    //MPI_Allreduce(&local_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //std::cout<<"rank "<< rank<< " lnorm " << local_norm << " norm " << norm << std::endl;
    ParallelJacobi(local_sim_mat, offsets, ncols, 1e-5, rank);
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
      printf( "n = %d, simulation time = %g seconds \n", nrows, simulation_time);
    } 
    MPI_Finalize( ); 
    return 0;
}
