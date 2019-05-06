#include <mpi.h>
#include <vector>
#include <cmath>
#include <iostream>
#include "jacobi_parallel.h"
#include <algorithm>
#include <iterator>

typedef struct dbl_twoindex_struct {
    double val;
    int    rank;
    int    k;
    int    j;
} dbl_twoindex;


void maxloc_dbl_twoindex(void *in, void *inout, int *len, MPI_Datatype *type){
    /* ignore type, just trust that it's our dbl_twoindex type */
    dbl_twoindex *invals    = static_cast<dbl_twoindex *>(in);
    dbl_twoindex *inoutvals = static_cast<dbl_twoindex *>(inout);

    for (int i=0; i<*len; i++) {
        if (std::abs(invals[i].val) > std::abs(inoutvals[i].val)) {
            inoutvals[i].val  = invals[i].val;
            inoutvals[i].rank = invals[i].rank;
            inoutvals[i].k = invals[i].k;
            inoutvals[i].j = invals[i].j;
        }

        else if (std::abs(invals[i].val) == std::abs(inoutvals[i].val)){
            if (invals[i].rank > inoutvals[i].rank){

                inoutvals[i].val  = invals[i].val;
                inoutvals[i].rank = invals[i].rank;
                inoutvals[i].k = invals[i].k;
                inoutvals[i].j = invals[i].j;

              }

        }
    }

    return;
}

double NormMatrix(Matrix m, const int nrows, const int ncols, const int row_offset) {
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

void find_abs_max(Matrix m, dbl_twoindex &max_ind, const int nrows, const int ncols, const int row_offset) {
  max_ind.val = m[0][row_offset + 1];//make sure last rank doesn't just have one row
  max_ind.k = row_offset + 1;
  max_ind.j = 0;
  //for (int j = 0; j < ncols; ++j) {
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      if (j != i+row_offset && std::abs(m[i][j]) > std::abs(max_ind.val)) {
        max_ind.val = m[i][j];
        max_ind.j = i;
        max_ind.k = j;
    }
  }
}
}

void ParallelJacobiRotate(Matrix m, Matrix v, dbl_twoindex max, std::vector<int> &offsets, std::vector<int> &nrows_vec, int rank) {

  double* mj = new double[nrows_vec[rank]];
  double* mk = new double[nrows_vec[rank]];

  std::vector<int>::iterator high;
  high = std::upper_bound (offsets.begin(), offsets.end(), max.k);
  int rank_k = high - offsets.begin() - 1;

  int local_k = max.k - offsets[rank_k];

  int rank_j = max.rank;
  int local_j = max.j;

  int jc = offsets[rank_j] + max.j;
  int kc = max.k;

  double mjj;
  double mjk;
  if(rank == rank_j)
  {
  mjj = m[local_j][jc];
  mjk = m[local_j][kc];
  }

  MPI_Bcast(&mjj, 1, MPI_DOUBLE, rank_j, MPI_COMM_WORLD);
  MPI_Bcast(&mjk, 1, MPI_DOUBLE, rank_j, MPI_COMM_WORLD);

  double mkk;
  if(rank == rank_k)
  {
  mkk = m[local_k][kc];
  }

  MPI_Bcast(&mkk, 1, MPI_DOUBLE, rank_k, MPI_COMM_WORLD);
  MPI_Scatterv(m[local_j],&nrows_vec[0],&offsets[0],MPI_DOUBLE, mj, nrows_vec[rank], MPI_DOUBLE, rank_j, MPI_COMM_WORLD);
  MPI_Scatterv(m[local_k],&nrows_vec[0],&offsets[0],MPI_DOUBLE, mk, nrows_vec[rank], MPI_DOUBLE, rank_k, MPI_COMM_WORLD);

  float c, s;
  if (mjj == mkk) {
    c = std::cos(M_PI / 4);
    s = std::sin(M_PI / 4);
  }
  else {
    float tau = (mjj - mkk) / (2 * mjk);
    float t = ((tau > 0) ? 1 : -1) / (std::abs(tau) + std::sqrt(1 + tau * tau));
    c = 1 / std::sqrt(1 + t * t);
    s = c * t;
  }

  if(rank == rank_k){  
      mj[local_k] = (c * c - s * s) * mjk + s * c * (mkk - mjj);
      mk[local_k] = s * s * mjj - 2 * s * c * mjk + c * c * mkk;
   }

  if(rank == rank_j)
  {
      mk[local_j] = (c * c - s * s) * mjk + s * c * (mkk - mjj);
      mj[local_j] = c * c * mjj + 2 * s * c * mjk + s * s * mkk;
  }

  double tmp_mjl;
  double tmp_vlj;

  int row_start = offsets[rank];
  int row_end = row_start + nrows_vec[rank];
  int row = 0;
  for (int l = row_start; l < row_end; ++l) {
    if (l != jc && l != kc) {
      tmp_mjl = mj[row];
      mj[row] = c * tmp_mjl + s * mk[row];
      mk[row] = s * tmp_mjl - c * mk[row];
      m[row][jc] = mj[row];
      m[row][kc] = mk[row];
    }
    tmp_vlj = v[row][jc];
    v[row][jc] =  c * tmp_vlj + s * v[row][kc];
    v[row][kc] = s * tmp_vlj - c * v[row][kc];
    row++;
  }

  MPI_Gatherv(mj, nrows_vec[rank],MPI_DOUBLE,m[local_j],&nrows_vec[0],&offsets[0],MPI_DOUBLE, rank_j, MPI_COMM_WORLD);
  MPI_Gatherv(mk, nrows_vec[rank],MPI_DOUBLE,m[local_k],&nrows_vec[0],&offsets[0],MPI_DOUBLE, rank_k, MPI_COMM_WORLD);
}

void ParallelJacobi(Matrix mat, Matrix v, std::vector<int> &offsets, std::vector<int> &nrows_vec, const int ncols, const float eps, const int rank) {

    int nrows = nrows_vec[rank];
    dbl_twoindex local_max, global_max;
    local_max.rank = rank;
    find_abs_max(mat, local_max, nrows, ncols, offsets[rank]);

    /* create our new data type */
    MPI_Datatype mpi_dbl_twoindex;
    MPI_Datatype types[4] = { MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT };
    MPI_Aint disps[4] = { offsetof(dbl_twoindex, val),
                     offsetof(dbl_twoindex, rank),
                     offsetof(dbl_twoindex, k),
                     offsetof(dbl_twoindex, j),  };
    int lens[4] = {1,1,1,1};
    MPI_Type_create_struct(4, lens, disps, types, &mpi_dbl_twoindex);
    MPI_Type_commit(&mpi_dbl_twoindex);

   /* create our operator */
    MPI_Op mpi_maxloc_dbl_twoindex;
    MPI_Op_create(maxloc_dbl_twoindex, 1, &mpi_maxloc_dbl_twoindex);
 
    MPI_Allreduce(&local_max, &global_max, 1, mpi_dbl_twoindex, mpi_maxloc_dbl_twoindex, MPI_COMM_WORLD);
    
    double local_norm = NormMatrix(mat, nrows, ncols, offsets[rank]);
    double global_norm;
    MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double tol = eps * global_norm;

    while (global_norm > tol) {
        
        ParallelJacobiRotate(mat, v, global_max, offsets, nrows_vec, rank);
        local_norm = NormMatrix(mat, nrows, ncols, offsets[rank]);
        MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        find_abs_max(mat, local_max, nrows, ncols, offsets[rank]);
        MPI_Allreduce(&local_max, &global_max, 1, mpi_dbl_twoindex, mpi_maxloc_dbl_twoindex, MPI_COMM_WORLD);
  }


}


