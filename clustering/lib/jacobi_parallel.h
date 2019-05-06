#include "common.h"
#ifndef JACOBIPARALLEL
#define JACOBIPARALLEL

void ParallelJacobi(Matrix mat, std::vector<int> &offsets, std::vector<int> &nrows_vec, const int ncols, const float eps, const int rank);
#endif

