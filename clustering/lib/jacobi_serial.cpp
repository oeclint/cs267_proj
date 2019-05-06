#include "jacobi_serial.h"
#include <cstdio>
#include <cmath>

struct Ind2D {
  Ind2D(int _i=0, int _j=0) {
    k = _i;
    j = _j;
  }
  int k;
  int j;
};

Ind2D find_abs_max(Matrix m, const int n) {
  double max = m[1][0];
  Ind2D max_ind(1, 0);
  for (int k = 0; k < n; ++k) {
    for (int j = 0; j < n; ++j) {
      if (j != k && std::abs(m[j][k]) > std::abs(max)) {
        max = m[j][k];
        max_ind.j = j;
        max_ind.k = k;
      }
    }
  }
  return max_ind;
}

double NormMatrix(Matrix m, const int n) {
  double norm = 0;
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < n; ++k) {
      if (j != k) {
        norm += m[j][k] * m[j][k];
      }
    }
  }
  return std::sqrt(norm);
}

void SerialJacobiRotate(Matrix m, Matrix v, const int j, const int k, const int n) {
  double c, s;
  if (m[j][j] == m[k][k]) {
    c = std::cos(M_PI / 4);
    s = std::sin(M_PI / 4);
  }
  else {
    double tau = (m[j][j] - m[k][k]) / (2 * m[j][k]);
    double t = ((tau > 0) ? 1 : -1) / (std::abs(tau) + std::sqrt(1 + tau * tau));
    c = 1 / std::sqrt(1 + t * t);
    s = c * t;
  }

  double tmp_mjk = m[j][k];
  double tmp_mjj = m[j][j];
  m[j][k] = (c * c - s * s) * tmp_mjk + s * c * (m[k][k] - m[j][j]);
  m[k][j] = m[j][k];
  m[j][j] = c * c * tmp_mjj + 2 * s * c * tmp_mjk + s * s * m[k][k];
  m[k][k] = s * s * tmp_mjj - 2 * s * c * tmp_mjk + c * c * m[k][k];

  double tmp_mjl;
  double tmp_vlj;

  for (int l = 0; l < n; ++l) {
    if (l != j && l != k) {
      tmp_mjl = m[j][l];
      m[j][l] = c * tmp_mjl + s * m[k][l];
      m[k][l] = s * tmp_mjl - c * m[k][l];
      m[l][j] = m[j][l];
      m[l][k] = m[k][l];}

      tmp_vlj = v[l][j];
      v[l][j] =  c * tmp_vlj + s * v[l][k];
      v[l][k] = s * tmp_vlj - c * v[l][k];

  }
}


void SerialJacobi(Matrix mat, Matrix v, const int n, const double eps) {
  Ind2D ind_max;
  ind_max = find_abs_max(mat, n);
  double norm = NormMatrix(mat, n);
  double tol = eps * norm;
//  printf("eps = %f, norm = %f, tol = %f\n",eps, norm, tol);
  while (norm > tol) {
    //printf("%f ", norm);
    SerialJacobiRotate(mat, v, ind_max.j, ind_max.k, n);
    norm = NormMatrix(mat, n);
    //PrintMatrix(mat, n);
    ind_max = find_abs_max(mat, n);
//    printf("eps = %f, norm = %f, tol = %f\n",eps, norm, tol);
  }
}

