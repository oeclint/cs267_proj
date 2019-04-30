#ifndef JACOBI
#define JACOBI

typedef double** Matrix;

void SerialJacobi(Matrix mat, const int n, const double eps);
void SerialJacobiV(Matrix mat, Matrix v, const int n, const double eps);

void ParallelJacobi(Matrix mat, const int n, const double eps);
void ParallelJacobiV(Matrix mat, Matrix v, const int n, const double eps);

void PrintMatrix(Matrix mat, const int i);
#endif // JACOBI

