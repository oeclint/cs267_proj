#include "mpi.h"
#include "stdio.h"
#include "time.h"
#include <cstdlib>
#include <cstring>
#include "jacobi.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include "common.h"
#include <Eigen/Eigenvalues>

/*
int main(int argc, char** argv){
  int n = std::atoi(argv[1]);
  std::cout << n << std::endl;
  double start, end;
  start = 0;
  end = 0;
  double elapsed_time = 0;

  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  elapsed_time -= MPI_Wtime();
  Matrix m = SquareMatrix(n);
  PrintMatrix(m, n);
  printf("\n");
  ParallelJacobi(m, n, 1e-5);
  elapsed_time += MPI_Wtime();
  *if (rank == 0) {
    for (int i = 0; i < n; ++i) {
      printf("%f ", m[i][i]);
    }
    printf("Dimension %i, time elapsed %f\n", n, elapsed_time);
    MPI_Abort(MPI_COMM_WORLD, MPI_SUCCESS);
  } *
  MPI_Finalize();
  for (int i = 0; i < n; ++i) {
    delete[] m[i];
  }
  delete[] m;
  return 0;
}
*/

/* Тестирование последовательного алгоритма Якоби*/
/*int main(int argc, char** argv){
  int dimesions[8] = {8, 16, 32, 64, 128, 256, 512, 1024};
  clock_t t;
  double elapsed_secs;
  for (int i = 0; i < 8; ++i) {
    t = clock();
    SerialJacobi(SquareMatrix(dimesions[i]), dimesions[i], 1e-5);
    elapsed_secs = double(clock() - t) / CLOCKS_PER_SEC;
    printf("Dimension %i, elapsed time %f\n", dimesions[i], elapsed_secs);
  }
  return 0;
}*/

int main(int argc, char **argv) {
    std::vector<std::pair<double, double> > points;
    char *filename = read_string(argc, argv, "-f", NULL);
    std::ifstream file(filename);
    std::string line;


    while (std::getline(file, line)) {
        std::stringstream linestream(line);
        std::string x;
        std::string y;

        linestream >> x >> y;
        points.push_back(std::make_pair(std::stod(x), std::stod(y)));

    }

    unsigned int size = points.size();
    double similarity;

    Matrix m = EmpMatrix(size, size);

    double simulation_time = read_timer( );
    for (unsigned int i=0; i < size; i++) {
        for (unsigned int j=0; j < size; j++) {
            similarity = calculate_similarity(points[i], points[j], 2.236);
            m[i][j] = similarity;
            m[j][i] = similarity;
        }

    }

    Matrix v = IdMatrix(size, size);
    SerialJacobiV(m, v, size, 1e-5);
    simulation_time = read_timer( ) - simulation_time;
    printf( "n = %d, simulation time = %g seconds \n", size, simulation_time);

    Eigen::MatrixXd A(size,size);

    for (unsigned int i=0; i < size; i++) {
        for (unsigned int j=0; j < size; j++) {
            similarity = calculate_similarity(points[i], points[j], 2.236);
            A(i,j) = similarity;
            A(j,i) = similarity;
        }

    }

    Eigen::EigenSolver<Eigen::MatrixXd> es(A, true);
    std::cout <<es.eigenvectors()<<std::endl;
    PrintMatrix(v, size);
    std::cout <<es.eigenvalues()<<std::endl;

    for(int i =0; i < size; i++){
        std::cout << m[i][i] << std::endl;

    }
  return 0;
}

