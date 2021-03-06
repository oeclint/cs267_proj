#include "stdio.h"
#include "time.h"
#include <cstdlib>
#include <cstring>
#include "jacobi_serial.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include "common.h"
#include <Eigen/Eigenvalues>

int main(int argc, char **argv) {
    std::vector<std::pair<double, double> > points;
    char *filename = read_string(argc, argv, "-f", NULL);
    int num_vec = read_int(argc, argv, "-k", 1);
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
    SerialJacobi(m, v, size, 1e-5);
    simulation_time = read_timer( ) - simulation_time;
    printf( "n = %d, simulation time = %g seconds \n", size, simulation_time);

    std::vector<double> eigs;
    for(int i =0; i < size; i++){
        eigs.push_back(m[i][i]);
    }
    
    std::vector<size_t> argsort = sort_indexes(eigs);
    Matrix vk = EmpMatrix(size, num_vec);

    for (int i =0; i < size; i++){
        for (int j =0; j < num_vec; j++) {
            vk[i][j] = v[i][argsort[j]];
        }
    }

    //Sanity checks

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
    PrintMatrix(vk, size, num_vec);
    std::cout <<es.eigenvalues()<<std::endl;

  return 0;
}

