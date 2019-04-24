#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <Eigen/Core>
#include "SpectralClustering.h"
#include <iostream>
#include <fstream>
#include "common.h"

double calculate_similarity(std::pair<double, double> a, std::pair<double, double> b) {
    return exp(-1 * (pow(a.first - b.first, 2) + pow(a.second - b.second, 2)) / 1000);
}

int main() {
    std::vector<std::pair<double, double> > points;
<<<<<<< HEAD
    std::ifstream file("large_points.txt");
=======
    std::ifstream file("points.txt");
>>>>>>> d3c6215a772efd8b0412c34bf9131a3d4dca7fac
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream linestream(line);
        std::string x;
        std::string y;

        linestream >> x >> y;
        points.push_back(std::make_pair(std::stod(x), std::stod(y)));

    }

    std::cout << points.size() << std::endl; 
    //std::vector<int> items = {1,2,3,4,5,16,17,18,19,20};

    // generate similarity matrix
 
    double simulation_time = read_timer( );
    unsigned int size = points.size();
    unsigned int numClusters = 2;
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(size,size);
    double similarity;
    for (unsigned int i=0; i < size; i++) {
        for (unsigned int j=0; j < size; j++) {
            //if (i == j) similarity = 0;
            //else {
            // generate similarity
            similarity = calculate_similarity(points[i], points[j]);
            //}
            //std::cout << similarity << std::endl;
            m(i,j) = similarity;
            m(j,i) = similarity;
        }
    }

    simulation_time = read_timer( ) - simulation_time;
  
    printf( "simulation time = %g seconds \n",  simulation_time);
    
    // the number of eigenvectors to consider. This should be near (but greater) than the number of clusters you expect. Fewer dimensions will speed up the clustering
    /*int numDims = 3;
    // do eigenvalue decomposition
    SpectralClustering* c = new SpectralClustering(m, numDims);

    // whether to use auto-tuning spectral clustering or kmeans spectral clustering
    bool autotune = false;

    std::vector<std::vector<int> > clusters;
    if (autotune) {
        // auto-tuning clustering
        clusters = c->clusterRotate();
    } else {
        // how many clusters you want
        clusters = c->clusterKmeans(numClusters);
    }

    // output clustered items
    // items are ordered according to distance from cluster centre
    /*for (unsigned int i=0; i < clusters.size(); i++) {
        std::cout << "Cluster " << i << ": " << "Item ";
        std::copy(clusters[i].begin(), clusters[i].end(), std::ostream_iterator<int>(std::cout, ", "));
        std::cout << std::endl;
    } */
}
