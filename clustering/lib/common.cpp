#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort

double calculate_similarity(std::pair<double, double> a, std::pair<double, double> b, double sig) {
    return exp(-1 * (pow(a.first - b.first, 2) + pow(a.second - b.second, 2)) / (2 * sig * sig));
}

Matrix EmpMatrix(int rows, int cols){
  Matrix matrix = new double*[rows];
  for (int i = 0; i < rows; i++) {
    matrix[i] = new double[cols];
  }
  return matrix;
  }

Matrix IdMatrix(int rows, int cols){
  Matrix matrix = EmpMatrix(rows, cols);
  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++){
      if(i==j){
        matrix[i][j] = 1;
      }
      else{
      matrix[i][j] = 0;
      }
    }
  }
  return matrix;
}

void PrintMatrix(Matrix m, int rows, int cols) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      printf("%2.2f ", m[i][j]);
    }
    printf("\n");
  }
}

void PrintMatrixFlat(double *m, int rows, int cols) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      printf("%2.2f ", m[i * cols + j]);
    }
    printf("\n");
  }
}

std::vector<size_t> sort_indexes(const std::vector<double> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}


//
//  I/O routines
//
/*
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}
*/

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
