#include <iostream>
#include <vector>
#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  timing routines
//
double read_timer( );

//
//  clustering routines
//
double calculate_similarity(std::pair<double, double> a, std::pair<double, double> b, double sig);

//
// matrix routines
//
typedef double** Matrix;
Matrix EmpMatrix(int rows, int cols);
Matrix IdMatrix(int rows, int cols);
void PrintMatrix(Matrix m, int rows, int cols);
std::vector<size_t> sort_indexes(const std::vector<double> &v);
//
//  I/O routines
//
/*
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );
*/
//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
