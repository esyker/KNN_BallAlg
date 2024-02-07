#ifndef _MEDIAN_
#define _MEDIAN_

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

int cmpfunc (const void * a, const void * b);
void printArray(double *a, int n_items);
void arraycpy(double* dest, double* sourc, int n_items);
double sorted_median(double *vector, int n_items);
double buildSet(double *setL, double *setR, int* counterL, int* counterR, double *vector, int n_items, double median);
double median(double *vector, int n_items);
double getSmallest(double *vector, int n_items);
double getKsmallest(double* vector, long k, long n_items);
long find_idx_from_value(double * projections, long n_points, double value);
void compare_with_median (double* projections, double** pts, double median, long n_points);
long getLowerNeighborIdx(double* vector, long n_items, double result);
void linerize(double** matrix, int nx, int ny, double* vector);
double getIntervalSet(double *interval, double *vector, int n_items, double comparatorL, double comparatorR);
int appendArray(double* dest, int begin_idx, double* source, int n_elem, double ignore);
double PSRS(double* vector, long n_items);
void memset_double(double* vector, double to_set, int n_elem);
int getVectorSum(int* vector, int n_elem);
int find_K_idx(int *vector, int k, int n_elem, int* median_holder);
double medianSort(double* vector, int n_items);
void checkMalloc(void* pointer, int tag);


#endif