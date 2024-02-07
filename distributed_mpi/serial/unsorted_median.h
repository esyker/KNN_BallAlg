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

#endif