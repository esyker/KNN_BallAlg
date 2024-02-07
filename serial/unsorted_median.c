#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "unsorted_median.h"

/**
 * Compare function used inside quicksort
 * @param   a   First number to compare
 * @param   b   Second number to compare
 * @return  Returns 1 if a>b, -1 if a<b and 0 if a=b  
 */
int cmpfunc (const void * a, const void * b)
{
  if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;
}

/**
 * Copies n_items from the source array to the destination array
 * @param   dest    Destination array
 * @param   sourc   Source array
 * @param   n_items Number of items to copy  
 */
void arraycpy(double* dest, double* sourc, int n_items)
{
    for(int i=0; i<n_items; i++)
    {
        dest[i] = sourc[i];
    }
}


/**
 * Compute the median of a vector using quicksort
 * @param   vector    Vector to find the median of
 * @param   n_items   Number of items of the vector
 * @return  Median of the vector  
 */
double sorted_median(double *vector, int n_items)
{
    double result;
    
    if(n_items==1)
    {
        result = (double)(vector[0]);
        return result;
    }

    if(n_items==2)
    {
        result = (double)(vector[0] + vector[1])/2;
        return result;
    }

    qsort(vector, n_items, sizeof(double), cmpfunc);

    if(n_items % 2 != 0)
    {
        return vector[n_items/2];
    }
    else
    {
        //printf("Sorted: %lf | %lf\n", vector[n_items/2], vector[n_items/2 - 1]);
        return 0.5 * ((double)vector[n_items/2] + vector[n_items/2 - 1]);
    }
}


/**
 * Builds two sets from the given vector, the left set has all the points lower than the median
 * and the right set has all the points higher than the median.
 * [Warning: The median is a prediction so the sets can be larger than half of the vector] 
 * @param   setL    Left set 
 * @param   setR    Rigth set
 * @param   counterL    Return by reference the number of elements on the left set
 * @param   counterR    Return by reference the number of elements on the rigth set
 * @param   vector  Vector to separate in two sets
 * @param   n_items     Number of elements on the vector
 * @param   median  Median to use to separate the vector
 * @return  Number of elements on the left set which is equal to the indice of the media
 */
double buildSet(double *setL, double *setR, int* counterL, int* counterR, double *vector, int n_items, double median)
{
    *counterL=0;
    *counterR=0;

    for(int i=0; i<n_items; i++)
    {
        //printArray(setL, *counterL);
        if(vector[i]<median)
        {
            (*counterL)++;
            setL[(*counterL)-1] = vector[i];
        }
        else
        {
            (*counterR)++;
            setR[(*counterR)-1] = vector[i];
        }
    }

    //printArray(setL, *counterL);
    //printf("Sets criados: setL %d items, setR %d items, idx da mediana: %d\n\n", *counterL, *counterR, (*counterL));
    return (*counterL);
}


/**
 * Compute the median of a vector using quicksort
 * @param   vector    Vector to find the median of
 * @param   n_items   Number of items of the vector
 * @return  Median of the vector
 */
double median(double *vector, int n_items)
{
    int full_splits = n_items/5;
    int semi_splits = n_items % 5;
    int i=0;
    double result;
    int n_medians = 0;
    double *medians = (double *)malloc(n_items/2 * sizeof(double));
    
    //printArray(vector, 17);

    for(i=0; i<full_splits; i++)
    {
        medians[i] = sorted_median(vector + 5*i, 5);
        //printf("Median do grupo %d: %lf\n", i, medians[i]);
    }

    if(semi_splits != 0)
    {
        medians[i] = sorted_median(vector + 5*full_splits, semi_splits);
        //printf("Median do semi grupo %d: %lf\n", i, medians[i]);
    }

    n_medians = full_splits + (semi_splits!=0 ? 1 : 0);

    if(n_medians <= 5)
    {
        result = sorted_median(medians, n_medians);
        //printf("Mediana das medianas: %lf\n", result);
        //result = buildSet(setL, setR, vector, n_items, result);
        free(medians);
        return result;
    }

    result = median(medians, full_splits + (semi_splits!=0 ? 1 : 0));
    free(medians);

    //printf("Mediana das medianas: %lf\n", result);

    //result = buildSet(setL, setR, vector, n_items, result);

    return result;
}


/**
 * Finds the smallest element of a vector
 * @param   vector    Vector to find the median of
 * @param   n_items   Number of items of the vector
 * @return  Values of the smallest element
 */
double getSmallest(double *vector, int n_items)
{
    double min = DBL_MAX;

    for(int i=0; i<n_items; i++)
    {
        if(vector[i] < min)
            min = vector[i];
    }
    return min;
}


/**
 * Finds the k smallest element of a vector
 * @param   vector  Vector to find the median of
 * @param   k   Element indice (counting from 0) to be found
 * @param   n_items Number of items of the array
 * @return  Value of the k smallest element
 */
double getKsmallest(double* vector, long k, long n_items)
{
    int L_items;
    int R_items;
    double result;
    int result_idx;
    double *setL = (double *)malloc(n_items * sizeof(double));
    double *setR = (double *)malloc(n_items * sizeof(double));
    double *vector_cpy = (double *)malloc(n_items * sizeof(double));
    
    arraycpy(vector_cpy, vector, n_items);
    
    //printArray(vector, n_items);

    result = median(vector_cpy, n_items);

    arraycpy(vector_cpy, vector, n_items);

    result_idx = buildSet(setL, setR, &L_items, &R_items, vector_cpy, n_items, result);

    //printf("Target idx: %d, result idx: %d\n", k, result_idx);

    while(k > result_idx || k < result_idx)
    {
        if(k < result_idx)
        {
            arraycpy(vector_cpy, setL, L_items);

            result = median(vector_cpy, L_items);
            result_idx = buildSet(setL, setR, &L_items, &R_items, vector_cpy, L_items, result);
            //printf("Target idx: %d, result idx: %d\n", k, result_idx);
        }

        if(k > result_idx)
        {
            arraycpy(vector_cpy, setR, R_items);
            {
                k = k - result_idx;
                result = median(vector_cpy, R_items);
                result_idx = buildSet(setL, setR, &L_items, &R_items, vector_cpy, R_items, result);
                //printf("Target idx: %d, result idx: %d\n", k, result_idx);
            }
        }
    }

    result = getSmallest(setR, R_items);
    free(vector_cpy);
    free(setL);
    free(setR);
    
    return result;
}


/**
 * Finds the index of the closer and lower element to result on the vector
 * @param   vector  Vector to find the element
 * @param   n_items Number of items of the array
 * @return  Index of the lower neightbor
 */
long getLowerNeighborIdx(double* vector, long n_items, double result)
{
    long neighbor_idx;
    double min_diff = -1;
    double diff;
    int flag = 0;

    for (int i = 0; i < n_items; i++)
    {
        if (vector[i] > result)
        {
            continue;
        } 
        else if ((vector[i] < result))
        {
            diff = result-vector[i];
            if ( diff < min_diff || min_diff == -1 )
            {
                min_diff = diff;
                neighbor_idx = i;
            }
        }
        else if (vector[i] == result)
        {
            if (flag == 0)
            {
                flag = 1;
            }
            else if (flag == 1)
            {
                neighbor_idx = i;
                break;
            }
        }
    }

    return neighbor_idx;
}


/**
 * Reorders the vectors projections and *pts. The left set has all the points with values of the projections
 * lower than the median and the right set has all the points higher than the median. The *pts array is swapped
 * acordingly to the projections array
 * @param   projections Reference array to reorder the arrays
 * @param   pts         Array to be swapped acordingly
 * @param   median      Median to use to separate the projections vector
 * @param   n_points    Number of elements on the array
 */
void compare_with_median (double* projections, double** pts, double median, long n_points)
{
    long left_pivot = 0;
    long right_pivot = n_points -1;
    long points_remaining;
    double aux;
    double* aux_pt;

    while (left_pivot < right_pivot) {
       if (projections[left_pivot] >= median) {
           if (projections[right_pivot] < median)
           {
                aux = projections[right_pivot];
                projections[right_pivot] = projections[left_pivot];
                projections[left_pivot] = aux;

                aux_pt = pts[right_pivot];
                pts[right_pivot] = pts[left_pivot];
                pts[left_pivot] = aux_pt;
            } else {
               right_pivot--;
           }
       } else {
           left_pivot++;
       }
    }

    points_remaining = n_points/2 - left_pivot;
    if (points_remaining == 0)
        return;

    for(int i=left_pivot; i<n_points; i++)
    {
        if(projections[i] == median)
        {
            aux = projections[left_pivot];
            projections[left_pivot] = projections[i];
            projections[i] = aux;

            aux_pt = pts[right_pivot];
            pts[right_pivot] = pts[left_pivot];
            pts[left_pivot] = aux_pt;
            
            left_pivot++;
            points_remaining--;
            if (points_remaining == 0)
                break;
        }
    }
}

/**
 * Searchs one value in an array and returns its idx
 * @param   projections    Array to search in
 * @param   n_points   number of points in the array 
 * @param   value   value to find
 * @return  Returns the index of value from the projection vector
 */
long find_idx_from_value(double *projections, long n_points, double value)
{
    for(int i=0; i<n_points; i++)
    {
        if(projections[i] == value)
            return i;
    }
    return -1;
}