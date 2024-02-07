#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <mpi.h>
#include <string.h>
#include "unsorted_median.h"
#include "gen_points.h"

//#define DEBUG_Lv1 //Without print of array (suitable for higher number of points)
//#define DEBUG_Lv2 //WARNING!! This level prints arrays (Use on smaller number of points)

#define IGNORE_NUM -1010101010
#define SAMPLE 3

int rank; //id of the current compute
int n_procs; //number of porcessors/processors

/**
 * Copies n_items from the source array to the destination array
 * @param   dest    Destination array
 * @param   sourc   Source array
 * @param   n_items Number of items to copy  
 */
void arraycpy(double* dest, double* sourc, int n_items)
{
    for(int i=0; i<n_items; i++)
        dest[i] = sourc[i];
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
        return 0.5 * ((double)vector[n_items/2] + vector[n_items/2 - 1]);
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
    double *medians = (double *)malloc(n_items * sizeof(double));

    for(i=0; i<full_splits; i++)
        medians[i] = sorted_median(vector + 5*i, 5);

    if(semi_splits != 0)
        medians[i] = sorted_median(vector + 5*full_splits, semi_splits);

    n_medians = full_splits + (semi_splits!=0 ? 1 : 0);

    if(n_medians <= 5)
    {
        result = sorted_median(medians, n_medians);
        free(medians);
 
        return result;
    }

    result = median(medians, full_splits + (semi_splits!=0 ? 1 : 0));
    free(medians);

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


int vectorSum(int* vector, int n)
{
    int count=0;

    for(int i=0; i<n; i++)
        count += vector[i];

    return count;
}


/**
 * find median and exchange points between computers
 * so that only have points bigger or smaller than 
 * the median
 * @param vector receives each computer points
 * @param n_items # of items
*/
double PSRS(double* vector, long n_items)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    int n_sample = n_procs+1;
    double* sample = (double *)malloc(n_sample * sizeof(double)); //Local vector for storing the samples before sending to root
    checkMalloc(sample, 0);

    double* recvValues; // memory position where rank 0 stores the received values
    
    if(!rank)//Only rank 0 needs the memory for storing the values
    {
        recvValues = (double*)malloc(n_procs*n_sample*sizeof(double));
        checkMalloc(recvValues, 1);
    }
    
    int n_total; // Total number of points in the problem (sum of all computers n_items)

    double** partitions = (double**)create_array_pts(n_items, n_procs); // partitions to send
    double** partitions2recv = (double**)create_array_pts(2*n_items, n_procs); // partitions to receive 
    //(2*n_items should guarantee that any computer can receive the amount of points in another computer)

    memset_double(partitions2recv[0], IGNORE_NUM, 2*n_items*n_procs); // Inicializes the vector with the a value 
                                                                      // that is going to be ignored, in case the 
                                                                      // recv vector has more space than needed

    double* finalVector = (double*)malloc(sizeof(double)*5*n_items); // Composition of the partitions belonging to the curr computer
    checkMalloc(finalVector, 2);
    
    int elems_recvd=0; // # of elements in the finalVector
    
    MPI_Request* requests = (MPI_Request*)malloc(sizeof(MPI_Request)*n_procs); //To pass as parameter into a non-blocking send call
    checkMalloc(requests, 3);

    double* partitions_limits = (double*)malloc(sizeof(double)*n_procs); //values where vector of the curr computer is splitted
    checkMalloc(partitions_limits, 4);
    //int n_elem; // Number of elements of the partition being worked on 
    int* n_elem = (int*)malloc(sizeof(int)*n_procs);
    checkMalloc(n_elem, 5);
    
    int n_to_recv; // Number of elements of the partition being received
    
    int i; //iteratorTempo

    double aux_limit; // saves the smallest/biggest value in the computer

    double finalMedians[2];

    double aux1_db;
    double aux2_db;
    
    int* aux_vector = (int*)malloc(n_procs*2*sizeof(int));
    checkMalloc(aux_vector, 6);
    
    #ifdef DEBUG_Lv1
        MPI_Barrier(MPI_COMM_WORLD);
        printf("Init vector, rank: %d with %ld elems\n", rank, n_items);
        #ifdef DEBUG_Lv2
            printArray(vector, n_items);
            fflush(stdout);
        #endif
    #endif

    //Each processor finds its samples to send to root(rank = 0)
    for(i=0; i<n_sample; i++)
    {
        sample[i] = getKsmallest(vector, (n_items/n_sample)*i, n_items);
        #ifdef DEBUG_Lv1
            printf("Sample %d = %lf, rank %d\n", i, sample[i], rank);
            fflush(stdout);
        #endif
    }
    //Each processor sends to root its samples
    MPI_Gather(&(sample[0]), n_sample, MPI_DOUBLE, &(recvValues[0]), n_sample, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    #ifdef DEBUG_Lv1
        MPI_Barrier(MPI_COMM_WORLD);
        if(!rank)
        {
            printf("Samples recebidas \n");
            
            #ifdef DEBUG_Lv2
                printArray(recvValues, n_sample*n_procs);
            #endif
            fflush(stdout);
        }
    #endif

    //Root(rank 0) computes the partitions limits from the samples received
    if(!rank)
    {
        //getKsmallest find the k value from zero
        for(i=1; i<n_procs; i++)
            partitions_limits[i-1] = getKsmallest(recvValues, n_sample*i-1, n_procs*n_sample);
            //partitions_limits[i] = getKsmallest(recvValues, n_sample*i-1, n_procs*n_sample);
    }

    //root sends partitions limits to every computer
    MPI_Bcast(&(partitions_limits[0]), n_procs-1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    #ifdef DEBUG_Lv1
        MPI_Barrier(MPI_COMM_WORLD);
        printf("Limits broadcasted\n");

        if(!rank)
        {
            printf("Limites\n");
            printArray(partitions_limits, n_procs-1);
        }
        fflush(stdout);
    #endif

    //Each computer generates their partitions (according with partition_limits)
    for(i=0; i<n_procs; i++)
    {
        //first and last partition lack limits on partition limit array
        if(i==0)
        {
            aux_limit = getKsmallest(vector, 0, n_items)-1;
            n_elem[i] = getIntervalSet(partitions[i], vector, n_items, aux_limit, partitions_limits[i]);
        }
        else if (i==n_procs-1)
        {
            aux_limit = getKsmallest(vector, n_items-1, n_items)+1;
            n_elem[i] = getIntervalSet(partitions[i], vector, n_items, partitions_limits[i-1], aux_limit);
        }
        else
        {
            n_elem[i] = getIntervalSet(partitions[i], vector, n_items, partitions_limits[i-1], partitions_limits[i]);
        }

        #ifdef DEBUG_Lv1
            MPI_Barrier(MPI_COMM_WORLD);
            if (n_elem[i] == 0)
            {
                printf("Partition empty!!!!\n");
                fflush(stdout);
            }
            printf("Particao %d do rank %d-----------------\n", i, rank);
            
            #ifdef DEBUG_Lv2
                printArray(partitions[i], n_elem[i]);
            #endif
            fflush(stdout);
        #endif

        if(i==rank) //This covers the case that the processor is building its own partition
        {
            elems_recvd = appendArray(finalVector, 0, partitions[i], n_elem[i], IGNORE_NUM);
        }
        else
        {
            MPI_Isend(&(n_elem[i]), 1, MPI_INT, i, 0, MPI_COMM_WORLD, &(requests[i]));
            if(n_elem[i])
                MPI_Isend(&(partitions[i][0]), n_elem[i], MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &(requests[i]));
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    #ifdef DEBUG_Lv1
        MPI_Barrier(MPI_COMM_WORLD);

        printf("Rank %d enviou as suas partições-----------------\n", rank);
        
        #ifdef DEBUG_Lv2
            printArray(partitions[i], n_elem[i]);
        #endif
        fflush(stdout);
    #endif
    
    //wait that all computers construct and send their partitions

    for(i=0; i<n_procs; i++)
    {
        if(i==rank)
            continue;
        else
        {
            MPI_Recv(&n_to_recv, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if(n_to_recv)
            {
                MPI_Recv(&(partitions2recv[i][0]), n_to_recv , MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                #ifdef DEBUG_Lv1
                    printf("--------Rank %d recebeu do rank %d--------\n", rank, i);
                    
                    #ifdef DEBUG_Lv2
                        printArray(partitions2recv[i], n_to_recv);
                    #endif
                    fflush(stdout);
                #endif
            }

            elems_recvd = appendArray(finalVector, elems_recvd, partitions2recv[i], n_to_recv, IGNORE_NUM);
        }
    }
    
    #ifdef DEBUG_Lv1
        //MPI_Barrier(MPI_COMM_WORLD);
        printf("Vetor de %d valores finais do rank %d completo\n", elems_recvd, rank);
        #ifdef DEBUG_Lv2
            printArray(finalVector, elems_recvd);
        #endif
        fflush(stdout);
    #endif

    //update the number of points in the computer, from now on the finalVector
    //will be used to refer to the points from this computer
    //n_items = elems_recvd;
    // shares with root the number of points in the curr computer
    
    MPI_Gather(&elems_recvd, 1, MPI_INT, aux_vector, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int median_holder1; // rank of the computer holding the median
    int median_holder2; // rank of the computer holding the second val to compute median (pair number of points)
    int median_idx1;    // index of the median in the computer thet holds it
    int median_idx2;    // index of the the second val to compute median (pair number of points)

    // root finds out the median index and which computer(s) holds it
    if(!rank)
    {
        //sum the number of points in each computer
        n_total = getVectorSum(aux_vector, n_procs);
        
        #ifdef DEBUG_Lv1
            printf("N total: %d\n", n_total);
            fflush(stdout);
        #endif

        // odd n_total -> median is one of the points
        if(n_total%2 != 0)
        {
            median_idx1 = find_K_idx(aux_vector, n_total/2, n_total, &median_holder1);
            
            // initializes vector which tells where the median is
            for(i = 0; i<n_procs*2; i++)
                aux_vector[i] = -1;

            //Fill the vector with the median_idx for the computer that holds it
            for(i=0; i<n_procs; i++)
            {
                if(i == median_holder1)
                    aux_vector[2*i] = median_idx1;
                MPI_Isend(&(aux_vector[2*i]), 2, MPI_INT, i, 2, MPI_COMM_WORLD, &(requests[i]));
            }
        }
        else
        {
            median_idx1 = find_K_idx(aux_vector, n_total/2, n_total, &median_holder1);
            median_idx2 = find_K_idx(aux_vector, (n_total/2)-1, n_total, &median_holder2);

            #ifdef DEBUG_Lv1
                printf("Median1: %d\n", median_idx1);
                printf("Median2: %d\n", median_idx2);
                printf("Holder1: %d\n", median_holder1);
                printf("Holder2: %d\n", median_holder1);
                fflush(stdout);
            #endif 
            
            for(i = 0; i<n_procs*2; i++)
                aux_vector[i] = -1;
            
            for(i = 0; i<n_procs; i++)
            {
                if(i == median_holder1)
                    aux_vector[2*i] = median_idx1;

                if(i == median_holder2)
                    aux_vector[2*i+1] = median_idx2;

                MPI_Isend(&(aux_vector[2*i]), 2, MPI_INT, i, 2, MPI_COMM_WORLD, &(requests[i]));
            }
        }
    }
    
    MPI_Recv(&(aux_vector[0]), 2, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    #ifdef DEBUG_Lv1
        MPI_Barrier(MPI_COMM_WORLD);
        printf("Rank %d recebeu: %d %d\n", rank, aux_vector[0], aux_vector[1]);
        fflush(stdout);
    #endif

    for(i=0; i<2; i++)
    {
        if(aux_vector[i] != -1)//aux_limit is reused
        {
            finalMedians[i] = getKsmallest(finalVector, aux_vector[i], elems_recvd);

            //sends median value to root
            MPI_Isend(&(finalMedians[i]), 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &(requests[i]));
        }
    }

    if(!rank)
    {
        if(n_total%2 == 0)
        {
            MPI_Recv(&aux1_db, 1, MPI_DOUBLE, median_holder1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&aux2_db, 1, MPI_DOUBLE, median_holder2, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //#ifdef DEBUG
                printf("#------------------Mediana final: %lf\n", (aux1_db + aux2_db)/2);
                fflush(stdout);
            //#endif
            return (aux1_db + aux2_db)/2;
        }
        else
        {
            MPI_Recv(&aux1_db, 1, MPI_DOUBLE, median_holder1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //#ifdef DEBUG
                printf("------------------Mediana final: %lf\n", aux1_db);
                fflush(stdout);
            //#endif
            return aux1_db;
        }
    }
    
    //Todos menos o rank 0 vão fazer return de 0. O rank 0 faz return da mediana.
    //Facilmente se faz um broadcast e todos retornam a mediana
    return 0;
}

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

double getIntervalSet(double *interval, double *vector, int n_items, double comparatorL, double comparatorR)
{
    int counter=0;

    for(int i=0; i<n_items; i++)
    {
        //printArray(setL, *counterL);
        if(vector[i]>comparatorL && vector[i]<=comparatorR)
        {
            counter++;
            interval[counter-1] = vector[i];
        }
    }
    //printArray(setL, *counterL);
    //printf("Sets criados: setL %d items, setR %d items, idx da mediana: %d\n\n", *counterL, *counterR, (*counterL));
    return counter;
}

int appendArray(double* dest, int begin_idx, double* source, int n_elem, double ignore)
{

    for(int i=0; i<n_elem; i++)
    {
        /*if(source[i] == (double)ignore)
        {
            //printf("rank aqui ----- %d", rank);
            //if(rank == 1)
            //    printf("entrou\n");
            continue;
        }*/

        //dest[j+begin_idx] = source[i];
        dest[i+begin_idx] = source[i];
        //j++;
    }

    return begin_idx+n_elem;
}

void printArray(double* vector, int n_elem)
{
    for(int i=0; i<n_elem; i++)
    {
        printf("%lf ", vector[i]);
    }
    printf("\n");
}

/**
 * Sets n_elems starting from vector to to_set double value
 */
void memset_double(double* vector, double to_set, int n_elem)
{
    for(int i=0; i<n_elem; i++)
    {
        vector[i] = to_set;
    }
}

/**
 * Finds who has the median and what is the median index inside the holder's vector
 */
int find_K_idx(int *vector, int k, int n_elem, int* median_holder)
{
    int counter = -1;

    for(int i=0; i<n_elem; i++)
    {
        counter += vector[i];

        if(counter>=k)
        {
            *median_holder = i;
            return k - (counter - vector[i])-1;
        }
    }
    return -1;
}

int getVectorSum(int* vector, int n_elem)
{
    int counter = 0;

    for(int i=0; i<n_elem; i++)
    {
        counter += vector[i];
    }
    return counter;
}

double medianSort(double* vector, int n_items)
{
    qsort(vector, n_items, sizeof(double), cmpfunc);

    if(n_items % 2 != 0)
    {
        return vector[n_items/2];
    }
    else
    {
        return 0.5 * (vector[n_items/2] + vector[n_items/2 - 1]);
    }
}

void checkMalloc(void* pointer, int tag)
{
    if(pointer==NULL)
    {
        printf("Error on malloc %d\n", tag);
        fflush(stdout);
        exit(-1);
    }
}