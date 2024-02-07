#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <mpi.h>
#include "gen_points.h"
#include "unsorted_median.h"

#define DEBUG

int main(int argc, char *argv[])
{
    double* pts;
    double* pts_test;
    long n_points = atol(argv[2]);;
    long n_points_test = atol(argv[2]);
    int rank;
    int n_procs;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    pts = get_points_1dim_mpi(argc, argv, MPI_COMM_WORLD, &n_points);


    #ifdef DEBUG
        if(!rank)
        {
            pts_test = get_points_1dim(argc, argv);
            printf("#-------------------------Mediana real: %lf\n", medianSort(pts_test, n_points_test));
            fflush(stdout);
        }
    #endif

    MPI_Barrier(MPI_COMM_WORLD);
    PSRS(pts, n_points);

    MPI_Finalize();
}