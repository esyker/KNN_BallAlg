#ifndef _GEN_POINTS_H_
#define _GEN_POINTS_H_
double **get_points_mpi(int argc, char *argv[], MPI_Comm, long* n_local_points);
double **create_array_pts(int n_dims, long np);
double* get_points_1dim_mpi(int argc, char *argv[], MPI_Comm comm, long* n_local_points);
#endif
