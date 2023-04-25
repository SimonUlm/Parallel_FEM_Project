#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

struct Matrix {
    size_t m;
    size_t n;
    ptrdiff_t incRow;
    ptrdiff_t incCol;
    double *data;
};

void
dgemv(size_t m, size_t n,
      double alpha,
      const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
      const double *x, ptrdiff_t incX,
      double beta,
      double *y, ptrdiff_t incY) {
    // Scale y <- beta * y
    if (beta)
        if (beta != 1)
            for (size_t i = 0; i < m; ++i)
                y[i * incY] *= beta;
        else
            for (size_t i = 0; i < m; ++i)
                y[i * incY] = 0; //set potential NaN to 0

    // Calculate y <- y + alpha * Ax
    if (alpha == 0)
        return;
    if (incColA < incRowA) /* row major */
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                y[i * incY] += alpha * A[i * incRowA + j * incColA] * x[j * incX];
    else /* col major */
        for (size_t j = 0; j < n; ++j)
            for (size_t i = 0; i < m; ++i)
                y[i * incY] += alpha * A[i * incRowA + j * incColA] * x[j * incX];
}

void
dgemv_MPI(int rank, size_t* size_y, size_t* offset_y, size_t* size_A, size_t* offset_A,
          double alpha, double beta,
          struct Matrix A, double *x, double *y) {
    // Define and allocate local share of vector y
    double *y_share = malloc(size_y[rank] * sizeof(double));
    assert(y_share);

    // Define and allocate local share of matrix A (called A_share)
    struct Matrix A_share;
    A_share.m = size_y[rank];
    A_share.n = A.n;
    A_share.incRow = A.incRow;
    A_share.incCol = 1;
    A_share.data = malloc(A_share.m * A_share.n * sizeof(double));
    assert(A_share.data);

    // Scatter
    MPI_Scatterv(A.data, size_A, offset_A, MPI_DOUBLE,
                 A_share.data, size_A[rank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    MPI_Scatterv(y, size_y, offset_y, MPI_DOUBLE,
                 y_share, size_y[rank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    MPI_Bcast(x, A.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculate GEMV of share
    dgemv(size_y[rank], A_share.n, alpha,
          A_share.data, A_share.incRow, A_share.incCol,
          x, 1,
          beta,
          y_share, 1);

    // Gather
    MPI_Gatherv(y_share, size_y[rank], MPI_DOUBLE,
                y, size_y, offset_y, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
}

void create_slices(size_t problem_size, int number_of_processes,
                   size_t* size, size_t* offset) {

    size_t block_size = problem_size / number_of_processes;
    size_t rest = problem_size % number_of_processes;
    size_t cum_off = 0;

    for (int i = 0; i < number_of_processes; ++i) {
        if (i < rest)
            size[i] = block_size + 1;
        else
            size[i] = block_size;
        offset[i] = cum_off;
        cum_off += size[i];
    }
}

int main(int argc, char **argv) {
    // Initialise MPI
    MPI_Init(&argc, &argv);
    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Define scalar values
    double alpha = 1;
    double beta = 1;

    // Define dimensions of matrix A
    struct Matrix A;
    A.m = 10;
    A.n = 10;
    A.incRow = 10;
    A.incCol = 1;

    // Calculate share of each process, seperate for y and A
    size_t size_y[nof_processes];
    size_t offset_y[nof_processes];
    create_slices(A.m, nof_processes, size_y, offset_y);
    size_t size_A[nof_processes];
    size_t offset_A[nof_processes];
    create_slices(A.m, nof_processes, size_A, offset_A);
    for (size_t i = 0; i < nof_processes; ++i) {
        size_A[i] *= A.n;
        offset_A[i] *= A.n;
    }

    // Allocate vectors
    double *x = malloc(A.n * sizeof(double));
    double *y = malloc(A.m * sizeof(double));
    assert(x);
    assert(y);

    // Initialise with sample values (root process only)
    if (rank == 0) {
        // Allocate matrix A
        A.data = malloc(A.m * A.n * sizeof(double));
        assert(A.data);
        // Initialise x
        for (size_t i = 0; i < A.n; ++i)
            x[i] = i;
        // Initialise y
        for (size_t i = 0; i < A.m; ++i)
            y[i] = 2 * i;
        // Initialise A
        for (size_t i = 0; i < A.m; ++i)
            for (size_t j = 0; j < A.n; ++j)
                A.data[i * A.incRow + j * A.incCol] = i * A.n + j;
    }

    // Parallalised GEMV
    dgemv_MPI(rank, size_y, offset_y, size_A, offset_A,
              alpha, beta, A, x, y);

    // Print result
    printf("y = \n");
    for (size_t i = 0; i < A.m; ++i) {
        printf("%g\n", y[i]);
    }

    MPI_Finalize();
}

