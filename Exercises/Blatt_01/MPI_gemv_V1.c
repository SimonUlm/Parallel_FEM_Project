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
dgemv_MPI(size_t share, double alpha, double beta,
          struct Matrix A, double *x, double *y) {
    // Define and allocate local share of vector y
    double *y_share = malloc(share * sizeof(double));
    assert(y_share);

    // Define and allocate local share of matrix A (called A_share)
    struct Matrix A_share;
    A_share.m = share;
    A_share.n = A.n;
    A_share.incRow = A.incRow;
    A_share.incCol = 1;
    A_share.data = malloc(A_share.m * A_share.n * sizeof(double));
    assert(A_share.data);

    // Scatter
    MPI_Scatter(A.data, share * A_share.n, MPI_DOUBLE,
                A_share.data, share * A_share.n, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Scatter(y, share, MPI_DOUBLE,
                y_share, share, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Bcast(x, A.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculate GEMV of share
    dgemv(share, A_share.n, alpha,
          A_share.data, A_share.incRow, A_share.incCol,
          x, 1,
          beta,
          y_share, 1);

    // Gather
    MPI_Gather(y_share, share, MPI_DOUBLE,
               y, share, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
}

int main(int argc, char **argv) {
    // Initialise MPI
    MPI_Init(&argc, &argv);
    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Define scalar values
    double alpha = 3;
    double beta = 5;

    // Define dimensions of matrix A
    struct Matrix A;
    A.m = 10;
    A.n = 10;
    A.incRow = 10;
    A.incCol = 1;

    // Calculate share of each process
    assert(A.m % nof_processes == 0);
    size_t share = A.m / nof_processes;

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
    dgemv_MPI(share, alpha, beta, A, x, y);

    // Print result
    printf("y = \n");
    for (size_t i = 0; i < A.m; ++i) {
        printf("%g\n", y[i]);
    }

    MPI_Finalize();
}

