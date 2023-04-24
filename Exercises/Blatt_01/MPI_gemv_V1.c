#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

void
dgemv(size_t m, size_t n,
        double alpha,
        const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
        const double *x, ptrdiff_t incX,
        double beta,
        double *y, ptrdiff_t incY) {
    if (beta) {
        if (beta != 1) {
            for (size_t i = 0; i < m; ++i) {
                y[i * incY] *= beta;
            }
        }
    } else {
        for (size_t i = 0; i < m; ++i) {
            y[i * incY] = 0; //set potential NaN to 0
        }
    }	
    if (alpha) {
        if (alpha == 1) {
            if (incColA < incRowA) /* row major */ {
                for (size_t i = 0; i < m; ++i) {
                    for (size_t j = 0; j < n; ++j) {
                        y[i * incY] += A[i * incRowA + j * incColA] * x[j * incX];
                    }
                }	    
            } else /* col major */ {
                for (size_t j = 0; j < n; ++j) {
                    for (size_t i = 0; i < m; ++i) {
                        y[i * incY] += A[i * incRowA + j * incColA] * x[j * incX];
                    }
                }
            }
        } else /* alpha is neither 0 nor 1*/ {
            if (incColA < incRowA) /* row major */ {
                for (size_t i = 0; i < m; ++i) {
                    for (size_t j = 0; j < n; ++j) {
                        y[i * incY] += alpha * A[i * incRowA + j * incColA] * x[j * incX];
                    }
                }	    
            } else /* col major */ {
                for (size_t j = 0; j < n; ++j) {
                    for (size_t i = 0; i < m; ++i) {
                        y[i * incY] += alpha * A[i * incRowA + j * incColA] * x[j * incX];
                    }
                }
            }
        }
    }
}

struct Matrix {
    size_t m;
    size_t n;
    ptrdiff_t incRow;
    ptrdiff_t incCol;
    double* data;
};

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    struct Matrix A;
    A.m = 10;
    A.n = 10;
    A.incRow = 10;
    A.incCol = 1;

    if(!rank) {
        A.data = malloc(A.m * A.n * sizeof(double));
        assert(A.data);
    }
    assert(A.m % nof_processes == 0);
    int share = A.m/nof_processes;

    double* x = malloc(A.n * sizeof(double));
    double* y = malloc(A.m * sizeof(double));
    double* z = malloc(share * sizeof(double));
    assert(x); assert(y); assert(z);


    struct Matrix B;
    B.m = share;
    B.n = A.n;
    B.incRow = A.incRow;
    B.incCol = 1;
    B.data = malloc(B.m * B.n * sizeof(double));
    assert(B.data);

    if (rank == 0) {
        // initializing
        for (size_t i = 0; i < A.n; ++i) {
            x[i] = i;
        }
        for (size_t i = 0; i < A.m; ++i) {
            y[i] = 2*i;
        }
        for (size_t i = 0; i < A.m; ++i) {
            for (size_t j = 0; j < A.n; ++j) {
                A.data[i*A.incRow + j*A.incCol] = i*A.n+j;
            }
        }


        MPI_Scatter(A.data, share*B.n, MPI_DOUBLE,
                B.data, share*B.n, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        MPI_Scatter(y, share, MPI_DOUBLE,
                z, share, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        MPI_Bcast(x, A.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        dgemv(share, B.n, 1,
                B.data, B.incRow, B.incCol,
                x, 1,
                1,
                z, 1);

        MPI_Gather(z, share, MPI_DOUBLE,
                y, share, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        printf("y = \n");
        for (size_t i = 0; i < A.m; ++i) {
            printf("%g\n", y[i]);
        }
    } else {
        MPI_Scatter(NULL, 0, NULL,
                B.data, share*B.n, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        MPI_Scatter(NULL, 0, NULL,
                z, share, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        MPI_Bcast(x, A.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        dgemv(share, B.n, 1,
                B.data, B.incRow, B.incCol,
                x, 1,
                1,
                z, 1);

        MPI_Gather(z, share, MPI_DOUBLE,
                NULL, 0, NULL,
                0, MPI_COMM_WORLD);
    }

    MPI_Finalize(); 
}

