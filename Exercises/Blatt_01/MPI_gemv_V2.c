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

void slices(size_t problem_size, int number_of_processes, 
        int* size, int* offset) {

    size_t block_size = problem_size / number_of_processes;
    size_t rest = problem_size % number_of_processes;
    size_t cum_off = 0;

    for (int i = 0; i < number_of_processes; ++i) {
        if(i < rest) {
            size[i] = block_size + 1;
        }else {
            size[i] = block_size;
        }
        offset[i] = cum_off;
        cum_off += size[i];
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


    int sizeA[nof_processes];
    int offsetA[nof_processes];
    slices(A.m, nof_processes, sizeA, offsetA);

    int sizey[nof_processes];
    int offsety[nof_processes];
    slices(A.m, nof_processes, sizey, offsety);

    double* x = malloc(A.n * sizeof(double));
    double* y = malloc(A.m * sizeof(double));
    double* z = malloc(sizey[rank] * sizeof(double));
    assert(x); assert(y); assert(z);

    struct Matrix B;
    B.m = sizey[rank];
    B.n = A.n;
    B.incRow = A.incRow;
    B.incCol = 1;
    B.data = malloc(B.m * B.n * sizeof(double));
    assert(B.data);

    for (size_t i = 0; i < nof_processes; ++i) {
        sizeA[i] *= A.n;
        offsetA[i] *= A.n;
    }


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


        MPI_Scatterv(A.data, sizeA, offsetA, MPI_DOUBLE,
                B.data, sizeA[rank], MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        MPI_Scatterv(y, sizey, offsety, MPI_DOUBLE,
                z, sizey[rank], MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        MPI_Bcast(x, A.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        dgemv(sizey[rank], B.n, 1,
                B.data, B.incRow, B.incCol,
                x, 1,
                1,
                z, 1);

        MPI_Gatherv(z, sizey[rank], MPI_DOUBLE,
                y, sizey, offsety, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        printf("y = \n");
        for (size_t i = 0; i < A.m; ++i) {
            printf("%g\n", y[i]);
        }
    } else {
        MPI_Scatterv(NULL, NULL, NULL, NULL,
                B.data, sizeA[rank], MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        MPI_Scatterv(NULL, NULL, NULL, NULL,
                z, sizey[rank], MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        MPI_Bcast(x, A.n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        dgemv(sizey[rank], B.n, 1,
                B.data, B.incRow, B.incCol,
                x, 1,
                1,
                z, 1);

        MPI_Gatherv(z, sizey[rank], MPI_DOUBLE,
                NULL, NULL, NULL, NULL,
                0, MPI_COMM_WORLD);
    }

    MPI_Finalize(); 
}

