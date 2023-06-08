# On the Mac, you need to first install Xcode.  Then
# type this command in the Terminal:
#       xcode-select --install
# to ensure that you have the command-line tools installed.
#
# To just compile the library, do:
#       make
# resp.  
#       make VERBOSE=1
#       make CC=gcc OPENMP=1
#       make CC=/opt/local/bin/clang-mp-3.7 OPT=1 OPENMP=1
#
# To run a demo using the library 
#       cd ../Demo ; make
#

CC = mpiCC

LIBRARY = libhpc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG -Wpedantic  # -mfpmath=sse #-march=native #
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -fPIC
I = -I../Include
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f

AR_TARGET = $(LIBRARY).a

all: install

HPC = \
     mesh_collect.o \
     mesh_create.o \
     mesh_create_rhs.o \
     mesh_create_stiffness.o \
     mesh_objects_print.o \
     mesh_objects_util.o \
     mesh_print.o \
     mesh_refine.o \
     mesh_scatter.o \
     mesh_serial_scatter.o \
     mpi_print_serial.o \
     skeleton_create.o \
     skeleton_create_local.o \
     skeleton_print.o \
     skeleton_scatter.o \
     solver_cg.o \
     solver_cg_parallel.o \
     solver_jacobi.o \
     solver_jacobi_parallel.o \
     solver_parallel_dot.o \
     util_clock.o \
     util_gematrix_from_sed.o \
     util_sedmatrix_add_val.o \
     util_matrix_blas.o \
     util_matrix_print.o \
     util_sedmatrix_zero.o \
     util_sedmatrix_add_val.o \
     util_vector_blas.o \
     skeleton_vector_converter.o
 
$(HPC): ../Include/hpc.hpp Makefile

ifneq ($(MPI), true)
    %.o: ../Source/%.cpp ../Include/hpc.hpp
	$(CC) $(CF) $(I) -D_MPI -c $<
    %.o: ../MPI/%.cpp ../Include/hpc.hpp
	$(CC) $(CF) $(I) -D_MPI -c $<
endif

static: $(AR_TARGET)

$(AR_TARGET): $(HPC)
	$(ARCHIVE)  $@ $^
	- $(RANLIB) $@

# install archive in this directory
install: $(AR_TARGET)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*

