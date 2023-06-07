FILE = main_test.cpp
NP = 12
SHELL=/bin/bash

CC = g++
ifeq ($(MPI), true)
    CC = mpiCC
endif

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

ifneq ($(MPI),true)
all: lib hpc_demo0
	- ./hpc_demo0 5
	- ./hpc_demo0 6
	- ./hpc_demo0 7
	- ./hpc_demo0 8
endif

ifeq ($(MPI),true)
all: lib hpc_demo0
	- mpirun -np $(NP) ./hpc_demo0 5
	- mpirun -np $(NP) ./hpc_demo0 6
	- mpirun -np $(NP) ./hpc_demo0 7
	- mpirun -np $(NP) ./hpc_demo0 8
endif
lib:
	( cd ../Lib ; $(MAKE) )

ifneq ($(MPI), true)
    hpc_demo0: lib $(FILE) Makefile
	    $(CC) $(CF) $(I) -o hpc_demo0 $(FILE) $(CS)
endif

ifeq ($(MPI), true)
    hpc_demo0: lib $(FILE) Makefile
	    $(CC) $(CF) $(I) -D_MPI -o hpc_demo0 $(FILE) $(CS)
endif


clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo0 \
    *.a *.dSYM *.obj *.dll
