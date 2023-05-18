
CC = g++
ifeq ($(MPI), true)
    CC = mpigCC
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

all: lib hpc_demo0
	- ./hpc_demo0


lib:
	( cd ../Lib ; $(MAKE) )

ifneq ($(MPI), true)
    hpc_demo0: lib test_create_mesh.cpp Makefile
	    $(CC) $(CF) $(I) -o hpc_demo0 test_create_mesh.cpp $(CS)
endif

ifeq ($(MPI), true)
    hpc_demo0: lib test_create_mesh.cpp Makefile
	    $(CC) $(CF) $(I) -o -D _MPI hpc_demo0 test_create_mesh.cpp $(CS)
endif


clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo0 \
    *.a *.dSYM *.obj *.dll
