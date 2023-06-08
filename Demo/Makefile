CASE=4

ifeq ($(CASE),1)
    FILE=testcase1.cpp
    MPI=true
    MIN_REFINES=1
endif

ifeq ($(CASE),2)
    FILE=testcase2.cpp
    MPI=true
    MIN_REFINES=1
    NP=6
endif

ifeq ($(CASE),3)
    FILE=testcase3.cpp
    MPI=true
    MIN_REFINES=1
    NP=24
    N_HOSTS = 4	    # must divide NP
    NP_SINGLE = $( $(NP) / $(N_HOSTS))
endif

ifeq ($(CASE),4)
    FILE=main_parallel.cpp
    MPI=true
    MIN_REFINES=1
    NP=6
endif

R_LEVEL1 = $(MIN_REFINES)
R_LEVEL2 = $$(( $(MIN_REFINES) + 1))
R_LEVEL3 = $$(( $(MIN_REFINES) + 2))
R_LEVEL4 = $$(( $(MIN_REFINES) + 3))
R_LEVEL5 = $$(( $(MIN_REFINES) + 4))
R_LEVEL6 = $$(( $(MIN_REFINES) + 5))
R_LEVEL7 = $$(( $(MIN_REFINES) + 6))
R_LEVEL8 = $$(( $(MIN_REFINES) + 7))
R_LEVEL9 = $$(( $(MIN_REFINES) + 8))

SHELL=/bin/bash

CC = mpiCC

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

ifeq ($(CASE),1)
all: lib hpc_demo0
	- ./hpc_demo0 $(R_LEVEL1)
	- ./hpc_demo0 $(R_LEVEL2)
	- ./hpc_demo0 $(R_LEVEL3)
	- ./hpc_demo0 $(R_LEVEL4)
	- ./hpc_demo0 $(R_LEVEL5)
	- ./hpc_demo0 $(R_LEVEL6)
	- ./hpc_demo0 $(R_LEVEL7)
endif

ifeq ($(CASE),2)
all: lib hpc_demo0
	- mpirun -np $(NP) ./hpc_demo0 $(R_LEVEL1)
	- mpirun -np $(NP) ./hpc_demo0 $(R_LEVEL2)
	- mpirun -np $(NP) ./hpc_demo0 $(R_LEVEL3)
	- mpirun -np $(NP) ./hpc_demo0 $(R_LEVEL4)
	- mpirun -np $(NP) ./hpc_demo0 $(R_LEVEL5)
	- mpirun -np $(NP) ./hpc_demo0 $(R_LEVEL6)
	- mpirun -np $(NP) ./hpc_demo0 $(R_LEVEL7)
	- mpirun -np $(NP) ./hpc_demo0 $(R_LEVEL8)
	- mpirun -np $(NP) ./hpc_demo0 $(R_LEVEL9)
endif


ifeq ($(CASE),3)
all: lib hpc_demo0
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL1)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL2)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL3)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL4)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL5)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL6)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL7)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL8)
endif

ifeq ($(CASE),4)
all: lib hpc_demo0
	- mpirun -np $(NP) ./hpc_demo0
endif

lib:
	( cd ../Lib ; $(MAKE) )

ifneq ($(MPI), true)
    hpc_demo0: lib $(FILE) Makefile
	    $(CC) $(CF) $(I) -o hpc_demo0 $(FILE) $(CS)
endif

ifeq ($(MPI), true)
    hpc_demo0: lib $(FILE) Makefile
	    $(CC) $(CF) $(I) -D_MPI -DMIN_REFINES=$(MIN_REFINES) -o hpc_demo0 $(FILE) $(CS)
endif


clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo0 \
    *.a *.dSYM *.obj *.dll
