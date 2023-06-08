FILE = main_test.cpp
NP = 48		    # number of processes
N_HOSTS = 8	    # must divide NP
NP_SINGLE = $( $(NP) / $(N_HOSTS))
MIN_REFINES = 1

R_LEVEL1 = $(MIN_REFINES)
R_LEVEL2 = $$(( $(MIN_REFINES) + 1))
R_LEVEL3 = $$(( $(MIN_REFINES) + 2))
R_LEVEL4 = $$(( $(MIN_REFINES) + 3))
R_LEVEL5 = $$(( $(MIN_REFINES) + 4))
R_LEVEL6 = $$(( $(MIN_REFINES) + 5))
R_LEVEL7 = $$(( $(MIN_REFINES) + 6))
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
	- ./hpc_demo0 $(MIN_REFINES)
endif

ifeq ($(MPI),true)
all: lib hpc_demo0
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE),krafft:$(NP_SINGLE),faulhaber:$(NP_SINGLE),syrlin:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL1)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE),krafft:$(NP_SINGLE),faulhaber:$(NP_SINGLE),syrlin:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL2)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE),krafft:$(NP_SINGLE),faulhaber:$(NP_SINGLE),syrlin:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL3)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE),krafft:$(NP_SINGLE),faulhaber:$(NP_SINGLE),syrlin:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL4)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE),krafft:$(NP_SINGLE),faulhaber:$(NP_SINGLE),syrlin:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL5)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE),krafft:$(NP_SINGLE),faulhaber:$(NP_SINGLE),syrlin:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL6)
	- mpirun -np $(NP) --hostfile hostfileE44 --host heim:$(NP_SINGLE),ensinger:$(NP_SINGLE),multscher:$(NP_SINGLE),unseld:$(NP_SINGLE),magirus:$(NP_SINGLE),krafft:$(NP_SINGLE),faulhaber:$(NP_SINGLE),syrlin:$(NP_SINGLE) ./hpc_demo0 $(R_LEVEL7)
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
