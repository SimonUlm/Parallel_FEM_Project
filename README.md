# Parallel FEM Project

This repository provides an implementation of two FEM solvers to solve the Poisson's problem on a
rectangular mesh with Neumann and Dirichlet boundary conditions. The problem is solved with the
conjugate gradient method and the \omega-Jacobi method. Both can be either be used serial or
parallized with MPI. This project was developed during the High Performance Computing II course
at the University of Ulm.

## Table of Contents

- [Summary](#summary)
- [Usage](#usage)
- [Results](#results)
- [References](#references)

## Summary


## Usage

The project is compiled with Make. The user can specify the the file with `FILE`, number of processes 
with `NP` and whether it should be compiled with MPI with `MPI=true`. For testing the following files
provided:

### Usage with Serial solver
To use it with serial solvers select a test file like so which directly executes the program (missing which refinements)

```
make FILE=test_cg_serial.cpp
```



### Usage with Parallel solver
To use the software with parallized algorithms the number of processes and 


## Results



## References




## Contributors
Dominik Aigner, Simon Deutrich, Oliver Kratz, Lukas Ramsperger
