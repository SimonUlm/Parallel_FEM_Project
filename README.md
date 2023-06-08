# Parallel FEM Project

This repository provides an implementation of two FEM solvers to solve the Poisson's problem on a
rectangular mesh with Neumann and Dirichlet boundary conditions. The problem is solved with the
conjugate gradient method and the w-Jacobi method. Both can be either be used serial or
parallized with MPI. This project was developed during the High Performance Computing II course
at the University of Ulm.

## Table of Contents

- [Summary](#summary)
- [Usage](#usage)
- [Results](#results)
- [References](#references)

## Summary


## Usage

The project is compiled with Make. Three pre made testfiles are provided to test
the general functionality of the software. The testcase is selected with e.g. 'CASE=1' when
calling make. The three test cases are:

1. Serial CG solver for 1 - 7x refinement
2. Parallel CG solver with 48 processors for 1 - 7x refinement
3. Parallel CG solver with 24 processors for 1 - 7x refinement

For compilation call make like so:

```console
make CASE=1
```

### Usage with Parallel solver
The parallel test cases are compiled so the 24/48 are evenly distributed between the
hosts:
heim, ensinger, multscher, unseld, magirus, krafft and faulhaber

## Contributors
Dominik Aigner, Simon Deutrich, Oliver Kratz, Lukas Ramsperger
