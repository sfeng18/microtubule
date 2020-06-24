# microtubule

This repository holds the simulation source code described in the paper [eLife 2020;9:e54077](https://doi.org/10.7554/eLife.54077). The code in this repository requires C compiler and MPICH to build and run. Please cite this work appropriately (Michaels, Thomas CT, Shuo Feng, Haiyi Liang, and L. Mahadevan. "Mechanics and kinetics of dynamic instability." Elife 9 (2020): e54077).

## Build and Run

Simply use `mpicc` and `mpirun`:

```bash
$ mpicc -o microtubule.out microtubule.c -lm -O3
$ mpirun -n 4 ./microtubule.out
```

## Description

This program builds a 3D structure of micritubule and does the molecular dynamics or monte-carlo (MD/MC) simulation. Macro


