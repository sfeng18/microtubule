# microtubule

This repository holds the simulation source code described in the paper [eLife 2020;9:e54077](https://doi.org/10.7554/eLife.54077). The code in this repository requires C compiler and MPICH to build and run. Please cite this work appropriately (Michaels, Thomas CT, Shuo Feng, Haiyi Liang, and L. Mahadevan. "Mechanics and kinetics of dynamic instability." Elife 9 (2020): e54077).

## Build and Run

Simply use `mpicc` and `mpirun`:

```bash
$ mpicc -o microtubule.out microtubule.c -lm -O3
$ mpirun -n 4 ./microtubule.out
```

## Description

This program builds a 3D structure of micritubule (MT) and does the molecular dynamics or monte-carlo (MD/MC) simulation.

All parameters of MT (including the curvature of tubulin, the distribution of GTP- and GDP- tubulins, the shape of initial microtubule, the growth rate, the hydrolyzation constant and) are controlled by macro definitions.

### Method

* `MD_ON` :         1 for MD simulation, 0 for MC simulation
* `Times` :         Output one frame after how many time steps
* `MTimes` :        Total simulation time (how many output steps)
* `NCellX` :        Box size along the axis `x`; Length of each cell is `CellL = 7`
* `NCellY` :        Box size along the axis `y`
* `NCellZ` :        Box size along the axis `z`

### Inital Structure

* `Whole_Struc` :   Geometric structure, 0 for sheet-like, 1 for tube-like 
* `My_Type` :       Initial type of tubulins, 0 for GTP, 1 for GDP
* `HELIX_ON` :      0 for cylindrical MT, 1 for helical MT
* `NInit` :         How many tubulin layers on the longitudinal direction
* `Lat_Length` :    How many tubulins on the lateral direction, works only when `Whole_Struc = 0`
* `InitFila` :      How many tubulins on the lateral direction, works only when `Whole_Struc = 1`

### Potential Function

The length unit in this code is 2nm, which means that the $r_0$ = 4nm in the paper is equivalent to $r_0$ = 2 here.
Similarly, the energy unit is $\epsilon_0$, same as that in the paper.

* `P_LJC` :         Parameter for longitudinal stretching stiffness, $a_long$
* `P_LJS` :         Parameter for lateral stretching stiffness, $a_lat$
* `KLong` :         Parameter for longitudinal bending stiffness, $b_long$
* `KLat` :          Parameter for lateral bending stiffness, $b_long$
* `TLong` :         Parameter for longitudinal twisting stiffness, $c_long$
* `GLat` :          Parameter for lateral twisting stiffness, $c_long$

### Dynamics

A layer of virtual tubulins is added on the MT at each growth step, the growth probability of each virtual tubulin is determined by its energy. At most, one virtual tubulin is chosen to grow.

* `GROW_ON` :       Let the MT grow or not, 3 for grow, 0 for not grow
* `HD_ON` :         Let tubulin hydrolyze or not
* `BREAK_ON` :      Let filament break or not
* `GSteps` :        Growth rate, grow one tubulin every `GSteps` steps
* `G_rate` :        Parameter for calculate the growth probability of each filament
* `EBase` :         Parameter for calculate the growth probability of each filament
* `MustGrow` :      Force to grow one tubulin at each growth step
* `HD_rate` :       Hydrolyzation constant
* `C_rate` :        Parameter for calculate the breaking probability of each filament

### Distribution of GTP-tubulins

* `RFLayers` :      How many layers of tubulin are reinforced
* `RFUnits` :       How many GTP-tubulin are reinforced (They will not hydrolyze); This value should be smaller than `InitFila * RFLayers`
* `QDType` :        Type of quenched disorder. 8-bit value, the higher 4-bit for longitudinal springs and lower 4-bit for lateral ones. 1 means distributing in this direction (see Eq. 14 in the paper), 2 means uniform distribution, 4 means same value for each filament, 8 means same value for each layer
* `QDDelta` :       1/k

### Other Parameters

* Parameters for the shape of MT tip: `TipType` `TipSinA` `TipSinL` `TipSinH` `TipSinPhase` `TipARatio`
* Parameters for the mechanical tests: `BCType` `ForceXYZ` `DeltaL` `DeltaA` `DeltaT` `DeltaXi` `Tst` `Ted`



