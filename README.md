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


>BoxOut		0x00	//High: 1:PatchyPoint 2:nst Low: 1:fp Box 2:vmd Box 4:fp connect 8:vmd connect
>NCellX			80
>NCellY			40
>NCellZ			40
>ELong			1.0
>P_LJC			10
>P_LJS			10//P_LJC
>KLat			4
>KLong			8
>TLong			2e1//(30.0*TLat)
>GLat			2e0
>DeltaL			0e-2	//UniformForce: Strain Other: Distance
>DeltaA			0e-3
>BCType			0x10	//High:-end Low:+end 0:Free 1:Hinged 2:Slide 4:Fix
>ForceXYZ		0x00	//Low: Direction:xyz(124) Rotate(8) High: 0:None 1:Uniform 2:+end 4:-end 8:ShearEdge
;>RFUnits		0//(RFLayers*InitFila)//25
>FixLink		1
>cutoff_s		2.7
>Max_Add_Link	12
>QDType			0x05	//High: chain Low: side 1:distribution 2: uniform 4: uni-axis 8: uni-Layer
>TipType		2		//0:Sin 1:Mix 2:Power 3:Power Shift 4:Double Peak 5:Linear Mix 6:Decay
>PhaseFile		0		//0:Given by paras 1:Given by File 'PhaseData'	2:Given by File 'PhaseData07'
>GROW_ON		3
>MustGrow		1
>BREAK_ON	0
>HD_ON		1
>G_rate		9.0
>EBase		0.71
>GSteps		5000
>G_Slow		500
>HD_Delay	100000
>C_rate		1e-12
>HD_rate		2e-7
>S_rate		1e7
>B_rate		5e-3

!Loop
>InitFila	13
>XI_0		5.0
>RFLayers	0
>RFUnits	24
>QDDelta	0.0
>TipSinA	0
>TipSinL	3
>TipSinH		1
>TipSinPhase	0.0
>TipARatio		0.0
;ELatAA			2.0
;ELatBA			1.64
>Tst				0
>Ted				100
>DeltaT			10
>DeltaXi			0


