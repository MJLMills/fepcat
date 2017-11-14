# Fepcat

Fepcat is a set of programs for computing Potentials of Mean Force (PMFs) via the Free Energy Perturbation/Umbrella Sampling (FEP/US) method.
The method is of [Warshel](http://laetro.usc.edu/), and has been implemented previously in the mapping program of [Molaris](http://laetro.usc.edu/software.html) and the QFep program of [Q](http://xray.bmc.uu.se/~aqwww/q/). 
The procedural forms these programs take make them difficult to modify and extend.
The Fepcat code is written in a modular fashion in order that the user may define additional new analyses with improved flexibility.
Fepcat also provides a modular reimplementation of Qfep for purposes of comparison.

Once an EVB model has been defined and FEP simulations performed, the first question following FEP or FEP/US analysis is "Why is the answer not what I expected?"
There are many parameters and assumptions built into EVB/FEP/US free energy calculations, and it is usually not obvious where the error lies.
For example, an insufficient number of FEP steps may have been chosen, or the length of each simulation may not be sufficient for relevant mean values to have converged.
Alongside its use as a diagnostic tool, Fepcat produces comma-separated value output files which can be easily read into plotting or analysis programs.
This is in contrast to earlier tools, which tend to write formatted text output files from which data has to be parsed, and allows the fast location of errors.
The program also offers optimization of EVB parameters to match experimental target free energy changes.

The initial release version is sparingly documented and is presented in the hope that researchers with Fortran skills can find something useful within.
Future refinement and a complete description of the code is intended - time permitting - to follow this release.

#### Compiling/Building the Program

A makefile for GNU make is provided, which assumes use of the gfortran compiler and features of Fortran95 (it has been tested on GNU Fortran (GCC) 4.9.2).
No other compilers have been tested.
To compile, issue 'make all' in the root directory.
This will produce the following executables:

*Fepcat* - The main program, perform FEP and FEP/US analyses of simulation trajectories including Flyvbjerg/Petersen statistical tests.

*Qfep* - Reproduce the results and output of a Qfep calculation. Included for sanity checking purposes.

*FepMovie* - Perform a FEP/US calculation, writing out the necessary input files to produce a movie of the data accumulation process (molecules and plots).

*Fep2D* - Perform a FEP/US calculation with a 2-dimensional reaction coordinate in place of the typical 1-dimensional energy gap.

*Fep* - Perform a FEP analysis only (no umbrella sampling corrections).

*AveGeom* - Compute average geometries over simulation trajectories.

#### Citing the Program

If you use this code in whole or part please cite:

Fepcat v1.0.0, M J L Mills, 2017, github.com/MJLMills/fepcat

#### Acknowledgements

This project was initially developed while the author was an employee of Sandia National Laboratories, in the Enzyme Optimization group at the DOE's [Joint BioEnergy Institute](https://www.jbei.org/) under the direction of Dr. Kenneth Sale.

#### License

This project is licensed under a modified BSD license - see LICENSE.txt for details.

#### Copyright Message

Fepcat Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.
 
If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation and Partnerships department at IPO@lbl.gov referring to " Fepcat (2017-166)."
 
NOTICE.  This software was developed under funding from the U.S. Department of Energy.  As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, and perform publicly and display publicly. The U.S. Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
