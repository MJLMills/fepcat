# Fepcat

Fepcat is a set of programs for computating Potentials of Mean Force (PMF) via the Free Energy Perturbation/Umbrella Sampling (FEP/US) method.
The method is of Warshel, and has been implemented previously in the mapping program of Molaris and the QFep program of Q. 
Their procedural forms make them difficult to modify and extend.
The Fepcat code is written in a modular fashion in order that the user may define additional new analyses with great flexibility.
Fepcat also provides a reimplementation of Qfep for purposes of comparison.

Free energy calculations reliant on the ergodic hypothesis involve computation of mean values of fluctuating quantities along a temporal trajectory.

Once an EVB model has been defined and FEP simulations performed, the first question following FEP or FEP/US analysis is "Why is the answer not what I expected?"
There are many parameters and assumptions built into EVB/FEP/US free energy calculations, and it is usually not obvious where the error lies.
For example, an insufficient number of FEP steps may have been chosen, or the length of each simulation may not be sufficient for relevant mean values to have converged.
Alongside its use as a diagnostic tool, Fepcat produces comma-separated value output files which can be easily read into plotting or analysis programs.
This is in contrast to earlier tools, which tend to write formatted text output files from which data has to be parsed.


Fepcat Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.
 
If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation and Partnerships department at IPO@lbl.gov referring to " Fepcat (2017-166)."
 
NOTICE.  This software was developed under funding from the U.S. Department of Energy.  As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, and perform publicly and display publicly. The U.S. Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
