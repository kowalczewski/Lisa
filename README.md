# Lisa

Lisa is a code to model the performance of solar cells using an extended Hovel model.

The idea is to solve the drift-diffusion equations analytically. One of the strengths of Lisa is the ability to perform batch calculations, and to do it fast. Therefore, Lisa can be very helpful in obtaining basic trends, say, as a function of the absorber thickness. The calculations are fast because we consider one-dimensional structures. If the structure is more complicated (more layers, vertical features of the structure cannot be neglected etc.), one should consider more sophisticated device simulators (like Finite-Element solvers).

## REFERENCE PAPER

The model implemented in Lisa is described and validated here:
 
Piotr Kowalczewski, Lisa Redorici, Angelo Bozzola, Lucio Claudio Andreani,
"Silicon solar cells reaching the efficiency limits: from simple to complex modelling,"
Journal of Optics 18, 054001 (2016)

If you use this code for your research, please cite the above paper.

Regarding the original Hovel model, please consult:

Hovel H J 1975 Semiconductors and Semimetals. Volume 11
Solar Cells (New York: Academic Press) 

## LICENSE

License: GNU General Public License v3.0 (see: LICENSE file)

This file is part of Lisa.

Lisa is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

If you find any errors or have any suggestions related to the code, please contact Piotr Kowalczewski.

## LIST OF FILES

- lisa.py -- solver
- parameters.py -- file with the parameters
- loop.py -- file to create a bash file to perform batch calculations
- loop.sh -- file resulting from loop.py 
- ./Data -- folder with the data files
- ./Results -- folder with the results
- plots.py -- file with functions to make plots

## REQUIRED PYTHON LIBRARIES

- numpy
- matplotlib
- scipy

## DATA FILES

The default optical functions of silicon are taken from:
Green, Sol. Energy Mat. Sol. Cells 92, 1305-10 (2008)

The default photon flux density (am15g_en.dat) is calculated based on the data taken from http://rredc.nrel.gov/solar/spectra/am1.5/
The data are recalculated to give the irradiance as a function of energy (see: parameters.py file).
It corresponds to the standard ASTM G173.
