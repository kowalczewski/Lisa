'''
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

'''

# ====================================== DATA FILES ====================================
# The files should be in Data/ directory (see example files).
# Solar photon flux
# Format: Energy (eV) AM1.5G (Wm-2eV-1)
# The default photon flux density (am15g_en.dat) is calculated based on the data taken from http://rredc.nrel.gov/solar/spectra/am1.5/
# It corresponds to the standard ASTM G173.
# The data are recalculated to give the irradiance as a function of energy (see the units above).
solar_f = "am15g_en.dat"
# Optical functions of the active layer
# Format: Energy (eV) Wavelength (nm) eps1	eps2	n	k	Absorption coefficient (cm-1)
# The default optical functions of silicon are taken from:
# Green, Sol. Energy Mat. Sol. Cells 92, 1305-10 (2008)
si_f = "cSi_Green.dat"
# ====================================== OUTPUT PARAMETERS =====================================
# If you want to do a loop over parameters, everything should be set to "false"
# Print parameters of the structures
params_flag = True
# Calculate the full JV characteristic and save it to the file (can take much longer than standard simulations).
jv_flag = False
# Show full results
results_flag = True

# ====================================== PARAMETERS ======================================
# temperature (K)
T = 300
# intrinsic carrier concentration (cm-3)
ni_0 = 9.6541E9                                  

# ====================================== STRUCTURE ======================================
# structure dimensions
# in cm, to be consisten with calculations
# thickness (cm)
# this value can be also given in the command line (as the 1st argument)
# if this value is given in the command line, the value below is overwritten
th = 1*1E-4
# emitter thickness (cm)
th_emitter = 0.005E-4

# ====================================== BASE ======================================
# base doping (cm-3)
# n-type or p-type doping, the other should be 0
N_a_B = 1E16
N_d_B = 0
# diffusion coefficient (cm2/s)
D_B = 25
# surface recombination velocities (cm/s)
# this value can be also given in the command line (as the 2nd argument)
# if this value is given in the command line, the value below is overwritten
S_B = 0

# ====================================== EMITTER ======================================
# emitter doping (cm-3)
# n-type or p-type doping, the other should be 0
N_a_E = 0
N_d_E = 1.5E18
# diffusion coefficient (cm2/s)
D_E = 25/2	
# surface recombination velocities (cm/s)
S_E = 0

# ====================================== PARAMETERS OF THE ENERGY VECTOR ======================================
en_start = 1.0
en_stop = 4.4
en_points = 1E3
