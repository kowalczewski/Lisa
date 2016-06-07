'''
This file includes functions to prepare plots.

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

import numpy as np
import matplotlib.pyplot as plt
import math
import sys

if __name__ == "__main__":
	
	"""
        Usage: python plots.py [type of the plot]
        Possible types: jv, panel, panel_comparison
	"""

	def myround(x, base=5.0):
		return int(base * math.ceil(x/base))

	def plot_JV():

		jv_filename = sys.argv[2]

		# read JV characteristic
		v_vec = []
		j_vec = []	
		with open("./Results/"+jv_filename) as jv_f:
			# read Jsc and Voc
			voc = float(jv_f.readline())
			jsc = float(jv_f.readline())
			for line in jv_f:
				line = line.split()
				v_vec.append(float(line[0]))
				j_vec.append(float(line[1]))

		plt.plot(v_vec, j_vec, linewidth=2)
		plt.ylim(0,myround(jsc))
		plt.xlabel('Voltage (V)')
		plt.ylabel('Current density (mA/cm2)')

		plt.show()
		
	def plot_panel():
		th_vec = []
		voc_vec = []
		jsc_vec = []
		ff_vec = []
		eff_vec = []
		
		panel_filename = sys.argv[2]
		
		with open("./Results/"+panel_filename) as panel_f:
			for line in panel_f:
				line = line.split()
				th_vec.append(float(line[0]))
				voc_vec.append(float(line[2]))	
				jsc_vec.append(float(line[3]))	
				ff_vec.append(float(line[4]))	
				eff_vec.append(float(line[5]))	
		
		plt.figure(1)
		
		# Voc
		# arguments: no of rows, number of columns, fig number
		plt.subplot(221)
		plt.xscale('log')
		plt.plot(th_vec, voc_vec, linewidth=2)
		plt.xlim(min(th_vec),max(th_vec))
		plt.ylabel('Open-circuit voltage (V)')
		plt.xlabel('Thickness (um)')
		
		# Jsc
		plt.subplot(222)
		plt.xscale('log')
		plt.plot(th_vec, jsc_vec, linewidth=2)
		plt.xlim(min(th_vec),max(th_vec))
		plt.ylabel('Short-circuit current (mA/cm2)')
		plt.xlabel('Thickness (um)')
		
		# FF
		plt.subplot(223)
		plt.xscale('log')
		plt.plot(th_vec, ff_vec, linewidth=2)
		plt.xlim(min(th_vec),max(th_vec))
		plt.ylim(0.86,0.92)
		plt.ylabel('Fill Factor')
		plt.xlabel('Thickness (um)')
		
		# Jsc
		plt.subplot(224)
		plt.xscale('log')
		plt.plot(th_vec, eff_vec, linewidth=2)
		plt.xlim(min(th_vec),max(th_vec))
		plt.ylabel('Efficiency (%)')
		plt.xlabel('Thickness (um)')
		
		plt.show()
		
	def plot_panel_comparison():
		th_vec = []
		voc_vec = []
		jsc_vec = []
		ff_vec = []
		eff_vec = []
		
		panel_filename_1 = sys.argv[2]
		panel_filename_2 = sys.argv[3]
		
		with open("./Results/"+panel_filename_1) as jv_f:
			for line in jv_f:
				line = line.split()
				th_vec.append(float(line[0]))
				voc_vec.append(float(line[2]))	
				jsc_vec.append(float(line[3]))	
				ff_vec.append(float(line[4]))	
				eff_vec.append(float(line[5]))	
		
		# results from the diode equation
		
		thD_vec = []
		vocD_vec = []
		jscD_vec = []
		ffD_vec = []
		effD_vec = []
		
		with open("./Results/"+panel_filename_2) as jv_f:
			for line in jv_f:
				line = line.split()
				thD_vec.append(float(line[0]))
				vocD_vec.append(float(line[2]))	
				jscD_vec.append(float(line[3]))	
				ffD_vec.append(float(line[4]))	
				effD_vec.append(float(line[5]))	
		
		plt.figure(2)
		
		# Voc
		# arguments: no of rows, number of columns, fig number
		plt.subplot(221)
		plt.xscale('log')
		plt.plot(th_vec, voc_vec, linewidth=2, label=panel_filename_1)
		plt.plot(thD_vec, vocD_vec, linewidth=2, label=panel_filename_2, linestyle='--', color='red')
		plt.xlim(min(th_vec),max(th_vec))
		plt.ylabel('Open-circuit voltage (V)')
		plt.xlabel('Thickness (um)')
		plt.legend(loc='best') 
		
		# Jsc
		plt.subplot(222)
		plt.xscale('log')
		plt.plot(th_vec, jsc_vec, linewidth=2, label=panel_filename_1)
		plt.plot(thD_vec, jscD_vec, linewidth=2, label=panel_filename_2, linestyle='--', color='red')
		plt.xlim(min(th_vec),max(th_vec))
		plt.ylabel('Short-circuit current (mA/cm2)')
		plt.xlabel('Thickness (um)')
		plt.legend(loc='best') 
		
		# FF
		plt.subplot(223)
		plt.xscale('log')
		plt.plot(th_vec, ff_vec, linewidth=2, label=panel_filename_1)
		plt.plot(thD_vec, ffD_vec, linewidth=2, label=panel_filename_2, linestyle='--', color='red')
		plt.xlim(min(th_vec),max(th_vec))
		plt.ylim(0.86,0.92)
		plt.ylabel('Fill Factor')
		plt.xlabel('Thickness (um)')
		plt.legend(loc='best') 
		
		# Jsc
		plt.subplot(224)
		plt.xscale('log')
		plt.plot(th_vec, eff_vec, linewidth=2, label=panel_filename_1)
		plt.plot(thD_vec, effD_vec, linewidth=2, label=panel_filename_2, linestyle='--', color='red')
		plt.xlim(min(th_vec),max(th_vec))
		plt.ylabel('Efficiency (%)')
		plt.xlabel('Thickness (um)')
		plt.legend(loc='best') 
		
		plt.show()

	# ============================= PLOTS =============================

	if (sys.argv[1] == "jv"): plot_JV()
	elif (sys.argv[1] == "panel"): plot_panel() 
	elif (sys.argv[1] == "panel_comparison"): plot_panel_comparison() 
	