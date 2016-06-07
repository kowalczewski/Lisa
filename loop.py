'''
This script prepares data for batch calculations.

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

from numpy import *

if __name__ == "__main__":
	
	"""
        Usage: python loop.py
    """

	th_vec = logspace(0,2.7,100)
	# srv_vec = logspace(-1,2,50)
	# for no SR
	srv_vec = [0]

	sh_f = open('loop.sh','w')

	for th in th_vec:
		for srv in srv_vec:
			print>>sh_f, 'python lisa.py ', th, srv
		
	sh_f.close()
