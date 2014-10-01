#!/usr/bin/python
# Initial author: Solaiappan Manimaran
# General utility functions

#	Pathoscope - Predicts strains of genomes in Nextgen seq alignment file (sam/bl8)
#	Copyright (C) 2013  Johnson Lab - Boston University
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os

# ===========================================================
def file_len(fname):
	with open(fname) as f:
		for i, _ in enumerate(f):
			pass
	return i + 1
# ===========================================================
def ensure_dir(d):
	if not os.path.exists(d):
		os.makedirs(d)
# ===========================================================
