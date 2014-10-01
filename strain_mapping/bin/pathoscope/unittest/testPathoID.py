#!/usr/bin/python
# Initial author: Solaiappan Manimaran
# Unit Test Functions to check the Pathoscope EM algorithm functions

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

import os,sys
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,parentdir) 

import PathoID
import unittest

class TestPathoscopeFunctions(unittest.TestCase):

	def setUp(self):
		self.maxIter = 10
		self.emEpsilon = 1e-7
		self.scoreCutoff = 0.01
		self.verbose = False

	## ex1.bl8 ##
	# 2 unique reads to genome1 and three reads to genomes 1,2,3 
	def test_conv_align2GRmat_bl8_1(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex1.bl8'
		aliFormat = 2
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {0: 0, 1: 0}
		self.assertEquals(expectedU, U, "Failed bl8 Example 1 Unique Reads Assertion")
		expectedNU = {2: [[1, 2, 0], [703681821, 703681821, 703681821], [0.3333, 0.3333, 0.3333]], 
					3: [[1, 2, 0], [703681821, 703681821, 703681821], [0.3333, 0.3333, 0.3333]], 
					4: [[1, 2, 0], [703681821, 703681821, 703681821], [0.3333, 0.3333, 0.3333]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed bl8 Example 1 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed bl8 Example 1 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed bl8 Example 1 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome1', 'genome3', 'genome2']
		self.assertEquals(expectedGenomes, genomes, "Failed bl8 Example 1 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_123', 'read4_123', 'read5_123']
		self.assertEquals(expectedReads, reads, "Failed bl8 Example 1 Reads Assertion")

	## ex2.bl8 ##
	# 2 reads to 1,2; 2 reads to 1,3 ; 3 reads to 1,2,3 
	def test_conv_align2GRmat_bl8_2(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex2.bl8'
		aliFormat = 2
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {}
		self.assertEquals(expectedU, U, "Failed bl8 Example 2 Unique Reads Assertion")
		expectedNU = {0: [[0, 1], [703681821, 703681821], [0.5, 0.5]], 
			1: [[0, 1], [703681821, 703681821], [0.5, 0.5]], 
			2: [[2, 1], [703681821, 703681821], [0.5, 0.5]], 
			3: [[2, 1], [703681821, 703681821], [0.5, 0.5]], 
			4: [[2, 0, 1], [703681821, 703681821, 703681821], [0.3333, 0.3333, 0.3333]], 
			5: [[2, 0, 1], [703681821, 703681821, 703681821], [0.3333, 0.3333, 0.3333]], 
			6: [[2, 0, 1], [703681821, 703681821, 703681821], [0.3333, 0.3333, 0.3333]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed bl8 Example 2 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed bl8 Example 2 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed bl8 Example 2 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome2', 'genome1', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed bl8 Example 2 Genomes Assertion")
		expectedReads = ['read1_12', 'read2_12', 'read3_13', 'read4_13', 'read5_123', 'read6_123', 'read7_123']
		self.assertEquals(expectedReads, reads, "Failed bl8 Example 2 Reads Assertion")

	## ex3.bl8 ##
	# 3 reads to 1; 2 reads to 2; 3 reads to 1,2,3 
	def test_conv_align2GRmat_bl8_3(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex3.bl8'
		aliFormat = 2
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {0: 0, 1: 0, 2: 0, 3: 1, 4: 1}
		self.assertEquals(expectedU, U, "Failed bl8 Example 3 Unique Reads Assertion")
		expectedNU = {5: [[2, 1, 0], [703681821, 703681821, 703681821], [0.3333, 0.3333, 0.3333]], 
			6: [[2, 1, 0], [703681821, 703681821, 703681821], [0.3333, 0.3333, 0.3333]], 
			7: [[2, 1, 0], [703681821, 703681821, 703681821], [0.3333, 0.3333, 0.3333]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed bl8 Example 3 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed bl8 Example 3 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed bl8 Example 3 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome1', 'genome2', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed bl8 Example 3 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_1', 'read4_3', 'read5_2', 'read6_123', 'read7_123', 'read8_123']
		self.assertEquals(expectedReads, reads, "Failed bl8 Example 3 Reads Assertion")

	## ex4.bl8 ##
	# 2 reads to 3,4; 1 read to 4; 1 read 1,2; 1 read to 1,3 
	def test_conv_align2GRmat_bl8_4(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex4.bl8'
		aliFormat = 2
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {2: 0}
		self.assertEquals(expectedU, U, "Failed bl8 Example 4 Unique Reads Assertion")
		expectedNU = {0: [[0, 1], [703681821, 703681821], [0.5, 0.5]], 
			1: [[0, 1], [703681821, 703681821], [0.5, 0.5]], 
			3: [[2, 3], [703681821, 703681821], [0.5, 0.5]], 
			4: [[1, 3], [703681821, 703681821], [0.5, 0.5]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed bl8 Example 4 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed bl8 Example 4 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed bl8 Example 4 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome4', 'genome3', 'genome2', 'genome1']
		self.assertEquals(expectedGenomes, genomes, "Failed bl8 Example 4 Genomes Assertion")
		expectedReads = ['read1_34', 'read2_34', 'read3_4', 'read4_12', 'read5_13']
		self.assertEquals(expectedReads, reads, "Failed bl8 Example 4 Reads Assertion")

	## ex1.g.sam ##
	# 2 unique reads to genome1 and three reads to genomes 1,2,3 
	def test_conv_align2GRmat_gsam_1(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex1.g.sam'
		aliFormat = 0
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {0: 0, 1: 0}
		self.assertEquals(expectedU, U, "Failed gnusam Example 1 Unique Reads Assertion")
		expectedNU = {2: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
					3: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
					4: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed gnusam Example 1 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed gnusam Example 1 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed gnusam Example 1 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome1', 'genome2', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed gnusam Example 1 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_123', 'read4_123', 'read5_123']
		self.assertEquals(expectedReads, reads, "Failed gnusam Example 1 Reads Assertion")

	## ex2.g.sam ##
	# 2 reads to 1,2; 2 reads to 1,3 ; 3 reads to 1,2,3 
	def test_conv_align2GRmat_gsam_2(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex2.g.sam'
		aliFormat = 0
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {}
		self.assertEquals(expectedU, U, "Failed gnusam Example 2 Unique Reads Assertion")
		expectedNU = {0: [[0, 1], [50, 50], [0.5, 0.5]], 1: [[0, 1], [50, 50], [0.5, 0.5]], 
			2: [[0, 2], [50, 50], [0.5, 0.5]], 3: [[0, 2], [50, 50], [0.5, 0.5]], 
			4: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
			5: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
			6: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed gnusam Example 2 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed gnusam Example 2 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed gnusam Example 2 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome1', 'genome2', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed gnusam Example 2 Genomes Assertion")
		expectedReads = ['read1_12', 'read2_12', 'read3_13', 'read4_13', 'read5_123', 'read6_123', 'read7_123']
		self.assertEquals(expectedReads, reads, "Failed gnusam Example 2 Reads Assertion")

	## ex3.g.sam ##
	# 3 reads to 1; 2 reads to 2; 3 reads to 1,2,3 
	def test_conv_align2GRmat_gsam_3(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex3.g.sam'
		aliFormat = 0
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {0: 0, 1: 0, 2: 0, 3: 1, 4: 1}
		self.assertEquals(expectedU, U, "Failed gnusam Example 3 Unique Reads Assertion")
		expectedNU = {5: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
			6: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
			7: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed gnusam Example 3 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed gnusam Example 3 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed gnusam Example 3 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome1', 'genome2', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed gnusam Example 3 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_1', 'read4_3', 'read5_2', 'read6_123', 'read7_123', 'read8_123']
		self.assertEquals(expectedReads, reads, "Failed gnusam Example 3 Reads Assertion")

	## ex4.g.sam ##
	# 2 reads to 3,4; 1 read to 4; 1 read 1,2; 1 read to 1,3 
	def test_conv_align2GRmat_gsam_4(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex4.g.sam'
		aliFormat = 0
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {2: 1}
		self.assertEquals(expectedU, U, "Failed gnusam Example 4 Unique Reads Assertion")
		expectedNU = {0: [[0, 1], [50, 50], [0.5, 0.5]], 1: [[0, 1], [50, 50], [0.5, 0.5]], 
			3: [[2, 3], [50, 50], [0.5, 0.5]], 4: [[2, 0], [50, 50], [0.5, 0.5]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed gnusam Example 4 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed gnusam Example 4 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed gnusam Example 4 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome3', 'genome4', 'genome1', 'genome2']
		self.assertEquals(expectedGenomes, genomes, "Failed gnusam Example 4 Genomes Assertion")
		expectedReads = ['read1_34', 'read2_34', 'read3_4', 'read4_12', 'read5_13']
		self.assertEquals(expectedReads, reads, "Failed gnusam Example 4 Reads Assertion")

	## ex1.sam ##
	# 2 unique reads to genome1 and three reads to genomes 1,2,3 
	def test_conv_align2GRmat_sam_1(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex1.sam'
		aliFormat = 1
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {0: 0, 1: 0}
		self.assertEquals(expectedU, U, "Failed sam Example 1 Unique Reads Assertion")
		expectedNU = {2: [[1, 0, 2], [21, 21, 21], [0.3333, 0.3333, 0.3333]], 
					3: [[1, 0, 2], [21, 21, 21], [0.3333, 0.3333, 0.3333]], 
					4: [[2, 1, 0], [21, 21, 21], [0.3333, 0.3333, 0.3333]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed sam Example 1 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed sam Example 1 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed sam Example 1 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome1', 'genome3', 'genome2']
		self.assertEquals(expectedGenomes, genomes, "Failed sam Example 1 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_123', 'read4_123', 'read5_123']
		self.assertEquals(expectedReads, reads, "Failed sam Example 1 Reads Assertion")

	## ex2.sam ##
	# 2 reads to 1,2; 2 reads to 1,3 ; 3 reads to 1,2,3 
	def test_conv_align2GRmat_sam_2(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex2.sam'
		aliFormat = 1
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {}
		self.assertEquals(expectedU, U, "Failed sam Example 2 Unique Reads Assertion")
		expectedNU = {0: [[0, 1], [21, 21], [0.5, 0.5]], 1: [[0, 1], [21, 21], [0.5, 0.5]], 
			2: [[2, 1], [21, 21], [0.5, 0.5]], 3: [[2, 1], [21, 21], [0.5, 0.5]], 
			4: [[0, 2, 1], [21, 21, 21], [0.3333, 0.3333, 0.3333]], 
			5: [[2, 1, 0], [21, 21, 21], [0.3333, 0.3333, 0.3333]], 
			6: [[0, 1, 2], [21, 21, 21], [0.3333, 0.3333, 0.3333]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed sam Example 2 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed sam Example 2 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed sam Example 2 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome2', 'genome1', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed sam Example 2 Genomes Assertion")
		expectedReads = ['read1_12', 'read2_12', 'read3_13', 'read4_13', 'read5_123', 'read6_123', 'read7_123']
		self.assertEquals(expectedReads, reads, "Failed sam Example 2 Reads Assertion")

	## ex3.sam ##
	# 3 reads to 1; 2 reads to 2; 3 reads to 1,2,3 
	def test_conv_align2GRmat_sam_3(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex3.sam'
		aliFormat = 1
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {0: 0, 1: 0, 2: 0, 3: 1, 4: 1}
		self.assertEquals(expectedU, U, "Failed sam Example 3 Unique Reads Assertion")
		expectedNU = {5: [[2, 0, 1], [21, 21, 21], [0.3333, 0.3333, 0.3333]], 
			6: [[1, 0, 2], [21, 21, 21], [0.3333, 0.3333, 0.3333]], 
			7: [[2, 1, 0], [21, 21, 21], [0.3333, 0.3333, 0.3333]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed sam Example 3 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed sam Example 3 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed sam Example 3 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome1', 'genome2', 'genome3']
		self.assertEquals(expectedGenomes, genomes, "Failed sam Example 3 Genomes Assertion")
		expectedReads = ['read1_1', 'read2_1', 'read3_1', 'read4_3', 'read5_2', 'read6_123', 'read7_123', 'read8_123']
		self.assertEquals(expectedReads, reads, "Failed sam Example 3 Reads Assertion")

	## ex4.sam ##
	# 2 reads to 3,4; 1 read to 4; 1 read 1,2; 1 read to 1,3 
	def test_conv_align2GRmat_sam_4(self):
		currentdir = os.path.dirname(os.path.abspath(__file__))
		ali_file = currentdir + os.sep + 'data' + os.sep + 'ex4.sam'
		aliFormat = 1
		(U, NU, genomes, reads) = PathoID.conv_align2GRmat(ali_file,self.scoreCutoff,aliFormat)
		print U, NU, genomes, reads
		expectedU = {2: 0}
		self.assertEquals(expectedU, U, "Failed sam Example 4 Unique Reads Assertion")
		expectedNU = {0: [[0, 1], [21, 21], [0.5, 0.5]], 1: [[1, 0], [21, 21], [0.5, 0.5]], 
			3: [[2, 3], [21, 21], [0.5, 0.5]], 4: [[2, 1], [21, 21], [0.5, 0.5]]}
		for read in expectedNU:
			self.assertEquals(expectedNU[read][0], NU[read][0], 
				"Failed sam Example 4 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(expectedNU[read][1], NU[read][1], 
				"Failed sam Example 4 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed sam Example 4 Non-Unique Read %d proportion Assertion" %read)
		expectedGenomes = ['genome4', 'genome3', 'genome1', 'genome2']
		self.assertEquals(expectedGenomes, genomes, "Failed sam Example 4 Genomes Assertion")
		expectedReads = ['read1_34', 'read2_34', 'read3_4', 'read4_12', 'read5_13']
		self.assertEquals(expectedReads, reads, "Failed sam Example 4 Reads Assertion")


	## Example 1 ##
	# 2 unique reads to genome1 and three reads to genomes 1,2,3 
	# Answer: all to genome1
	def test_pathoscope_em_1(self):
		U = {0: 0, 1: 0}
		NU = {2: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
			3: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
			4: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]]}
		### Genome hash
		genomes = {0:"genome1", 1:"genome2", 2:"genome3"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, self.maxIter, self.emEpsilon, self.verbose)
		print initPi, pi, theta, NU
		expectedInitPi = [0.6, 0.2, 0.2]
		for j in range(len(expectedInitPi)):
			self.assertAlmostEquals(initPi[j], expectedInitPi[j], 4, "Failed EM Example 1 Initial PI Assertion")
		expectedPi = [1.0, 0.0, 0.0]
		for j in range(len(expectedPi)):
			self.assertAlmostEquals(pi[j], expectedPi[j], 4, "Failed EM Example 1 PI Assertion")
		expectedTheta = [1.0, 0.0, 0.0]
		for j in range(len(expectedTheta)):
			self.assertAlmostEquals(theta[j], expectedTheta[j], 4, "Failed EM Example 1 Theta Assertion")
		expectedNU = {2: [[0, 1, 2], [33, 33, 33], [1.0, 0.0, 0.0]], 
			3: [[0, 1, 2], [33, 33, 33], [1.0, 0.0, 0.0]], 
			4: [[0, 1, 2], [33, 33, 33], [1.0, 0.0, 0.0]]}
		for read in expectedNU:
			self.assertEquals(NU[read][0], expectedNU[read][0], 
				"Failed EM Example 1 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(NU[read][1], expectedNU[read][1], 
				"Failed EM Example 1 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed EM Example 1 Non-Unique Read %d proportion Assertion" %read)

	## Example 2 ##
	# 2 reads to 1,2; 2 reads to 1,3 ; 3 reads to 1,2,3 
	# Answer: all to genome1
	def test_pathoscope_em_2(self):
		U = {}
		NU = {0: [[0, 1], [50, 50], [0.5, 0.5]], 1: [[0, 1], [50, 50], [0.5, 0.5]], 
			2: [[0, 2], [50, 50], [0.5, 0.5]], 3: [[0, 2], [50, 50], [0.5, 0.5]], 
			4: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
			5: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
			6: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]]}
		### Genome hash
		genomes = {0:"genome1", 1:"genome2", 2:"genome3"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, self.maxIter, self.emEpsilon, self.verbose)
		print initPi, pi, theta, NU
		expectedInitPi = [0.4286, 0.2857, 0.2857]
		for j in range(len(expectedInitPi)):
			self.assertAlmostEquals(initPi[j], expectedInitPi[j], 4, "Failed EM Example 2 Initial PI Assertion j=%d" %j)
		expectedPi = [1.0, 0.0, 0.0]
		for j in range(len(expectedPi)):
			self.assertAlmostEquals(pi[j], expectedPi[j], 4, "Failed EM Example 2 PI Assertion")
		expectedTheta = [1.0, 0.0, 0.0]
		for j in range(len(expectedTheta)):
			self.assertAlmostEquals(theta[j], expectedTheta[j], 4, "Failed EM Example 2 Theta Assertion")
		expectedNU = {0: [[0, 1], [50, 50], [1.0, 0.0]], 1: [[0, 1], [50, 50], [1.0, 0.0]], 
			2: [[0, 2], [50, 50], [1.0, 0.0]], 3: [[0, 2], [50, 50], [1.0, 0.0]], 
			4: [[0, 1, 2], [33, 33, 33], [1.0, 0.0, 0.0]], 5: [[0, 1, 2], [33, 33, 33], [1.0, 0.0, 0.0]], 
			6: [[0, 1, 2], [33, 33, 33], [1.0, 0.0, 0.0]]}
		for read in expectedNU:
			self.assertEquals(NU[read][0], expectedNU[read][0], 
				"Failed EM Example 2 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(NU[read][1], expectedNU[read][1], 
				"Failed EM Example 2 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 4, 
					"Failed EM Example 2 Non-Unique Read %d proportion Assertion" %read)

	## Example 3 ##
	# 3 reads to 1; 2 reads to 2; 3 reads to 1,2,3 
	# Answer: 6 to genome1; 2 to genome2
	def test_pathoscope_em_3(self):
		U = {0: 0, 1: 0, 2: 0, 3: 1, 4: 1}
		NU = {5: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
			6: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]], 
			7: [[0, 1, 2], [33, 33, 33], [0.3333, 0.3333, 0.3333]]}
		### Genome hash
		genomes = {0:"genome1", 1:"genome2", 2:"genome3"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, self.maxIter, self.emEpsilon, self.verbose)
		print initPi, pi, theta, NU
		expectedInitPi = [0.5, 0.375, 0.125]
		for j in range(len(expectedInitPi)):
			self.assertAlmostEquals(initPi[j], expectedInitPi[j], 4, "Failed EM Example 3 Initial PI Assertion j=%d" %j)
		expectedPi = [0.75, 0.25, 0.0]
		for j in range(len(expectedPi)):
			self.assertAlmostEquals(pi[j], expectedPi[j], 2, "Failed EM Example 3 PI Assertion")
		expectedTheta = [1.0, 0.0, 0.0]
		for j in range(len(expectedTheta)):
			self.assertAlmostEquals(theta[j], expectedTheta[j], 2, "Failed EM Example 3 Theta Assertion")
		expectedNU = {5: [[0, 1, 2], [33, 33, 33], [1.0, 0.0, 0.0]], 
			6: [[0, 1, 2], [33, 33, 33], [1.0, 0.0, 0.0]], 
			7: [[0, 1, 2], [33, 33, 33], [1.0, 0.0, 0.0]]}
		for read in expectedNU:
			self.assertEquals(NU[read][0], expectedNU[read][0], 
				"Failed EM Example 3 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(NU[read][1], expectedNU[read][1], 
				"Failed EM Example 3 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 2, 
					"Failed EM Example 3 Non-Unique Read %d proportion Assertion" %read)

	## Example 4 ##
	# 2 reads to 3,4; 1 read to 4; 1 read 1,2; 1 read to 1,3 
	# Answer: 3 reads to genome4; 2 reads to genome1
	def test_pathoscope_em_4(self):
		U = {2: 1}
		NU = {0: [[0, 1], [50, 50], [0.5, 0.5]], 1: [[0, 1], [50, 50], [0.5, 0.5]], 
			3: [[2, 3], [50, 50], [0.5, 0.5]], 4: [[2, 0], [50, 50], [0.5, 0.5]]}
		### Genome hash
		genomes = {0:"genome3", 1:"genome4", 2:"genome1", 3:"genome2"}
		(initPi, pi, theta, NU) = PathoID.pathoscope_em(U, NU, genomes, self.maxIter, self.emEpsilon, self.verbose)
		print initPi, pi, theta, NU
		expectedInitPi = [0.3, 0.4, 0.2, 0.1]
		for j in range(len(expectedInitPi)):
			self.assertAlmostEquals(initPi[j], expectedInitPi[j], 4, "Failed EM Example 4 Initial PI Assertion j=%d" %j)
		expectedPi = [0.57, 0.20, 0.23, 0.0]
		for j in range(len(expectedPi)):
			self.assertAlmostEquals(pi[j], expectedPi[j], 2, "Failed EM Example 4 PI Assertion")
		expectedTheta = [0.72, 0.0, 0.28, 0.0]
		for j in range(len(expectedTheta)):
			self.assertAlmostEquals(theta[j], expectedTheta[j], 2, "Failed EM Example 4 Theta Assertion")
		expectedNU = {0: [[0, 1], [50, 50], [1.0, 0.0]], 1: [[0, 1], [50, 50], [1.0, 0.0]], 
			3: [[2, 3], [50, 50], [1.0, 0.0]], 4: [[2, 0], [50, 50], [0.14, 0.86]]}
		for read in expectedNU:
			self.assertEquals(NU[read][0], expectedNU[read][0], 
				"Failed EM Example 4 Non-Unique Read %d genome mapping Assertion" %read)
			self.assertEquals(NU[read][1], expectedNU[read][1], 
				"Failed EM Example 4 Non-Unique Read %d score Assertion" %read)
			for j in range(len(expectedNU[read][2])):
				self.assertAlmostEquals(NU[read][2][j], expectedNU[read][2][j], 2, 
					"Failed EM Example 4 Non-Unique Read %d proportion Assertion" %read)


if __name__ == '__main__':
	unittest.main()

