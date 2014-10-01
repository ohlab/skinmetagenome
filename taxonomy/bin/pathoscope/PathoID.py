#!/usr/bin/python
# Initial author: Solaiappan Manimaran
# Functions to read alignment file (sam/gnu-sam or bl8), run EM algorithm
# and output report file that can be opened in Excel and also 
# output updated alignment file (sam/gnu-sam or bl8)  

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

import pathoscope_util, os, math, csv
# ===========================================================
def conv_align2GRmat(aliDfile,pScoreCutoff,aliFormat):
	in1 = open(aliDfile,'r')
	U = {}
	NU = {}
	h_readId = {}
	h_refId = {}
	genomes = []
	read =[]
	gCnt = 0
	rCnt = 0

	mxBitSc = 700
	sigma2 = 3
	for ln in in1:
		if (ln[0] == '@' or ln[0] == '#'):
			continue

		l = ln.split('\t')
		
		readId=l[0]
		if (aliFormat == 0 or aliFormat == 1): # gnu-sam or sam
			refId=l[2]
		elif (aliFormat == 2): # bl8
			refId=l[1]
		
		if refId == '*':
			continue

		if (aliFormat == 0): # gnu-sam
			if len(l) < 13:
				print "line has no score entry: %s" %(ln)
				continue
			scoreStringList = l[12].split(':')
			if len(scoreStringList) < 3:
				print "Number format error in line: %s" %(ln)
				continue
			pScore = float(scoreStringList[2])
			if (pScore == float('Inf')):
				pScore = 1.0
			if (pScore < pScoreCutoff):
				continue
		elif (aliFormat == 1): # sam
			mapq = float(l[4])
			mapq2 = mapq/(-10.0)
			pScore = 1.0 - pow(10,mapq2)
			if (pScore < pScoreCutoff):
				continue
		elif (aliFormat == 2): # bl8
			eVal = float(l[10])
			if (eVal > pScoreCutoff):
				continue
			bitSc = float(l[11])/sigma2
			if bitSc > mxBitSc:
				bitSc = mxBitSc
			pScore = math.exp(bitSc)
		pScore = int(round(pScore*100)) # Converting to integer to conserve memory space
		if pScore < 1:
			continue
		
		gIdx = h_refId.get(refId,-1)
		if gIdx == -1:
			gIdx = gCnt
			h_refId[refId] = gIdx
			genomes.append(refId)
			gCnt += 1

		rIdx = h_readId.get(readId,-1)
		if rIdx == -1:
			#hold on this new read
			#first, wrap previous read profile and see if any previous read has a same profile with that!
			rIdx = rCnt
			h_readId[readId] = rIdx
			read.append(readId)
			rCnt += 1
			U[rIdx] = [[gIdx], [pScore], [float(pScore)]]
		else:
			if (rIdx in U):
				if gIdx in U[rIdx][0]:
					continue
				NU[rIdx] = U[rIdx]
				del U[rIdx]
			if gIdx in NU[rIdx][0]:
				continue
			NU[rIdx][0].append(gIdx)
			NU[rIdx][1].append(pScore)
			NU[rIdx][2][0] += pScore
#			length = len(NU[rIdx][1])
#			NU[rIdx][2] = [1.0/length]*length

	in1.close()

	del h_refId, h_readId
	for rIdx in U:
		U[rIdx] = U[rIdx][0][0] #keep gIdx only
	for rIdx in NU:
		pScoreSum = NU[rIdx][2][0]
		NU[rIdx][2] = [k/pScoreSum for k in NU[rIdx][1]] #Normalizing pScore

	return U, NU, genomes, read

# ===========================================================
def pathoscope_reassign(out_matrix, verbose, scoreCutoff, expTag, ali_format, ali_file, outdir, emEpsilon, maxIter, upalign):
	
	if ali_format == 'gnu-sam':
		aliFormat = 0
		if verbose:
			print "parsing gnu-sam file/likelihood score/reads and mapped genomes..."
	elif ali_format == 'sam': #standard sam
		aliFormat = 1
		if verbose:
			print "parsing sam file/likelihood score/reads and mapped genomes..."
	elif ali_format == 'bl8': #blat m8 format
		aliFormat = 2
		if verbose:
			print "parsing bl8 file/likelihood score/reads and mapped genomes..."
	else:
		print "unknown alignment format file..."
		return
	(U, NU, genomes, read) = conv_align2GRmat(ali_file,scoreCutoff,aliFormat)
	
	nG = len(genomes)
	nR = len(read)
	if verbose:
		print "EM iteration..."
		print "(G,R)=%dx%d" % (nG, nR)
	
	if out_matrix:
		if verbose:
			print "writing initial alignment ..."
		out_initial_align_matrix(genomes, read, U, NU, expTag, ali_file, outdir)	

	(bestHitInitialReads, bestHitInitial, level1Initial, level2Initial) = \
		computeBestHit(U, NU, genomes, read)
	
	(initPi, pi, _, NU) = pathoscope_em(U, NU, genomes, maxIter, emEpsilon, verbose)
	tmp = zip(initPi,genomes)
	tmp = sorted(tmp,reverse=True) #similar to sort row
	
	if out_matrix:
		initialGuess = outdir + os.sep + expTag + '-initGuess.txt'
		oFp = open(initialGuess,'wb')
		csv_writer = csv.writer(oFp, delimiter='\t')
		csv_writer.writerows(tmp)
		oFp.close()
	
	del tmp
	
	(bestHitFinalReads, bestHitFinal, level1Final, level2Final) = \
		computeBestHit(U, NU, genomes, read)

	if out_matrix:
		finalGuess = outdir + os.sep + expTag + '-finGuess.txt'
		oFp = open(finalGuess,'wb')
		tmp = zip(pi,genomes)
		tmp = sorted(tmp,reverse=True)
		csv_writer = csv.writer(oFp, delimiter='\t')
		csv_writer.writerows(tmp)
		oFp.close()

	finalReport = outdir + os.sep + expTag +'-'+ ali_format + '-report.tsv'
	oFp = open(finalReport,'wb')
	tmp = zip(pi,genomes, initPi, bestHitInitial, bestHitInitialReads, bestHitFinal, bestHitFinalReads, \
		level1Initial, level2Initial, level1Final, level2Final)
	tmp = sorted(tmp,reverse=True) # Sorting based on Final Guess
	x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11 = zip(*tmp)
	for i in range(len(x10)):
		if (x10[i] <= 0 and x11[i] <= 0):
			break
		if i == (len(x10)-1):
			i += 1
	tmp = zip (x2[:i], x1[:i], x6[:i], x7[:i], x10[:i], x11[:i], x3[:i], x4[:i], x5[:i], x8[:i], x9[:i]) # Changing the column order here
	csv_writer = csv.writer(oFp, delimiter='\t')
	header = ['Total Number of Aligned Reads:', nR, 'Total Number of Mapped Genomes:', nG]
	csv_writer.writerow(header)
	header = ['Genome', 'Final Guess', 'Final Best Hit', 'Final Best Hit Read Numbers', 'Final High Confidence Hits', \
		'Final Low Confidence Hits', 'Initial Guess', 'Initial Best Hit', 'Initial Best Hit Read Numbers', \
		'Initial High Confidence Hits', 'Initial Low Confidence Hits']
	csv_writer.writerow(header)
	csv_writer.writerows(tmp)
	oFp.close()
	
	reAlignfile = ali_file
	if upalign:
		reAlignfile = rewrite_align(U, NU, ali_file,scoreCutoff, aliFormat,outdir)

	return (finalReport, x2, x3, x4, x5, x1, x6, x7, x8, x9, x10, x11, reAlignfile)

# ===========================================================
# This is the main EM algorithm
# ===========================================================
def pathoscope_em(U, NU, genomes, maxIter, emEpsilon, verbose):
	G = len(genomes)

	### Initial values
	pi = [1./G for _ in genomes]
	initPi = pi
	theta = [1./G for _ in genomes]
	
	pisum0=[0 for i in genomes]
	for i in U: 
		pisum0[U[i]]+=1
	
	### data structure
	### 3 unique reads: 2 reads to genome1 and 1 read to genome4
	# U = {0: 0, 1: 0, 2: 3}
	### non-unique reads: 3 total reads  1:[[genomes],[qij],[xij]]
	#NU = {0: [[0, 2, 3], [0.4, 0.2, 0.4] , [0.33, 0.33, 0.33]], 1: [[0, 1], [0.6, 0.4] , [0.5,0.5]], 2: [[1, 3], [0.5, 0.5] , [0.5,0.5]]}
	### Genome hash
	#genomes = {0:"ecoli", 1:"strep", 2:"anthrax", 3:"plague"}
	lenNU=len(NU)
	if lenNU==0:
		lenNU=1

	for i in range(maxIter):  ## EM iterations--change to convergence 
		pi_old = pi
		thetasum=[0 for k in genomes]
		
		# E Step 

		for j in NU: #for each non-uniq read, j
			z = NU[j] 
			ind = z[0] #a set of any genome mapping with j
			pitmp = [pi[k] for k in ind]			### get relevant pis for the read
			thetatmp = [theta[k] for k in ind]  	### get relevant thetas for the read
			xtmp = [1.*pitmp[k]*thetatmp[k]*z[1][k] for k in range(len(ind))]  ### Calculate unormalized xs
			xsum = sum(xtmp)
			xnorm = [1.*k/xsum for k in xtmp]  		### Normalize new xs
			
			NU[j][2] = xnorm     					## Update x in NU 
			
			for k in range(len(ind)):
				thetasum[ind[k]] += xnorm[k]   		### Keep running tally for theta
		# M step	
		pisum = [thetasum[k]+pisum0[k] for k in range(len(thetasum))]   ### calculate tally for pi
		pip = 0 # pi prior - may be updated later
		pi = [(1.*k+pip)/(len(U)+len(NU)+pip*len(pisum)) for k in pisum]  		## update pi
		
		#pi = [1.*k/G for k in pisum]  		## update pi
		if (i == 0):
			initPi = pi
		
		thetap=0 # theta prior - may be updated later
		theta = [(1.*k+thetap)/(lenNU+thetap*len(thetasum)) for k in thetasum]
		
		cutoff = 0.0
		for k in range(len(pi)):
			cutoff += abs(pi_old[k]-pi[k])
		if verbose:
			print "[%d]%g" % (i,cutoff)
		if (cutoff <= emEpsilon or lenNU==1):
			break

	return initPi, pi, theta, NU

def out_initial_align_matrix(ref, read, U, NU, expTag, ali_file, outdir):
	genomeId = outdir + os.sep + expTag + '-genomeId.txt'
	oFp = open(genomeId,'wb')
	csv_writer = csv.writer(oFp, delimiter='\n')
	csv_writer.writerows([ref])
	oFp.close()

	readId = outdir + os.sep + expTag + '-readId.txt'
	oFp = open(readId,'wb')
	csv_writer = csv.writer(oFp, delimiter='\n')
	csv_writer.writerows([read])
	oFp.close()
	
def computeBestHit(U, NU, genomes, read):
	bestHitReads=[0.0 for _ in genomes]
	level1Reads=[0.0 for _ in genomes]
	level2Reads=[0.0 for _ in genomes]
	for i in U: 
		bestHitReads[U[i]]+=1
		level1Reads[U[i]] += 1
	for j in NU:
		z = NU[j]
		ind = z[0]
		xnorm = z[2]
		bestGenome = max(xnorm)
		numBestGenome = 0
		for i in range(len(xnorm)):
			if (xnorm[i] == bestGenome):
				numBestGenome += 1
		if (numBestGenome == 0):
			numBestGenome = 1
		for i in range(len(xnorm)):
			if (xnorm[i] == bestGenome):
				bestHitReads[ind[i]] += 1.0/numBestGenome
				if (xnorm[i] >= 0.5):
					level1Reads[ind[i]] += 1
				elif (xnorm[i] >= 0.01):
					level2Reads[ind[i]] += 1
		
	nG = len(genomes)
	nR = len(read)
	bestHit = [bestHitReads[k]/nR for k in range(nG)]
	level1 = [level1Reads[k]/nR for k in range(nG)]
	level2 = [level2Reads[k]/nR for k in range(nG)]
	return bestHitReads, bestHit, level1, level2

# ===========================================================
def rewrite_align(U, NU, aliDfile, pScoreCutoff, aliFormat,outdir):
	pathoscope_util.ensure_dir(outdir)
	f = os.path.basename(aliDfile)
	reAlignfile = outdir + os.sep + 'updated_' + f

	with open(reAlignfile,'w') as of:
		with open(aliDfile,'r') as in1:
			h_readId = {}
			h_refId = {}
			genomes = []
			read =[]
			gCnt = 0
			rCnt = 0
		
			mxBitSc = 700
			sigma2 = 3
			for ln in in1:
				if (ln[0] == '@' or ln[0] == '#'):
					of.write(ln)
					continue
		
				l = ln.split('\t')
				
				readId=l[0]
				if (aliFormat == 0 or aliFormat == 1): # gnu-sam or sam
					refId=l[2]
				elif (aliFormat == 2): # bl8
					refId=l[1]
				
				if refId == '*':
					continue
		
				if (aliFormat == 0): # gnu-sam
					pScore = float(l[12].split(':')[2])
					if (pScore < pScoreCutoff):
						continue
				elif (aliFormat == 1): # sam
					mapq = float(l[4])
					mapq2 = mapq/(-10.0)
					pScore = 1.0 - pow(10,mapq2)
					if (pScore < pScoreCutoff):
						continue
				elif (aliFormat == 2): # bl8
					eVal = float(l[10])
					if (eVal > pScoreCutoff):
						continue
					bitSc = float(l[11])/sigma2
					if bitSc > mxBitSc:
						bitSc = mxBitSc
					pScore = math.exp(bitSc)
				pScore = int(round(pScore*100)) # Converting to integer to conserve memory space
				if pScore < 1:
					continue
				
				gIdx = h_refId.get(refId,-1)
				if gIdx == -1:
					gIdx = gCnt
					h_refId[refId] = gIdx
					genomes.append(refId)
					gCnt += 1
		
				rIdx = h_readId.get(readId,-1)
				if rIdx == -1:
					#hold on this new read
					#first, wrap previous read profile and see if any previous read has a same profile with that!
					rIdx = rCnt
					h_readId[readId] = rIdx
					read.append(readId)
					rCnt += 1
					if rIdx in U:
						of.write(ln)
						continue
							
				if rIdx in NU:
					if (aliFormat == 0): # gnu-sam
						scoreComponents = l[12].split(':')
						pScore = float(scoreComponents[2])
						if (pScore < pScoreCutoff):
							continue
						(upPscore, pscoreSum) = find_updated_score(NU, rIdx, gIdx)
						scoreComponents[2] = str(upPscore*pscoreSum)
						if (scoreComponents[2] < pScoreCutoff):
							continue
						l[12] = ':'.join(scoreComponents)
						ln = '\t'.join(l)
						of.write(ln)
					elif (aliFormat == 1): # sam
						mapq = float(l[4])
						mapq2 = mapq/(-10.0)
						pScore = 1.0 - pow(10,mapq2)
						if (pScore < pScoreCutoff):
							continue
						(upPscore, pscoreSum) = find_updated_score(NU, rIdx, gIdx)
						if (upPscore < pScoreCutoff):
							continue
						if (upPscore >= 1.0):
							upPscore = 0.999999
						mapq2 = math.log10(1 - upPscore)
						l[4] = str(int(round(-10.0*mapq2)))
						ln = '\t'.join(l)
						of.write(ln)
					elif (aliFormat == 2): # bl8
						eVal = float(l[10])
						if (eVal > pScoreCutoff):
							continue
						(upPscore, pscoreSum) = find_updated_score(NU, rIdx, gIdx)
						score = upPscore*pscoreSum
						if score <= 0.0:
							continue
						bitSc = math.log(score)
						if bitSc > mxBitSc:
							bitSc = mxBitSc
						l[10] = str(bitSc*sigma2)
						ln = '\t'.join(l)
						of.write(ln)

	return reAlignfile

def find_updated_score(NU, rIdx, gIdx):
	index = NU[rIdx][0].index(gIdx);
	pscoreSum = 0.0
	for pscore in NU[rIdx][1]:
		pscoreSum += pscore
	pscoreSum /= 100
	upPscore = NU[rIdx][2][index]
	return (upPscore, pscoreSum)

