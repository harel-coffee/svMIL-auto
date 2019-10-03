"""
	Do some quick checks to see what nice things we can get from the method. Can be cleaned up into proper scripts later. 

"""

import sys
import numpy as np
from inputParser import InputParser

#1. Find how many of the genes are recurrent across nc samples, how many have coding evidence in other samples, and how many of these are DEG.

ncScores = np.loadtxt(sys.argv[1], dtype='object')

#compute how many pairs also have coding evidence in other samples
totalPairs = []
pairsWithCodingEvidence = []
for gene in ncScores:
	
	samples = gene[31]
	if samples != "None":
		splitSamples = samples.split(",")
	
		for sample in splitSamples:
			totalPairs.append(gene[0] + "_" + sample)
			codingSamples = gene[32]
			
			if codingSamples != "None":
				splitCodingSamples = codingSamples.split(",")
				
				pairsWithCodingEvidence.append(gene[0] + "_" + sample)
				

print(len(totalPairs))
print(len(pairsWithCodingEvidence))

#compute how many of these pairs have DEG evidence. 
#degPairs = np.load(sys.argv[2], allow_pickle='True')

#For non-coding and coding pairs, check which are overlapping TAD boundaries. Also, disrupting TAD boundaries (smaller SVs).

ncPairs = np.loadtxt(sys.argv[3], dtype='object')
codingPairs = np.loadtxt(sys.argv[4], dtype='object')

ncSVs = []
for pair in ncPairs[:,0]:
	
	splitPair = pair.split("_")
	ncSVs.append([splitPair[1], int(splitPair[2]), int(splitPair[3]), splitPair[4], int(splitPair[5]), int(splitPair[6])])

ncSVs = np.array(ncSVs, dtype='object')

#Get the non-coding pairs from the separate pairs file
#Get the TADs
#Compute how many TAD boundaries are overlapped in the nc case
tads = InputParser().getTADsFromFile(sys.argv[5])

#First just for the general overlap. How many TAD boundaries are covered by an SV?
overlapNc = []
for tad in tads:
	
	#check if there is any SV overlapping this boundary
	#we can exclude translocations here, so no need to look at chr2
	#if the same SV overlaps both the start and the end, count it only once

	svChrSubset = ncSVs[ncSVs[:,0] == tad[0],:]
	
	startMatches = (tad[1] >= svChrSubset[:,1]) * (tad[1] <= svChrSubset[:,5])
	endMatches = (tad[2] >= svChrSubset[:,1]) * (tad[2] <= svChrSubset[:,5])
	
	#Use sum to include either, and both only once
	matches = startMatches + endMatches
	
	matchingSVs = svChrSubset[matches]
	overlapNc.append(len(matchingSVs))
	

#Then get the coding pairs from the same set
#repeat the analysis

codingSVs = []
for pair in codingPairs:
	
	splitPair = pair.split("_")
	codingSVs.append([splitPair[1], int(splitPair[2]), int(splitPair[3]), splitPair[4], int(splitPair[5]), int(splitPair[6])])

codingSVs = np.array(codingSVs, dtype='object')


#First just for the general overlap. How many TAD boundaries are covered by an SV?
overlapC = []
for tad in tads:
	
	#check if there is any SV overlapping this boundary
	#we can exclude translocations here, so no need to look at chr2
	#if the same SV overlaps both the start and the end, count it only once

	svChrSubset = codingSVs[codingSVs[:,0] == tad[0],:]
	
	startMatches = (tad[1] >= svChrSubset[:,1]) * (tad[1] <= svChrSubset[:,5])
	endMatches = (tad[2] >= svChrSubset[:,1]) * (tad[2] <= svChrSubset[:,5])
	
	#Use sum to include either, and both only once
	matches = startMatches + endMatches
	
	matchingSVs = svChrSubset[matches]
	overlapC.append(len(matchingSVs))

print(np.mean(overlapNc))
print(np.mean(overlapC))

#To make it more specific, look at TAD disrupting SVs, meaning that they actually start and end within a TAD. How many of these are really overlapping a boundary? 

