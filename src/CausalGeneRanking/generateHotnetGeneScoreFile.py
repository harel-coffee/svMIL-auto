"""
	Process the tool output to get input scores for hotnet. 

"""

import sys
import numpy as np

#Get the non-coding SV-gene pairs
nonCodingPairs = np.loadtxt(sys.argv[1], dtype='object')

#Get the coding SV-gene pairs
codingPairs = np.loadtxt(sys.argv[2], dtype='object')


#Count how many unique samples have an SV linked to that gene and use it as a score.
nonCodingGeneSamples = dict()
for pair in nonCodingPairs[:,0]:
	splitPair = pair.split("_")
	gene = splitPair[0]
	sample = splitPair[len(splitPair)-1]
	
	if gene not in nonCodingGeneSamples:
		nonCodingGeneSamples[gene] = []
	if sample not in nonCodingGeneSamples:
		nonCodingGeneSamples[gene].append(sample)
	
codingGeneSamples = dict()
for pair in codingPairs:
	splitPair = pair.split("_")
	gene = splitPair[0]
	sample = splitPair[len(splitPair)-1]
	
	if gene not in codingGeneSamples:
		codingGeneSamples[gene] = []
	if sample not in codingGeneSamples[gene]:
		codingGeneSamples[gene].append(sample)	
	

nonCodingGeneScores = []
for gene in nonCodingGeneSamples:
	nonCodingGeneScores.append([gene, len(nonCodingGeneSamples[gene])])

nonCodingGeneScores = np.array(nonCodingGeneScores)

codingGeneScores = []
for gene in codingGeneSamples:
	codingGeneScores.append([gene, len(codingGeneSamples[gene])])

codingGeneScores = np.array(codingGeneScores)

#generate the score files
np.savetxt('Output/hotnet_noncoding_scores.txt', nonCodingGeneScores, fmt='%s', header='geneName\tgeneScore')
np.savetxt('Output/hotnet_coding_scores.txt', codingGeneScores, fmt='%s', header='geneName\tgeneScore')

