"""
	For every gene, compute how often the features are higher than in the random case
"""

from __future__ import absolute_import
import sys
import numpy as np
import os
from os import listdir
from os.path import isfile, join
from six.moves import range

inFile = sys.argv[1]
permutationFolder = sys.argv[2]
outFile = sys.argv[3]

#Read the real gene scores
geneScores = np.loadtxt(inFile, dtype="object")

delta = np.zeros(geneScores.shape)
delta = np.array(delta, dtype="object")
delta[:,0] = geneScores[:,0]

permutationFiles = [f for f in listdir(permutationFolder) if isfile(join(permutationFolder, f))]

for permutationFile in permutationFiles:
	if permutationFile == "realSVs_geneScores_chr.txt":
		continue
	
	permutationScores = np.loadtxt(permutationFolder + permutationFile, dtype="object")
	
	#Compute the delta per gene.
	for row in range(0, geneScores.shape[0]):
		for col in range(4, geneScores.shape[1]-1):
			if float(geneScores[row,col]) > float(permutationScores[row,col]):
				delta[row,col] += 1

#Compute total
for row in range(0, delta.shape[0]):
	
	rowSum = np.sum(delta[row,1:delta.shape[1]-1])
	delta[row,delta.shape[1]-1] = rowSum
	
delta = delta[delta[:,30].argsort()[::-1]] #Select the column  to rank by	
		
header = "geneName\tgeneScore\teQTLGains\teQTLLosses\tenhancerGains\tenhancerLosses\tpromoterGains\tpromoterLosses\tcpgGains\tcpgLosses\ttfGains\ttfLosses\thicGains\thicLosses\th3k9me3Gains\th3k9me3Losses\th3k4me3Gains\th3k4me3Losses\th3k27acGains\th3k27acLosses\th3k27me3Gains\th3k27me3Losses\th3k4me1Gains\th3k4me1Losses\th3k36me3Gains\th3k36me3Losses\tdnaseIGains\tdnaseILosses\ttotal\tsamples"
				
#Write to numpy output file	
np.savetxt(outFile, delta, delimiter='\t', fmt='%s', header=header)	

