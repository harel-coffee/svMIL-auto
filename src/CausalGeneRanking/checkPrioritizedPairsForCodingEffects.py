"""
	for every SV-gene pair from the windowed, or TAD, approach, check if any of the linked genes also have a coding SV in the same sample that is DEG. 

"""

import sys
import numpy as np

prioritizedPairs = np.load(sys.argv[1], allow_pickle=True, encoding='latin1')
codingDegPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')


print("Starting number of pairs: ", len(prioritizedPairs))
print("Number of coding DEG pairs: ", len(codingDegPairs))

#First split the coding pairs into gene-sample, we do not check for the same SV
splitCodingDegPairs = []
for pair in codingDegPairs[:,0]:
	splitPair = pair.split("_")
	newPair = splitPair[0] + "_" + splitPair[len(splitPair)-1]
	splitCodingDegPairs.append(newPair)

pairsWithCodingEvidence = 0
for pair in prioritizedPairs:
	
	splitPair = pair[0].split("_")
	newPair = splitPair[0] + "_" + splitPair[len(splitPair)-1]
	
	if newPair in splitCodingDegPairs: #check if there 
		pairsWithCodingEvidence += 1
		
print("Explained with coding evidence for: ", pairsWithCodingEvidence)		
	

