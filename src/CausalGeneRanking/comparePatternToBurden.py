"""
	Compare the genes in the ranking between pattern and burden. 
"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
from six.moves import range

burdenFile = sys.argv[1]
patternFile = sys.argv[2]

burden = np.loadtxt(burdenFile, dtype="object")
pattern = np.loadtxt(patternFile, dtype="object")

burdenGenes = []
for row in range(0,burden.shape[0]):
	
	if float(burden[row,28]) > 0:
		burdenGenes.append(burden[row,0])

patternGenes = []
for row in range(0,pattern.shape[0]):
	
	if float(pattern[row,28]) > 0:
		patternGenes.append(pattern[row,0])


#compute set difference
print("everything in burden not in pattern: ", np.setdiff1d(burdenGenes, patternGenes))
print("everything in pattern not in burden: ", len(np.setdiff1d(patternGenes, burdenGenes)))

print("overlap: ", len(np.intersect1d(burdenGenes, patternGenes)))
print("total burden & pattern: ", len(burdenGenes), len(patternGenes))


