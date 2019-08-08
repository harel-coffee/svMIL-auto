"""
	In a ranking file, check how many genes have no non-coding SVs, but do have coding SVs.
	For each gene, check if there are samples in the coding column that are not in the non-coding column. These are coding only.
	The same for non-coding, and for the mixed case. 

"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np

geneRanks = np.loadtxt(sys.argv[1], dtype="object")

codingOnly = 0
nonCodingOnly = 0
mixedNC = 0
mixedC = 0
both = 0
for rank in geneRanks:
	
	ncSamples = rank[31]
	codingSamples = rank[32]
	
	if ncSamples != "None":
		ncSamples = ncSamples.split(",")
	else:
		ncSamples = []
	if codingSamples != "None":
		codingSamples = codingSamples.split(",")
	else:
		codingSamples = []
	
	for codingSample in codingSamples:
		if codingSample not in ncSamples:
			codingOnly += 1
			break
	for nonCodingSample in ncSamples:
		if nonCodingSample not in codingSamples:
			nonCodingOnly += 1
			break
	
	if len(ncSamples) > 0 and len(codingSamples) > 0:
		both += 1
	if (len(ncSamples) > 0 and len(codingSamples) == 0):
		mixedNC +=1
	
	if (len(ncSamples) == 0 and len(codingSamples) > 0):
		mixedC += 1
		
print("Number of genes with non coding effects only in at least 1 sample: ", nonCodingOnly)
print("Number of genes with coding effects only in at least 1 sample: ", codingOnly)
print("Number of genes found with both coding and non-coding: ", both)
print("Number of genes found only when looking at coding: ", mixedC)
print("Number of genes found only when looking at non-coding: ", mixedNC)
			
	
