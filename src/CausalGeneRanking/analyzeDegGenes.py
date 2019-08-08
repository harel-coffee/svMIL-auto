"""
	Do some analysis on the DEG per patient significant gene set. 

"""

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import sys
from six.moves import range

genePatientPairs = np.loadtxt(sys.argv[1], dtype="object")


#1. Check how many genes are in the set multiple times
#2. Check how many samples are in the set multiple times, are some patients overrepresented?

geneCounts = dict()
patientCounts = dict()
for pair in genePatientPairs[:,0]:
	
	splitPair = pair.split("_")
	
	if splitPair[0] not in geneCounts:
		geneCounts[splitPair[0]] = 0
	geneCounts[splitPair[0]] += 1
	
	if splitPair[1] not in patientCounts:
		patientCounts[splitPair[1]] = 0
	patientCounts[splitPair[1]] += 1
	
print("Unique number of genes: ", len(geneCounts))
print("Unique number of patients: ", len(patientCounts))

geneCountsArray = np.empty([len(geneCounts), 2], dtype="object")
patientCountsArray = np.empty([len(patientCounts), 2], dtype="object")

for geneCountInd in range(0, len(list(geneCounts.keys()))):
	
	geneCountsArray[geneCountInd,0] = list(geneCounts.keys())[geneCountInd]
	geneCountsArray[geneCountInd,1] = list(geneCounts.values())[geneCountInd]
	
sortedGeneCounts = geneCountsArray[geneCountsArray[:,1].argsort()][::-1]
print(sortedGeneCounts[0:20,:])

for patientCountInd in range(0, len(list(patientCounts.keys()))):
	
	patientCountsArray[patientCountInd,0] = list(patientCounts.keys())[patientCountInd]
	patientCountsArray[patientCountInd,1] = list(patientCounts.values())[patientCountInd]
	
sortedPatientCounts = patientCountsArray[patientCountsArray[:,1].argsort()][::-1]
print(sortedPatientCounts[0:20,:])


# #What is the distribution of the number of times that we see specific genes and patients DEG?
# 
# import matplotlib.pyplot as plt
# 
# # plt.hist(geneCounts.values())
# # plt.show()
# 
# plt.hist(patientCounts.values())
# plt.show()
# 
# exit()
# 





