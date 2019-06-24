import sys
import numpy as np

geneScores = np.loadtxt(sys.argv[1], dtype="object")

rnkFormat = np.empty([geneScores.shape[0],2],dtype="object")
rnkFormat[:,0] = geneScores[:,0]
#rnkFormat[:,1] = geneScores[:,30]

noOfSamples = []
for rank in range(0, len(geneScores[:,31])):
	if geneScores[rank,31] == "None":
		noOfSamples.append(0)
	else:
		noOfSamples.append(len(geneScores[rank,31].split(",")))

rnkFormat[:,1] = np.array(noOfSamples)
 
np.savetxt(sys.argv[2], rnkFormat, delimiter='\t', fmt='%s')





