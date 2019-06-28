"""
	Get the total score file, a subset with names of genes, and output the scores of these genes to a new file for further analyses

"""

import sys
import numpy as np

scores = np.loadtxt(sys.argv[1], dtype='object')
subset = np.loadtxt(sys.argv[2], dtype='object')
print subset

subsetScores = []

for gene in scores[:,0]:
	if gene in subset:
		continue
	geneScores = scores[np.where(scores[:,0] == gene)[0],:][0]
	#subsetScores.append([geneScores[0], geneScores[30]])
	subsetScores.append(geneScores)



# subsetScores = np.array(subsetScores, dtype="object")
# subsetScores[:,1] = subsetScores[:,1].astype(float)
# subsetScores = subsetScores[subsetScores[:,1].argsort()[::-1]] #Select the column  to rank by

subsetScores = np.array(subsetScores, dtype="object")
subsetScores[:,30] = subsetScores[:,30].astype(float)
subsetScores = subsetScores[subsetScores[:,30].argsort()[::-1]] #Select the column  to rank by

np.savetxt(sys.argv[3], subsetScores, delimiter='\t', fmt='%s')
print subsetScores
exit()

#Count the number of samples that these subset genes have
sampleLengths = []
rankings = []
for gene in subsetScores:
	
	samples = gene[31]
	splitSamples = samples.split(",")
	sampleLengths.append(len(splitSamples))
	rankings.append(gene[30])

print sampleLengths

import matplotlib.pyplot as plt

# plt.hist(sampleLengths)
# plt.show()

plt.hist(rankings)
plt.show()
