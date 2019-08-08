

from __future__ import absolute_import
import sys
import numpy as np

scores = np.loadtxt(sys.argv[1], dtype="object")
subset = np.loadtxt(sys.argv[2], dtype="object")

#get the number of patients per gene
patientDistribution = dict()
for gene in subset[:,0]:
	
	#get the score entry
	geneScore = scores[np.where(scores[:,0] == gene)[0]]
	if len(geneScore) < 1:
		continue
	geneScore = geneScore[0]
	
	splitPatients = geneScore[31].split(",")
	
	
	numberOfPatients = len(splitPatients)
	
	if numberOfPatients not in patientDistribution:
		patientDistribution[numberOfPatients] = 0
		
	patientDistribution[numberOfPatients] += 1	


import matplotlib.pyplot as plt

plt.bar(list(patientDistribution.keys()), list(patientDistribution.values()))
plt.show()



