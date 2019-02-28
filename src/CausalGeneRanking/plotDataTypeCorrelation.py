"""
	Make a plot of the correlation between the scores of all the different data types.
	
	On the y-axis, show the genes.
	On the x-axis, show the scores of the data types. Each data type should get a line with a different color. 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt

scoresFile = sys.argv[1]
header = []

with open(scoresFile, 'r') as inF:
	
	lineCount = 0
	for line in inF:
		line = line.strip()
		splitLine = line.split("\t")
		if lineCount < 1:
			for ind in range(1, len(splitLine)-1):
				header.append(splitLine[ind])
			break
header = np.array(header)
print header
scores = np.loadtxt(scoresFile, dtype="object")


#Correlate histone marks with promoters & elements
enhancerGains = [float(score) for score in scores[:,4]]
enhancerLosses = [float(score) for score in scores[:,5]]
promoterGains = [float(score) for score in scores[:,6]]
promoterLosses = [float(score) for score in scores[:,7]]

h3k4me3Gains = [float(score) for score in scores[:,22]]
h3k4me3Losses = [float(score) for score in scores[:,23]]

h3k27acGains = [float(score) for score in scores[:,18]]
h3k27acLosses = [float(score) for score in scores[:,19]]

print "Correlation between enhancer gains and H3K27ac gains: ", np.corrcoef(enhancerGains, h3k27acGains)[0,1]
print "Correlation between enhancer losses and H3K27ac losses: ", np.corrcoef(enhancerLosses, h3k27acLosses)[0,1]
print "Correlation between promoter gains and h3k4me3 losses: ", np.corrcoef(promoterGains, h3k4me3Gains)[0,1]
print "Correlation between promoter losses and h3k4me3 losses: ", np.corrcoef(promoterLosses, h3k4me3Losses)[0,1]

exit()


print scores

scores = scores[1:100, :]

 

for col in range(1, scores.shape[1]-1):
	
	plt.plot(range(0, scores.shape[0]), list(scores[:,col].astype("float")), linestyle='-')
	if col % 10 == 0: #we can only show 10 lines at the same time 

		plt.legend(header[col-10:col])	
		plt.show()
		plt.clf()
		

plt.plot(range(0, scores.shape[0]), list(scores[:,17].astype("float")), linestyle='-')
plt.plot(range(0, scores.shape[0]), list(scores[:,18].astype("float")), linestyle='-')
plt.plot(range(0, scores.shape[0]), list(scores[:,19].astype("float")), linestyle='-')
plt.plot(range(0, scores.shape[0]), list(scores[:,20].astype("float")), linestyle='-')
plt.plot(range(0, scores.shape[0]), list(scores[:,21].astype("float")), linestyle='-')
plt.plot(range(0, scores.shape[0]), list(scores[:,22].astype("float")), linestyle='-')
plt.plot(range(0, scores.shape[0]), list(scores[:,23].astype("float")), linestyle='-')
plt.plot(range(0, scores.shape[0]), list(scores[:,24].astype("float")), linestyle='-')
plt.plot(range(0, scores.shape[0]), list(scores[:,25].astype("float")), linestyle='-')
plt.plot(range(0, scores.shape[0]), list(scores[:,26].astype("float")), linestyle='-')
plt.plot(range(0, scores.shape[0]), list(scores[:,27].astype("float")), linestyle='-')
print header[17:27]
plt.legend(header[17:27])	
plt.show()
plt.clf()