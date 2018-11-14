"""
	Get an input file with SVs
	Plot the size distribution

"""

import sys
import matplotlib.pyplot as plt
import numpy as np

deletionsFile = sys.argv[1]

deletionsLengths = []
deletionsLengthsDict = dict()
with open(deletionsFile, 'r') as delF:
	lineCount = 0
	for line in delF:
		if lineCount < 1:
			lineCount += 1
			continue
		
		line = line.strip()
		splitLine = line.split("\t")
		
		#Get the length of the deletion. The deletions are intrachromosomal and in this file have the same start and end for both chromosomes. 
		start = int(splitLine[1])
		end = int(splitLine[2])
		length = end - start
		deletionsLengths.append(length)
		if length not in deletionsLengthsDict:
			deletionsLengthsDict[length] = 0
		deletionsLengthsDict[length] += 1
		
		
#Make a boxplot
plt.boxplot(deletionsLengths)
#plt.show()

#Some statistics
print np.min(deletionsLengths)
print np.max(deletionsLengths)
print np.median(deletionsLengths)

print "number of dels with length 0: ", deletionsLengthsDict[0]
