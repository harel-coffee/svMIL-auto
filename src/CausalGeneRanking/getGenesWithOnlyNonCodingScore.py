from __future__ import absolute_import
from __future__ import print_function
import sys

rankedGenes = sys.argv[1]
geneCount = 0 
with open(rankedGenes, 'r') as inF:
	
	lineCount = 0
	
	for line in inF:
		
		if lineCount < 1:
			lineCount += 1
			continue
	

		line = line.strip()
		splitLine = line.split("\t")
		
		geneScore = splitLine[1]
		total = splitLine[len(splitLine)-1]
		
		if float(total) > 0 and float(geneScore) == 0:
			#print "Gene ", splitLine[0], " has non-coding score but no gene score"
			print(splitLine[0], total)
			geneCount += 1
		
print(geneCount)		
