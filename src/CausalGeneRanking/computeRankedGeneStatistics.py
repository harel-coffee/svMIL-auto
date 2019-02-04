"""

	Investigate how many of the top ranking genes are is COSMIC. 
"""

import sys


#1. How many genes with a score >0 are in COSMIC?
#2. How mant genes with a score == 0 are in COSMIC?

#3. What is the maximum length of string of COSMIC genes starting from the top? 

rankedGenesFile = sys.argv[1]
cosmicGenesFile = sys.argv[2]
breastCancerGenesFile = sys.argv[3]

cosmicGenes = []
with open(cosmicGenesFile, 'rb') as f:
	lineCount = 0
	for line in f:
		if lineCount == 0:
			lineCount += 1
			continue
		
		splitLine = line.split("\t")
		
		geneName = splitLine[0]
		cosmicGenes.append(geneName)

print "total number of cosmic genes: ", len(cosmicGenes)		
	
#Also read the breast cancer genes specifically

breastCancerGenes = []
with open(breastCancerGenesFile, 'r') as f:
	
	for line in f:
		
		line = line.strip()
		
		breastCancerGenes.append(line)

		
cosmicCountGoodScore = 0
cosmicCountBadScore = 0
allGenesGoodScore = 0
allGenesBadScore = 0
bcCountGoodScore = 0
bcCountBadScore = 0
genes = dict()
with open(rankedGenesFile, 'rb') as f:
	
	for line in f:
		line = line.strip()
		splitLine = line.split("\t")
		
		
		if float(splitLine[3]) == 0:
			if splitLine[0] in cosmicGenes:
				print "COSMIC gene: ", splitLine[0]
				cosmicCountGoodScore += 1
			if splitLine[0] in breastCancerGenes:
				print "BC gene: ", splitLine[0]
				bcCountGoodScore += 1
			allGenesGoodScore += 1
		else:
			if splitLine[0] in cosmicGenes:
				if splitLine[0] not in genes:
					genes[splitLine[0]] = 0

				genes[splitLine[0]] += 1
				# if genes[splitLine[0]] > 1:
				# 	print splitLine[0]
				cosmicCountBadScore += 1
			if splitLine[0] not in breastCancerGenes:
				bcCountBadScore += 1
			allGenesBadScore += 1	
				
print cosmicCountGoodScore, " out of ", allGenesGoodScore, " are in COSMIC and score 0"
print cosmicCountBadScore, " out of ", allGenesBadScore, " are in COSMIC and score > 0"

print bcCountGoodScore, " out of ", allGenesGoodScore, " are known breast cancer genes and score 0"
print bcCountBadScore, " out of ", allGenesBadScore, " are known breast cancer genes and score > 0"
