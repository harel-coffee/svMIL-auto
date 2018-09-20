"""
	Given a list of breast cancer genes and the ranking outcome of our tool, compute how many of these are significant? And what are the ranks of these genes? 
"""

import re

#1. Compute how many breast cancer genes are significant in the results

knownBreastCancerGenesFile = "../../data/Genes/breastCancerCausalGenes.txt"
genes = []

with open(knownBreastCancerGenesFile, 'r') as f:
	
	for line in f:
		
		line = line.strip()
		
		genes.append(line)
		
			
toolRankingFile = "RankedGenes/breastCancerRanking_diff0.8.txt"

geneRank = dict()
significance = dict()

totalSign = 0
for gene in genes:
	
	lineCount = 0
	with open(toolRankingFile, 'r') as f:
		
		for line in f:
			
			line = line.strip()
			
			
			
			
			if re.search(gene, line):
				
				geneRank[gene] = lineCount
				
			
				
				#obtain the p-value and see if it is significant
				splitLine = line.split(" ")
				interactionSignificance = float(splitLine[6])
				#print "sign: ", interactionSignificance
				if interactionSignificance < 0.01:
				#	print "True"
					significance[gene] = 'True'
					totalSign += 1
				else:
				#	print "False"
					significance[gene] = 'False'
			
			lineCount += 1
			
			
print geneRank
print significance

print totalSign / float(len(genes))

#Get the scores of the breast cancer genes

finalScoresFile = "RankedGenes/0/breast/realSVs_geneScores.txt"

breastCancerGeneScores = dict()

for gene in genes:
	
	lineCount = 0
	with open(finalScoresFile, 'r') as f:
		
		for line in f:
			
			line = line.strip()

			if re.search(gene, line):
				
				#obtain the p-value and see if it is significant
				
				splitLine = line.split("\t")
				geneScore = float(splitLine[1])
				tadScore = float(splitLine[2])
				eQTLScore = float(splitLine[3])
				totalScore = geneScore + tadScore + eQTLScore
				breastCancerGeneScores[gene] = totalScore	
			
			lineCount += 1
			
#print breastCancerGeneScores
print sum(breastCancerGeneScores.values())

#Repeat but then for the cosmic genes

cosmicGenesFile = "../../data/Genes/Census_allTue Apr 10 14_56_44 2018.tsv"
cosmicGeneNames = []
with open(cosmicGenesFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		if lineCount < 1:
			lineCount += 1
			continue
		
		splitLine = line.split("\t")
		
		geneName = splitLine[0]
		if geneName not in genes: #don't include the breast cancer genes here
			cosmicGeneNames.append(geneName)

cosmicGeneScores = dict()

for gene in cosmicGeneNames:
	
	lineCount = 0
	with open(finalScoresFile, 'r') as f:
		
		for line in f:
			
			line = line.strip()

			if re.search(gene, line):
				
				#obtain the p-value and see if it is significant
				
				splitLine = line.split("\t")
				geneScore = float(splitLine[1])
				tadScore = float(splitLine[2])
				eQTLScore = float(splitLine[3])
				totalScore = geneScore + tadScore + eQTLScore
				cosmicGeneScores[gene] = totalScore	
			
			lineCount += 1

#print cosmicGeneScores
print sum(cosmicGeneScores.values())

#Repeat for the non-cosmic genes

nonCosmicGeneScores = dict()
with open(finalScoresFile, 'r') as f:
	
	for line in f:
		
		line = line.strip()
	
		splitLine = line.split("\t")
		geneName = splitLine[0]
		
		if geneName not in cosmicGeneNames and geneName not in genes: #check if it is not in cosmic or in the set of breast cancer genes
			
			geneScore = float(splitLine[1])
			tadScore = float(splitLine[2])
			eQTLScore = float(splitLine[3])
			totalScore = geneScore + tadScore + eQTLScore
			nonCosmicGeneScores[geneName] = totalScore
			

#print nonCosmicGeneScores
print sum(nonCosmicGeneScores.values())
