"""
	Given a list of breast cancer genes and the ranking outcome of our tool, compute how many of these are significant? And what are the ranks of these genes? 
"""

import re

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
				if gene == 'MAP3K1':
					print "MAP3K1"
					print lineCount
				if gene == "MAP3K13":
					print "MAP3K13"
					print lineCount
				
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