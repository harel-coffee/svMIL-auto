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




			
# toolRankingFile = "rankedGenes_test.txt"
# 
# geneRank = dict()
# significance = dict()
# 
# totalSign = 0
# for gene in genes:
# 	
# 	lineCount = 0
# 	with open(toolRankingFile, 'r') as f:
# 		
# 		for line in f:
# 			
# 			line = line.strip()
# 			
# 			
# 			
# 			
# 			if re.search(gene, line):
# 				
# 				geneRank[gene] = lineCount
# 				
# 			
# 				
# 				#obtain the p-value and see if it is significant
# 				splitLine = line.split("\t")
# 				signCount = float(splitLine[3])
# 				#print "sign: ", interactionSignificance
# 				if signCount < 0.05:
# 				#	print "True"
# 					significance[gene] = 'True'
# 					totalSign += 1
# 				else:
# 				#	print "False"
# 					significance[gene] = 'False'
# 			
# 			lineCount += 1
# 			
# 			
# print geneRank
# print significance
# 
# print totalSign / float(len(genes))

#Get the scores of the breast cancer genes

finalScoresFile = "RankedGenes/0/breast/realSVs_geneScores.txt"
#finalScoresFile = "rankedGenes_test.txt"

breastCancerGeneScores = dict()
bcGeneScores = dict()
bcTadScores = dict()
bceQTLScores = dict()
bceQTLGainScores = dict()

bcNonZeroAlpha = 0
bcNonZeroBeta = 0
bcNonZeroAlphaZeroBeta = 0
bcNonZeroBetaZeroAlpha = 0

for gene in genes:
	
	lineCount = 0
	with open(finalScoresFile, 'r') as f:
		
		for line in f:
			
			line = line.strip()
			splitLine = line.split("\t")

			if splitLine[0] == gene:
				
				#obtain the p-value and see if it is significant

				eQTLScore = float(splitLine[3])
				print "gene: ", gene
				print "current eQTLScore: ", eQTLScore
			
				bceQTLScores[gene] = eQTLScore
				
				breastCancerGeneScores[gene] = eQTLScore
			
			lineCount += 1
			
#print breastCancerGeneScores
print "known BC scores:"
#print sum(breastCancerGeneScores.values()) / float(len(breastCancerGeneScores))
#print sum(bcGeneScores.values()) / float(len(bcGeneScores))
#print sum(bcTadScores.values()) / float(len(bcTadScores))
#print sum(bceQTLScores.values()) / float(len(bceQTLScores))

#Count how many breast cancer genes have a score higher than 0.

print bceQTLScores

higherThanZeroCount = 0
for gene in bceQTLScores:
	geneScore = bceQTLScores[gene]
	if geneScore > 0:
		higherThanZeroCount += 1
		
print "Number of known breast cancer genes with a score higher than 0: ", higherThanZeroCount

print "Number of non-zero alphas for known breast cancer genes: ", bcNonZeroAlpha
print "Number of non-zero betas for known breast cancer genes: ", bcNonZeroBeta

print "Number of breast cancer genes with alpha = 0 and beta > 0: ", bcNonZeroBetaZeroAlpha
print "Number of breast cancer genes with alpha > 0 and beta == 0: ", bcNonZeroAlphaZeroBeta


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

cosmicEQtlScores = dict()
cosmicGeneScores = dict()

cosmicNonZeroAlpha = 0
cosmicNonZeroBeta = 0
cosmicNonZeroAlphaZeroBeta = 0
cosmicNonZeroBetaZeroAlpha = 0

higherThanZeroCount = 0
for gene in cosmicGeneNames:
	
	lineCount = 0
	with open(finalScoresFile, 'r') as f:
		
		for line in f:
			
			line = line.strip()
			splitLine = line.split("\t")

			if splitLine[0] == gene:

				
				#obtain the p-value and see if it is significant
				
				eQTLScore = float(splitLine[3])
				
				cosmicGeneScores[gene] = eQTLScore
				
			
			lineCount += 1




#Count how many cosmic cancer genes have a score higher than 0.

print "Number of cosmic cancer genes with a score higher than 0: ", higherThanZeroCount

print "Number of non-zero alphas for cosmic cancer genes: ", cosmicNonZeroAlpha
print "Number of non-zero betas for cosmic cancer genes: ", cosmicNonZeroBeta
print "Number of cosmic cancer genes with alpha > 0 and beta == 0: ", cosmicNonZeroAlphaZeroBeta
print "Number of cosmic cancer genes with alpha == 0 and beta > 0: ", cosmicNonZeroBetaZeroAlpha

#Repeat for the non-cosmic genes

nonCosmicGeneScores = dict()
nceQTLScores = dict()

nonCosmicNonZeroAlpha = 0
nonCosmicNonZeroBeta = 0
nonCosmicNonZeroAlphaZeroBeta = 0
nonCosmicNonZeroBetaZeroAlpha = 0

higherThanZeroCount = 0
with open(finalScoresFile, 'r') as f:
	
	for line in f:
		
		line = line.strip()
	
		splitLine = line.split("\t")
		geneName = splitLine[0]
		
		if geneName not in cosmicGeneNames and geneName not in genes: #check if it is not in cosmic or in the set of breast cancer genes
			
			eQTLScore = float(splitLine[3])
			
			nceQTLScores[geneName] = eQTLScore
			
			nonCosmicGeneScores[geneName] = eQTLScore
			

#print nonCosmicGeneScores
print "Non cosmic gene sscores"


#Count how many non-cosmic  genes have a score higher than 0.

		
print "Number of non-cosmic  genes with a score higher than 0: ", higherThanZeroCount

print "Number of non-zero alphas for non cancer genes: ", nonCosmicNonZeroAlpha
print "Number of non-zero betas for non cancer genes: ", nonCosmicNonZeroBeta
print "Number of non cosmic genes with alpha >0 and beta == 0: ", nonCosmicNonZeroAlphaZeroBeta
print "Number of non cosmic genes with alpha == 0 and beta > 0: ", nonCosmicNonZeroBetaZeroAlpha


import numpy as np
import matplotlib.pyplot as plt
import numpy as np



tadScores = [np.mean(breastCancerGeneScores.values()), np.mean(cosmicGeneScores.values()), np.mean(nonCosmicGeneScores.values())]
tadScoresStdLower = [np.percentile(breastCancerGeneScores.values(), 25), np.percentile(cosmicGeneScores.values(), 25), np.percentile(nonCosmicGeneScores.values(), 25)]
tadScoresStdUpper = [np.percentile(breastCancerGeneScores.values(), 75), np.percentile(cosmicGeneScores.values(), 75), np.percentile(nonCosmicGeneScores.values(), 75)]

plt.errorbar([1,1.5,2], tadScores, [tadScoresStdLower, tadScoresStdUpper], linestyle='None', marker='*')
plt.xlim(0.5,2.5)
plt.xticks([1,1.5,2], ('Known breast cancer genes', 'COSMIC genes', 'Non-COSMIC genes'), rotation='vertical')
plt.ylabel("Gene causality score")
plt.tight_layout()
#plt.show()
plt.savefig("eQTLScores.svg")

exit()
#print A
#print B
#print C



fig = figure()
ax = axes()
hold(True)

# first boxplot pair
bp = boxplot(A, positions = [1], widths = 0.6)
#setBoxColors(bp)

# second boxplot pair
bp = boxplot(B, positions = [3], widths = 0.6)
#setBoxColors(bp)

# thrid boxplot pair
bp = boxplot(C, positions = [5], widths = 0.6)
#setBoxColors(bp)

xlim(0,7)
# #ylim(0,9)
ax.set_xticklabels(['A', 'B', 'C'])
ax.set_xticks([1,3,5])

show()