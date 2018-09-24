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
		
			
toolRankingFile = "rankedGenes_test.txt"

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
				splitLine = line.split("\t")
				signCount = splitLine[4]
				#print "sign: ", interactionSignificance
				if signCount > 0:
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
#finalScoresFile = "rankedGenes_test.txt"

breastCancerGeneScores = dict()
bcGeneScores = dict()
bcTadScores = dict()
bceQTLScores = dict()

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
				
				
				bcGeneScores[gene] = geneScore
				
				
				bcTadScores[gene] = tadScore
				
				bceQTLScores[gene] = eQTLScore
				
				totalScore = geneScore + tadScore + eQTLScore
				
				
				breastCancerGeneScores[gene] = totalScore
			
			lineCount += 1
			
#print breastCancerGeneScores
print "known BC scores:"
#print sum(breastCancerGeneScores.values()) / float(len(breastCancerGeneScores))
#print sum(bcGeneScores.values()) / float(len(bcGeneScores))
#print sum(bcTadScores.values()) / float(len(bcTadScores))
#print sum(bceQTLScores.values()) / float(len(bceQTLScores))

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
cosmicGeneLayerScores = dict()
cosmicTadScores = dict()
cosmicEQtlScores = dict()

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
				
				
				cosmicGeneLayerScores[gene] = geneScore
				
				cosmicTadScores[gene] = tadScore
				
				cosmicEQtlScores[gene] = eQTLScore
				
				totalScore = geneScore + tadScore + eQTLScore
				
				cosmicGeneScores[gene] = totalScore
			
			lineCount += 1

#print cosmicGeneScores
print "Cosmic gene scores: "
#print sum(cosmicGeneScores.values()) / float(len(cosmicGeneScores))
#print sum(cosmicGeneLayerScores.values()) / float(len(cosmicGeneLayerScores))
#print sum(cosmicTadScores.values()) / float(len(cosmicTadScores))
#print sum(cosmicEQtlScores.values()) / float(len(cosmicEQtlScores))

#Repeat for the non-cosmic genes

nonCosmicGeneScores = dict()

ncGeneLayerScores = dict()
ncTadScores = dict()
nceQTLScores = dict()

with open(finalScoresFile, 'r') as f:
	
	for line in f:
		
		line = line.strip()
	
		splitLine = line.split("\t")
		geneName = splitLine[0]
		
		if geneName not in cosmicGeneNames and geneName not in genes: #check if it is not in cosmic or in the set of breast cancer genes
			
			geneScore = float(splitLine[1])
			tadScore = float(splitLine[2])
			eQTLScore = float(splitLine[3])
			
			
			
			ncGeneLayerScores[geneName] = geneScore
			
			
			ncTadScores[geneName] = tadScore
			
			nceQTLScores[geneName] = eQTLScore
			
			
			
			totalScore = geneScore + tadScore + eQTLScore
			
			
			nonCosmicGeneScores[geneName] = totalScore
			

#print nonCosmicGeneScores
print "Non cosmic gene sscores:"
#print sum(nonCosmicGeneScores.values()) / float(len(nonCosmicGeneScores))
#print sum(ncGeneLayerScores.values()) / float(len(ncGeneLayerScores))
#print sum(ncTadScores.values()) / float(len(ncTadScores))
#print sum(nceQTLScores.values()) / float(len(nceQTLScores))

#Make the boxplots, one with the average scores and one for the different layers

# 
# from pylab import plot, show, savefig, xlim, figure, \
#                 hold, ylim, legend, boxplot, setp, axes, hist
# # 
# # # function for setting the colors of the box plots pairs
# # def setBoxColors(bp):
# # 	setp(bp['boxes'][0], color='blue')
# # 	setp(bp['caps'][0], color='blue')
# # 	setp(bp['caps'][1], color='blue')
# # 	setp(bp['whiskers'][0], color='blue')
# # 	setp(bp['whiskers'][1], color='blue')
# # 	setp(bp['fliers'][0], color='blue')
# # 	setp(bp['fliers'][1], color='blue')
# # 	setp(bp['medians'][0], color='blue')
# # 
# # 	setp(bp['boxes'][1], color='red')
# # 	setp(bp['caps'][2], color='red')
# # 	setp(bp['caps'][3], color='red')
# # 	setp(bp['whiskers'][2], color='red')
# # 	setp(bp['whiskers'][3], color='red')
# # 	setp(bp['fliers'][2], color='red')
# # 	setp(bp['fliers'][3], color='red')
# # 	setp(bp['medians'][1], color='red')
# # 	
# # 	setp(bp['boxes'][1], color='orange')
# # 	setp(bp['caps'][2], color='orange')
# # 	setp(bp['caps'][3], color='orange')
# # 	setp(bp['whiskers'][2], color='orange')
# # 	setp(bp['whiskers'][3], color='orange')
# # 	setp(bp['fliers'][2], color='orange')
# # 	setp(bp['fliers'][3], color='orange')
# # 	setp(bp['medians'][1], color='orange')
# 
# # Some fake data to plot
# A= [bcGeneScores.values(), bcTadScores.values(), bceQTLScores.values()]
# B = [cosmicGeneLayerScores.values(), cosmicTadScores.values(), cosmicEQtlScores.values()]
# C = [ncGeneLayerScores.values(), ncTadScores.values(), nceQTLScores.values()]
# D = [breastCancerGeneScores.values(), cosmicGeneScores.values(), nonCosmicGeneScores.values()]
# 
# fig = figure()
# ax = axes()
# hold(True)
# 
# # first boxplot pair
# bp = boxplot(D, positions = [1, 2, 3], widths = 0.6)
# #setBoxColors(bp)
# 
# # second boxplot pair
# bp = boxplot(A, positions = [5, 6, 7], widths = 0.6)
# #setBoxColors(bp)
# 
# # thrid boxplot pair
# bp = boxplot(B, positions = [9, 10, 11], widths = 0.6)
# #setBoxColors(bp)
# 
# bp = boxplot(C, positions = [13, 14, 15], widths = 0.6)
# 
# # set axes limits and labels
# xlim(0,16)
# #ylim(0,9)
# ax.set_xticklabels(['A', 'B', 'C'])
# ax.set_xticks([1.5, 4.5, 7.5])
# 
# # draw temporary red and blue lines and use them to create a legend
# # hB, = plot([1,1],'b-')
# # hR, = plot([1,1],'r-')
# # legend((hB, hR),('Apples', 'Oranges'))
# # hB.set_visible(False)
# # hR.set_visible(False)
# 
# #savefig('boxcompare.png')
# show()
# exit()
# # 
# #Make the boxplots separately for the different layers.
# A= [bcGeneScores.values()]
# B = [cosmicGeneLayerScores.values()]
# C = [ncGeneLayerScores.values()]
# 
# 
# fig = figure()
# ax = axes()
# hold(True)
# 
# # first boxplot pair
# bp = boxplot(A, positions = [1], widths = 0.6)
# #setBoxColors(bp)
# 
# # second boxplot pair
# bp = boxplot(B, positions = [3], widths = 0.6)
# #setBoxColors(bp)
# 
# # thrid boxplot pair
# bp = boxplot(C, positions = [5], widths = 0.6)
# #setBoxColors(bp)
# 
# xlim(0,7)
# # #ylim(0,9)
# ax.set_xticklabels(['A', 'B', 'C'])
# ax.set_xticks([1,3,5])
# 
# show()

import numpy as np
import matplotlib.pyplot as plt
import numpy as np


#plt.show()



bcTotalScores = breastCancerGeneScores.values()
cosmicTotalScores = cosmicGeneScores.values()
nonCosmicTotalScores = nonCosmicGeneScores.values()

tadScores = [np.mean(bcTotalScores), np.mean(cosmicTotalScores), np.mean(nonCosmicTotalScores)]
tadScoresStdLower = [np.percentile(bcTotalScores, 25), np.percentile(cosmicTotalScores, 25), np.percentile(nonCosmicTotalScores, 25)]
tadScoresStdUpper = [np.percentile(bcTotalScores, 75), np.percentile(cosmicTotalScores, 75), np.percentile(nonCosmicTotalScores, 75)]

plt.errorbar([1,1.5,2], tadScores, [tadScoresStdLower, tadScoresStdUpper], linestyle='None', marker='*')
plt.xlim(0.5,2.5)
plt.xticks([1,1.5,2], ('Known breast cancer genes', 'COSMIC genes', 'Non-COSMIC genes'), rotation='vertical')
plt.ylabel("Gene causality score")
#plt.show()
plt.tight_layout()
plt.savefig("total_scores.svg")

plt.clf()
# 
tadScores = [np.mean(bcTadScores.values()), np.mean(cosmicTadScores.values()), np.mean(ncTadScores.values())]
tadScoresStdLower = [np.percentile(bcTadScores.values(), 25), np.percentile(cosmicTadScores.values(), 25), np.percentile(ncTadScores.values(), 25)]
tadScoresStdUpper = [np.percentile(bcTadScores.values(), 75), np.percentile(cosmicTadScores.values(), 75), np.percentile(ncTadScores.values(), 75)]

plt.errorbar([1,1.5,2], tadScores, [tadScoresStdLower, tadScoresStdUpper], linestyle='None', marker='^')
plt.xlim(0.5,2.5)
plt.xticks([1,1.5,2], ('Known breast cancer genes', 'COSMIC genes', 'Non-COSMIC genes'), rotation='vertical')
plt.ylabel("Gene causality score")
#plt.show()
plt.tight_layout()
plt.savefig("tadScores.svg")
plt.clf()
# 
tadScores = [np.mean(bceQTLScores.values()), np.mean(cosmicEQtlScores.values()), np.mean(nceQTLScores.values())]
tadScoresStdLower = [np.percentile(bceQTLScores.values(), 25), np.percentile(cosmicEQtlScores.values(), 25), np.percentile(nceQTLScores.values(), 25)]
tadScoresStdUpper = [np.percentile(bceQTLScores.values(), 75), np.percentile(cosmicEQtlScores.values(), 75), np.percentile(nceQTLScores.values(), 75)]

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