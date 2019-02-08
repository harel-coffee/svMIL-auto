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
degFile = sys.argv[4]
snvFile = sys.argv[5] #file with genes that have SNVs

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
		
		
		if float(splitLine[1]) + float(splitLine[4]) > 0:
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
			if splitLine[0] in breastCancerGenes:
				bcCountBadScore += 1
			allGenesBadScore += 1	
				
print cosmicCountGoodScore, " out of ", allGenesGoodScore, " are in COSMIC and score > 0"
print cosmicCountBadScore, " out of ", allGenesBadScore, " are in COSMIC and score 0"

print bcCountGoodScore, " out of ", allGenesGoodScore, " are known breast cancer genes and score > 0"
print bcCountBadScore, " out of ", allGenesBadScore, " are known breast cancer genes and score 0"

#In addition, test how many of the genes also have SNVs affecting them in at least 1 sample
snvGenes = []
with open(snvFile, 'r') as inf:
	for line in inf:
		line = line.strip()
		snvGenes.append(line)

#Positive for the genes with score > 0, negative for others. 		
snvCountPos = 0
nonSnvCountPos = 0
snvCountNeg = 0
nonSnvCountNeg = 0

with open(rankedGenesFile, 'rb') as f:
	
	for line in f:
		line = line.strip()
		splitLine = line.split("\t")
		
		if float(splitLine[1]) + float(splitLine[4]) > 0:
			if splitLine[0] in snvGenes:
				snvCountPos += 1
			else:
				nonSnvCountPos += 1
		else:
			if splitLine[0] in snvGenes:
				snvCountNeg += 1
			else:
				nonSnvCountNeg += 1

			
print snvCountPos, " out of ", snvCountPos + nonSnvCountPos, " score > 0 and have SNVs"
print snvCountNeg, " out of ", snvCountNeg + nonSnvCountNeg, " score 0 and have SNVs"


#Then, determine how many of the genes are differentially expressed
degGenes = []
with open(degFile, 'r') as inf:
	for line in inf:
		line = line.strip()
		degGenes.append(line)

#Positive for the genes with score > 0, negative for others. 		
degCountPos = 0
nonDegCountPos = 0
degCountNeg = 0
nonDegCountNeg = 0

with open(rankedGenesFile, 'rb') as f:
	
	for line in f:
		line = line.strip()
		splitLine = line.split("\t")
		
		
		if float(splitLine[1]) + float(splitLine[4]) > 0:
			if splitLine[0] in degGenes:
				degCountPos += 1
			else:
				nonDegCountPos += 1
		else:
			if splitLine[0] in degGenes:
				degCountNeg += 1
			else:
				nonDegCountNeg += 1

			
print degCountPos, " out of ", degCountPos + nonDegCountPos, " score > 0 and are DEG"
print degCountNeg, " out of ", degCountNeg + nonDegCountNeg, " score 0 and are DEG"


#Finally, show how many of the genes have all of the above, these are likely the most interesting genes. 

degGenesPos = []
degGenesNeg = []
cosmicGenesPos = []
cosmicGenesNeg = []
snvGenesPos = []
snvGenesNeg = []

with open(rankedGenesFile, 'rb') as f:
	
	for line in f:
		line = line.strip()
		splitLine = line.split("\t")
		
		if float(splitLine[1]) + float(splitLine[4]) > 0:
			if splitLine[0] in degGenes:
				degGenesPos.append(splitLine[0])
			if splitLine[0] in cosmicGenes:
				cosmicGenesPos.append(splitLine[0])
			if splitLine[0] in snvGenes:
				snvGenesPos.append(splitLine[0])
		else:
			if splitLine[0] in degGenes:
				degGenesNeg.append(splitLine[0])
			if splitLine[0] in cosmicGenes:
				cosmicGenesNeg.append(splitLine[0])
			if splitLine[0] in snvGenes:
				snvGenesNeg.append(splitLine[0])
				
#Do some intersect things here

allCriteriaIntersect = list(set(degGenesPos) & set(cosmicGenesPos) & set(snvGenesPos))
cosmicSNVsIntersect = list(set(cosmicGenesPos) & set(snvGenesPos))
cosmicDEGsIntersect = list(set(cosmicGenesPos) & set(degGenesPos))
snvDEGsIntersect = list(set(snvGenesPos) & set(degGenesPos))
print "Number of genes that are in COSMIC, have SNVs and are DEG: ", len(allCriteriaIntersect)
print "Number of genes that are in COSMIC and have SNVs: ", len(cosmicSNVsIntersect)
print "Number of genes that are in COSMIC and are DEG: ", len(cosmicDEGsIntersect)
print "Number of genes that have SNV and are DEG: ", len(snvDEGsIntersect)

import pylab as plt
from matplotlib_venn import venn3, venn3_circles

#3 categories, genes in COSMIC, genes with SNVs, genes that are DEG

v = venn3(subsets=(len(cosmicGenesPos),len(snvGenesPos),len(cosmicSNVsIntersect), len(degGenesPos),len(cosmicDEGsIntersect), len(snvDEGsIntersect),len(allCriteriaIntersect)),
		  set_labels=('COSMIC', 'SNVs', 'DEGs'))
v.get_label_by_id('100').set_text(len(cosmicGenesPos))
v.get_label_by_id('010').set_text(len(snvGenesPos))
v.get_label_by_id('001').set_text(len(degGenesPos))
plt.title("")
plt.show()


allCriteriaIntersect = list(set(degGenesNeg) & set(cosmicGenesNeg) & set(snvGenesNeg))
cosmicSNVsIntersect = list(set(cosmicGenesNeg) & set(snvGenesNeg))
cosmicDEGsIntersect = list(set(cosmicGenesNeg) & set(degGenesNeg))
snvDEGsIntersect = list(set(snvGenesNeg) & set(degGenesNeg))
print "Negative set:"
print "Number of genes that are in COSMIC, have SNVs and are DEG: ", len(allCriteriaIntersect)
print "Number of genes that are in COSMIC and have SNVs: ", len(cosmicSNVsIntersect)
print "Number of genes that are in COSMIC and are DEG: ", len(cosmicDEGsIntersect)
print "Number of genes that have SNV and are DEG: ", len(snvDEGsIntersect)


v = venn3(subsets=(len(cosmicGenesNeg),len(snvGenesNeg),len(cosmicSNVsIntersect), len(degGenesNeg),len(cosmicDEGsIntersect), len(snvDEGsIntersect),len(allCriteriaIntersect)),
		  set_labels=('COSMIC', 'SNVs', 'DEGs'))
v.get_label_by_id('100').set_text(len(cosmicGenesNeg))
v.get_label_by_id('010').set_text(len(snvGenesNeg))
v.get_label_by_id('001').set_text(len(degGenesNeg))
plt.title("")
plt.show()
