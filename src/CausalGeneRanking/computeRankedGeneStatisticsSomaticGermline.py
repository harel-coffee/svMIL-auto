"""

	Compute ranked gene statistics but then for 2 groups, also showing the overlap between somatic and germline and their overlap with the groups. 
"""

import sys


somaticFile = sys.argv[1]
germlineFile = sys.argv[2]
cosmicGenesFile = sys.argv[3]
breastCancerGenesFile = sys.argv[4]
degFile = sys.argv[5]
snvFile = sys.argv[6] #file with genes that have SNVs
degGt3File = sys.argv[7]

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

#In addition, test how many of the genes also have SNVs affecting them in at least 1 sample
snvGenes = []
with open(snvFile, 'r') as inf:
	for line in inf:
		line = line.strip()
		snvGenes.append(line)
	
#Then, determine how many of the genes are differentially expressed
degGenes = []
with open(degFile, 'r') as inf:
	for line in inf:
		line = line.strip()
		degGenes.append(line)

degGt3Genes = []
with open(degGt3File, 'r') as inf:
	for line in inf:
		line = line.strip()
		splitLine = line.split("\t")
		degGt3Genes.append(splitLine[0])

#Make re-usable functions that give the same result for somatic and germline

#Compute the category count for a given set of scores

def getSomaticSVSubset(rankedGenesFile, sampleThreshold):
	
	somaticSubset = []
	svCount = 0
	with open(rankedGenesFile, 'r') as f:
		lineCount = 0
		for line in f:
			line = line.strip()
			splitLine = line.split("\t")
			
			if lineCount < 1:
				lineCount += 1
				continue

			if float(splitLine[30]) > 0:
				
				samples = splitLine[31]
				splitSamples = samples.split(",")
				if len(splitSamples) < sampleThreshold:
					continue
				svCount += len(splitSamples)
				somaticSubset.append(splitLine)
	return somaticSubset, svCount

def getGermlineSVSubset(rankedGenesFile, svCount):
	
	print "Somatic SV count: ", svCount
	
	germlineSubset = []
	
	with open(rankedGenesFile, 'r') as f:
		lineCount = 0
		for line in f:
			line = line.strip()
			splitLine = line.split("\t")
			
			if lineCount < 1:
				lineCount += 1
				continue

			if float(splitLine[30]) > 0:
				# samples = splitLine[31]
				# splitSamples = samples.split(",")
				# if len(splitSamples) < sampleThreshold:
				# 	continue
				
				germlineSubset.append(splitLine)
				
	#Select random SVs in the same size as the somatic set size
	#Randomly go through the genes in the germline subset and keep adding genes until we have exactly the same number of SVs (or approximately)
	import random
	random.shuffle(germlineSubset)
	
	germlineSVCount = 0
	germlineSVSubset = []
	for score in germlineSubset:
	
		samples = score[31]
		splitSamples = samples.split(",")
	
		if germlineSVCount + len(splitSamples) > svCount:
			continue
		
		
		germlineSVSubset.append(score)
		germlineSVCount += len(splitSamples)
	
	print "Germline SV count: ", svCount
	return germlineSVSubset					 
	#return germlineSubset
	

def computeCategoryOverlap(scores):
		
	
	#Finally, show how many of the genes have all of the above, these are likely the most interesting genes. 
	
	degGenesPos = []
	degGenesNeg = []
	cosmicGenesPos = []
	cosmicGenesNeg = []
	snvGenesPos = []
	snvGenesNeg = []
	allGenesPos = []
	allGenesNeg = []
	degGt3GenesPos = []
	degGt3GenesNeg = []
	
	for splitLine in scores:

		allGenesPos.append(splitLine[0])

		if splitLine[0] in degGenes:
			#print "deg: ", splitLine[0]
			degGenesPos.append(splitLine[0])
		if splitLine[0] in cosmicGenes:
			cosmicGenesPos.append(splitLine[0])
		if splitLine[0] in snvGenes:
			snvGenesPos.append(splitLine[0])
		if splitLine[0] in degGt3Genes:
			degGt3GenesPos.append(splitLine[0])
	else:
		allGenesNeg.append(splitLine[0])
		if splitLine[0] in degGenes:
			degGenesNeg.append(splitLine[0])
		if splitLine[0] in cosmicGenes:
			cosmicGenesNeg.append(splitLine[0])
		if splitLine[0] in snvGenes:
			snvGenesNeg.append(splitLine[0])
		if splitLine[0] in degGt3Genes:
			degGt3GenesNeg.append(splitLine[0])
				
	return allGenesPos, allGenesNeg





# #Make a venn diagram for the intersect between the genes of somatic & germline, and their overlaps with the deg, cosmic and snv genes.  
# import matplotlib.pyplot as plt
# #matplotlib.use('Agg')
# 
# import venn
# 
# 
# labels = venn.get_labels([allGenesPosS, allGenesPosG, cosmicGenes, snvGenes, degGenes], fill=['number'])
# fig, ax = venn.venn5(labels, names=['Somatic', 'Germline', 'COSMIC', 'DEG all patients', 'SNV'])
# fig.show()
# plt.show()

#Do chi-square tests for each category combination
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
from scipy.stats import chi2
import numpy as np
# contingency table

#Compute how many genes intersect between somatic & cosmic and germline & cosmic
#make this a function
def getCategoryOverlapCounts(somaticGenes, germlineGenes, categoryGenes):
	
	#1. How many genes are in somatic & in the category
	somaticCategoryIntersect = list(set(somaticGenes) & set(categoryGenes))
	
	#2. How many genes are in germline & in the category
	germlineCategoryIntersect = list(set(germlineGenes) & set(categoryGenes))
	
	#3. How many genes are in somatic but not in the category
	somaticCategoryDifference = list(np.setdiff1d(somaticGenes, categoryGenes))
	
	#4. How many genes are in germline but not in the category
	germlineCategoryDifference = list(np.setdiff1d(germlineGenes, categoryGenes))
	
	#return contingency table
	table = [[len(somaticCategoryIntersect), len(germlineCategoryIntersect)], [len(somaticCategoryDifference), len(germlineCategoryDifference)]]
	table = [[len(somaticCategoryIntersect), len(somaticCategoryDifference)], [len(germlineCategoryIntersect), len(germlineCategoryDifference)]]
	return table

def computeSignificance(somaticGenes, germlineGenes, categoryGenes):
	table = getCategoryOverlapCounts(somaticGenes, germlineGenes, categoryGenes)
	print(table)
	#stat, p, dof, expected = chi2_contingency(table)
	odds, p = fisher_exact(table)
	print p
	# print('dof=%d' % dof)
	# print(expected)
	# # interpret p-value
	# alpha = 0.05
	# print('significance=%.3f, p=%.3f' % (alpha, p))
	# if p <= alpha:
	# 	print('Significant')
	# else:
	# 	print('Not significant')
	return p

#Repeat these tests at different thresholds. We need to make different sets based on the total number of samples
pValues = np.zeros([18,7])
maxThreshold = 18 #there are no genes with more samples than this
for threshold in range(0, maxThreshold):
	print "Threshold: ", threshold
	
	somaticSVs, svCount = getSomaticSVSubset(somaticFile, threshold)
	germlineSVs = getGermlineSVSubset(germlineFile, svCount)
	
	
	allGenesPosS, allGenesNegS = computeCategoryOverlap(somaticSVs)
	allGenesPosG, allGenesNegG = computeCategoryOverlap(germlineSVs)
	
	print "No of somatic genes: ", len(allGenesPosS)
	print "No of germline genes: ", len(allGenesPosG)
		
	print "COSMIC:"
	cosmic = computeSignificance(allGenesPosS, allGenesPosG, cosmicGenes)
	print "DEG:"
	deg = computeSignificance(allGenesPosS, allGenesPosG, degGenes)
	print "SNV:"
	snv = computeSignificance(allGenesPosS, allGenesPosG, snvGenes)
	
	print "COSMIC+DEG:"
	cosmicDeg = computeSignificance(allGenesPosS, allGenesPosG, list(set(cosmicGenes) & set(degGenes)))
	print "COSMIC+SNVs:"
	cosmicSnv = computeSignificance(allGenesPosS, allGenesPosG, list(set(cosmicGenes) & set(snvGenes)))
	print "SNVs+DEGs:"
	snvDeg = computeSignificance(allGenesPosS, allGenesPosG, list(set(degGenes) & set(snvGenes)))
	
	print "COSMIC+DEGs+SNVs:"
	cosmicSnvDeg = computeSignificance(allGenesPosS, allGenesPosG, list(set(cosmicGenes) & set(degGenes) & set(snvGenes)))											
	
	pValues[threshold, 0] = cosmic
	pValues[threshold, 1] = deg
	pValues[threshold, 2] = snv
	pValues[threshold, 3] = cosmicDeg
	pValues[threshold, 4] = cosmicSnv
	pValues[threshold, 5] = snvDeg
	pValues[threshold, 6] = cosmicSnvDeg
	
	table = getCategoryOverlapCounts(allGenesPosS, allGenesNegS, degGt3Genes)
	print "DEGs for current threshold: "
	print table
	

print pValues


print pValues.flatten()
flattenedPValues = list(pValues.flatten())

from statsmodels.sandbox.stats.multicomp import multipletests
reject, pAdjusted, _, _ = multipletests(flattenedPValues, method='fdr_bh') #fdr_bh

print pAdjusted
print reject
#reshape the pvalues
pAdjustedReshaped = np.reshape(pAdjusted, pValues.shape)
rejectReshaped = np.reshape(reject, pValues.shape)

print pAdjustedReshaped[rejectReshaped]

#Output a table to a file
np.savetxt("categoryOverlap_sign.txt", pAdjustedReshaped, delimiter='\t', fmt='%s')



#Make a venn diagram for the intersect between the genes of somatic & germline, and their overlaps with the deg, cosmic and snv genes.  
import matplotlib.pyplot as plt
#matplotlib.use('Agg')

import venn

somaticSVs, svCount = getSomaticSVSubset(somaticFile, 0)
germlineSVs = getGermlineSVSubset(germlineFile, svCount)

allGenesPosS, allGenesNegS = computeCategoryOverlap(somaticSVs)
allGenesPosG, allGenesNegG = computeCategoryOverlap(germlineSVs)

labels = venn.get_labels([allGenesPosS, allGenesPosG], fill=['number'])
fig, ax = venn.venn2(labels, names=['Somatic', 'Germline'])
fig.show()
plt.show()


# labels = venn.get_labels([allGenesPosS, allGenesPosG, cosmicGenes, snvGenes, degGenes], fill=['number'])
# fig, ax = venn.venn5(labels, names=['Somatic', 'Germline', 'COSMIC', 'DEG all patients', 'SNV'])
# fig.show()
# plt.show()

exit()



#Do some intersect things here

allCriteriaIntersect = list(set(degGenesPos) & set(cosmicGenesPos) & set(snvGenesPos))
cosmicSNVsIntersect = list(set(cosmicGenesPos) & set(snvGenesPos))
cosmicDEGsIntersect = list(set(cosmicGenesPos) & set(degGenesPos))
snvDEGsIntersect = list(set(snvGenesPos) & set(degGenesPos))
print "Number of genes that are in COSMIC, have SNVs and are DEG: ", len(allCriteriaIntersect)
print "Number of genes that are in COSMIC and have SNVs: ", len(cosmicSNVsIntersect)
print "Number of genes that are in COSMIC and are DEG: ", len(cosmicDEGsIntersect)
print "Number of genes that have SNV and are DEG: ", len(snvDEGsIntersect)


print "genes in the total intersect: "
print allCriteriaIntersect
intersectScores = []
with open(rankedGenesFile, 'rb') as f:
	lineCount = 0
	for line in f:
		line = line.strip()
		splitLine = line.split("\t")
		
		if lineCount < 1:
			lineCount += 1
			continue
		
		if float(splitLine[30]) > 0:
		#if float(splitLine[1]) < 0.05:
			if splitLine[0] in allCriteriaIntersect:
				intersectScores.append(float(splitLine[30]))
import matplotlib.pyplot as plt
plt.hist(intersectScores)
plt.show()
#exit()

#exit()

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
