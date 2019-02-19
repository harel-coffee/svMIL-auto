
#1. Read the data from the folder with real and permuted intersect values
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join

permutationDataFolder = sys.argv[1]
realScoreCountsCosmic = dict() #keep per threshold what the scores are
realScoreCountsSNVs = dict()
realScoreCountsDEGS = dict()
realScoreCountsCosmicSNVs = dict()
realScoreCountsCosmicDEGs = dict()
realScoreCountsSNVDEGs = dict()
realScoreCountsAll = dict()
maxScore = 0
geneScoreFiles = [f for f in listdir(permutationDataFolder) if isfile(join(permutationDataFolder, f))]
print geneScoreFiles
for geneScoreFile in geneScoreFiles:
	
	if geneScoreFile == 'realSVs_geneScores.txt':
		#Read the file
		geneScores = np.loadtxt(permutationDataFolder + "/" + geneScoreFile, dtype="object")
		
		#Go through the same thresholds
		for row in range(0, geneScores.shape[0]):
			
			threshold = int(geneScores[row][0])
			cosmic = int(geneScores[row][1])
			snvs = int(geneScores[row][2])
			cosmicSNVs = int(geneScores[row][3])
			degs = int(geneScores[row][4])
			cosmicDEGs = int(geneScores[row][5])
			snvDEGs = int(geneScores[row][6])
			allCriteria = int(geneScores[row][7])
			#[cosmic, snvs, cosmicSNVs, degs, cosmicDEGs, snvDEGs, allCriteria] = computeCategoryMatches(geneScores, threshold)
			if threshold > maxScore:
				maxScore = threshold
	
			realScoreCountsCosmic[threshold] = cosmic
			realScoreCountsSNVs[threshold] = snvs
			realScoreCountsDEGS[threshold] = degs
			realScoreCountsCosmicSNVs[threshold] = cosmicSNVs
			realScoreCountsCosmicDEGs[threshold] = cosmicDEGs
			realScoreCountsSNVDEGs[threshold] = snvDEGs
			realScoreCountsAll[threshold] = allCriteria
	


permutedScoreCountsCosmic = dict() #keep per threshold what the scores are
permutedScoreCountsSNVs = dict()
permutedScoreCountsDEGS = dict()
permutedScoreCountsCosmicSNVs = dict()
permutedScoreCountsCosmicDEGs = dict()
permutedScoreCountsSNVDEGs = dict()
permutedScoreCountsAll = dict()
for threshold in range(0, maxScore): #prepare dictionaries
	permutedScoreCountsCosmic[threshold] = []
	permutedScoreCountsSNVs[threshold] = []
	permutedScoreCountsDEGS[threshold] = []
	permutedScoreCountsCosmicSNVs[threshold] = []
	permutedScoreCountsCosmicDEGs[threshold] = []
	permutedScoreCountsSNVDEGs[threshold] = []
	permutedScoreCountsAll[threshold] = []



for geneScoreFile in geneScoreFiles:
	
	if geneScoreFile == 'realSVs_geneScores.txt':
		continue
	else:
		#Read the file
		geneScores = np.loadtxt(permutationDataFolder + "/" + geneScoreFile, dtype="object")
		
		#Go through the same thresholds
		for row in range(0, geneScores.shape[0]):
			print row
			exit()
			
			if len(geneScores[row]) > 8: #if there is something wird in the file
				continue
			
			threshold = int(geneScores[row][0])
			if threshold >= maxScore:
				continue
			cosmic = int(geneScores[row][1])
			snvs = int(geneScores[row][2])
			cosmicSNVs = int(geneScores[row][3])
			degs = int(geneScores[row][4])
			cosmicDEGs = int(geneScores[row][5])
			snvDEGs = int(geneScores[row][6])
			allCriteria = int(geneScores[row][7])
			#[cosmic, snvs, cosmicSNVs, degs, cosmicDEGs, snvDEGs, allCriteria] = computeCategoryMatches(geneScores, threshold)
			
	
			permutedScoreCountsCosmic[threshold].append(cosmic)
			permutedScoreCountsSNVs[threshold].append(snvs)
			permutedScoreCountsDEGS[threshold].append(degs)
			permutedScoreCountsCosmicSNVs[threshold].append(cosmicSNVs)
			permutedScoreCountsCosmicDEGs[threshold].append(cosmicDEGs)
			permutedScoreCountsSNVDEGs[threshold].append(snvDEGs)
			permutedScoreCountsAll[threshold].append(allCriteria)
			#if threshold not in permutedScoreCounts:
			#	permutedScoreCounts[threshold] = []	
			#permutedScoreCounts[threshold].append([cosmic, snvs, cosmicSNVs, degs, cosmicDEGs, snvDEGs, allCriteria])

#2. Make the plots

	
#Add the points to the plots, make seprate plots
#The y axis is the number in the overlap, the x axis the threshold
#All intersection values are individual points in the plot
def plotData(realScores, permutedScores, maxScore):
	plt.clf()
	ax = plt.subplot(1,1,1)
	for threshold in range(0, maxScore):
		#Plot the real score
		ax.plot(threshold, realScores[threshold], 'bo')
		
		#Take mean and std of permuted scores and plot
		#ax.plot(threshold, np.mean(permutedScores[threshold]), 'ko')
		thrSum = np.sum(permutedScores[threshold])
		
		thrMean = 0
		if thrSum > 0:
			thrMean = np.mean(permutedScores[threshold])
			
		ax.errorbar(threshold, thrMean, np.std(permutedScores[threshold]), marker='o', mfc='black', mec='black')
print "plotting data"
plotData(realScoreCountsCosmic, permutedScoreCountsCosmic, maxScore)
plt.savefig("cosmic.svg")

plotData(realScoreCountsSNVs, permutedScoreCountsSNVs, maxScore)
plt.savefig("snvs.svg")

plotData(realScoreCountsDEGS, permutedScoreCountsDEGS, maxScore)
plt.savefig("degs.svg")

plotData(realScoreCountsCosmicSNVs, permutedScoreCountsCosmicSNVs, maxScore)
plt.savefig("cosmicSNVs.svg")

plotData(realScoreCountsCosmicDEGs, permutedScoreCountsCosmicDEGs, maxScore)
plt.savefig("cosmicDEGs.svg")

plotData(realScoreCountsSNVDEGs, permutedScoreCountsSNVDEGs, maxScore)
plt.savefig("snvDEGs.svg")

plotData(realScoreCountsAll, permutedScoreCountsAll, maxScore)
plt.savefig("allCriteria.svg")


#Make the z-score plots.
#For every threshold, compute the z-score (value - mean / std) and plot.

def plotThresholdEnrichment(realScores, permutedScores, maxScore):
	plt.clf()
	ax = plt.subplot(1,1,1)
	for threshold in range(0, maxScore):
		zScore = (realScores[threshold] - np.mean(permutedScores[threshold])) / np.std(permutedScores[threshold])
		ax.plot(threshold, zScore, 'bo')

# plotThresholdEnrichment(realScoreCountsCosmic, permutedScoreCountsCosmic, 3)
# plt.savefig("cosmic_enrichtment.svg")
print "plot enrichment"
plotThresholdEnrichment(realScoreCountsCosmic, permutedScoreCountsCosmic, maxScore)
plt.savefig("cosmic_enrichment.svg")

plotThresholdEnrichment(realScoreCountsSNVs, permutedScoreCountsSNVs, maxScore)
plt.savefig("snvs_enrichment.svg")

plotThresholdEnrichment(realScoreCountsDEGS, permutedScoreCountsDEGS, maxScore)
plt.savefig("degs_enrichment.svg")

plotThresholdEnrichment(realScoreCountsCosmicSNVs, permutedScoreCountsCosmicSNVs, maxScore)
plt.savefig("cosmicSNVs_enrichment.svg")

plotThresholdEnrichment(realScoreCountsCosmicDEGs, permutedScoreCountsCosmicDEGs, maxScore)
plt.savefig("cosmicDEGs_enrichment.svg")

plotThresholdEnrichment(realScoreCountsSNVDEGs, permutedScoreCountsSNVDEGs, maxScore)
plt.savefig("snvDEGs_enrichment.svg")

plotThresholdEnrichment(realScoreCountsAll, permutedScoreCountsAll, maxScore)
plt.savefig("allCriteria_enrichment.svg")
