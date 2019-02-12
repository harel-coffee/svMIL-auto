	
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
		ax.errorbar(threshold, np.mean(permutedScores[threshold]), np.std(permutedScores[threshold]), marker='o', mfc='black', mec='black')

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
		ax.plot(threshold. zScore, 'bo')

# plotThresholdEnrichment(realScoreCountsCosmic, permutedScoreCountsCosmic, 3)
# plt.savefig("cosmic_enrichtment.svg")

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
