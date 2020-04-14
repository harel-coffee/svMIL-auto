"""
	Get the feature values of the coding & non-coding SVs, and see if these can be separated using PCA.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy import interp
from scipy.stats import chi2_contingency
import random
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from statsmodels.sandbox.stats.multicomp import multipletests

#re-write for having all SVs in the same plot
### 4 separate plots, but gains and losses in one

def plotSignificances(significances, label, title, dataType, minScale, maxScale, finalOutDir):


	#plot heatmap instead.
	#make a matrix of features by SV types. First test with gains only
	significanceMatrix = np.zeros([len(significances['DEL'][1])*2, len(significances)])
	plotInd = [0]*len(significances['DEL'][1])
	for svTypeInd in range(0, len(significances)):
		svType = list(significances.keys())[svTypeInd]

		lossSignificances = significances[svType][0]
		gainSignificances = significances[svType][1]

		lossInd = 0
		gainInd = 1
		for sign in range(0, len(lossSignificances)):

			if lossSignificances[sign] < 0.05:
				#significanceMatrix[lossInd, svTypeInd] = -1
				significanceMatrix[lossInd, svTypeInd] = (np.log(lossSignificances[sign]) - minScale) / (maxScale - minScale)

				#significanceMatrix[lossInd, svTypeInd] = np.log(lossSignificances[sign]) #keep losses negative
			else:
				significanceMatrix[lossInd, svTypeInd] = 0

			plotInd[sign] = lossInd+1
			#plotInd[sign] = lossInd+1
			lossInd += 2


			if gainSignificances[sign] < 0.05:
				#significanceMatrix[gainInd, svTypeInd] = 1
				significanceMatrix[gainInd, svTypeInd] = (-1) * (np.log(gainSignificances[sign]) - minScale) / (maxScale - minScale)
				# print((np.log(gainSignificances[sign]) - minScale) / (maxScale - minScale))
				# print(maxScale - minScale)
				# print(gainSignificances[sign])
				# print(np.log(gainSignificances[sign]))
				# print(minScale)
				# print(maxScale)
				# exit()

				#significanceMatrix[gainInd, svTypeInd] = (-1) * np.log(gainSignificances[sign]) #gains positive
			else:
				significanceMatrix[gainInd, svTypeInd] = 0

			gainInd += 2

	print(significanceMatrix)

	print(plotInd)

	fig =plt.figure(figsize=(15,10))
	if dataType != 'germline':
		data = pd.DataFrame(significanceMatrix)
		#plot heat map
		g=sns.heatmap(data.T,annot=False,square=True, linewidths=0.5,
					  cmap='coolwarm', vmin=-1, vmax=1,
					  #cmap=ListedColormap(['#0055d4ff', '#b3b3b3ff', '#c83737ff']),
					  yticklabels=['Deletions', 'Duplications', 'Inversions', 'Translocations'])


	else:
		data = pd.DataFrame(significanceMatrix[:,0:3]) #exclude translocations, these are not there for germline.
		g=sns.heatmap(data.T,annot=False,square=True, linewidths=0.5,
					  #cmap=ListedColormap(['#0055d4ff', '#b3b3b3ff', '#c83737ff']),
					  cmap='coolwarm', vmin=-1, vmax=1,
					  yticklabels=['Deletions', 'Duplications', 'Inversions'])

	g.set_yticklabels(g.get_yticklabels(), horizontalalignment='right',fontsize='small')
	plt.xticks(plotInd, ['eQTLs', 'Enhancers', 'Promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
													  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'RNA pol II', 'chromHMM CTCF', 'chromHMM CTCF+enhancer',
													  'chromHMM CTCF+promoter', 'chromHMM enhancer', 'chromHMM heterochromatin', 'chromHMM poised promoter',
													  'chromHMM promoter', 'chromHMM repeat', 'chromHMM repressed', 'chromHMM transcribed', 'Super enhancer', 'CTCF'], rotation=45, horizontalalignment='right')
	#plt.xticks([0.5,1.5,2.5,3.5], ['Deletions', 'Duplications', 'Inversions', 'Translocations'], rotation=45)

	plt.tight_layout()

	plt.savefig(finalOutDir + '/' + title + '.svg')
	plt.show()


def getSignificance(features, shuffledFeatures, rangeA, rangeB):

	#for each loss feature, check if more/less than by random chance.
	lossSignificances = []
	for i in range(rangeA, rangeB):

		degLosses = len(np.where(features[:,i] == '1.0')[0])
		degNoLosses = len(np.where(features[:,i] == '0.0')[0])

		shuffledLosses = len(np.where(shuffledFeatures[:,i] == '1.0')[0])
		shuffledNoLosses = len(np.where(shuffledFeatures[:,i] == '0.0')[0])


		if degLosses == 0 or degNoLosses == 0 or shuffledLosses == 0 or shuffledNoLosses == 0:
			lossSignificances.append(1)
			continue

		obs = np.array([[degLosses, degNoLosses], [shuffledLosses, shuffledNoLosses]])

		result = chi2_contingency(obs)
		p = result[1]

		lossSignificances.append(p)

	#do MTC
	reject, pAdjusted, _, _ = multipletests(lossSignificances, method='bonferroni')


	return pAdjusted

outDir = sys.argv[1]
finalOutDir = outDir + '/featureComparisonHeatmap/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

#The truth data, DEG pairs
degData = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_pathogenicPairsFeatures.txt', dtype='object')

nonDegData = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_nonPathogenicPairsFeatures.txt', dtype='object')
germlineData = np.loadtxt(outDir + '/linkedSVGenePairs/germline/nonCoding_geneSVPairs.txt_', dtype='object')
shuffledData = np.loadtxt(outDir + '/linkedSVGenePairs/random/nonCoding_geneSVPairs.txt_0', dtype='object')

svTypes = ['DEL', 'DUP', 'INV', 'ITX']
perTypeResults = dict()
perTypeResultsGL = dict()
perTypeResultsS = dict()
for svType in svTypes:

	nonDEGs = []

	for degPair in nonDegData:

		sv = degPair[0].split("_")
		if svType != '':
			typeMatch = re.search(svType, sv[12], re.IGNORECASE)
			if typeMatch is None:
				continue
		nonDEGs.append(degPair)

	germlinePairs = []
	for pair in germlineData:

		sv = pair[0].split("_")

		if svType != '':
			typeMatch = re.match(svType, sv[12], re.IGNORECASE)
			if typeMatch is None:
				continue
		germlinePairs.append(pair)

	shuffledPairs = []
	for pair in shuffledData:

		sv = pair[0].split("_")
		if svType != '':
			typeMatch = re.match(svType, sv[12], re.IGNORECASE)
			if typeMatch is None:
				continue
		shuffledPairs.append(pair)

	degs = []
	for degPair in degData:

		sv = degPair[0].split("_")

		if svType != '':
			typeMatch = re.match(svType, sv[12], re.IGNORECASE)
			if typeMatch is None:
				continue



		degs.append(degPair)

	degs = np.array(degs)
	nonDEGs = np.array(nonDEGs)
	germlinePairs = np.array(germlinePairs)
	shuffledPairs = np.array(shuffledPairs)

	unAnnotatedDEGs = []
	for degPair in degs:

		features = degPair[1:]

		unAnnotatedDEGs.append(features)

	unAnnotatedDEGs = np.array(unAnnotatedDEGs)

	lossSignificancesTmp = [2]*26

	if svType == 'INV':
		lossSignificances = getSignificance(unAnnotatedDEGs, nonDEGs[:,1:], 0, 26)
		lossSignificancesGL = getSignificance(unAnnotatedDEGs, germlinePairs[:,1:], 0, 26)
		lossSignificancesS = getSignificance(unAnnotatedDEGs, shuffledPairs[:,1:], 0, 26)
	elif svType == 'ITX':
		lossSignificances = getSignificance(unAnnotatedDEGs, nonDEGs[:,1:], 0, 26)
		lossSignificancesGL = lossSignificancesTmp
		lossSignificancesS = getSignificance(unAnnotatedDEGs, shuffledPairs[:,1:], 0, 26)
	else:
		lossSignificances = lossSignificancesTmp
		lossSignificancesGL = lossSignificancesTmp
		lossSignificancesS = lossSignificancesTmp


	gainSignificancesTmp = [2]*26

	if svType == 'ITX':
		gainSignificances = getSignificance(unAnnotatedDEGs, nonDEGs[:,1:], 26, 52)
		gainSignificancesGL = gainSignificancesTmp
		gainSignificancesS = getSignificance(unAnnotatedDEGs, shuffledPairs[:,1:], 26, 52)
	elif svType == 'DUP':
		gainSignificances = getSignificance(unAnnotatedDEGs, nonDEGs[:,1:], 26, 52)
		gainSignificancesGL = gainSignificancesTmp
		gainSignificancesS = getSignificance(unAnnotatedDEGs, shuffledPairs[:,1:], 26, 52)
	else:

		gainSignificances = getSignificance(unAnnotatedDEGs, nonDEGs[:,1:], 26, 52)
		gainSignificancesGL = getSignificance(unAnnotatedDEGs, germlinePairs[:,1:], 26, 52)
		gainSignificancesS = getSignificance(unAnnotatedDEGs, shuffledPairs[:,1:], 26, 52)

	#MTC on results


	#collect results per SV type
	perTypeResults[svType] = [lossSignificances, gainSignificances]
	perTypeResultsGL[svType] = [lossSignificancesGL, gainSignificancesGL]
	perTypeResultsS[svType] = [lossSignificancesS, gainSignificancesS]

#Normalize the p-values so that the scale is the same between plots
def getMinMaxScale(results, svTypes):
	minScale = float('inf') #the lowest log(p-value) for the losses
	maxScale = float('inf') #the lowest log(p-value) for the gains

	for svType in svTypes:

		losses, gains = results[svType]

		for loss in losses:
			if loss < 0.05:
				if np.log(loss) < minScale:
					minScale = np.log(loss)

		for gain in gains:
			if gain < 0.05:
				if np.log(gain) < maxScale:
					maxScale = np.log(gain)

	return minScale, maxScale

pnMin, pnMax = getMinMaxScale(perTypeResults, svTypes)
pglMin, pglMax = getMinMaxScale(perTypeResultsGL, svTypes)
psMin, psMax = getMinMaxScale(perTypeResultsS, svTypes)

minScale = min(pnMin, pglMin, psMin)
maxScale = min(pnMax, pglMax, psMax)

#anchor around 0 so that non-significant is shown as grey.
# if maxScale < minScale:
# 	minScale = maxScale
# else:
# 	maxScale = minScale

#plot for each SV type
plotSignificances(perTypeResults, 'Positive pairs vs. negative pairs', 'gains_losses_pos_neg', 'somatic', minScale, maxScale, finalOutDir)
plotSignificances(perTypeResultsGL, 'Positive pairs vs. negative pairs', 'gains_losses_pos_gl', 'germline', minScale, maxScale, finalOutDir)
plotSignificances(perTypeResultsS, 'Positive pairs vs. negative pairs', 'gains_losses_pos_s', 'random', minScale, maxScale, finalOutDir)
