"""
	Make figure 2e

	Gets all the gains and losses for the pathogenic SVs, non-pathogenic SVs, germline SVs, and random SVs, and plots which are significant in the
	pathogenic group compared to the others.
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

import matplotlib
matplotlib.use('Agg')

def plotSignificances(significances, label, title, dataType, finalOutDir):
	"""
		Plot the heatmaps of fig 2e.

		significances (dictionary): keys are the SV types, contains a list where ind 0 are the losses, 1 the gains.
		label (str): unused
		title (str): title that the plot will get in the output file.
		dataType (str): if germline, we will not plot anything for translocations, since these are not present in germline.
		finalOutDir (str): the final output directory that the plot will be written to.

	"""

	#plot heatmap
	#make a matrix of features by SV types.
	significanceMatrix = np.zeros([len(significances['DEL'][1])*2, len(significances)])
	plotInd = [0]*len(significances['DEL'][1])

	#below this we call it 'very' significant.
	signCutoff = 1e-5

	for svTypeInd in range(0, len(significances)):
		svType = list(significances.keys())[svTypeInd]

		lossSignificances = significances[svType][0]
		gainSignificances = significances[svType][1]

		lossInd = 0
		gainInd = 1
		for sign in range(0, len(lossSignificances)):

			if lossSignificances[sign] < 0.05:

				if lossSignificances[sign] < signCutoff:
					significanceMatrix[lossInd, svTypeInd] = -2
				else:
					significanceMatrix[lossInd, svTypeInd] = -1
				
			else:
				significanceMatrix[lossInd, svTypeInd] = 0

			plotInd[sign] = lossInd+1
			lossInd += 2


			if gainSignificances[sign] < 0.05:

				if gainSignificances[sign] < signCutoff:
					significanceMatrix[gainInd, svTypeInd] = 2
				else:
					significanceMatrix[gainInd, svTypeInd] = 1

				
			else:
				significanceMatrix[gainInd, svTypeInd] = 0

			gainInd += 2
	

	fig =plt.figure(figsize=(15,10))
	if dataType != 'germline':
		data = pd.DataFrame(significanceMatrix)
		#plot heat map
		g=sns.heatmap(data.T,annot=False,square=True, linewidths=0.5,
					  cmap=ListedColormap(['#0055d4ff', '#0055d47d', '#f7f6f6ff', '#c8373780', '#c83737ff']),
					  yticklabels=['Deletions', 'Duplications', 'Inversions', 'Translocations'])


	else:
		data = pd.DataFrame(significanceMatrix[:,0:3]) #exclude translocations, these are not there for germline.
		g=sns.heatmap(data.T,annot=False,square=True, linewidths=0.5,
					  cmap=ListedColormap(['#0055d4ff', '#0055d47d', '#f7f6f6ff', '#c8373780', '#c83737ff']),
					  yticklabels=['Deletions', 'Duplications', 'Inversions'])

	g.set_yticklabels(g.get_yticklabels(), horizontalalignment='right',fontsize='small')
	plt.xticks(plotInd, ['eQTLs', 'Enhancers', 'Promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
													  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'RNA pol II', 'chromHMM CTCF', 'chromHMM CTCF+enhancer',
													  'chromHMM CTCF+promoter', 'chromHMM enhancer', 'chromHMM heterochromatin', 'chromHMM poised promoter',
													  'chromHMM promoter', 'chromHMM repeat', 'chromHMM repressed', 'chromHMM transcribed', 'Super enhancer', 'CTCF'], rotation=45, horizontalalignment='right')

	plt.tight_layout()

	plt.savefig(finalOutDir + '/' + title + '.svg')
	


def getSignificance(features, negativeFeatures, rangeA, rangeB):
	"""
		For each feature (gain/loss), get the significance compared to the negative group.

		features (numpy array): all the features for the positive group
		negativeFeatures (numpy array): all the features for the negative group
		rangeA (int): start from this feature
		rangeB (int): until this feature

		return:
		pAdjusted (numpy array): bonferroni adjusted p-values.

	"""

	#for each loss feature, check if more/less than by random chance.

	significances = []
	for i in range(rangeA, rangeB):

		pathogenic = len(np.where(features[:,i] == '1.0')[0])
		pathogenicNo = len(np.where(features[:,i] == '0.0')[0])

		negative = len(np.where(negativeFeatures[:,i] == '1.0')[0])
		negativeNo = len(np.where(negativeFeatures[:,i] == '0.0')[0])

		if pathogenic == 0 or pathogenicNo == 0 or negative == 0 or negativeNo == 0:
			significances.append(1)
			continue

		obs = np.array([[pathogenic, pathogenicNo], [negative, negativeNo]])

		result = chi2_contingency(obs)
		p = result[1]

		significances.append(p)

	#do MTC
	reject, pAdjusted, _, _ = multipletests(significances, method='bonferroni')


	return pAdjusted

outDir = sys.argv[1]
finalOutDir = outDir + '/figure2/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

#get the gains and losses for each group
pathogenicData = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_pathogenicPairsFeatures.txt', dtype='object')
nonPathogenicData = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_nonPathogenicPairsFeatures.txt', dtype='object')
germlineData = np.loadtxt(outDir + '/linkedSVGenePairs/germline/nonCoding_geneSVPairs.txt_', dtype='object')
shuffledData = np.loadtxt(outDir + '/linkedSVGenePairs/random/nonCoding_geneSVPairs.txt_0', dtype='object')

svTypes = ['DEL', 'DUP', 'INV', 'ITX']
perTypeResults = dict() #pathogenic to non-pathogenic
perTypeResultsGL = dict() #pathogenic to germline
perTypeResultsS = dict() #pathogenic to shuffled
for svType in svTypes:

	nonPathogenics = []

	for pair in nonPathogenicData:

		sv = pair[0].split("_")
		if svType != '':
			typeMatch = re.search(svType, sv[12], re.IGNORECASE)
			if typeMatch is None:
				continue
		nonPathogenics.append(pair)

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

	pathogenics = []
	for pair in pathogenicData:

		sv = pair[0].split("_")

		if svType != '':
			typeMatch = re.match(svType, sv[12], re.IGNORECASE)
			if typeMatch is None:
				continue



		pathogenics.append(pair)

	pathogenics = np.array(pathogenics)
	nonPathogenics = np.array(nonPathogenics)
	germlinePairs = np.array(germlinePairs)
	shuffledPairs = np.array(shuffledPairs)

	unAnnotatedPathogenics = []
	for pair in pathogenics:

		features = pair[1:]

		unAnnotatedPathogenics.append(features)

	unAnnotatedPathogenics = np.array(unAnnotatedPathogenics)

	#Get all the significances.
	#for losses, these are between 0 and 26, while gains are 26 to 52.
	#set values to 2 by default, so that when there are no losses (e.g. duplications and deletions),
	#we can easily distinguish between not-significant and simply not checked.
	lossSignificancesTmp = [2]*26

	if svType == 'INV':
		lossSignificances = getSignificance(unAnnotatedPathogenics, nonPathogenics[:,1:], 0, 26)
		lossSignificancesGL = getSignificance(unAnnotatedPathogenics, germlinePairs[:,1:], 0, 26)
		lossSignificancesS = getSignificance(unAnnotatedPathogenics, shuffledPairs[:,1:], 0, 26)
	elif svType == 'ITX':
		lossSignificances = getSignificance(unAnnotatedPathogenics, nonPathogenics[:,1:], 0, 26)
		lossSignificancesGL = lossSignificancesTmp
		lossSignificancesS = getSignificance(unAnnotatedPathogenics, shuffledPairs[:,1:], 0, 26)
	else:
		lossSignificances = lossSignificancesTmp
		lossSignificancesGL = lossSignificancesTmp
		lossSignificancesS = lossSignificancesTmp


	gainSignificancesTmp = [2]*26

	if svType == 'ITX':
		gainSignificances = getSignificance(unAnnotatedPathogenics, nonPathogenics[:,1:], 26, 52)
		gainSignificancesGL = gainSignificancesTmp
		gainSignificancesS = getSignificance(unAnnotatedPathogenics, shuffledPairs[:,1:], 26, 52)
	elif svType == 'DUP':
		gainSignificances = getSignificance(unAnnotatedPathogenics, nonPathogenics[:,1:], 26, 52)
		gainSignificancesGL = gainSignificancesTmp
		gainSignificancesS = getSignificance(unAnnotatedPathogenics, shuffledPairs[:,1:], 26, 52)
	else:

		gainSignificances = getSignificance(unAnnotatedPathogenics, nonPathogenics[:,1:], 26, 52)
		gainSignificancesGL = getSignificance(unAnnotatedPathogenics, germlinePairs[:,1:], 26, 52)
		gainSignificancesS = getSignificance(unAnnotatedPathogenics, shuffledPairs[:,1:], 26, 52)

	#collect results per SV type
	perTypeResults[svType] = [lossSignificances, gainSignificances]
	perTypeResultsGL[svType] = [lossSignificancesGL, gainSignificancesGL]
	perTypeResultsS[svType] = [lossSignificancesS, gainSignificancesS]


#plot for each SV type
plotSignificances(perTypeResults, 'Positive pairs vs. negative pairs', 'fig2e_gains_losses_pos_neg', 'somatic',finalOutDir)
plotSignificances(perTypeResultsGL, 'Positive pairs vs. negative pairs', 'fig2e_gains_losses_pos_gl', 'germline', finalOutDir)
plotSignificances(perTypeResultsS, 'Positive pairs vs. negative pairs', 'fig2e_gains_losses_pos_s', 'random', finalOutDir)
