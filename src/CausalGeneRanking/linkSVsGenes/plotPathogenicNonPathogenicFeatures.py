"""
	Get the feature values of the coding & non-coding SVs, and see if these can be separated using PCA.
"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from six.moves import range
import re
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from scipy import interp
from scipy.stats import chi2_contingency
import random
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from statsmodels.sandbox.stats.multicomp import multipletests

#re-write for having all SVs in the same plot
### 4 separate plots, but gains and losses in one

def plotSignificances(significances, label, title, dataType):


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
				significanceMatrix[lossInd, svTypeInd] = -1
			else:
				significanceMatrix[lossInd, svTypeInd] = 0

			plotInd[sign] = lossInd+1
			#plotInd[sign] = lossInd+1
			lossInd += 2


			if gainSignificances[sign] < 0.05:
				significanceMatrix[gainInd, svTypeInd] = 1
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
					  cmap=ListedColormap(['#0055d4ff', '#b3b3b3ff', '#c83737ff']),
					  yticklabels=['Deletions', 'Duplications', 'Inversions', 'Translocations'])


	else:
		data = pd.DataFrame(significanceMatrix[:,0:3]) #exclude translocations, these are not there for germline.
		g=sns.heatmap(data.T,annot=False,square=True, linewidths=0.5,
					  cmap=ListedColormap(['#0055d4ff', '#b3b3b3ff', '#c83737ff']),
					  yticklabels=['Deletions', 'Duplications', 'Inversions'])

	g.set_yticklabels(g.get_yticklabels(), horizontalalignment='right',fontsize='small')
	plt.xticks(plotInd, ['eQTLs', 'Enhancers', 'Promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
													  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'RNA pol II', 'chromHMM CTCF', 'chromHMM CTCF+enhancer',
													  'chromHMM CTCF+promoter', 'chromHMM enhancer', 'chromHMM heterochromatin', 'chromHMM poised promoter',
													  'chromHMM promoter', 'chromHMM repeat', 'chromHMM repressed', 'chromHMM transcribed', 'Super enhancer', 'CTCF'], rotation=45, horizontalalignment='right')
	#plt.xticks([0.5,1.5,2.5,3.5], ['Deletions', 'Duplications', 'Inversions', 'Translocations'], rotation=45)

	plt.tight_layout()

	plt.savefig(title + '.svg')
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

#The truth data, DEG pairs
degData = np.loadtxt(sys.argv[1], dtype='object')

#pathwayAnnotation = np.loadtxt(sys.argv[1] + '_pathwayAnnotation.txt', dtype='object')

#Other comparisons: non-DEG pairs, coding SV-gene pairs, germline pairs, shuffled pairs.
nonDegData = np.loadtxt(sys.argv[2], dtype='object')
#codingData = np.loadtxt(sys.argv[3], dtype='object')
germlineData = np.loadtxt(sys.argv[3], dtype='object')
shuffledData = np.loadtxt(sys.argv[4], dtype='object')

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

	print(degs)
	print(nonDEGs)
	print(germlinePairs)
	print(shuffledPairs)

	print(degs.shape)
	print(nonDEGs.shape)
	print(germlinePairs.shape)
	print(shuffledPairs.shape)

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

#plot for each SV type
plotSignificances(perTypeResults, 'Positive pairs vs. negative pairs', 'gains_losses_pos_neg', 'somatic')
plotSignificances(perTypeResultsGL, 'Positive pairs vs. negative pairs', 'gains_losses_pos_gl', 'germline')
plotSignificances(perTypeResultsS, 'Positive pairs vs. negative pairs', 'gains_losses_pos_s', 'random')
exit()
#plotGainsLossesSamePlot(unAnnotatedLossesNormND, unAnnotatedGainsNormND, lossSignificances, gainSignificances, 'Positive pairs vs. negative pairs', typeLabel, 'log(% of positive pairs / % of negative pairs)', 1)
#plotGainsLossesSamePlot(unAnnotatedLossesNormGL, unAnnotatedGainsNormGL, lossSignificancesGL, gainSignificancesGL, 'Positive pairs vs. germline pairs', typeLabel, 'log(% of positive pairs / % of germline pairs)', 2)
#plotGainsLossesSamePlot(unAnnotatedLossesNormC, unAnnotatedGainsNormC, 'DEG pairs vs. coding pairs', typeLabel, 'log(% of DEG pairs / % of coding pairs)', 3)
#plotGainsLossesSamePlot(unAnnotatedLossesNormS, unAnnotatedGainsNormS, lossSignificancesS, gainSignificancesS, 'Positive pairs vs. shuffled pairs', typeLabel, 'log(% of positive pairs / % of shuffled pairs)', 3)
#plt.tight_layout()
#plt.savefig('gains_losses_' + svType + '.svg')
#plt.show()

np.random.seed(0)
positive = degs[:,1:52]

indices = np.arange(nonDEGs.shape[0])
rnd_indices = np.random.choice(indices, size=positive.shape[0])

negative = nonDEGs[rnd_indices][:,1:52]

allFeatures = np.concatenate((positive, negative), axis=0).astype(float)

#normalize
normalized = np.zeros(allFeatures.shape)
for col in range(0, allFeatures.shape[1]):

	if np.min(allFeatures[:,col]) == np.max(allFeatures[:,col]):
		continue

	print(col)
	print(np.min(allFeatures[:,col]))
	print(np.max(allFeatures[:,col]))

	normalized[:,col] = (allFeatures[:,col] - np.min(allFeatures[:,col])) / (np.max(allFeatures[:,col] - np.min(allFeatures[:,col])))

#normalized = (allFeatures-np.min(allFeatures))/(np.max(allFeatures)-np.min(allFeatures))
allFeatures = normalized
print(allFeatures)



#remove some features
# allFeaturesFiltered = []
# for row in range(0, allFeatures.shape[0]):
# 	rowFeatures = []
# 	for featureInd in range(0, allFeatures.shape[1]):
#
# 		if featureInd > 70 and featureInd < 75:
# 			continue
# 		else:
# 			rowFeatures.append(allFeatures[row][featureInd])
#
# 	allFeaturesFiltered.append(rowFeatures)
#
# allFeatures = np.array(allFeaturesFiltered)
#Make a PCA plot for the left/right set and see if these are really different

from sklearn.decomposition import PCA

#Get subset of PCA

pca = PCA(n_components=2)

projected = pca.fit_transform(allFeatures)
# projectedWithOffset = projected
#
# jitter = [0.01, -0.01]
# for row in range(0, projected.shape[0]):
# 	for col in range(0, projected.shape[1]):
# 		projectedWithOffset[row][col] += np.random.normal(-1, 1) * 0.1
#
# projected = projectedWithOffset

colorLabels = ['r']*positive.shape[0] + ['b']*negative.shape[0]
labels = [1]*positive.shape[0] + [0]*negative.shape[0]

fig,ax=plt.subplots(figsize=(7,5))
plt.scatter(projected[:, 0], projected[:, 1], edgecolors=colorLabels, facecolors='none')
plt.show()

#rasterize the PCA plot and make a density heatmap
# import math
#
# #
# colorLabels = np.array(colorLabels)
#
# #Get the minimum and maximum to determine the bounds of the plot.
# xmin = np.min(projected[:,0])
# xmax = np.max(projected[:,0])
# ymin = np.min(projected[:,1])
# ymax = np.max(projected[:,1])
#
# #Define the box size and how many boxes we should make
# print(xmin, xmax, ymin, ymax)
#
# #round the values to get covering boxes
# xmin = round(xmin)
# xmax = round(xmax)
# ymin = round(ymin)
# ymax = round(ymax)
#
# boxWidth = 0.2
# #Take the ceil to get the maximum possible without leaving out points
# xBoxNum = int(math.ceil((xmax - xmin) / boxWidth))
# yBoxNum = int(math.ceil((ymax - ymin) / boxWidth))
#
# #Placeholder for smoothed data
# plotGrid = np.zeros([xBoxNum, yBoxNum])
#
# #Loop through the data and show the data in the boxes
# yBoxStart = ymin
# yBoxEnd = ymin + boxWidth
# xBoxStart = xmin
# xBoxEnd = xmin + boxWidth
# for yInd in range(0, yBoxNum):
# 	for xInd in range(0, xBoxNum):
#
# 		#Find all data points that are within the current box
# 		xStartMatches = projected[:,0] >= xBoxStart
# 		xEndMatches = projected[:,0] <= xBoxEnd
#
# 		xMatches = xStartMatches * xEndMatches
#
# 		yStartMatches = projected[:,1] >= yBoxStart
# 		yEndMatches = projected[:,1] <= yBoxEnd
#
# 		yMatches = yStartMatches * yEndMatches
#
# 		dataInBox = projected[xMatches * yMatches]
# 		boxLabels = colorLabels[xMatches * yMatches]
#
# 		if len(dataInBox) > 0:
# 			#print dataInBox
#
# 			posCount = len(np.where(boxLabels == 'r')[0]) + 0.01
# 			negCount = len(np.where(boxLabels == 'b')[0]) + 0.01
#
# 			#Normalize for the total count of that label
# 			posCount = posCount / len(np.where(colorLabels == 'r')[0])
# 			negCount = negCount / len(np.where(colorLabels == 'b')[0])
#
# 			if negCount > 0:
# 				plotGrid[xInd,yInd] = np.log(posCount / float(negCount))
#
#
# 		#Move the box along x
# 		xBoxStart += boxWidth
# 		xBoxEnd += boxWidth
#
# 	yBoxStart += boxWidth
# 	yBoxEnd += boxWidth
# 	#Reset the box on x
# 	xBoxStart = xmin
# 	xBoxEnd = xmin + boxWidth
#
# plotGrid = np.ma.masked_where(plotGrid == 0, plotGrid)
# cmap = plt.cm.seismic
# cmap.set_bad(color='white')
# print(plotGrid)
# plt.imshow(plotGrid, cmap=cmap, interpolation='nearest')
# plt.show()

#very simple ml test
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import auc, precision_recall_curve
from sklearn import model_selection
from sklearn.metrics import average_precision_score, plot_roc_curve
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

from random import shuffle

def cvClassification(similarityMatrix, bagLabels, clf):


	scoring = {'accuracy' : make_scorer(accuracy_score),
			   'precision' : make_scorer(precision_score),
			   'recall' : make_scorer(recall_score),
			   'f1_score' : make_scorer(f1_score),
			   'average_precision' : make_scorer(average_precision_score),
			   'auc' : make_scorer(roc_auc_score)}

	kfold = model_selection.StratifiedKFold(n_splits=10, random_state=10, shuffle=True)

	results = model_selection.cross_validate(estimator=clf,
											  X=similarityMatrix,
											  y=bagLabels,
											  cv=kfold,
											  scoring=scoring)

	#print('accuracy: ', np.mean(results['test_accuracy']), np.std(results['test_accuracy']))
	#print('precision: ', np.mean(results['test_precision']), np.std(results['test_precision']))
	#print('recall: ', np.mean(results['test_recall']), np.std(results['test_recall']))
	print('F1 score: ', np.mean(results['test_f1_score']), np.std(results['test_f1_score']))
	print('AP: ', np.mean(results['test_average_precision']), np.std(results['test_average_precision']))
	print('AUC: ', np.mean(results['test_auc']), np.std(results['test_auc']))

	bagLabels = np.array(bagLabels)
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)

	fig, ax = plt.subplots()
	for i, (train, test) in enumerate(kfold.split(similarityMatrix, bagLabels)):
		clf.fit(similarityMatrix[train], bagLabels[train])
		viz = plot_roc_curve(clf, similarityMatrix[test], bagLabels[test],
							 name='ROC fold {}'.format(i),
							 alpha=0.3, lw=1, ax=ax)
		interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(viz.roc_auc)
	print(aucs)
	print(np.mean(aucs))
	print(np.std(aucs))

	ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
			label='Chance', alpha=.8)

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	# ax.plot(mean_fpr, mean_tpr, color='b',
	# 		label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
	# 		lw=2, alpha=.8)
	ax.plot(mean_fpr, mean_tpr, color='b',
			label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (np.mean(aucs), np.std(aucs)),
			lw=2, alpha=.8)

	std_tpr = np.std(tprs, axis=0)
	tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
	ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
					label=r'$\pm$ 1 std. dev.')

	ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
		   title="Receiver operating characteristic: deletions")
	ax.legend(loc="lower right")
	plt.show()
	exit()

	return np.mean(results['test_f1_score']), np.mean(results['test_average_precision'])

from sklearn.ensemble import RandomForestClassifier

names = ["Random Forest",]

classifiers = [
    #KNeighborsClassifier(3),
    #SVC(kernel="linear", C=0.025),
    #SVC(gamma=2, C=1),
    #GaussianProcessClassifier(1.0 * RBF(1.0)),
    #DecisionTreeClassifier(max_depth=100),
    RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)]
    #MLPClassifier(alpha=1, max_iter=1000),
    #AdaBoostClassifier(),
    #GaussianNB(),
    #QuadraticDiscriminantAnalysis()]

for classifierInd in range(0, len(classifiers)):
	print(names[classifierInd])
	rfClassifier = classifiers[classifierInd]

	featureCount = allFeatures.shape[1]-1
	#featureCount = 1
	f1s = []
	aps = []
	for featureInd in range(featureCount, allFeatures.shape[1]):

		f1, ap = cvClassification(allFeatures[:,0:featureInd], labels, rfClassifier)
		f1s.append(f1)
		aps.append(ap)

print('shuffled: ')

for classifierInd in range(0, len(classifiers)):
	print(names[classifierInd])
	rfClassifier = classifiers[classifierInd]
	f1Shuffled = []
	apShuffled = []
	for featureInd in range(featureCount, allFeatures.shape[1]):


		shuffle(labels)
		f1, ap = cvClassification(allFeatures[:,0:featureInd], labels, rfClassifier)
		f1Shuffled.append(f1)
		apShuffled.append(ap)