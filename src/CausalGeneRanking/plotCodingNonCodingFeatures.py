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

degData = np.loadtxt(sys.argv[1], dtype='object')
nonDegData = np.loadtxt(sys.argv[2], dtype='object')

#Split into features for deg pairs and non-deg pairs.
#allow option to split this per SV type

svTypeDeg = "transl_inter"
svTypeNonDeg = 'inversion'

if svTypeDeg != "":
	
	degFeatures = []
	for sv in degData:
		splitSV = sv[0].split("_")
		if len(splitSV) == 9: #for del, dup, inv
			if re.search(splitSV[8], svTypeDeg, re.IGNORECASE):
				degFeatures.append(sv[1:])
		else:
			if re.search("_".join([splitSV[8], splitSV[9]]), svTypeDeg, re.IGNORECASE):
				degFeatures.append(sv[1:])
			
	nonDegFeatures = []
	for sv in nonDegData:
		
		splitSV = sv[0].split("_")
		if len(splitSV) == 9: #for del, dup, inv
			if re.search(splitSV[8], svTypeNonDeg, re.IGNORECASE):
				nonDegFeatures.append(sv[1:])
		else:
			if re.search("_".join([splitSV[8], splitSV[9]]), svTypeNonDeg, re.IGNORECASE):
				nonDegFeatures.append(sv[1:])
	
	degFeatures = np.array(degFeatures)
	nonDegFeatures = np.array(nonDegFeatures)
else:
	
	degFeatures = degData[:,1:]
	nonDegFeatures = nonDegData[:,1:]

print(degFeatures)
print(nonDegFeatures)

#Make the bar plots
leftFeatures = nonDegFeatures
rightFeatures = degFeatures

eQTLLosses = leftFeatures[:,0].astype(float)
enhancerLosses = leftFeatures[:,1].astype(float)
promoterLosses = leftFeatures[:,2].astype(float)
cpgLosses = leftFeatures[:,3].astype(float)
tfLosses = leftFeatures[:,4].astype(float)
hicLosses = leftFeatures[:,5].astype(float)
h3k9me3Losses = leftFeatures[:,6].astype(float)
h3k4me3Losses = leftFeatures[:,7].astype(float)
h3k27acLosses = leftFeatures[:,8].astype(float)
h3k27me3Losses = leftFeatures[:,9].astype(float)
h3k4me1Losses = leftFeatures[:,10].astype(float)
h3k36me3Losses = leftFeatures[:,11].astype(float)
dnaseLosses = leftFeatures[:,12].astype(float)

ctcfLosses = leftFeatures[:,13].astype(float)
ctcfEnhancerLosses = leftFeatures[:,14].astype(float)
ctcfPromoterLosses = leftFeatures[:,15].astype(float)
chromHmmEnhancerLosses = leftFeatures[:,16].astype(float)
heterochromatinLosses = leftFeatures[:,17].astype(float)
poisedPromoterLosses = leftFeatures[:,18].astype(float)
chromHmmPromoterLosses = leftFeatures[:,19].astype(float)
repeatLosses = leftFeatures[:,20].astype(float)
repressedLosses = leftFeatures[:,21].astype(float)
transcribedLosses = leftFeatures[:,22].astype(float)

deletionsLeft = leftFeatures[:,46].astype(float)
duplicationsLeft = leftFeatures[:,47].astype(float)
inversionsLeft = leftFeatures[:,48].astype(float)
translocationsLeft = leftFeatures[:,49].astype(float)
cosmicLeft = leftFeatures[:,50].astype(float)

lossData = [np.sum(eQTLLosses), np.sum(enhancerLosses), np.sum(promoterLosses), np.sum(cpgLosses),
			np.sum(tfLosses), np.sum(hicLosses), np.sum(h3k9me3Losses), np.sum(h3k4me3Losses), np.sum(h3k27acLosses),
			np.sum(h3k27me3Losses), np.sum(h3k4me1Losses), np.sum(h3k36me3Losses), np.sum(dnaseLosses),
			np.sum(ctcfLosses), np.sum(ctcfEnhancerLosses), np.sum(ctcfPromoterLosses), np.sum(chromHmmEnhancerLosses),
			np.sum(heterochromatinLosses), np.sum(poisedPromoterLosses), np.sum(chromHmmPromoterLosses),
			np.sum(repeatLosses), np.sum(repressedLosses), np.sum(transcribedLosses),
			np.sum(deletionsLeft),np.sum(duplicationsLeft), np.sum(inversionsLeft), np.sum(translocationsLeft),
			np.sum(cosmicLeft)]
lossData = np.array(lossData)

lossData = (lossData / float(leftFeatures.shape[0]))
print(lossData)
#lossData = np.log(lossData)


width = 0.35

plt.barh(np.arange(len(lossData)), lossData, width, label='Germline pairs', color='blue')
# plt.xticks(range(0, len(lossData)),
	   # ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'], rotation=90)
# plt.show()

eQTLLosses = rightFeatures[:,0].astype(float)
enhancerLosses = rightFeatures[:,1].astype(float)
promoterLosses = rightFeatures[:,2].astype(float)
cpgLosses = rightFeatures[:,3].astype(float)
tfLosses = rightFeatures[:,4].astype(float)
hicLosses = rightFeatures[:,5].astype(float)
h3k9me3Losses = rightFeatures[:,6].astype(float)
h3k4me3Losses = rightFeatures[:,7].astype(float)
h3k27acLosses = rightFeatures[:,8].astype(float)
h3k27me3Losses = rightFeatures[:,9].astype(float)
h3k4me1Losses = rightFeatures[:,10].astype(float)
h3k36me3Losses = rightFeatures[:,11].astype(float)
dnaseLosses = rightFeatures[:,12].astype(float)


ctcfLosses = rightFeatures[:,13].astype(float)
ctcfEnhancerLosses = rightFeatures[:,14].astype(float)
ctcfPromoterLosses = rightFeatures[:,15].astype(float)
chromHmmEnhancerLosses = rightFeatures[:,16].astype(float)
heterochromatinLosses = rightFeatures[:,17].astype(float)
poisedPromoterLosses = rightFeatures[:,18].astype(float)
chromHmmPromoterLosses = rightFeatures[:,19].astype(float)
repeatLosses = rightFeatures[:,20].astype(float)
repressedLosses = rightFeatures[:,21].astype(float)
transcribedLosses = rightFeatures[:,22].astype(float)

deletionsRight = rightFeatures[:,46].astype(float)
duplicationsRight = rightFeatures[:,47].astype(float)
inversionsRight = rightFeatures[:,48].astype(float)
translocationsRight = rightFeatures[:,49].astype(float)
cosmicRight = rightFeatures[:,50].astype(float)

lossData = [np.sum(eQTLLosses), np.sum(enhancerLosses), np.sum(promoterLosses), np.sum(cpgLosses),
			np.sum(tfLosses), np.sum(hicLosses), np.sum(h3k9me3Losses), np.sum(h3k4me3Losses), np.sum(h3k27acLosses),
			np.sum(h3k27me3Losses), np.sum(h3k4me1Losses), np.sum(h3k36me3Losses), np.sum(dnaseLosses),
			np.sum(ctcfLosses), np.sum(ctcfEnhancerLosses), np.sum(ctcfPromoterLosses), np.sum(chromHmmEnhancerLosses),
			np.sum(heterochromatinLosses), np.sum(poisedPromoterLosses), np.sum(chromHmmPromoterLosses),
			np.sum(repeatLosses), np.sum(repressedLosses), np.sum(transcribedLosses),
			np.sum(deletionsRight),np.sum(duplicationsRight), np.sum(inversionsRight), np.sum(translocationsRight),
			np.sum(cosmicRight)]
lossData = np.array(lossData)
lossData = (lossData / float(rightFeatures.shape[0]))
print(lossData)
#lossData = np.log(lossData)


plt.barh(np.arange(len(lossData)) + width, lossData, width, label='Somatic pairs', color='red')
plt.yticks(np.arange(len(lossData) + width / 2), ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
												  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'CTCF', 'CTCF+Enhancer',
												  'CTCF+Promoter', 'chromHMM Enhancer', 'heterochromatin', 'poised promoter',
												  'chromHMM Promoter', 'repeat', 'repressed', 'transcribed',
												  'Deletions', 'Duplications', 'Inversions', 'Translocations',
												  'COSMIC'])
plt.xlim([0,1])
plt.legend(loc='best')
plt.tight_layout()
plt.show()
#plt.savefig('Output/degNonDeg_losses.svg')
plt.clf()
#Gains

eQTLGains = leftFeatures[:,23].astype(float)
enhancerGains = leftFeatures[:,24].astype(float)
promoterGains = leftFeatures[:,25].astype(float)
cpgGains = leftFeatures[:,26].astype(float)
tfGains = leftFeatures[:,27].astype(float)
hicGains = leftFeatures[:,28].astype(float)
h3k9me3Gains = leftFeatures[:,29].astype(float)
h3k4me3Gains = leftFeatures[:,30].astype(float)
h3k27acGains = leftFeatures[:,31].astype(float)
h3k27me3Gains = leftFeatures[:,32].astype(float)
h3k4me1Gains = leftFeatures[:,33].astype(float)
h3k36me3Gains = leftFeatures[:,34].astype(float)
dnaseGains = leftFeatures[:,35].astype(float)

ctcfGains = leftFeatures[:,36].astype(float)
ctcfEnhancerGains = leftFeatures[:,37].astype(float)
ctcfPromoterGains = leftFeatures[:,38].astype(float)
chromHmmEnhancerGains = leftFeatures[:,39].astype(float)
heterochromatinGains = leftFeatures[:,40].astype(float)
poisedPromoterGains = leftFeatures[:,41].astype(float)
chromHmmPromoterGains = leftFeatures[:,42].astype(float)
repeatGains = leftFeatures[:,43].astype(float)
repressedGains = leftFeatures[:,44].astype(float)
transcribedGains = leftFeatures[:,45].astype(float)

gainData = [np.sum(eQTLGains), np.sum(enhancerGains), np.sum(promoterGains), np.sum(cpgGains),
			np.sum(tfGains), np.sum(hicGains), np.sum(h3k9me3Gains), np.sum(h3k4me3Gains), np.sum(h3k27acGains),
			np.sum(h3k27me3Gains), np.sum(h3k4me1Gains), np.sum(h3k36me3Gains), np.sum(dnaseGains),
			np.sum(ctcfGains), np.sum(ctcfEnhancerGains), np.sum(ctcfPromoterGains), np.sum(chromHmmEnhancerGains),
			np.sum(heterochromatinGains), np.sum(poisedPromoterGains), np.sum(chromHmmPromoterGains),
			np.sum(repeatGains), np.sum(repressedGains), np.sum(transcribedGains),
			np.sum(deletionsLeft),np.sum(duplicationsLeft), np.sum(inversionsLeft), np.sum(translocationsLeft),
			np.sum(cosmicLeft)]

gainData = np.array(gainData)
gainData = (gainData / float(leftFeatures.shape[0]))
print(gainData)
#gainData = np.log(gainData)


plt.barh(np.arange(len(gainData)), gainData, width, label='Germline pairs', color='blue')
# plt.xticks(np.arange(len(gainData)),
# 		   ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'], rotation=90)
# plt.show()
	
eQTLGains = rightFeatures[:,23].astype(float)
enhancerGains = rightFeatures[:,24].astype(float)
promoterGains = rightFeatures[:,25].astype(float)
cpgGains = rightFeatures[:,26].astype(float)
tfGains = rightFeatures[:,27].astype(float)
hicGains = rightFeatures[:,28].astype(float)
h3k9me3Gains = rightFeatures[:,29].astype(float)
h3k4me3Gains = rightFeatures[:,30].astype(float)
h3k27acGains = rightFeatures[:,31].astype(float)
h3k27me3Gains = rightFeatures[:,32].astype(float)
h3k4me1Gains = rightFeatures[:,33].astype(float)
h3k36me3Gains = rightFeatures[:,34].astype(float)
dnaseGains = rightFeatures[:,35].astype(float)

ctcfGains = rightFeatures[:,36].astype(float)
ctcfEnhancerGains = rightFeatures[:,37].astype(float)
ctcfPromoterGains = rightFeatures[:,38].astype(float)
chromHmmEnhancerGains = rightFeatures[:,39].astype(float)
heterochromatinGains = rightFeatures[:,40].astype(float)
poisedPromoterGains = rightFeatures[:,41].astype(float)
chromHmmPromoterGains = rightFeatures[:,42].astype(float)
repeatGains = rightFeatures[:,43].astype(float)
repressedGains = rightFeatures[:,44].astype(float)
transcribedGains = rightFeatures[:,45].astype(float)

gainData = [np.sum(eQTLGains), np.sum(enhancerGains), np.sum(promoterGains), np.sum(cpgGains),
			np.sum(tfGains), np.sum(hicGains), np.sum(h3k9me3Gains), np.sum(h3k4me3Gains), np.sum(h3k27acGains),
			np.sum(h3k27me3Gains), np.sum(h3k4me1Gains), np.sum(h3k36me3Gains), np.sum(dnaseGains),
			np.sum(ctcfGains), np.sum(ctcfEnhancerGains), np.sum(ctcfPromoterGains), np.sum(chromHmmEnhancerGains),
			np.sum(heterochromatinGains), np.sum(poisedPromoterGains), np.sum(chromHmmPromoterGains),
			np.sum(repeatGains), np.sum(repressedGains), np.sum(transcribedGains),
			np.sum(deletionsRight),np.sum(duplicationsRight), np.sum(inversionsRight), np.sum(translocationsRight),
			np.sum(cosmicRight)]

gainData = np.array(gainData)
gainData = (gainData / float(rightFeatures.shape[0]))
print(gainData)
#gainData = np.log(gainData)



plt.barh(np.arange(len(gainData)) + width, gainData, width, label='Somatic pairs',color='red')
plt.yticks(np.arange(len(gainData) + width / 2),['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
												 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI', 'CTCF', 'CTCF+Enhancer',
												 'CTCF+Promoter', 'chromHMM Enhancer', 'heterochromatin', 'poised promoter',
												 'chromHMM Promoter', 'repeat', 'repressed', 'transcribed',
												 'Deletions', 'Duplications', 'Inversions', 'Translocations',
												 'COSMIC'])

plt.xlim([0,1])
plt.legend(loc='best')
plt.tight_layout()
plt.show()
#plt.savefig('Output/degNonDeg_gains.svg')





allFeatures = np.concatenate((leftFeatures, rightFeatures), axis=0)

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

colorLabels = []
labels = []

for i in range(0, allFeatures.shape[0]):
	
	if i < leftFeatures.shape[0]:
		colorLabels.append('b')
		labels.append(0)
	elif i >= leftFeatures.shape[0] and i < (leftFeatures.shape[0] + rightFeatures.shape[0]):
		colorLabels.append('r')
		labels.append(1)

fig,ax=plt.subplots(figsize=(7,5))
plt.scatter(projected[:, 0], projected[:, 1], edgecolors=colorLabels, facecolors='none')
plt.show()

#rasterize the PCA plot and make a density heatmap
import math

#
colorLabels = np.array(colorLabels)

#Get the minimum and maximum to determine the bounds of the plot.
xmin = np.min(projected[:,0])
xmax = np.max(projected[:,0])
ymin = np.min(projected[:,1])
ymax = np.max(projected[:,1])

#Define the box size and how many boxes we should make
print(xmin, xmax, ymin, ymax)

#round the values to get covering boxes
xmin = round(xmin)
xmax = round(xmax)
ymin = round(ymin)
ymax = round(ymax)

boxWidth = 0.2
#Take the ceil to get the maximum possible without leaving out points
xBoxNum = int(math.ceil((xmax - xmin) / boxWidth))
yBoxNum = int(math.ceil((ymax - ymin) / boxWidth))

#Placeholder for smoothed data
plotGrid = np.zeros([xBoxNum, yBoxNum])

#Loop through the data and show the data in the boxes
yBoxStart = ymin
yBoxEnd = ymin + boxWidth
xBoxStart = xmin
xBoxEnd = xmin + boxWidth
for yInd in range(0, yBoxNum):
	for xInd in range(0, xBoxNum):
		
		#Find all data points that are within the current box
		xStartMatches = projected[:,0] >= xBoxStart
		xEndMatches = projected[:,0] <= xBoxEnd
		
		xMatches = xStartMatches * xEndMatches
		
		yStartMatches = projected[:,1] >= yBoxStart
		yEndMatches = projected[:,1] <= yBoxEnd
		
		yMatches = yStartMatches * yEndMatches
		
		dataInBox = projected[xMatches * yMatches]
		boxLabels = colorLabels[xMatches * yMatches]
		
		if len(dataInBox) > 0:
			#print dataInBox
			
			posCount = len(np.where(boxLabels == 'r')[0]) + 0.01
			negCount = len(np.where(boxLabels == 'b')[0]) + 0.01
			
			#Normalize for the total count of that label
			posCount = posCount / len(np.where(colorLabels == 'r')[0])
			negCount = negCount / len(np.where(colorLabels == 'b')[0])
			
			if negCount > 0:
				plotGrid[xInd,yInd] = np.log(posCount / float(negCount))
			

		#Move the box along x
		xBoxStart += boxWidth
		xBoxEnd += boxWidth
	
	yBoxStart += boxWidth
	yBoxEnd += boxWidth
	#Reset the box on x
	xBoxStart = xmin
	xBoxEnd = xmin + boxWidth

plotGrid = np.ma.masked_where(plotGrid == 0, plotGrid)
cmap = plt.cm.seismic
cmap.set_bad(color='white')
print(plotGrid)
plt.imshow(plotGrid, cmap=cmap, interpolation='nearest')		
plt.show()

#very simple ml test
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import auc, precision_recall_curve

from random import shuffle


def performCV(featureMatrix, labels, clf):
	
	cv = StratifiedKFold(n_splits=10)
	
	#dirty
	X = featureMatrix
	y = np.array(labels)
	
	tprs = []
	aucs = []
	auprcs = []
	scores = []
	mean_fpr = np.linspace(0, 1, 100)
	
	i = 0
	for train, test in cv.split(X, y):
		clf.fit(X[train], y[train])
		
		predictions = clf.predict(X[test])
		score = np.average(y[test] == np.sign(predictions))
		scores.append(score)
		
		probas_ = clf.predict_proba(X[test])
		
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		roc_auc = auc(fpr, tpr)
		aucs.append(roc_auc)

		precision, recall, thresholds = precision_recall_curve(y[test], predictions)
		aucScore = auc(recall, precision)
		auprcs.append(aucScore)
		
		i += 1
	
	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	mean_score = np.mean(scores)
	
	
	
	print("Score: ", mean_score)
	print("AUC: ", mean_auc)
	print("AUPRC: ", np.mean(auprcs))

print("Random forest")
from sklearn.ensemble import RandomForestClassifier
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
performCV(allFeatures, labels, rfClassifier)

print("Shuffled:")
shuffle(labels)
performCV(allFeatures, labels, rfClassifier)
