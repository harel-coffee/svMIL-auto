"""
	Plot the ROC and PR curves for the permuted labels in the bag classifier

"""

from __future__ import absolute_import
import sys
import matplotlib.pyplot as plt
import numpy as np
import os
from os import listdir
from os.path import isfile, join

from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from scipy import interp
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.utils.fixes import signature

folder = "bagMIL/"
#Get all the files from this folder

scoreFiles = [f for f in listdir(folder) if isfile(join(folder, f))]

aucs = []
accuracies = []
tprs = []
trueLabels = []
truePredictions = []

for scoreFile in scoreFiles:
	if scoreFile != "true.txt":
		continue
	
	lineCount = 0
	with open(folder + scoreFile, 'r') as inF:
		
		for line in inF:
			
			line = line.strip()
			if lineCount == 0:
				aucs.append(float(line))
			if lineCount == 1:
				accuracies.append(float(line))
			if lineCount == 2:
				splitLine = line.split(",")
				tprs.append([float(i) for i in splitLine])
			if lineCount == 3:
				splitLine = line.split(",")
				trueLabels.append([float(i) for i in splitLine])
			if lineCount == 4:
				splitLine = line.split(",")
				truePredictions.append([float(i) for i in splitLine])
						
			lineCount += 1	
		

totalAccuracy = np.mean(accuracies)
#Make ROC curve

plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
mean_fpr = np.linspace(0, 1, 100)
mean_tpr = np.mean(tprs, axis=0)
#mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
plt.plot(mean_fpr, mean_tpr, color='b',
         label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
         lw=2, alpha=.8)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                 label=r'$\pm$ 1 std. dev.')

plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
title = 'Receiver operating characteristic, average accuracy: ' + str(totalAccuracy)
plt.title(title)
plt.legend(loc="lower right")
plt.show()
exit()

#Make PR curve

testReal = np.concatenate(trueLabels)
testProbs = np.concatenate(truePredictions)

precision, recall, _ = precision_recall_curve(testReal, testProbs)

# In matplotlib < 1.5, plt.fill_between does not have a 'step' argument
step_kwargs = ({'step': 'post'}
               if 'step' in signature(plt.fill_between).parameters
               else {})

plt.step(recall, precision, color='b', alpha=0.2,
         where='post')
plt.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
#plt.title('2-class Precision-Recall curve')
plt.show()






