## for eah cancer type, get the performance from all swap folders
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

cancerTypes = ['Breast', 'Colorectal', 'Esophagus', 'Lung', 'Kidney', 'NervousSystem', 'Ovary',
			   'Pancreas', 'Prostate', 'Skin', 'UrinaryTract', 'Uterus']
baselines = ['hmec', 'coad', 'esophagus', 'luad', 'kidney', 'nervousSystem', 'ov',
			 'pancreas', 'prostate', 'skin', 'urinaryTract', 'uterus']

allTypes = ['hmec', 'coad', 'esophagus', 'luad', 'kidney', 'nervousSystem', 'ov', 'pancreas',
			'prostate', 'skin', 'urinaryTract', 'uterus', 'gm12878']

svTypes = ['DEL', 'DUP', 'INV', 'ITX']

differencesAcrossCancerTypes = []
for cancerTypeInd in range(0, len(cancerTypes)):
	cancerType = cancerTypes[cancerTypeInd]
	baseline = baselines[cancerTypeInd]

	aucs = dict()
	for swapType in allTypes:

		swapFolder = 'HMF_' + cancerType + '_' + swapType
		aucs[swapFolder] = []

		for svType in svTypes:
			outFile =  'output/' + swapFolder + '/multipleInstanceLearning/leaveOnePatientOutCV/leaveOnePatientOutCV_' + svType + '_FINAL_AUC.txt'

			#skip runs for which there was no output due to e.g. no SVs
			if os.path.isfile(outFile) == True:
				aucData = np.loadtxt(outFile, dtype='object')
				#if cancerType == 'HMF_Uterus':
				aucData = aucData[0]
				svAuc = float(aucData[0])
				aucs[swapFolder].append(svAuc)
			else:
				aucs[swapFolder].append(0)

	#compare the performance to the baseline
	baseData = 'HMF_' + cancerType + '_' + baseline
	baseAUCs = aucs[baseData]

	allDifferences = []
	for swapType in allTypes:

		swapFolder = 'HMF_' + cancerType + '_' + swapType

		differenceFromBase = 0
		for svTypeInd in range(0, len(svTypes)):
			baseAUC = baseAUCs[svTypeInd]
			swapAUC = aucs[swapFolder][svTypeInd]

			if swapAUC == 0:

				continue

			differenceFromBase += (swapAUC - baseAUC)


		#allDifferences[swapFolder] = differenceFromBase

		allDifferences.append(differenceFromBase)

	zDifferences = []
	for difference in allDifferences:

		if np.mean(allDifferences) == 0 or np.std(allDifferences) == 0:
			zDifferences.append(0)
			continue

		z = (difference - np.mean(allDifferences)) / np.std(allDifferences)

		if z >= 1 and z < 2:
			zDifferences.append(1)
		elif z >= 2:
			zDifferences.append(2)
		elif z <= -1 and z > -2:
			zDifferences.append(-1)
		elif z <= -2:
			zDifferences.append(-2)
		else:
			zDifferences.append(0)

	differencesAcrossCancerTypes.append(zDifferences)

#make a plot showing the differences

fig =plt.figure(figsize=(5,5))

data = pd.DataFrame(differencesAcrossCancerTypes)
g=sns.heatmap(data,annot=False,square=True, linewidths=0.5,
			  xticklabels=allTypes, yticklabels=cancerTypes,
			  cmap="vlag", center=0, vmin=-2, vmax=2)
plt.tight_layout()
plt.savefig('output/figures/figure4a.svg')
plt.show()