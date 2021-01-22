"""
	Plotting performance (AUC)-related figures, Figure 1C and 4C.
	
"""
import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
import os.path
import pandas as pd

def plotAUC(cancerTypes, outFileName, svTypes = ['DEL', 'DUP', 'INV', 'ITX']):

	#auto-read the AUC from the output folders
	aucs = dict()
	stdevs = dict()
	for cancerType in cancerTypes:
		#split the cancer types for clearer names in the plots
		splitCancerType = cancerType.split('_')
		cancerType2 = '_'.join(splitCancerType[1:2])

		if outFileName == 'figure4B':
			cancerType2 += '_CTCF'

		aucs[cancerType2] = []
		stdevs[cancerType2] = []

		for svType in svTypes:
			#read the AUC file
			outFile = 'output/' + cancerType + '/multipleInstanceLearning/leaveOnePatientOutCV/leaveOnePatientOutCV_' + svType + '_FINAL_AUC.txt'

			#skip runs for which there was no output due to e.g. no SVs
			if os.path.isfile(outFile) == True:
				aucData = np.loadtxt(outFile, dtype='object')

				aucData = aucData[0]
				stdev = float(aucData[1])
				svAuc = float(aucData[0])
				aucs[cancerType2].append(svAuc)
				stdevs[cancerType2].append(stdev)
			else:
				aucs[cancerType2].append(0)
				stdevs[cancerType2].append(0)

	#show the auc as dots in a scatterplot
	cancerTypeInd = 0
	plottedCancerTypes = []
	svTypeColors = ['#e41a1cff', '#377eb8ff', '#ff6000ff', '#984ea3ff']
	#make sure that the points do not overlap
	jitter = [-0.15, -0.05, 0.05, 0.15]
	plotData = []
	cancerTypeNames = []
	for cancerType in aucs:
		for svTypeInd in range(0, len(svTypes)):
			plotData.append([cancerTypeInd+jitter[svTypeInd], aucs[cancerType][svTypeInd], svTypeColors[svTypeInd], stdevs[cancerType][svTypeInd]])
		cancerTypeInd += 1


	data = pd.DataFrame(plotData)
	data.columns = ['cancer type', 'AUC', 'color', 'stdev']
	data = data.drop_duplicates()

	print(data)

	fig, ax = plt.subplots(1,1)
	plt.axhline(y=0.5, color='k', linestyle='--', linewidth=0.5)
	sns.scatterplot(data=data, x='cancer type', y='AUC', hue=data.color,
					cmap=ListedColormap(svTypeColors), legend=False,
					s = 60, edgecolor = 'k')

	#set separators
	ax.set_xticks(np.arange(0, len(aucs)-1)+0.5, minor=True)
	ax.grid(b=True, which='minor', linewidth=0.5, linestyle='--')

	plt.ylim([0.4,1])

	plt.xticks(np.arange(0, len(aucs)), list(aucs.keys()), rotation='vertical')
	plt.tight_layout()
	plt.savefig('output/figures/' + outFileName + '.svg')

	plt.show()

#1. Make figure 1C for the TAD-based results
cancerTypes = ['HMF_Breast_hmec', 'HMF_Ovary_ov', 'HMF_Lung_luad', 'HMF_Colorectal_coad',
			   'HMF_UrinaryTract_urinaryTract', 'HMF_Prostate_prostate', 'HMF_Esophagus_esophagus', 'HMF_Skin_skin',
			   'HMF_Pancreas_pancreas', 'HMF_Uterus_uterus', 'HMF_Kidney_kidney', 'HMF_NervousSystem_nervousSystem']

plotAUC(cancerTypes, 'figure1C')

cancerTypesCTCF = ['HMF_Breast_CTCF']
plotAUC(cancerTypesCTCF, 'figure4B')