"""

	Here the MIL classifier performance will be tested on the input similarity matrices
	This can work with a single similarity matrix, but can also be run multiple times
	in the feature elimination way.

"""

import sys
import numpy as np

featureElimination = False
svTypes = ['ITX']
outDir = sys.argv[3]

def cvClassification(similarityMatrix, bagLabels, clf, svType, title, plot):

	#get the kfold model
	kfold = model_selection.StratifiedKFold(n_splits=10, shuffle=True, random_state=10)

	#make the ROC curve and compute the AUC
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)

	fig, ax = plt.subplots()
	importances = []
	for i, (train, test) in enumerate(kfold.split(similarityMatrix, bagLabels)):
		clf.fit(similarityMatrix[train], bagLabels[train])
		viz = plot_roc_curve(clf, similarityMatrix[test], bagLabels[test],
							 name='ROC fold {}'.format(i),
							 alpha=0.3, lw=1, ax=ax)
		interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(viz.roc_auc)
		importances.append(clf.feature_importances_)
	print('aucs: ')
	print(aucs)
	print('mean auc: ', np.mean(aucs))
	print('std of auc: ', np.std(aucs))

	if plot == True:

		ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
				label='Chance', alpha=.8)

		mean_tpr = np.mean(tprs, axis=0)
		mean_tpr[-1] = 1.0
		mean_auc = auc(mean_fpr, mean_tpr)
		std_auc = np.std(aucs)

		ax.plot(mean_fpr, mean_tpr, color='b',
				label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (np.mean(aucs), np.std(aucs)),
				lw=2, alpha=.8)

		std_tpr = np.std(tprs, axis=0)
		tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
		tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
		ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
						label=r'$\pm$ 1 std. dev.')

		ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
			   title="Receiver operating characteristic: " + title)
		ax.legend(loc="lower right")
		plt.tight_layout()
		plt.savefig('miles_' + svType + '.svg')
		plt.show()

for svType in svTypes:

	if svType == 'DEL':
		classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'deletions'
	elif svType == 'DUP':
		#classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=2, min_samples_leaf=2, max_features='sqrt', max_depth=110, bootstrap=False)
		classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		#classifier = RandomForestClassifier(n_estimators= 1200, min_samples_split=2, min_samples_leaf=4, max_features='sqrt', max_depth=70, bootstrap=False)
		title = 'duplications'
	elif svType == 'INV':
		classifier = RandomForestClassifier(n_estimators= 200, min_samples_split=5, min_samples_leaf=4, max_features='auto', max_depth=10, bootstrap=True)
		title = 'inversions'
	elif svType == 'ITX':
		classifier = RandomForestClassifier(n_estimators= 1000, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'translocations'
	else:
		classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'All SV types'

	#obtain the right similarity matrix and bag labels
	
	dataPath = outDir + '/multipleInstanceLearning/similarityMatrices/'
	similarityMatrix = np.load(dataPath + '/similarityMatrix_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	bagLabels = np.load(dataPath + '/bagLabels' + svType + '.npy', encoding='latin1', allow_pickle=True)
	
	plot = True
	if featureSelection == True:
		plot = False

	cvClassification(similarityMatrix, bagLabels, classifier, svType, title, plot)

	#repeat, but then with random labels.
	#mind here, if multiple iterations, the baglabels are permanently shuffled!
	if featureSelection == False:
		shuffle(bagLabels)

		cvClassification(similarityMatrix, bagLabels, classifier, svType, title, plot)

