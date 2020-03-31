"""
	Interface for Multiple Instance Learning on SV-gene pairs to classify these as having transcriptional consequences or not. 
	
	To do:
	- move parts to separate functions
	- move things around to settings
	- feature elimination
	- feature importance to learn more about biology

"""

import numpy as np
import sys
from inputParser import InputParser
from featureMatrixDefiner import FeatureMatrixDefiner
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import math
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import Lasso
from sklearn.metrics import auc, precision_recall_curve
from sklearn.svm import LinearSVC
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from scipy import interp
import pickle as pkl
from random import shuffle
import random
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from copy import deepcopy

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold
from inspect import signature
from sklearn import model_selection
from sklearn.metrics import average_precision_score
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import train_test_split

import pandas as pd

"""

	TO DO: update this code into functions!!

"""

#settings for running in different scenarios
svTypes = [sys.argv[3]]
#svTypes = ['DEL', 'DUP', 'INV', 'ITX']
#svTypes = ['DEL', 'INV', 'ITX']
svTypes = ['ITX']
normalize = False #re-normalize, or use the saved file for speed? 
optimize = False #optimize classifier? 
test = True #test classifier performance with CV?
featureSelection = True #randomize features 1 by 1? 
featureImportance = False
featureLoad = False #re-create feature importances or just load from file? 

adjustedPValues = dict()
allFeatureZScores = dict()


if featureLoad == False:

			#do RF optimization with random parameter search
			if optimize == True:
				#Number of trees in random forest
				n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
				# Number of features to consider at every split
				max_features = ['auto', 'sqrt']
				# Maximum number of levels in tree
				max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
				max_depth.append(None)
				# Minimum number of samples required to split a node
				min_samples_split = [2, 5, 10]
				# Minimum number of samples required at each leaf node
				min_samples_leaf = [1, 2, 4]
				# Method of selecting samples for training each tree
				bootstrap = [True, False]# Create the random grid
				random_grid = {'n_estimators': n_estimators,
							   'max_features': max_features,
							   'max_depth': max_depth,
							   'min_samples_split': min_samples_split,
							   'min_samples_leaf': min_samples_leaf,
							   'bootstrap': bootstrap}
				
				#initial classifier to optimize from
				rfClassifier = RandomForestClassifier(n_estimators= 500)
				rf_random = RandomizedSearchCV(estimator = rfClassifier, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)
				
				X_train, X_test, y_train, y_test = train_test_split(similarityMatrices[featureStart], bagLabels, test_size=0.33, random_state=42)
				rf_random.fit(X_train, y_train)
				print('best params; ')
				print(rf_random.best_params_)
				print('new score: ')
				print(rf_random.score(X_test, y_test))
				
				print('base score :')
				rfClassifier.fit(X_train, y_train)
				print(rfClassifier.score(X_test, y_test))
				exit()
				
				
			#some results of optimization: 
			#best DEL: n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True
			#best ITX: n_estimators= 1000, min_samples_split=2, min_samples_leaf=1, max_features='sqrt', max_depth=110, bootstrap=True
			#best inv: n_estimators= 400, min_samples_split=10, min_samples_leaf=2, max_features='auto', max_depth=40, bootstrap=False
			#best DUP: n_estimators= 1200, min_samples_split=2, min_samples_leaf=4, max_features='sqrt', max_depth=70, bootstrap=False
			
			#run2
			
			#best itx: n_estimators= 1000, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True
			#best inv: {'n_estimators': 200, 'min_samples_split': 5, 'min_samples_leaf': 4, 'max_features': 'auto', 'max_depth': 10, 'bootstrap': True}
			#best dup: 'n_estimators': 600, 'min_samples_split': 2, 'min_samples_leaf': 2, 'max_features': 'sqrt', 'max_depth': 110, 'bootstrap': False
			#best del:
			