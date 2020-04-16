import numpy as np
import sys
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

outDir = sys.argv[1]

svTypes = ['DEL', 'DUP', 'INV', 'ITX']

for svType in svTypes:
	
	#load the full similarity matrix for this SV type, and also the labels
	dataPath = outDir + '/multipleInstanceLearning/similarityMatrices/'
	similarityMatrix = np.load(dataPath + '/similarityMatrix_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	bagLabels = np.load(dataPath + '/bagLabels_' + svType + '.npy', encoding='latin1', allow_pickle=True)

	
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
	
	X_train, X_test, y_train, y_test = train_test_split(similarityMatrix, bagLabels, test_size=0.33, random_state=42)
	rf_random.fit(X_train, y_train)
	print('best params; ')
	print(rf_random.best_params_)
	print('new score: ')
	print(rf_random.score(X_test, y_test))
	
	print('base score :')
	rfClassifier.fit(X_train, y_train)
	print(rfClassifier.score(X_test, y_test))