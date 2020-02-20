"""
	Given our annotations for causal and non-causal SVs, check which features are important.

"""

import numpy as np
import sys
import re
import matplotlib.pyplot as plt

#1. load the data

degData = np.loadtxt(sys.argv[1], dtype='object')
nonDegData = np.loadtxt(sys.argv[2], dtype='object')

#Filter by SV type if required.  
svType = ''
nonDEGs = []
for degPair in nonDegData:
	
	sv = degPair[0].split("_")
	if svType != '':
		typeMatch = re.search(svType, sv[12], re.IGNORECASE)
		if typeMatch is None:
			continue
	nonDEGs.append(degPair)

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

positive = degs[:,1:]
negative = nonDEGs[:,1:]

allFeatures = np.concatenate((positive, negative), axis=0).astype(float)

labels = [1]*positive.shape[0] + [0]*negative.shape[0]

featuresWithLabel = []

for ind in range(0, allFeatures.shape[0]):
	
	newFeatures = list(allFeatures[ind])
	newFeatures += [labels[ind]]
	featuresWithLabel.append(newFeatures)

featuresWithLabel = np.array(featuresWithLabel)

for feature in range(0, featuresWithLabel.shape[1]):
	
	#print(feature, np.var(featuresWithLabel[:,feature]))
	print(feature, np.max(featuresWithLabel[:,feature]))

exit()

# 
# maxPlots = 10
# pos = 0
# for ind in range(0, featuresWithLabel.shape[1]):
# 	print(featuresWithLabel[:,ind])
# 	
# 	print('feature: ', ind)
# 	print('min: ', np.min(featuresWithLabel[:,ind]))
# 	print('max: ', np.max(featuresWithLabel[:,ind]))
# 	print('avg: ', np.mean(featuresWithLabel[:,ind]))
# 	
# 	if np.mean(featuresWithLabel[:,ind]) == 0:
# 		continue #skip this one, not applicable for the SV type.
# 	
# 	plt.boxplot(featuresWithLabel[:,ind], positions = pos)
# 	
# 	pos += 1
# 	
# 	if pos == maxPlots:
# 		plt.show()
# 		plt.clf()
# 		pos = 0
# exit()

#correlation maps
import pandas as pd
import numpy as np
import seaborn as sns


#get correlations of each features in dataset

data = pd.DataFrame(featuresWithLabel)
corrmat = data.corr()
top_corr_features = corrmat.index
plt.figure(figsize=(20,20))
#plot heat map
g=sns.heatmap(data[top_corr_features].corr(),annot=False,cmap="RdYlGn")
plt.show()