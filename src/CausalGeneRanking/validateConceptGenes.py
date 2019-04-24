"""
	
	The purpose of this script is to check if the genes in the concept are selected correctly.
	
	- What is the distance from the concept genes to all genes that are not in the concept? Is this distance larger to the negative genes than to the positive genes? 

"""

import numpy as np
import matplotlib.pyplot as plt


bags = np.load("bags.txt.npy")
pairNames = np.load("pairNames.txt.npy") #the sv-gene pair names of each bag entry
labels = np.load("labels.txt.npy")
similarityMatrix = np.load("similarityMatrix.txt.npy")
conceptIndices = np.load("conceptIndices.txt.npy")



#For every concept gene index

#Get the similarity matrix distances from every concept to all the positive bags & the negative bags

conceptSimilarities = similarityMatrix[:,conceptIndices]

#Split by positive and negative bags
positiveIndices = np.where(np.array(labels) == 1)[0]
negativeIndices = np.where(np.array(labels) == 0)[0]

positiveConceptSimilarities = conceptSimilarities[positiveIndices,:]
negativeConceptSimilarities = conceptSimilarities[negativeIndices,:]

#Average the distance across the concept genes, show the distribution across the bags.

distanceToBagsPositive = np.mean(positiveConceptSimilarities, axis=1)
distanceToBagsNegative = np.mean(negativeConceptSimilarities, axis=1)

#Show in a boxplot
plt.figure()
plt.boxplot([distanceToBagsPositive, distanceToBagsNegative])
plt.show()


#Make a second plot where we show the distribution across the concept genes, and average across the bags. 

distanceToBagsPositive = np.mean(positiveConceptSimilarities, axis=0)
distanceToBagsNegative = np.mean(negativeConceptSimilarities, axis=0)

#Show in a boxplot
plt.figure()
plt.boxplot([distanceToBagsPositive, distanceToBagsNegative])
plt.show()
