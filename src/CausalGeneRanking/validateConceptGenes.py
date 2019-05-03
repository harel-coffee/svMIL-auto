"""
	
	The purpose of this script is to check if the genes in the concept are selected correctly.
	
	- What is the distance from the concept genes to all genes that are not in the concept? Is this distance larger to the negative genes than to the positive genes? 

"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# 
bags = np.load("SomaticGermline/bags.txt.npy")
pairNames = np.load("SomaticGermline/pairNames.txt.npy") #the sv-gene pair names of each bag entry
labels = np.load("SomaticGermline/labels.txt.npy")
similarityMatrix = np.load("SomaticGermline/similarityMatrix.txt.npy")
conceptIndices = np.load("SomaticGermline/conceptIndices.txt.npy")


#Perform some clustering on the similarity matrix

# #from scipy.cluster.hierarchy import dendrogram, linkage
# import fastcluster
# from matplotlib import pyplot as plt
# 
# linked = fastcluster.linkage(similarityMatrix, 'single')
# 
# print similarityMatrix.shape
# print len(pairNames)
# exit()
# 
# 
# plt.figure(figsize=(10, 7))  
# dendrogram(linked,  
#             orientation='top',
#             labels=labelList,
#             distance_sort='descending',
#             show_leaf_counts=True)
# plt.show() 


# Compute the variance in the similarity matrix
# 
# variance = np.var(similarityMatrix, axis=0) #get the variance for every instance
# 
# print variance
# print len(variance)
# print similarityMatrix.shape
# print np.sort(variance)

positiveIndices = np.where(np.array(labels) == 1)[0]
negativeIndices = np.where(np.array(labels) == -1)[0]

positiveBagVariance = np.var(similarityMatrix[positiveIndices,:], axis=0)
negativeBagVariance = np.var(similarityMatrix[negativeIndices,:], axis=0)

print np.sort(positiveBagVariance)
print np.sort(negativeBagVariance)

#compute the variability index
varIndexPos = positiveBagVariance / np.mean(similarityMatrix[positiveIndices,:], axis=0)
varIndexNeg = negativeBagVariance / np.mean(similarityMatrix[negativeIndices,:], axis=0)


plt.figure()
plt.boxplot([positiveBagVariance, negativeBagVariance])
plt.show()

plt.figure()
plt.boxplot([varIndexPos, varIndexNeg])
plt.show()

exit()


#Make a plot where we show the order of the concept genes in the ranking by nc score and in the DEGs
# 
# conceptGenes = np.loadtxt(sys.argv[1], dtype="object")
# rankedGenes = np.loadtxt(sys.argv[2], dtype="object")
# genePatientPairs = np.loadtxt(sys.argv[3], dtype="object")
# 
# #Process the deg pairs
# geneCounts = dict()
# for pair in genePatientPairs[:,0]:
# 	
# 	splitPair = pair.split("_")
# 	
# 	if splitPair[0] not in geneCounts:
# 		geneCounts[splitPair[0]] = 0
# 	geneCounts[splitPair[0]] += 1
# 	
# geneCountsArray = np.empty([len(geneCounts), 2], dtype="object")
# 
# for geneCountInd in range(0, len(geneCounts.keys())):
# 	
# 	geneCountsArray[geneCountInd,0] = geneCounts.keys()[geneCountInd]
# 	geneCountsArray[geneCountInd,1] = geneCounts.values()[geneCountInd]
# 	
# sortedGeneCounts = geneCountsArray[geneCountsArray[:,1].argsort()][::-1]
# 
# conceptRankX = []
# conceptDegX = []
# for gene in conceptGenes[:,0]:
# 	geneInd = np.where(rankedGenes[:,0] == gene)[0][0]
# 	conceptRankX.append(geneInd)
# 	
# 	degInd = np.where(sortedGeneCounts[:,0] == gene)[0]
# 	
# 	if len(degInd) == 0: #if the gene is not a DEG
# 		print gene
# 		continue
# 	
# 	conceptDegX.append(degInd[0])
# 	
# print conceptRankX
# print conceptDegX
# # 
# # for ind in conceptRankX:
# # 	plt.axvline(ind)
# # 
# # plt.xlim(0,rankedGenes.shape[0])
# # plt.show()
# # exit()
# # plt.clf()
# 
# for ind in conceptDegX:
# 	plt.axvline(ind)
# 
# plt.xlim(0,sortedGeneCounts[:,0].shape[0])
# plt.show()
# 
# exit()

# 
# # 
from sklearn.decomposition import PCA

pca = PCA(n_components=2)

projected = pca.fit_transform(similarityMatrix)
projectedWithOffset = projected

jitter = [0.01, -0.01]
for row in range(0, projected.shape[0]):
	for col in range(0, projected.shape[1]):
		projectedWithOffset[row][col] += np.random.normal(-1, 1) * 25

colorLabels = []

for label in labels:
	
	if label == 1:
		colorLabels.append('r')
	else:
		colorLabels.append('b')

#plt.scatter(projected[:, 0], projected[:, 1], c=colorLabels)
# plt.scatter(projectedWithOffset[:, 0], projectedWithOffset[:, 1], c=colorLabels)
# plt.show()

alpha={0:.3, 1:.5}
cdict={0:'blue',1:'red'}

fig,ax=plt.subplots(figsize=(7,5))
fig.patch.set_facecolor('white')

for l in np.unique(labels):
 ix=np.where(labels==l)
 ax.scatter(projected[:,0][ix],projected[:,1][ix],c=cdict[l],s=40,
           alpha=alpha[l])

plt.show()
#rasterize the PCA plot and make a density heatmap

exit()




# # 
# from tsne import bh_sne
# 
# colorLabels = []
# posClass = []
# for label in labels:
# 	
# 	if label == 1:
# 		colorLabels.append('r')
# 	else:
# 		colorLabels.append('b')
# 
# 
# vis_data = bh_sne(similarityMatrix)
# 
# # plot the result
# vis_x = vis_data[:, 0]
# vis_y = vis_data[:, 1]
# 
# 
# 
# plt.scatter(vis_x, vis_y, c=colorLabels)
# plt.show()
# 
# exit()



#For every concept gene index


#Get the similarity matrix distances from every concept to all the positive bags & the negative bags

conceptSimilarities = similarityMatrix[:,conceptIndices]

#Split by positive and negative bags
positiveIndices = np.where(np.array(labels) == 1)[0]
negativeIndices = np.where(np.array(labels) == -1)[0]

print len(positiveIndices)
print len(negativeIndices)
exit()
positiveConceptSimilarities = conceptSimilarities[positiveIndices,:]
negativeConceptSimilarities = conceptSimilarities[negativeIndices,:]
# 
# #Average the distance across the concept genes, show the distribution across the bags.

distanceToBagsPositive = np.mean(positiveConceptSimilarities, axis=1)
distanceToBagsNegative = np.mean(negativeConceptSimilarities, axis=1)

#Show in a boxplot
plt.figure()
plt.boxplot([distanceToBagsPositive, distanceToBagsNegative])
plt.show()


#Make a second plot where we show the distribution across the concept genes, and average across the bags. 
# 
# distanceToBagsPositive = np.mean(positiveConceptSimilarities, axis=0)
# distanceToBagsNegative = np.mean(negativeConceptSimilarities, axis=0)
# 
# #Show in a boxplot
# plt.figure()
# plt.boxplot([distanceToBagsPositive, distanceToBagsNegative])
# plt.show()
