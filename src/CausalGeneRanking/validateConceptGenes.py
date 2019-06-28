"""
	
	The purpose of this script is to check if the genes in the concept are selected correctly.
	
	- What is the distance from the concept genes to all genes that are not in the concept? Is this distance larger to the negative genes than to the positive genes? 

"""

import numpy as np
import matplotlib.pyplot as plt
import sys
# 
# #Plot recurrence of concept genes
# conceptGenes = np.loadtxt(sys.argv[1], dtype="object")
# 
# counts = dict()
# 
# for gene in conceptGenes:
# 	
# 	if int(gene[1]) not in counts:
# 		counts[int(gene[1])] = 0
# 	counts[int(gene[1])] += 1
# 
# 
# plt.bar(counts.keys(), counts.values())
# plt.show()
# 
# # plt.bar(conceptGenes[:,0], conceptGenes[:,1])
# # plt.show()
# 
# 
# 
# exit()



# # # 
# bags = np.load("lassoPerPatient/bags.txt.npy")
# pairNames = np.load("lassoPerPatient/pairNames.txt.npy") #the sv-gene pair names of each bag entry
labels = np.load("Output/similarityMatrices/somaticGermlineLabels.txt.npy")
similarityMatrix = np.load("Output/similarityMatrices/somaticGermline.txt.npy")
pairNames = np.load("Output/similarityMatrices/somaticGermlinePairNames.txt.npy")

#non-coding only DEG pairs
degPairs = ['RPL7_brcaA07B', 'PSMB7_brcaA130', 'TTLL11_brcaA1LL','RABEPK_brcaA130', 'PCDHB6_brcaA0RE', 'MCM4_brcaA14X', 'VDAC3_brcaA2DF', 'STAU2_brcaA07B', 'MANEAL_brcaA04D', 'RIMS1_brcaA04X',
'UBE2V2_brcaA14X', 'GRK6_brcaA3XL', 'ZNF695_brcaA1FC', 'ITPKB_brcaA0RU', 'LAMC2_brcaA0BW', 'FAM3C_brcaA0D2', 'TBX19_brcaA0J4', 'TMED9_brcaA3XL', 'ZNF496_brcaA1FC', 'GPR20_brcaA14X',
'PI4KB_brcaA1P8', 'GNL2_brcaA0D0', 'SLC35F3_brcaA04D', 'PRELID1_brcaA3XL', 'SH3BGRL3_brcaA1LL', 'PRKDC_brcaA14X', 'CDC25A_brcaA1AY', 'EIF1B_brcaA27H', 'GYPA_brcaA3Y0', 'CDCA8_brcaA04D',
'ZNF124_brcaA1FC', 'ZNF727_brcaA0J4']




 # show the order of the concept genes in the ranking by nc score and in the DEGs
# 
# conceptGenes = np.loadtxt(sys.argv[1], dtype="object")
# rankedGenes = np.loadtxt(sys.argv[2], dtype="object")
# genePatientPairs = np.loadtxt(sys.argv[3], dtype="object")
# 
# # Process the deg pairs
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
# 	# degInd = np.where(sortedGeneCounts[:,0] == gene)[0]
# 	# 
# 	# if len(degInd) == 0: #if the gene is not a DEG
# 	# 	print gene
# 	# 	continue
# 	# 
# 	# conceptDegX.append(degInd[0])
# 	
# print conceptRankX
# # print conceptDegX
# 
# for ind in conceptRankX:
# 	plt.axvline(ind)
# 
# plt.xlim(0,rankedGenes.shape[0])
# plt.show()
# exit()
# plt.clf()

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

#Get subset of PCA

pca = PCA(n_components=2)

projected = pca.fit_transform(similarityMatrix)

rightPairs = []
leftPairs = []
thresholdX = 400
for x in range(0, projected.shape[0]):
	if projected[x,0] > 400:
		rightPairs.append(pairNames[x])
	else:
		leftPairs.append(pairNames[x])

# np.save('Output/similarityMatrices/uniqueSomaticGermlineRecurrentLeftPairs.txt', leftPairs)
# np.save('Output/similarityMatrices/uniqueSomaticGermlineRecurrentRightPairs.txt', rightPairs)
# exit()

#Do PCA
pca = PCA(n_components=2)

projected = pca.fit_transform(similarityMatrix)


# projectedWithOffset = projected
# 
# jitter = [0.01, -0.01]
# for row in range(0, projected.shape[0]):
# 	for col in range(0, projected.shape[1]):
# 		projectedWithOffset[row][col] += np.random.normal(-1, 1) * 25
# 
colorLabels = []

for label in labels:
	
	if label == 1:
		colorLabels.append('r')
	else:
		colorLabels.append('b')

fig,ax=plt.subplots(figsize=(7,5))
plt.scatter(projected[:, 0], projected[:, 1], c=colorLabels)

#plt.show()
colorLabels2 = []
pairInd = 0
added = 0
for pair in pairNames:
	splitPair = pair.split("_")
	gene = splitPair[0]
	sample = splitPair[len(splitPair)-1]
	shortPairName = gene + "_" + sample
	
	if shortPairName in degPairs:
		ax.scatter(projected[:,0][pairInd],projected[:,1][pairInd],c='yellow',s=40)
		#if projected[:,0][pairInd] > 2000:
		print pair
		added +=1
	
	pairInd += 1
	
print added

#plt.scatter(projected[:, 0], projected[:, 1], c=colorLabels)

#plt.show()

plt.show()

# for l in np.unique(labels):
#  ix=np.where(labels==l)
#  ax.scatter(projected[:,0][ix],projected[:,1][ix],c=cdict[l],s=40,
#            alpha=alpha[l])
# 
# plt.show()

exit()
# plt.scatter(projectedWithOffset[:, 0], projectedWithOffset[:, 1], c=colorLabels)
# plt.show()
# 
# alpha={0:.3, 1:.5}
# cdict={0:'blue',1:'red'}
# 
# fig,ax=plt.subplots(figsize=(7,5))
# fig.patch.set_facecolor('white')
# 
# for l in np.unique(labels):
#  ix=np.where(labels==l)
#  ax.scatter(projected[:,0][ix],projected[:,1][ix],c=cdict[l],s=40,
#            alpha=alpha[l])
# 
# plt.show()
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
print xmin, xmax, ymin, ymax

#round the values to get covering boxes
xmin = round(xmin)
xmax = round(xmax)
ymin = round(ymin)
ymax = round(ymax)

boxWidth = 100
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
print plotGrid
plt.imshow(plotGrid, cmap=cmap, interpolation='nearest')		
plt.show()
exit()


import scipy as sp
import scipy.ndimage
from matplotlib import pyplot as plt

sigma = [3, 3]
y = sp.ndimage.filters.gaussian_filter(projected, sigma, mode='constant')

print projected.shape
exit()

plt.scatter(y[:,0], y[:,1], c=colorLabels)
plt.show()

exit()

# Display filtered array
plt.imshow(y, cmap='hot', interpolation='nearest')
plt.show()
exit()



from matplotlib import cm as cm
from matplotlib import mlab as ml

#n = 1e5
#x = y = NP.linspace(-5, 5, 100)
x = projected[:, 0]
y = projected[:, 1]

#z = np.array([np.log(j/i) for j in y for i in x])

logScaled = []
for j in y:
	for i in x:
		
		#logScaled.append(-np.log(j/i))
		logScaled.append(j/float(i))
z = np.array(logScaled)

Z = z.reshape(x.shape[0], y.shape[0])
print Z
# print Z
plt.imshow(Z, cmap='hot', interpolation='nearest')
plt.show()
exit()

X, Y = np.meshgrid(x, y)
plt.imshow(Z,origin='lower',interpolation='bilinear')
plt.show()
exit()


# logScaled = np.log(y/x)
# z = np.meshgrid(logScaled)


# Z1 = ML.bivariate_normal(X, Y, 2, 2, 0, 0)
# Z2 = ML.bivariate_normal(X, Y, 4, 1, 1, 1)
# ZD = Z2 - Z1
#x = X.ravel()
#y = Y.ravel()
#z = Z.ravel()
# print x.shape
# print y.shape
# print z.shape
# print z

gridsize=30
PLT.subplot(111)

# if 'bins=None', then color of each hexagon corresponds directly to its count
# 'C' is optional--it maps values to x-y coordinates; if 'C' is None (default) then 
# the result is a pure 2D histogram 

PLT.hexbin(x, y, C=z, gridsize=gridsize, cmap=CM.jet, bins=None)
PLT.axis([x.min(), x.max(), y.min(), y.max()])

cb = PLT.colorbar()
cb.set_label('mean value')
PLT.show()   

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
