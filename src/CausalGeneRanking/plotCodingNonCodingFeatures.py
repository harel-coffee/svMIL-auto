"""
	Get the feature values of the coding & non-coding SVs, and see if these can be separated using PCA. 

"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from six.moves import range

nonCodingFeatureData = np.loadtxt(sys.argv[1], dtype='object')
codingFeatureData = np.loadtxt(sys.argv[2], dtype='object')

print(nonCodingFeatureData.shape)
print(codingFeatureData.shape)

nonCodingFeatures = nonCodingFeatureData[:,1:]
codingFeatures = codingFeatureData[:,1:]

#Make the bar plots
leftFeatures = nonCodingFeatures
rightFeatures = codingFeatures

eQTLLosses = leftFeatures[:,0].astype(float)
enhancerLosses = leftFeatures[:,1].astype(float)
promoterLosses = leftFeatures[:,2].astype(float)
cpgLosses = leftFeatures[:,3].astype(float)
tfLosses = leftFeatures[:,4].astype(float)
hicLosses = leftFeatures[:,5].astype(float)
h3k9me3Losses = leftFeatures[:,6].astype(float)
h3k4me3Losses = leftFeatures[:,7].astype(float)
h3k27acLosses = leftFeatures[:,8].astype(float)
h3k27me3Losses = leftFeatures[:,9].astype(float)
h3k4me1Losses = leftFeatures[:,10].astype(float)
h3k36me3Losses = leftFeatures[:,11].astype(float)
dnaseLosses = leftFeatures[:,12].astype(float)	
	

lossData = [np.sum(eQTLLosses), np.sum(enhancerLosses), np.sum(promoterLosses), np.sum(cpgLosses),
			np.sum(tfLosses), np.sum(hicLosses), np.sum(h3k9me3Losses), np.sum(h3k4me3Losses), np.sum(h3k27acLosses),
			np.sum(h3k27me3Losses), np.sum(h3k4me1Losses), np.sum(h3k36me3Losses), np.sum(dnaseLosses)]
lossData = np.array(lossData)
lossData = (lossData / float(leftFeatures.shape[0])) * 100
#lossData = -np.log(lossData)
print(lossData)

width = 0.35

plt.barh(np.arange(len(lossData)), lossData, width, label='Non-DEG pairs', color='blue')
# plt.xticks(range(0, len(lossData)),
	   # ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'], rotation=90)
# plt.show()

eQTLLosses = rightFeatures[:,0].astype(float)
enhancerLosses = rightFeatures[:,1].astype(float)
promoterLosses = rightFeatures[:,2].astype(float)
cpgLosses = rightFeatures[:,3].astype(float)
tfLosses = rightFeatures[:,4].astype(float)
hicLosses = rightFeatures[:,5].astype(float)
h3k9me3Losses = rightFeatures[:,6].astype(float)
h3k4me3Losses = rightFeatures[:,7].astype(float)
h3k27acLosses = rightFeatures[:,8].astype(float)
h3k27me3Losses = rightFeatures[:,9].astype(float)
h3k4me1Losses = rightFeatures[:,10].astype(float)
h3k36me3Losses = rightFeatures[:,11].astype(float)
dnaseLosses = rightFeatures[:,12].astype(float)

lossData = [np.sum(eQTLLosses), np.sum(enhancerLosses), np.sum(promoterLosses), np.sum(cpgLosses),
			np.sum(tfLosses), np.sum(hicLosses), np.sum(h3k9me3Losses), np.sum(h3k4me3Losses), np.sum(h3k27acLosses),
			np.sum(h3k27me3Losses), np.sum(h3k4me1Losses), np.sum(h3k36me3Losses), np.sum(dnaseLosses)]
lossData = np.array(lossData)
lossData = (lossData / float(rightFeatures.shape[0])) * 100
#lossData = -np.log(lossData)
print(lossData)

plt.barh(np.arange(len(lossData)) + width, lossData, width, label='DEG pairs', color='red')
plt.yticks(np.arange(len(lossData) + width / 2),
		   ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'])
plt.xlim([0,100])
plt.legend(loc='best')
plt.tight_layout()
#plt.show()
plt.savefig('Output/degNonDeg_losses.svg')
plt.clf()
#Gains

eQTLGains = leftFeatures[:,13].astype(float)
enhancerGains = leftFeatures[:,14].astype(float)
promoterGains = leftFeatures[:,15].astype(float)
cpgGains = leftFeatures[:,16].astype(float)
tfGains = leftFeatures[:,17].astype(float)
hicGains = leftFeatures[:,18].astype(float)
h3k9me3Gains = leftFeatures[:,19].astype(float)
h3k4me3Gains = leftFeatures[:,20].astype(float)
h3k27acGains = leftFeatures[:,21].astype(float)
h3k27me3Gains = leftFeatures[:,22].astype(float)
h3k4me1Gains = leftFeatures[:,23].astype(float)
h3k36me3Gains = leftFeatures[:,24].astype(float)
dnaseGains = leftFeatures[:,25].astype(float)

gainData = [np.sum(eQTLGains), np.sum(enhancerGains), np.sum(promoterGains), np.sum(cpgGains),
			np.sum(tfGains), np.sum(hicGains), np.sum(h3k9me3Gains), np.sum(h3k4me3Gains), np.sum(h3k27acGains),
			np.sum(h3k27me3Gains), np.sum(h3k4me1Gains), np.sum(h3k36me3Gains), np.sum(dnaseGains)]

gainData = np.array(gainData)
gainData = (gainData / float(leftFeatures.shape[0])) * 100
#gainData = -np.log(gainData)
print(gainData)

plt.barh(np.arange(len(gainData)), gainData, width, label='Non-DEG pairs', color='blue')
# plt.xticks(np.arange(len(gainData)),
# 		   ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'], rotation=90)
# plt.show()
	
eQTLGains = rightFeatures[:,13].astype(float)
enhancerGains = rightFeatures[:,14].astype(float)
promoterGains = rightFeatures[:,15].astype(float)
cpgGains = rightFeatures[:,16].astype(float)
tfGains = rightFeatures[:,17].astype(float)
hicGains = rightFeatures[:,18].astype(float)
h3k9me3Gains = rightFeatures[:,19].astype(float)
h3k4me3Gains = rightFeatures[:,20].astype(float)
h3k27acGains = rightFeatures[:,21].astype(float)
h3k27me3Gains = rightFeatures[:,22].astype(float)
h3k4me1Gains = rightFeatures[:,23].astype(float)
h3k36me3Gains = rightFeatures[:,24].astype(float)
dnaseGains = rightFeatures[:,25].astype(float)

gainData = [np.sum(eQTLGains), np.sum(enhancerGains), np.sum(promoterGains), np.sum(cpgGains),
			np.sum(tfGains), np.sum(hicGains), np.sum(h3k9me3Gains), np.sum(h3k4me3Gains), np.sum(h3k27acGains),
			np.sum(h3k27me3Gains), np.sum(h3k4me1Gains), np.sum(h3k36me3Gains), np.sum(dnaseGains)]

gainData = np.array(gainData)
gainData = (gainData / float(rightFeatures.shape[0])) * 100
#gainData = -np.log(gainData)

print(gainData)

plt.barh(np.arange(len(gainData)) + width, gainData, width, label='DEG pairs',color='red')
plt.yticks(np.arange(len(gainData) + width / 2),
		   ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'])

plt.xlim([0,100])
plt.legend(loc='best')
plt.tight_layout()
#plt.show()
plt.savefig('Output/degNonDeg_gains.svg')


exit()




allFeatures = np.concatenate((nonCodingFeatures, codingFeatures), axis=0)

#Make a PCA plot for the left/right set and see if these are really different

from sklearn.decomposition import PCA

#Get subset of PCA

pca = PCA(n_components=2)

projected = pca.fit_transform(allFeatures)
# projectedWithOffset = projected
# 
# jitter = [0.01, -0.01]
# for row in range(0, projected.shape[0]):
# 	for col in range(0, projected.shape[1]):
# 		projectedWithOffset[row][col] += np.random.normal(-1, 1) * 0.1
# 		
# projected = projectedWithOffset

colorLabels = []

for i in range(0, allFeatures.shape[0]):
	
	if i < nonCodingFeatures.shape[0]:
		colorLabels.append('b')
	elif i >= nonCodingFeatures.shape[0] and i < (nonCodingFeatures.shape[0] + codingFeatures.shape[0]):
		colorLabels.append('r')

fig,ax=plt.subplots(figsize=(7,5))
plt.scatter(projected[:, 0], projected[:, 1], c=colorLabels)
plt.show()

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
print(xmin, xmax, ymin, ymax)

#round the values to get covering boxes
xmin = round(xmin)
xmax = round(xmax)
ymin = round(ymin)
ymax = round(ymax)

boxWidth = 1
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
print(plotGrid)
plt.imshow(plotGrid, cmap=cmap, interpolation='nearest')		
plt.show()
