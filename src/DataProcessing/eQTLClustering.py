import sys
import numpy as np


#read eQTL file
eQTLData = np.loadtxt('breast_eQTLs.bed', dtype='object')

#group the eQTLs by genes to make it easier to access them.
perGeneEQTLs = dict()
for eQTL in eQTLData:

	gene = eQTL[3]
	if gene not in perGeneEQTLs:
		perGeneEQTLs[gene] = []
	perGeneEQTLs[gene].append([eQTL[0], int(eQTL[1]), int(eQTL[2]), eQTL[3]])

#convert to np arrays
for gene in perGeneEQTLs:

	perGeneEQTLs[gene] = np.array(perGeneEQTLs[gene], dtype='object')

#test with 1 chromsome first
#per gene, make a sliding window.

#if we find another eQTL for that same gene within 10kb, add it to the cluster.

windowSize = 1000
clusters = dict()
currentCluster = []
geneCount = 0
for gene in perGeneEQTLs:

	#get all other eQTLs from this gene.
	#print(gene)
	#if geneCount > 2:
	#	break

	geneClusters = []

	#if the next eQTL is within the window size, append it to the cluster.
	for eQTLInd in range(0, perGeneEQTLs[gene].shape[0]):

		eQTL = perGeneEQTLs[gene][eQTLInd]

		#if this is the first eQTL, it is always the start of a cluster.
		if eQTLInd == 0:
			currentCluster.append([eQTL[0], eQTL[1], eQTL[2], eQTL[3]])
		else:
			#otherwise, check if the distance to the previous is within the window size.
			previousEQTL = perGeneEQTLs[gene][eQTLInd-1]
			if eQTL[1] - previousEQTL[1] <= windowSize:
				currentCluster.append([eQTL[0], eQTL[1], eQTL[2], eQTL[3]])
			else:
				geneClusters.append(currentCluster)
				currentCluster = []
				currentCluster.append([eQTL[0], eQTL[1], eQTL[2], eQTL[3]])

	#remaining eQTLs
	geneClusters.append(currentCluster)
	geneClusters = np.array(geneClusters, dtype='object')
	currentCluster = []
	
	clusters[gene] = geneClusters
	
	geneCount += 1
	
#write clusters to a new file
clusterFile = 'clusters_test.txt'
with open(clusterFile, 'w') as outF:
	for gene in clusters:

		for cluster in clusters[gene]:

			clusterStart = 0
			clusterEnd = 0
			eQTLInd = 0
			for eQTL in cluster:
				if eQTLInd == 0:
					clusterStart = eQTL[1]
					eQTLInd += 1
				clusterEnd = eQTL[2]


			chrom = cluster[0][0]
			gene = cluster[0][3]
			
			clusterLine = chrom + '\t' + str(clusterStart) + '\t' + str(clusterEnd) + '\t' + gene + '\n'
			outF.write(clusterLine)



	

