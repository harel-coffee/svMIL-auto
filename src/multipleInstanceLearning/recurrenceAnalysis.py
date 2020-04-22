"""
	This script serves to do recurrence analysis on the sv-gene pairs identified
	
	We do the following things:
	1. From the top 100 of each SV type (so top 400), which genes are there? Which are the top 15 most recurrent?
	2. For these genes, also check which other mutations are found in these genes in different patients.
	3. Then also check which genes are recurrent if we ignore the top 100, and just look across all positive SV-gne pairs.

"""

import sys
import numpy as np
import random
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.pyplot as plt
import os

import matplotlib
matplotlib.use('Agg')

outDir = sys.argv[1]

finalOutDir = outDir + '/figure4/'
if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)
finalOutDirFullFigure = outDir + '/figS5/'
if not os.path.exists(finalOutDirFullFigure):
	os.makedirs(finalOutDirFullFigure)

#load the sv-gene pairs
positivePairs = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_pathogenicPairsFeatures.txt', dtype='object')
print(positivePairs.shape)

topPairs = dict()
topPairGenes = dict() #all unique genes
svTypes = ['DEL', 'DUP', 'INV', 'ITX']
for svType in svTypes:
	svPairs = np.loadtxt(outDir + '/featureImportance/pairLabels_top100_' + svType + '.txt', dtype='object')

	rankedPairs = []
	ind = len(svPairs)
	for pair in svPairs:
		splitPair = pair.split('_')
		topPairGenes[splitPair[0]] = 0
		rankedPairs.append([pair, ind])
		ind -= 1
	topPairs[svType] = rankedPairs

degPairs = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object')

#format: gene as key, as values, the first is the number of times we found the gene,
#the second how many were dels, then dups, invs, itx.
splitPairs = dict()
genes = dict()
for pair in positivePairs: #get stats for all pairs
	splitPair = pair[0].split('_')


	if splitPair[0] + '_' + splitPair[7] not in splitPairs:
		splitPairs[splitPair[0] + '_' + splitPair[7]] = []
	splitPairs[splitPair[0] + '_' + splitPair[7]].append(pair)

	if splitPair[0] not in genes:

		#count, cross-patient count, nc: del, dup, inv, itx, snv, cnv amp, cnv del, sv del, sv dup, sv inv, sv itx
		genes[splitPair[0]] = [0]*13
		genes[splitPair[0]].append([]) #add the patient names here
		genes[splitPair[0]].append(0) #negative dels
		genes[splitPair[0]].append(0) #negative dups
		genes[splitPair[0]].append(0) #negative invs
		genes[splitPair[0]].append(0) #negative itx

	genes[splitPair[0]][0] += 1

	if splitPair[12] == 'DEL':
		genes[splitPair[0]][1] += 1
	elif splitPair[12] == 'DUP':
		genes[splitPair[0]][2] += 1
	elif splitPair[12] == 'INV':
		genes[splitPair[0]][3] += 1
	elif splitPair[12] == 'ITX':
		genes[splitPair[0]][4] += 1

	patient = splitPair[7]
	genes[splitPair[0]][13].append(patient)

#convert to numpy array for easy ranking
recurrentGenes = []
for gene in genes:

	#count how many unique patients affect the gene for recurrence
	uniquePatients = np.unique(genes[gene][13])

	data = [gene] + [len(uniquePatients)] + genes[gene]

	recurrentGenes.append(data)

recurrentGenes = np.array(recurrentGenes, dtype='object')

#sort
sortedGenes = recurrentGenes[np.argsort(recurrentGenes[:,1])[::-1]]

sortedGenesTop = []
for gene in sortedGenes:

	if gene[0] not in topPairGenes:
		continue
	sortedGenesTop.append(gene)

sortedGenesTop = np.array(sortedGenesTop, dtype='object')

#make a matrix in which we show visually which genes are affected in which patients
#this matrix is genes x patients
uniquePatients = dict()
top = 15 #making the matrix only for the top X genes
ind = 0
for gene in sortedGenesTop:

	if ind >= top:
		continue

	patients = gene[15]
	for patient in patients:
		if patient not in uniquePatients:
			uniquePatients[patient] = 0
		uniquePatients[patient] += 1

	ind += 1

#make a matrix of genes by patients
recurrenceMatrix = np.zeros([top, len(uniquePatients)])
ind = 0
patientOrder = dict() #order of patients in the matrix

for patientInd in range(0, len(uniquePatients)):

	patient = list(uniquePatients.keys())[patientInd]

	patientOrder[patient] = patientInd

for gene in sortedGenesTop:

	if ind >= top:
		continue

	patients = gene[15]
	for patient in patients:

		patientInd = patientOrder[patient]
		recurrenceMatrix[ind, patientInd] += 1

	ind += 1

print(recurrenceMatrix)

#make a grid plot, showing the different SV types that the patients have
#color the genes with -/+ direction, see if it correlates with the SV types.
fig, ax = plt.subplots()
for row in range(0, recurrenceMatrix.shape[0]):

	if row < recurrenceMatrix.shape[0]-1:
		ax.axhline(row+0.5, linestyle='--', color='k', linewidth=0.5)

	for col in range(0, recurrenceMatrix.shape[1]):

		if col < recurrenceMatrix.shape[1]-1:
			ax.axvline(col+0.5, linestyle='--', color='k', linewidth=0.5)

		if recurrenceMatrix[row,col] > 0:

			#get the sv type to see which symbol to assign
			gene = sortedGenesTop[row, 0]
			patient = list(uniquePatients.keys())[col]

			pairs = splitPairs[gene + '_' + patient]

			#generate some random offsets to avoid overlapping data
			offsetsX = random.sample(range(-30,30), len(pairs))
			offsetsX = [i / float(100) for i in offsetsX]

			offsetsY = random.sample(range(-30,30), len(pairs))
			offsetsY = [i / float(100) for i in offsetsY]

			ind = 0
			for pair in pairs:

				splitPair = pair[0].split('_')
				svType = splitPair[12]

				markerType = '.'
				if svType == 'DEL':
					markerType = '.'
				elif svType == 'DUP':
					markerType = 's'
				elif svType == 'INV':
					markerType = '^'
				elif svType == 'ITX':
					markerType = '*'

				#also get up/down color
				if patient + '_' + gene in degPairs[:,0]:

					#get the z-score of the pair.
					degPairInfo = degPairs[degPairs[:,0] == patient + '_' + gene][0]

					color = 'red'
					if float(degPairInfo[5]) > 1.5:
						color = 'red'
					elif float(degPairInfo[5]) < -1.5:
						color = 'blue'
					else:
						color = 'grey'
				else:
					continue #this is a pair with likely coding mutations, skip it
				plt.scatter(col + offsetsY[ind], offsetsX[ind] + (recurrenceMatrix.shape[0] - row -1), marker=markerType, edgecolor=color,
							facecolor='none', s=35)
				ind += 1
#the genes are swapped around to show the most recurrent on top, so reverse thelabels as well
plt.yticks(range(0, recurrenceMatrix.shape[0]), sortedGenesTop[0:top,0][::-1])
plt.xticks(range(0, recurrenceMatrix.shape[1]), list(uniquePatients.keys()), rotation=90)
#plt.grid()
plt.tight_layout()
plt.savefig(finalOutDir + '/recurrence_top400.svg')
plt.clf()

#Next, we are interested in patients with alternative mutations.
#So here, for each gene, first show how many patients have an SNV, CNV, or SV
#keep in mind that a duplication could be non-coding if it is in the same patient
#this will later become obvious in the visualization

#load the patient-gene mutation pairs
mutationDir = outDir + '/patientGeneMutationPairs/'
snvPatients = np.load(mutationDir + 'snvPatients.npy', allow_pickle=True, encoding='latin1').item()

svPatientsDel = np.load(mutationDir + 'svPatientsDel.npy', allow_pickle=True, encoding='latin1').item()
svPatientsDup = np.load(mutationDir + 'svPatientsDup.npy', allow_pickle=True, encoding='latin1').item()
svPatientsInv = np.load(mutationDir + 'svPatientsInv.npy', allow_pickle=True, encoding='latin1').item()
svPatientsItx = np.load(mutationDir + 'svPatientsItx.npy', allow_pickle=True, encoding='latin1').item()

cnvPatientsDel = np.load(mutationDir + 'cnvPatientsDel.npy', allow_pickle=True, encoding='latin1').item()
cnvPatientsAmp = np.load(mutationDir + 'cnvPatientsAmp.npy', allow_pickle=True, encoding='latin1').item()

#also show the non-coding SVs that do not lead to expression changes
allPairs = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_', dtype='object')

for pair in allPairs:

	splitPair = pair[0].split('_')

	gene = splitPair[0]
	patient = splitPair[7]

	sortedGeneInd = np.where(sortedGenes[:,0] == gene)[0]

	if gene in snvPatients[patient]:
		sortedGenes[sortedGeneInd, 5] += 1
	if gene in cnvPatientsDel[patient]:
		sortedGenes[sortedGeneInd, 6] += 1
	if gene in cnvPatientsAmp[patient]:
		sortedGenes[sortedGeneInd, 7] += 1
	if gene in svPatientsDel[patient]:
		sortedGenes[sortedGeneInd, 8] += 1
	if gene in svPatientsDup[patient]:
		sortedGenes[sortedGeneInd, 9] += 1
	if gene in svPatientsInv[patient]:
		sortedGenes[sortedGeneInd, 10] += 1
	if gene in svPatientsItx[patient]:
		sortedGenes[sortedGeneInd, 11] += 1

	#for the current pair, only add it if it is not in the positive set.
	if pair[0] not in positivePairs[:,0]:
		#then check the type of SV, and add it to the right gene.

		svType = splitPair[12]
		if svType == 'DEL':
			sortedGenes[sortedGeneInd, 16] += 1
		elif svType == 'DUP':
			sortedGenes[sortedGeneInd, 17] += 1
		elif svType == 'INV':
			sortedGenes[sortedGeneInd, 18] += 1
		elif svType == 'ITX':
			sortedGenes[sortedGeneInd, 19] += 1

print(sortedGenesTop[0:15,:])

#show these data in a bar plot.
#for each type of mutation, add to the stacked bar chart.
#fig,ax = plt.subplots()
geneInd = 0
ymax = 0
for gene in sortedGenes:

	if gene[0] not in sortedGenesTop[0:15,0]:
		continue
	print(gene)


	plt.bar(geneInd, gene[5], color='#ffcc00ff')
	plt.bar(geneInd, gene[6], bottom=gene[5], color='#9955ffff')
	plt.bar(geneInd, gene[7], bottom=gene[5]+gene[6], color='#ff6600b5')
	plt.bar(geneInd, gene[8], bottom=gene[5]+gene[6]+gene[7], color='#0000ffb4')
	plt.bar(geneInd, gene[9], bottom=gene[5]+gene[6]+gene[7]+gene[8], color='#d40000c6')
	plt.bar(geneInd, gene[10], bottom=gene[5]+gene[6]+gene[7]+gene[8]+gene[9], color='#ff00ccb8')
	plt.bar(geneInd, gene[11], bottom=gene[5]+gene[6]+gene[7]+gene[8]+gene[9]+gene[10], color='#808080ff')

	if gene[5]+gene[6]+gene[7]+gene[8]+gene[9]+gene[10]+gene[11] > ymax:
		ymax = gene[5]+gene[6]+gene[7]+gene[8]+gene[9]+gene[10]+gene[11] + 1

	geneInd += 1

plt.ylim(0,ymax+1)
plt.tight_layout()
plt.savefig(finalOutDir + '/recurrence_bars.svg')
plt.clf()
exit()

###Also make the full recurrence plot for all patients.
#this is quick and dirty, should have been a re-usable function.
#load the sv-gene pairs

topPairs = dict()
topPairGenes = dict() #all unique genes
svTypes = ['DEL', 'DUP', 'INV', 'ITX']

degPairs = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object')

#format: gene as key, as values, the first is the number of times we found the gene,
#the second how many were dels, then dups, invs, itx.
splitPairs = dict()
genes = dict()
for pair in positivePairs: #get stats for all pairs
	splitPair = pair[0].split('_')


	if splitPair[0] + '_' + splitPair[7] not in splitPairs:
		splitPairs[splitPair[0] + '_' + splitPair[7]] = []
	splitPairs[splitPair[0] + '_' + splitPair[7]].append(pair)

	if splitPair[0] not in genes:

		#count, cross-patient count, nc: del, dup, inv, itx, snv, cnv amp, cnv del, sv del, sv dup, sv inv, sv itx
		genes[splitPair[0]] = [0]*13
		genes[splitPair[0]].append([]) #add the patient names here
		genes[splitPair[0]].append(0) #negative dels
		genes[splitPair[0]].append(0) #negative dups
		genes[splitPair[0]].append(0) #negative invs
		genes[splitPair[0]].append(0) #negative itx

	genes[splitPair[0]][0] += 1
	
	if splitPair[12] == 'DEL':
		genes[splitPair[0]][1] += 1
	elif splitPair[12] == 'DUP':
		genes[splitPair[0]][2] += 1
	elif splitPair[12] == 'INV':
		genes[splitPair[0]][3] += 1
	elif splitPair[12] == 'ITX':
		genes[splitPair[0]][4] += 1

	patient = splitPair[7]
	genes[splitPair[0]][13].append(patient)

#convert to numpy array for easy ranking
recurrentGenes = []
for gene in genes:

	#count how many unique patients affect the gene for recurrence
	uniquePatients = np.unique(genes[gene][13])
	
	data = [gene] + [len(uniquePatients)] + genes[gene]

	recurrentGenes.append(data)

recurrentGenes = np.array(recurrentGenes, dtype='object')

#sort
sortedGenes = recurrentGenes[np.argsort(recurrentGenes[:,1])[::-1]]

sortedGenesTop = []
for gene in sortedGenes:

	#if gene[0] not in topPairGenes:
	#	continue
	sortedGenesTop.append(gene)
	
sortedGenesTop = np.array(sortedGenesTop, dtype='object')


#make a matrix in which we show visually which genes are affected in which patients
#this matrix is genes x patients
uniquePatients = dict()
top = 50 #making the matrix only for the top X genes
ind = 0
for gene in sortedGenesTop:
	
	if ind >= top:
		continue

	patients = gene[15]
	for patient in patients:
		if patient not in uniquePatients:
			uniquePatients[patient] = 0
		uniquePatients[patient] += 1

	ind += 1

#make a matrix of genes by patients
recurrenceMatrix = np.zeros([top, len(uniquePatients)])
ind = 0
patientOrder = dict() #order of patients in the matrix

for patientInd in range(0, len(uniquePatients)):
	
	patient = list(uniquePatients.keys())[patientInd]
	
	patientOrder[patient] = patientInd

for gene in sortedGenesTop:
	
	if ind >= top:
		continue
	
	patients = gene[15]
	for patient in patients:

		patientInd = patientOrder[patient]
		recurrenceMatrix[ind, patientInd] += 1
		
	ind += 1	
		
print(recurrenceMatrix)

#make a grid plot, showing the different SV types that the patients have
#color the genes with -/+ direction, see if it correlates with the SV types.
fig, ax = plt.subplots(figsize=(20,10))
for row in range(0, recurrenceMatrix.shape[0]):

	if row < recurrenceMatrix.shape[0]-1:
		ax.axhline(row+0.5, linestyle='--', color='k', linewidth=0.5)

	for col in range(0, recurrenceMatrix.shape[1]):

		if col < recurrenceMatrix.shape[1]-1:
			ax.axvline(col+0.5, linestyle='--', color='k', linewidth=0.5)
		
		if recurrenceMatrix[row,col] > 0:
			
			#get the sv type to see which symbol to assign
			gene = sortedGenesTop[row, 0]
			patient = list(uniquePatients.keys())[col]
			
			pairs = splitPairs[gene + '_' + patient]
			
			#generate some random offsets to avoid overlapping data
			offsetsX = random.sample(range(-30,30), len(pairs))
			offsetsX = [i / float(100) for i in offsetsX]
			
			offsetsY = random.sample(range(-30,30), len(pairs))
			offsetsY = [i / float(100) for i in offsetsY]
							
			ind = 0
			for pair in pairs:
				
				splitPair = pair[0].split('_')
				svType = splitPair[12]
			
				markerType = '.'
				if svType == 'DEL':
					markerType = '.'
				elif svType == 'DUP':
					markerType = 's'
				elif svType == 'INV':
					markerType = '^'
				elif svType == 'ITX':
					markerType = '*'
					
				#also get up/down color
				if patient + '_' + gene in degPairs[:,0]: 
					
					#get the z-score of the pair. 
					degPairInfo = degPairs[degPairs[:,0] == patient + '_' + gene][0]
		
					color = 'red'
					if float(degPairInfo[5]) > 1.5:
						color = 'red'
					elif float(degPairInfo[5]) < -1.5:
						color = 'blue'
					else:
						color = 'grey'
				else:
					continue #this is a pair with likely coding mutations, skip it
				plt.scatter(col + offsetsY[ind], offsetsX[ind] + (recurrenceMatrix.shape[0] - row -1), marker=markerType, edgecolor=color,
							facecolor='none', s=35)
				ind += 1
#the genes are swapped around to show the most recurrent on top, so reverse thelabels as well
plt.yticks(range(0, recurrenceMatrix.shape[0]), sortedGenesTop[0:top,0][::-1])
plt.xticks(range(0, recurrenceMatrix.shape[1]), list(uniquePatients.keys()), rotation=90)
#plt.grid()
plt.tight_layout()
plt.savefig(finalOutDirFullFigure + '/recurrence_allPatients.svg')
plt.clf()

exit()
