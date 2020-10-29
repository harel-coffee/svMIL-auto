#read all hmf data and merge it into 1 matrix

import sys
import numpy as np
import os
from os import listdir
from os.path import isfile, join
import glob
import re
import gzip
path = sys.argv[1]
sys.path.insert(1, path)
sys.path.insert(1, 'linkSVsGenes/')

from inputParser import InputParser
import settings

metadataFile = settings.files['metadataHMF']

#save the IDs of the patients with this cancer type
cancerTypeIds = dict()
with open(metadataFile, 'rb') as inF:

	for line in inF:
		line = line.decode('ISO-8859-1')

		splitLine = line.split('\t')

		sampleId = splitLine[1]
		patientId = splitLine[0]

		#here already skip patients for which there is no SV/expression data.
		#we don't need to process these patients
		matchedFiles = glob.glob(settings.files['svDir'] + '/*_' + patientId + '/' + sampleId + '.purple.sv.ann.vcf.gz')

		#if we don't have SVs for this sample, skip it.
		if len(matchedFiles) < 1:
			continue

		#### check if we have expression for this sample.
		expressionDir = settings.files['expressionDir']
		matchedExpressionFiles = glob.glob(expressionDir + sampleId)

		if len(matchedExpressionFiles) < 1:
			continue

		cancerTypeIds[sampleId] = patientId

print('processing for ', len(cancerTypeIds), ' patients')

#then get the expression data for these.
expressionDir = settings.files['expressionDir']

#check 1 file to get the number of genes in total.
sampleId = list(cancerTypeIds.keys())[0]
data = np.loadtxt(expressionDir + sampleId + '/' + sampleId + '.isf.gene_data.csv', dtype='object', skiprows=1, delimiter=',')
genes = data[:,1]

expressionMatrix = np.empty([len(genes)+1, len(cancerTypeIds)+1], dtype='object')
expressionMatrix[0,:] = [''] + list(cancerTypeIds.keys())
sampleInd = -1
for sampleId in cancerTypeIds:
	sampleInd += 1

	if sampleInd % 50 == 0:
		print('sample: ', sampleInd)

	#read the file
	data = np.loadtxt(expressionDir + sampleId + '/' + sampleId + '.isf.gene_data.csv', skiprows=1, usecols=range(7,8), delimiter=',')

	#extract the expression
	#offset by 1 to account for column and row names
	expressionMatrix[1:,sampleInd+1] = data

geneNames = np.array([''] + list(genes), dtype='object')

expressionMatrix[:,0] = geneNames

np.savetxt('../data/expression/HMF_merged.txt', expressionMatrix, delimiter='\t', fmt='%s')