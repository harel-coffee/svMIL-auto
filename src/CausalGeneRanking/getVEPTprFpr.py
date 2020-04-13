"""
	Check the output of VEP for each SV, and assign a score.

"""

import sys
import numpy as np
import glob
import re

vepOutDir = sys.argv[1]

#get all files with VEP output

vepFiles = glob.glob(vepOutDir + '*.svTypes.passed_fixed.vcf')

predictions = dict()
bndImpact= dict()
for outFile in vepFiles:

	print(outFile)

	#read the file and get the consequences
	with open(outFile, 'r') as inF:

		sampleName = re.search('.*\/([A-Z\d]+)\.', outFile).group(1)
		print(sampleName)

		for line in inF:
			line = line.strip()
			if re.search('^#', line):
				continue

			splitLine = line.split('\t')

			svType = splitLine[2]
			if svType == 'deletion':
				svType = 'DEL'
			elif svType == 'inversion':
				svType = 'INV'
			elif svType == 'duplication':
				svType = 'DUP'
			else:
				#for BND there are no impact predictions other than modifier, so skip this
				#vep doesn't work with translocations.
				continue
			
			posInfo = splitLine[1]
			splitPosInfo = posInfo.split(':')
			chr1 = splitPosInfo[0]
			positions = splitPosInfo[1]

			pos1 = int(positions.split('-')[0])
			pos2 = positions.split('-')[1]

			#default positions
			s1 = str(pos1-1)
			e1 = str(pos1-1)
			s2 = pos2
			e2 = pos2

			finalChr1 = 'chr' + chr1
			finalChr2 = 'chr' + chr1

			#use a name for the SV so that we can easily match
			svStr = finalChr1 + '_' + s1 + '_' + e1 + '_' + finalChr2 + '_' + s2 + '_' + e2 + '_' + sampleName + '_' + svType

			#get the field with the impact
			infoField = splitLine[13]

			splitInfoField = infoField.split(';')
			for field in splitInfoField:

				splitField = field.split('=')

				if splitField[0] == 'IMPACT':
					predictions[svStr] = splitField[1]

print(predictions)

#load the true labels. How many are correct?
positivePairs = np.loadtxt(sys.argv[2], dtype='object')

fixedNamesPairs = dict() #change the names to match the VEP output
for pair in positivePairs:

	splitPair = pair[0].split('_')
	newName = splitPair[1] + '_' + splitPair[2] + '_' + splitPair[3] + '_' + splitPair[4] + '_' + splitPair[5] + '_' + splitPair[6] + '_' + splitPair[7] + '_' + splitPair[12]

	fixedNamesPairs[newName] = 0

#assign scores.
svTypes = ['DEL', 'DUP', 'INV']
for svType in svTypes:
	print(svType)

	TP = 0
	FP = 0
	FN = 0
	TN = 0
	for prediction in predictions:

		splitPred = prediction.split('_')
		if splitPred[7] != svType:
			continue

		if predictions[prediction] in ['HIGH', 'MODERATE'] and prediction in fixedNamesPairs:
			TP += 1
		elif predictions[prediction] in ['HIGH', 'MODERATE'] and prediction not in fixedNamesPairs:
			FP += 1
		elif predictions[prediction] in ['MODIFIER'] and prediction in fixedNamesPairs:
			FN += 1
		elif predictions[prediction] in ['MODIFIER'] and prediction not in fixedNamesPairs:
			TN += 1

	print(TP, FP, FN, TN)

	tpr = TP / float((TP + FP))
	fpr = FP / float((TN + FP))
	print('tpr', tpr)
	print('fpr', fpr)





