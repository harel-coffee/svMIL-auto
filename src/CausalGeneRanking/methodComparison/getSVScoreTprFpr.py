"""

	For SVScore, read all the output VCFs and get the scores of the SVs.

"""


import sys
import re
import numpy as np
import glob

path = sys.argv[1]
sys.path.insert(1, path)

import settings

#First read the vcfs and get the scores per SV

svScoreOutDir = settings.files['svDir']

svScoreFiles = glob.glob(svScoreOutDir + '/**/*.svScore*', recursive=True)

predictions = dict()
for outFile in svScoreFiles:

	#read the file and get the consequences
	with open(outFile, 'r') as inF:

		sampleName = re.search('.*\/([A-Z\d]+)\.', outFile).group(1)
		print(sampleName)

		for line in inF:
			line = line.strip()
			if re.search('^#', line):
				continue

			splitLine = line.split("\t")
			filterInfo = splitLine[6]
			if filterInfo != 'PASS':
				continue

			chr1 = splitLine[0]
			pos1 = int(splitLine[1])
			pos2Info = splitLine[4]

			#match the end position and orientation. if there is no orientation info, this is an insertion, which we can skip.
			if not re.search(':', pos2Info):
				continue

			if re.match('[A-Z]*\[.*\:\d+\[$', pos2Info):
				o1 = '+'
				o2 = '-'
			elif re.match('[A-Z]*\].*\:\d+\]$', pos2Info):
				o1 = '-'
				o2 = '+'
			elif re.match('^\].*\:\d+\][A-Z]*', pos2Info):
				o1 = '+'
				o2 = '+'
			elif re.match('^\[.*\:\d+\[[A-Z]*', pos2Info):
				o1 = '-'
				o2 = '-'
			else:
				print('unmatched: ', pos2Info)
				print(line)
				exit()

			#get the chr2 information
			chr2 = re.search('[\[\]]+(.*):(\d+).*', pos2Info).group(1)
			pos2 = int(re.search('.*\:(\d+).*', pos2Info).group(1))

			infoField = splitLine[7]
			splitInfoField = infoField.split(";")
			svType = ''
			svScore = 0
			for field in splitInfoField:
				splitField = field.split("=")
				if splitField[0] == 'SVTYPE':
					svType = splitField[1]
				if splitField[0] == 'SVSCORETOP10':
					svScore = float(splitField[1])

			if svType not in ['DEL', 'DUP', 'INV', 'BND']:
				continue

			#default positions
			s1 = pos1
			e1 = pos1
			s2 = pos2
			e2 = pos2
			orderedChr1 = chr1
			orderedChr2 = chr2

			#switch chromosomes if necessary
			if chr1 != chr2:
				if chr1 == 'Y' and chr2 == 'X':
					orderedChr1 = chr2
					orderedChr2 = chr1
				if (chr1 == 'X' or chr1 == 'Y' or chr1 == 'MT') and (chr2 != 'X' and chr2 != 'Y' and chr2 != 'MT'):
					orderedChr1 = chr2
					orderedChr2 = chr1
				if (chr1 != 'X' and chr1 != 'Y' and chr1 != 'MT') and (chr2 != 'X' and chr2 != 'Y' and chr2 != 'MT'):
					if int(chr1) > int(chr2):
						orderedChr1 = chr2
						orderedChr2 = chr1
				if (chr1 in ['X', 'Y', 'MT']) and (chr2 in ['X', 'Y', 'MT']): #order these as well
					if chr1 == 'Y' and chr2 == 'X':
						orderedChr1 = chr2
						orderedChr2 = chr1
					if chr1 == 'MT' and chr2 in ['X', 'Y']:
						orderedChr1 = chr2
						orderedChr2 = chr1


				#always switch the coordinates as well if chromosomes are switched.
				if orderedChr1 == chr2:
					s1 = pos2
					e1 = pos2
					s2  = pos1
					e2 = pos1

			else: #if the chr are the same but the positions are reversed, change these as well.
				if pos2 < pos1:
					s1 = pos2
					e1 = pos2
					s2  = pos1
					e2 = pos1

			s1 = str(s1)
			e1 = str(e1)
			s2 = str(s2)
			e2 = str(e2)

			finalChr1 = 'chr' + orderedChr1
			finalChr2 = 'chr' + orderedChr2
			svStr = finalChr1 + '_' + s1 + '_' + e1 + '_' + finalChr2 + '_' + s2 + '_' + e2 + '_' + sampleName + '_' + svType


			predictions[svStr] = svScore

#all positive pairs to compare to
positivePairs = np.loadtxt(sys.argv[2], dtype='object')

fixedNamesPairs = dict() #change the names to match the VEP output
for pair in positivePairs:

	splitPair = pair[0].split('_')
	newName = splitPair[1] + '_' + splitPair[2] + '_' + splitPair[3] + '_' + splitPair[4] + '_' + splitPair[5] + '_' + splitPair[6] + '_' + splitPair[7] + '_' + splitPair[12]

	fixedNamesPairs[newName] = 0


#Collect all the scores per SV type, and then select the ones that are likely pathogenic.
#in the paper they use the top 10% and 90 percentile.

svTypes = ['DEL', 'DUP', 'INV', 'BND']

for svType in svTypes:
	print(svType)

	TP = 0
	FP = 0
	FN = 0
	TN = 0
	allScores = []
	for prediction in predictions:

		splitPred = prediction.split('_')
		if splitPred[7] != svType:
			continue

		allScores.append(predictions[prediction])

	allScores = np.array(allScores)
	#get the 90 percentile
	filteredScores = allScores[allScores >= np.percentile(allScores, 90)]

	print(min(filteredScores), max(filteredScores))

	#only keep the matchng pathogenic SVs
	#check for these what the TPR/FPR is.
	for prediction in predictions:
		if predictions[prediction] in filteredScores: #pathogenic predicted

			splitPred = prediction.split('_')
			if splitPred[7] != svType:
				continue
			
			if prediction in fixedNamesPairs:
				TP += 1
			elif prediction not in fixedNamesPairs:
				FP += 1
		else: #not pathogenic predicted
			if prediction in fixedNamesPairs:
				FN += 1
			elif prediction not in fixedNamesPairs:
				TN += 1

	print(TP, FP, FN, TN)

	tpr = TP / float((TP + FP))
	fpr = FP / float((TN + FP))
	print('tpr', tpr)
	print('fpr', fpr)



#Then get the true positive SV pairs, and compute TPR/FPR