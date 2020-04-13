"""
	Check the output of VEP for each SV, and assign a score.

"""

import sys
import numpy as np
import glob
import re

vepOutDir = sys.argv[1]

#get all files with VEP output

vepFiles = glob.glob(vepOutDir + '*.svTypes.passed.txt', recursive=True)

predictions = dict()
for outFile in vepFiles:

	print(outFile)

	#read the file and get the consequences
	with open(outFile, 'r') as inF:

		for line in inF:
			line = line.strip()
			if re.search('^#', line):
				continue

			splitLine = line.split('\t')
			#get the field with the impact
			infoField = splitLine[13]

			splitInfoField = infoField.split(';')
			for field in splitInfoField:

				splitField = field.split('=')

				if splitField[0] == 'IMPACT':
					if splitField[1] not in predictions:
						predictions[splitField[1]] = 0
					predictions[splitField[1]] += 1

print(predictions)






