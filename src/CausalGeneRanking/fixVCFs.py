"""
	Make versions of the HMF vcfs that SVScore and VEP can read

	- remove non-pass, add end, replace svtype

"""

import sys
import glob
import re

vcfs = glob.glob(sys.argv[1] + '/**/*.svTypes.passed')

for vcf in vcfs:

	sampleName = re.search('.*\/([A-Z\d]+)\.', vcf).group(1)

	if sampleName != 'CPCT02020344T':
		continue

	newVCF = []
	with open(vcf, 'r') as inF:

		for line in inF:

			if re.search('^#', line): #skip header
				newVCF.append(line)
				continue


			#skip the SV if it did not pass.
			splitLine = line.split("\t")

			if splitLine[2] == 'gridss0_18445o':
				print(line)
				exit()

			filterInfo = splitLine[6]
			if filterInfo != 'PASS':
				continue

			#fix the SV type
			infoField = splitLine[7]
			splitInfoField = infoField.split(";")
			svType = ''
			oldSVType = ''
			for field in splitInfoField:

				splitField = field.split("=")
				if splitField[0] == 'SIMPLE_TYPE':
					svType = splitField[1]

					if svType == 'ITX':
						svType = 'BND'
				if splitField[0] ==  'SVTYPE':
					oldSVType = splitField[1]

			pos1 = int(splitLine[1])
			pos2Info = splitLine[4]
			if not re.search(':', pos2Info):
				continue
			pos2 = int(re.search('.*\:(\d+).*', pos2Info).group(1))

			start = pos1
			end = pos2
			if pos1 > pos2:
				start = pos2
				end = pos1


			#replace the old SV type with this new one
			line = line.replace('SVTYPE=' + oldSVType, 'SVTYPE=' + svType)
			#replace the simple type with END
			line = line.replace('SIMPLE_TYPE=' + svType, 'END=' + str(end))

			splitLine = line.split("\t")
			splitLine[1] = str(start) #make sure that the start and end are swapped correctly.

			line = '\t'.join(splitLine)


			newVCF.append(line)

	#with open(vcf + '_fixed.vcf', 'w') as outF:
	#	for line in newVCF:
	#		outF.write(line)
		