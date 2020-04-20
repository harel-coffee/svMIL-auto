"""

	Get the SVs from PCAWG for one cancer type in a single file that we can easily use.

"""

import sys
path = sys.argv[1]
sys.path.insert(1, path)

import settings
import glob
import gzip

inFolder = sys.argv[2]
metadataFile = sys.argv[3]
outFile = sys.argv[4]

#go through the metadata and map the WGS names to the RNA-seq names
nameMap = dict()
with open(metadataFile, 'r') as inF:

	header = dict()
	lineCount = 0
	for line in inF:
		line = line.strip()
		splitLine = line.split('\t')
		if lineCount < 1:

			for colInd in range(0, len(splitLine)):
				header[splitLine[colInd]] = colInd

			lineCount += 1
			continue

		#only get the ones with the right study that we selected in the settings
		if settings.general['cancerType'] == 'OV':
			if splitLine[header['study']] != 'Ovarian Cancer - AU':
				continue
		if settings.general['cancerType'] == 'LIVER':
			if splitLine[header['study']] != 'Liver Cancer - RIKEN, JP':
				continue


		if header['matched_wgs_aliquot_id'] < len(splitLine):
			wgsName = splitLine[header['matched_wgs_aliquot_id']]
			rnaName = splitLine[header['aliquot_id']]

			nameMap[wgsName] = rnaName


#Go through the files in the SV folder, and get the SVs for the right patients
allSVFiles = glob.glob(inFolder + '/*bedpe.gz')
with open(outFile, 'w') as outF:
	outF.write("chr1\ts1\te1\to1\tchr2\ts2\te2\to2\tsource\tsample_name\tsv_type\tcancer_type\n")
	for svFile in allSVFiles:

		splitFileName = svFile.split('/')
		suffix = splitFileName[6]
		rnaSampleName = suffix.split('.')[0]

		if rnaSampleName not in nameMap:
			continue

		#if this is the right sample
		#read this file and parse it into the same format we used for TCGA data


		with gzip.open(svFile, 'rb') as inF:
			lineCount = 0
			for line in inF:

				if lineCount < 1:
					lineCount += 1
					continue

				line = line.strip().decode('utf-8')
				splitLine = line.split('\t')

				#change to rna id
				sampleName = nameMap[rnaSampleName]

				#Parse to: chr1	s1	e1	o1	chr2	s2	e2	o2	source	sample_name	sv_type	cancer_type

				#depending on the cancer type, parse to a different sample name format
				#to match with the expression data file

				chr1 = splitLine[0]
				s1 = splitLine[1]
				e1 = splitLine[2]
				o1 = splitLine[8]

				chr2 = splitLine[3]
				s2 = splitLine[4]
				e2 = splitLine[5]
				o2 = splitLine[9]

				svType = splitLine[10]

				outF.write(chr1 + "\t" + s1 + "\t" + e1 + "\t" + o1 + "\t" + chr2 + "\t" + s2 + "\t" + e2 + "\t" + o2 + "\t" + "PCAWG" + "\t" + sampleName + "\t" + svType + "\t" + settings.general['cancerType'] + "\n")



outF.close()




