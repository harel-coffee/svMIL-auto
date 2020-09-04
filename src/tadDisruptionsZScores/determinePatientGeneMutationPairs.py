"""
	Generate files with dictionaries where each patient has a list of genes that are affected by a certain mutation type.
	We use this to filter out genes that are not affected just by the non-coding SV, but also by other stuff.

	Split by: SNVs, CNV amplifications (> 2.3), CNV deletions ( < 1.7), (coding) SV deletions, SV duplications, SV inversions, SV translocations.

	This works with HMF data, PCAWG data and TCGA data.

"""


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

outDir = sys.argv[2]


#Get the CNVs per gene
def getPatientsWithCNVGeneBased_hmf(cnvDir):

	tsvs = glob.glob(cnvDir + '/**/*.gene.tsv', recursive=True)

	cnvPatientsDel = dict()
	cnvPatientsAmp = dict()

	for tsv in tsvs:

		#get the samplename from the vcf
		sampleName = re.search('.*\/([A-Z\d]+)\.', tsv).group(1)

		if sampleName not in cnvPatientsAmp:
			cnvPatientsAmp[sampleName] = []
		if sampleName not in cnvPatientsDel:
			cnvPatientsDel[sampleName] = []

		#open the .gz file
		with open(tsv, 'r') as inF:

			lineCount = 0
			for line in inF:

				if lineCount < 1: #skip header
					lineCount += 1
					continue

				splitLine = line.split("\t")

				gene = splitLine[3]


				if float(splitLine[5]) > 1.7 and float(splitLine[5]) < 2.3: #these are not CNVs
					continue

				if float(splitLine[5]) > 2.3:

					cnvPatientsAmp[sampleName].append(gene)
				elif float(splitLine[5]) < 1.7:

					cnvPatientsDel[sampleName].append(gene)

	return cnvPatientsAmp, cnvPatientsDel

#Get the SNVs per gene
def getPatientsWithSNVs_hmf(snvDir):

	geneNameConversionMap = dict()
	geneNameConversionFile = settings.files['geneNameConversionFile']
	with open(geneNameConversionFile, 'r') as inF:

		for line in inF:
			line = line.strip()
			splitLine = line.split("\t")
			ensgId = splitLine[3]
			splitEnsgId = ensgId.split('.') #we only keep everything before the dot

			geneName = splitLine[4]
			geneNameConversionMap[splitEnsgId[0]] = geneName

	#search through the SNVs and link these to genes.
	vcfs = glob.glob(snvDir + '/**/*.somatic.vcf.gz', recursive=True)

	patientsWithSNVs = dict()
	for vcf in vcfs:

		#get the samplename from the vcf
		sampleName = re.search('.*\/([A-Z\d]+)\.', vcf).group(1)


		#open the .gz file
		with gzip.open(vcf, 'rb') as inF:

			for line in inF:
				line = line.strip().decode('utf-8')

				if re.search('^#', line): #skip header
					continue

				#skip the SV if it did not pass.
				splitLine = line.split("\t")
				filterInfo = splitLine[6]
				if filterInfo != 'PASS':
					continue

				#Check if this SNV has any affiliation with a gene. This means that in the info field, a gene is mentioned somewhere. That is, there is an ENSG identifier.
				infoField = splitLine[7]

				geneSearch = re.search('(ENSG\d+)', infoField)
				if geneSearch:
					geneMatch = re.search('(ENSG\d+)', infoField).group(1)
					#skip genes for which we do not know the name
					if geneMatch not in geneNameConversionMap:
						continue
					geneName = geneNameConversionMap[geneMatch]

					if sampleName not in patientsWithSNVs:
						patientsWithSNVs[sampleName] = []
					patientsWithSNVs[sampleName].append(geneName)

	return patientsWithSNVs

# #Get the SV per gene
def getPatientsWithSVs_hmf(svDir, allGenes):

	#Get all parsed and annotated SV type files from the main dir
	#use all genes because there is no gene in the file, so by overlap we determine which genes are affected by the SVs. 

	vcfs = glob.glob(svDir + '/**/*.svTypes.passed', recursive=True)

	svPatientsDel = dict()
	svPatientsDup = dict()
	svPatientsInv = dict()
	svPatientsItx = dict()

	for vcf in vcfs:

		#get the samplename from the vcf
		sampleName = re.search('.*\/([A-Z\d]+)\.', vcf).group(1)
		if sampleName not in svPatientsDel:
			svPatientsDel[sampleName] = []
		if sampleName not in svPatientsDup:
			svPatientsDup[sampleName] = []
		if sampleName not in svPatientsInv:
			svPatientsInv[sampleName] = []
		if sampleName not in svPatientsItx:
			svPatientsItx[sampleName] = []


		#open the .gz file
		with open(vcf, 'r') as inF:

			for line in inF:

				if re.search('^#', line): #skip header
					continue

				#skip the SV if it did not pass.
				splitLine = line.split("\t")
				filterInfo = splitLine[6]
				if filterInfo != 'PASS':
					continue

				#Check if the SV is a deletion
				infoField = splitLine[7]
				splitInfoField = infoField.split(";")
				svType = ''
				for field in splitInfoField:

					splitField = field.split("=")
					if splitField[0] == 'SIMPLE_TYPE':
						svType = splitField[1]

				#skip non-deletions
				if svType not in ['DEL', 'DUP', 'INV', 'ITX']:
				#if svType not in ['DUP']:
					continue

				chr1 = splitLine[0]
				pos1 = int(splitLine[1])
				pos2Info = splitLine[4]
				pos2 = int(re.search('.*\:(\d+).*', pos2Info).group(1))
				chr2 = re.search('[\[\]]+(.*):(\d+).*', pos2Info).group(1)

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

				chr1 = 'chr' + orderedChr1
				chr2 = 'chr' + orderedChr2

				#Check which genes are overlapped by this SV.
				#Keep track of the disrupted genes in the patient.

				#intrachromosomal SV
				if chr1 == chr2:

					geneChrSubset = allGenes[allGenes[:,0] == chr1]

					geneMatches = geneChrSubset[(geneChrSubset[:,1] <= e2) * (geneChrSubset[:,2] >= s1)]

					if svType == 'DEL':
						for match in geneMatches:
							svPatientsDel[sampleName].append(match[3].name)


					elif svType == 'DUP':
						for match in geneMatches:
							svPatientsDup[sampleName].append(match[3].name)
					elif svType == 'INV':
						for match in geneMatches:
							svPatientsInv[sampleName].append(match[3].name)

				else:

					#find breakpoints in the gene for each side of the SV
					geneChr1Subset = allGenes[allGenes[:,0] == chr1]
					geneChr2Subset = allGenes[allGenes[:,0] == chr2]

					#check if the bp start is within the gene.
					geneChr1Matches = geneChr1Subset[(s1 >= geneChr1Subset[:,1]) * (s1 <= geneChr1Subset[:,2])]
					geneChr2Matches = geneChr2Subset[(s2 >= geneChr2Subset[:,1]) * (s2 <= geneChr2Subset[:,2])]

					for match in geneChr1Matches:
						svPatientsItx[sampleName].append(match[3].name)



					for match in geneChr2Matches:
						svPatientsItx[sampleName].append(match[3].name)



	return svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx

def getPatientsWithSVs_tcga(svFile, allGenes):
	#use all genes because there is no gene in the file, so by overlap we determine which genes are affected by the SVs. 
	
	#load the svs
	svData = InputParser().getSVsFromFile(svFile, '')
	
	svPatientsDel = dict()
	svPatientsDup = dict()
	svPatientsInv = dict()
	svPatientsItx = dict()

	for sv in svData:
		
		sampleName = sv[7]
		if sampleName not in svPatientsDel:
			svPatientsDel[sampleName] = []
		if sampleName not in svPatientsDup:
			svPatientsDup[sampleName] = []
		if sampleName not in svPatientsInv:
			svPatientsInv[sampleName] = []
		if sampleName not in svPatientsItx:
			svPatientsItx[sampleName] = []

		chr1 = sv[0]
		s1 = sv[1]
		e1 = sv[2]
		chr2 = sv[3]
		s2 = sv[4]
		e2 = sv[5]
		svType = sv[8].svType
		#intrachromosomal SV
		if chr1 == chr2:

			geneChrSubset = allGenes[allGenes[:,0] == chr1]

			geneMatches = geneChrSubset[(geneChrSubset[:,1] <= e2) * (geneChrSubset[:,2] >= s1)]

			if svType == 'DEL':
				for match in geneMatches:
					svPatientsDel[sampleName].append(match[3].name)


			elif svType == 'DUP':
				for match in geneMatches:
					svPatientsDup[sampleName].append(match[3].name)
			elif svType == 'INV':
				for match in geneMatches:
					svPatientsInv[sampleName].append(match[3].name)

		else:

			#find breakpoints in the gene for each side of the SV
			geneChr1Subset = allGenes[allGenes[:,0] == chr1]
			geneChr2Subset = allGenes[allGenes[:,0] == chr2]

			#check if the bp start is within the gene.
			geneChr1Matches = geneChr1Subset[(s1 >= geneChr1Subset[:,1]) * (s1 <= geneChr1Subset[:,2])]
			geneChr2Matches = geneChr2Subset[(s2 >= geneChr2Subset[:,1]) * (s2 <= geneChr2Subset[:,2])]

			for match in geneChr1Matches:
				svPatientsItx[sampleName].append(match[3].name)

			for match in geneChr2Matches:
				svPatientsItx[sampleName].append(match[3].name)

	return svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx

def getPatientsWithSNVs_tcga(snvDir):

	allFiles = [f for f in listdir(snvDir) if isfile(join(snvDir, f))]

	snvPatients = dict()

	for currentFile in allFiles:

		if currentFile == "MANIFEST.txt":
			continue
		splitFileName = currentFile.split(".")
		patientID = splitFileName[0]
		splitPatientID = patientID.split("-")
		shortPatientID = settings.general['cancerType'] + splitPatientID[2]
		
		if shortPatientID not in snvPatients:
			snvPatients[shortPatientID] = []

		#Load the contents of the file
		with open(snvDir + "/" + currentFile, 'r') as inF:
			lineCount = 0
			for line in inF:
				line = line.strip() #remove newlines
				if lineCount < 1: #only read the line if it is not a header line
					lineCount += 1
					continue

				splitLine = line.split("\t")
				geneName = splitLine[0]

				if splitLine[8] == 'Silent': #skip the ones with no effect
					continue

				snvPatients[shortPatientID].append(geneName)

	return snvPatients

def getMetadataPCAWG(metadataFile):
	#get the metadata file to extract the mapping from wgs to rna-seq identifiers

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

			if header['matched_wgs_aliquot_id'] < len(splitLine):
				wgsName = splitLine[header['matched_wgs_aliquot_id']]
				rnaName = splitLine[header['aliquot_id']]

				nameMap[wgsName] = rnaName

	return nameMap

def getPatientsWithSNVs_pcawg(snvDir, allGenes):
	#use all genes because there is no gene in the file, so by overlap we determine which genes are affected by the SNVs.
	
	#get the cancer type from the settings
	cancerType = settings.general['cancerType']
	#read in the metadata file and get the right file identifiers
	metadataFile = settings.files['pcawgMetadata']
	
	#save the IDs of the patients with this cancer type
	cancerTypeIds = dict()
	with open(metadataFile, 'r') as inF:
		
		for line in inF:
			
			if re.search(cancerType, line):
				splitLine = line.split('\t')
				sampleId = splitLine[1]
				cancerTypeIds[sampleId] = 0
				
	
	#Then read the SNV files for the right cancer type
	patientsWithSNVs = dict()
	count = 0
	for sampleId in cancerTypeIds:
		
		#use glob to find the right file
		matchedFolder = glob.glob(snvDir + '/' + sampleId + '*.somatic.snv_mnv')
		
		#if we don't have SVs for this sample, skip it. 
		if len(matchedFolder) < 1:
			print(sampleId)
			continue
		
		#there should be just 1 file
		sampleSNVFolder = matchedFolder[0]
		sampleSNVFile = glob.glob(sampleSNVFolder + '/*PASS.vcf.gz')[0]
	
		
		with gzip.open(sampleSNVFile, 'rb') as inF:

			for line in inF:
				line = line.strip().decode('utf-8')

				if re.search('^#', line): #skip header
					continue
				
				
				splitLine = line.split("\t")

				#Check if this SNV has any affiliation with a gene. This means that in the info field, a gene is mentioned somewhere. That is, there is an ENSG identifier.
				infoField = splitLine[7]

				#check if the mutation is silent
				splitInfoField = infoField.split(';')
				skip = False
				for field in splitInfoField:

					splitField = field.split('=')
					if splitField[0] == 'Variant_Classification' and splitField[1] == 'Silent':
						skip = True

				if skip == True:
					continue

				#link it to the gene it is in
				geneChrSubset = allGenes[allGenes[:,0] == 'chr' + splitLine[0]]

				geneMatch = geneChrSubset[(geneChrSubset[:,1] <= int(splitLine[1])) * (geneChrSubset[:,2] >= int(splitLine[1]))]

				if len(geneMatch) < 1:
					continue

				geneName = geneMatch[0][3].name

				if sampleId not in patientsWithSNVs:
					patientsWithSNVs[sampleId] = []
				patientsWithSNVs[sampleId].append(geneName)

	return patientsWithSNVs

def getPatientsWithCNVGeneBased_pcawg(cnvDir, allGenes):
	#use nameMap to map the WGS identifiers to the rna-seq identifiers

	#get the cancer type from the settings
	cancerType = settings.general['cancerType']
	#read in the metadata file and get the right file identifiers
	metadataFile = settings.files['pcawgMetadata']

	#save the IDs of the patients with this cancer type
	cancerTypeIds = dict()
	cnvPatientsDel = dict()
	cnvPatientsAmp = dict()
	with open(metadataFile, 'r') as inF:

		for line in inF:

			if re.search(cancerType, line):
				splitLine = line.split('\t')
				sampleId = splitLine[1]
				cancerTypeIds[sampleId] = 0
				cnvPatientsAmp[sampleId] = []
				cnvPatientsDel[sampleId] = []


	#Then read the SNV files for the right cancer type
	count = 0
	for sampleId in cancerTypeIds:

		#use glob to find the right file
		matchedFile = glob.glob(cnvDir + '/' + sampleId + '*copyNumberEstimation*.somatic.cnv.vcf.gz')

		#if we don't have SVs for this sample, skip it.
		if len(matchedFile) < 1:
			print(sampleId)
			continue

		#there should be just 1 file
		cnvFile = matchedFile[0]
		print(cnvFile)

		with gzip.open(cnvFile, 'rb') as inF:

			for line in inF:
				line = line.strip().decode('utf-8')

				if re.search('^#', line): #skip header
					continue


				splitLine = line.split("\t")

				#check which gene is affected by this CNV
				formatField = splitLine[9]

				splitFormatField = formatField.split(':')
				cn = int(splitFormatField[0]) #get the total copy number (TCN)
				if cn == 2:
					continue

				#link it to the gene it is in
				geneChrSubset = allGenes[allGenes[:,0] == 'chr' + splitLine[0]]
				cnvStart = int(splitLine[1])
				#get the end positiion from the vcf
				cnvEnd = 0
				splitInfoField = splitLine[7].split(';')
				for field in splitInfoField:
					splitField = field.split('=')
					if splitField[0] == 'END':
						cnvEnd = int(splitField[1])

				geneMatch = geneChrSubset[(geneChrSubset[:,1] <= cnvEnd) * (geneChrSubset[:,2] >= cnvStart)]

				if len(geneMatch) < 1:
					continue
				
				for gene in geneMatch:
					geneName = gene[3].name
					
					if cn > 2:
						cnvPatientsAmp[sampleId].append(geneName)
					elif cn < 2:
						cnvPatientsDel[sampleId].append(geneName)

	return cnvPatientsAmp, cnvPatientsDel
		
def getPatientsWithSVs_pcawg(svDir, allGenes):

	#Get all parsed and annotated SV type files from the main dir
	#use all genes because there is no gene in the file, so by overlap we determine which genes are affected by the SVs.

	#get the cancer type from the settings
	cancerType = settings.general['cancerType']
	#read in the metadata file and get the right file identifiers
	metadataFile = settings.files['pcawgMetadata']

	#save the IDs of the patients with this cancer type
	cancerTypeIds = dict()
	svPatientsDel = dict()
	svPatientsDup = dict()
	svPatientsInv = dict()
	svPatientsItx = dict()
	with open(metadataFile, 'r') as inF:

		for line in inF:

			if re.search(cancerType, line):
				splitLine = line.split('\t')
				sampleId = splitLine[1]
				cancerTypeIds[sampleId] = 0


				svPatientsDel[sampleId] = []
				svPatientsDup[sampleId] = []
				svPatientsInv[sampleId] = []
				svPatientsItx[sampleId] = []


	#Then read the SV files for the right cancer type
	allSVs = []
	count = 0
	for sampleId in cancerTypeIds:

		#use glob to find the right file
		matchedFiles = glob.glob(svDir + '/' + sampleId + '*sv.vcf.gz')

		#if we don't have SVs for this sample, skip it.
		if len(matchedFiles) < 1:
			print(sampleId)
			continue


		#there should be just 1 file
		svFile = matchedFiles[0]

		with gzip.open(svFile, 'rb') as inF:

			for line in inF:
				line = line.strip().decode('utf-8')

				if re.search('^#', line): #skip header
					continue

				#skip the SV if it did not pass.
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
				for field in splitInfoField:

					splitField = field.split("=")
					if splitField[0] == 'SVCLASS':

						if re.search('DEL', splitField[1]):
							svType = 'DEL'
						elif re.search('DUP', splitField[1]):
							svType = 'DUP'
						elif re.search('INV', splitField[1]):
							svType = 'INV'
						elif re.search('TRA', splitField[1]):
							svType = 'ITX'
						else:
							print('unknown: ', field)
							exit()

				#skip SV types that we do not consider
				if svType not in ['DEL', 'DUP', 'ITX', 'INV']:
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

				finalChr1 = 'chr' + orderedChr1
				finalChr2 = 'chr' + orderedChr2
				
				#Check which genes are overlapped by this SV.
				#Keep track of the disrupted genes in the patient.

				#intrachromosomal SV
				if chr1 == chr2:

					geneChrSubset = allGenes[allGenes[:,0] == finalChr1]

					geneMatches = geneChrSubset[(geneChrSubset[:,1] <= e2) * (geneChrSubset[:,2] >= s1)]

					
					if svType == 'DEL':
						for match in geneMatches:
							svPatientsDel[sampleId].append(match[3].name)

					elif svType == 'DUP':
						for match in geneMatches:
							svPatientsDup[sampleId].append(match[3].name)
					elif svType == 'INV':
						for match in geneMatches:
							svPatientsInv[sampleId].append(match[3].name)

				else:

					#find breakpoints in the gene for each side of the SV
					geneChr1Subset = allGenes[allGenes[:,0] == finalChr1]
					geneChr2Subset = allGenes[allGenes[:,0] == finalChr2]

					#check if the bp start is within the gene.
					geneChr1Matches = geneChr1Subset[(s1 >= geneChr1Subset[:,1]) * (s1 <= geneChr1Subset[:,2])]
					geneChr2Matches = geneChr2Subset[(s2 >= geneChr2Subset[:,1]) * (s2 <= geneChr2Subset[:,2])]

					for match in geneChr1Matches:
						svPatientsItx[sampleId].append(match[3].name)



					for match in geneChr2Matches:
						svPatientsItx[sampleId].append(match[3].name)

	return svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx
			

causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set.
allGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#get the right data based on the data input source.
if settings.general['source'] == 'HMF':
	cnvPatientsAmp, cnvPatientsDel = getPatientsWithCNVGeneBased_hmf(settings.files['cnvDir'])
	snvPatients = getPatientsWithSNVs_hmf(settings.files['snvDir'])
	svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx = getPatientsWithSVs_hmf(settings.files['svDir'], allGenes)
elif settings.general['source'] == 'TCGA':
	svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx = getPatientsWithSVs_tcga(settings.files['svFile'], allGenes)
	cnvPatientsAmp = dict() #not available for TCGA. Only SNP6 arrays but there is no copy number or purity, so hard to interpret.
	cnvPatientsDel = dict()
	snvPatients = getPatientsWithSNVs_tcga(settings.files['snvDir'])
elif settings.general['source'] == 'PCAWG':
	#nameMap = getMetadataPCAWG(settings.files['metaDataFile'])
	#snvPatients = getPatientsWithSNVs_pcawg(settings.files['snvDir'], allGenes)
	#cnvPatientsAmp, cnvPatientsDel = getPatientsWithCNVGeneBased_pcawg(settings.files['cnvDir'], allGenes)
	svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx = getPatientsWithSVs_pcawg(settings.files['svDir'], allGenes)


finalOutDir = outDir + '/patientGeneMutationPairs/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

np.save(finalOutDir + 'svPatientsDel.npy', svPatientsDel)
np.save(finalOutDir + 'svPatientsDup.npy', svPatientsDup)
np.save(finalOutDir + 'svPatientsInv.npy', svPatientsInv)
np.save(finalOutDir + 'svPatientsItx.npy', svPatientsItx)
np.save(finalOutDir + 'cnvPatientsDel.npy', cnvPatientsDel)
np.save(finalOutDir + 'cnvPatientsAmp.npy', cnvPatientsAmp)
np.save(finalOutDir + 'snvPatients.npy', snvPatients)