import sys
import numpy as np
import os

outDir = sys.argv[2]


#Get a list of all genes * patients that have mutations.

#Get the CNVs per gene
def getPatientsWithCNVGeneBased(cnvDir):

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

cnvPatientsAmp, cnvPatientsDel = getPatientsWithCNVGeneBased(sys.argv[1])

#Get the SNVs per gene

def getPatientsWithSNVs(snvDir):
	import gzip
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

snvPatients = getPatientsWithSNVs(sys.argv[1])

# #Get the SV per gene
def getPatientsWithSVs(svDir, allGenes):

	#Get all parsed and annotated SV type files from the main dir

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

svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx = getPatientsWithSVs(sys.argv[1], allGenes)

finalOutDir = outDir + '/patientGeneMutationPairs/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

np.save(finalOutDir + 'svPatientsDel.npy', svPatientsDel)
np.save(finalOutDir + 'svPatientsDup.npy', svPatientsDup)
np.save(finalOutDir + 'svPatientsInv.npy', svPatientsInv)
np.save(finalOutDir + 'svPatientsItx.npy', svPatientsItx)
np.save(finalOutDir + 'cnvPatientsDel.npy', cnvPatientsDel)
np.save(finaOutDir + 'cnvPatientsAmp.npy', cnvPatientsAmp)
np.save(finalOutDir + 'snvPatients.npy', snvPatients)