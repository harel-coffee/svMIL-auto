import settings
from inputParser import InputParser
import numpy as np

class FeatureMatrixDefiner:
	enhancerData = []
	tadData = []
	geneData = []
	
	"""
		Make a feature matrix for deep learning
	
	"""
	
	
	def defineFeatureMatrix(self, sortedPairs):
		
		window = 2000000 #look in a 2 MB window
		binSize = 1000
		
		#Make a bin map, where each possible position within the window is assigned a bin.
		binMap = dict()
		currentBin = 0
		
		for pos in range(0, (window*2)+1):
			
			binMap[pos] = currentBin
			if pos % binSize == 0: #new bin
				currentBin += 1
		
		
		if settings.general['genes'] == True:
			windowFeaturesGenes = []
			
			print("Making features for pairs for genes")
			currentChr = None
			count = 0
			for pair in sortedPairs:
				
				pairFeatures = np.zeros(int(window*2/binSize+1))
			
				sv = pair[1:len(pair)]
				
				#get the gene by name
				gene = self.geneData[self.geneData[:,3] == pair[0]][0]
				
				#which window side is this gene on?
				#The SV cannot overlap the gene, so it must either be starting to the right of the SV, or ending before the left of the SV. 
				if gene[1] > sv[5]:
					#collect bins
					bins = dict()
					#opt: start and end, then infer everything in the middle.
					
					featurePos1 = window + abs(gene[1] - int(sv[5]))
					featurePos2 = window + abs(gene[2] - int(sv[5]))
					
					bin1 = binMap[featurePos1]
					
					
					#get a range of possible bins.
					if featurePos2 <= window*2:
						bin2 = pairFeatures.shape[0] #get the max bin
					
					binRange = range(bin1, bin2)
					for newBin in binRange:	
						pairFeatures[newBin] = 1
				else:
					#collect bins
					bins = dict()
					#opt: start and end, then infer everything in the middle.
					
					
					featurePos1 = window - abs(gene[1] - int(sv[1]))
					featurePos2 = window - abs(gene[2] - int(sv[1]))
					
					
					bin2 = binMap[featurePos2]
					
					#get a range of possible bins.
					if featurePos1 < 0: 
						bin1 = 0 #get the max bin
					else:
						bin1 = binMap[featurePos1]
					binRange = range(bin1, bin2)
					for newBin in binRange:	
						pairFeatures[newBin] = 1
				
				
				windowFeaturesGenes.append(pairFeatures)
				
				#Also do the reverse
		
		print(np.unique(np.where(np.array(windowFeaturesGenes) == 1)[0]).shape)
		
		if settings.general['tads'] == True:
			
			windowFeaturesTads = []
			
			print("Making features for pairs for TADs")
			currentChr = None
			count = 0
			for pair in sortedPairs:
				
				pairFeatures = np.zeros(int(window*2/binSize+1))
			
				sv = pair[1:len(pair)]
				
				#only get a new chr subset if the chr changes
				if sv[0] != currentChr:
					#Get the subset of eQTls on this chromosome
					chr1Subset = self.tadData[np.where(self.tadData[:,0] == sv[0])]
					print(chr1Subset.shape)
					currentChr = sv[0]
				
				#The only eQTLs we consider are the ones that are starting or ending within the window, and are not within the SV.
				#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
				
				#Match for TADs that start in the start window
				startStartMatches = (chr1Subset[:,1] <= sv[1]) * (chr1Subset[:,1] >= (sv[1] - window))
				#Match for TADs that end in the start window		
				startEndMatches = (chr1Subset[:,2] <= sv[1]) * (chr1Subset[:,2] >= (sv[1] - window))		
						
				#The reverse for the end of the SV.
				endStartMatches = (chr1Subset[:,1] >= sv[5]) * (chr1Subset[:,1] <= (sv[5] + window))
				endEndMatches = (chr1Subset[:,2] >= sv[5]) * (chr1Subset[:,2] <= (sv[5] + window))
				
				matchingStart = chr1Subset[startStartMatches + startEndMatches] # must be either matching on the left or right.
				matchingEnd = chr1Subset[endStartMatches + endEndMatches] # must be either matching on the left or right.
				
				#For every position, determine the position relative to the SV and add it at the right place in the feature vector.
				for match in chr1Subset[startStartMatches]:
					#the position in the feature vector is 2 MB - abs(the SV start - element pos)
					featurePos = window - abs(match[1] - int(sv[1]))
					
					correctBin = binMap[featurePos]
					pairFeatures[correctBin] = 1
				for match in chr1Subset[startEndMatches]:
					#the position in the feature vector is 2 MB - abs(the SV start - element pos)
					featurePos = window - abs(match[2] - int(sv[1]))
					
					correctBin = binMap[featurePos]
					pairFeatures[correctBin] = 1	
				
				for match in chr1Subset[endStartMatches]:
					#Here, we do + 2MB
					featurePos = window + abs(match[1] - int(sv[5]))
					correctBin = binMap[featurePos]
					pairFeatures[correctBin] = 1
				for match in chr1Subset[endEndMatches]:
					#Here, we do + 2MB
					featurePos = window + abs(match[2] - int(sv[5]))
					correctBin = binMap[featurePos]
					pairFeatures[correctBin] = 1
				
				windowFeaturesTads.append(pairFeatures)
		
		#print np.unique(np.where(np.array(windowFeaturesTads) == 1)[0]).shape
		
		#Add feature layer for enhancer data
		if settings.general['enhancers'] == True:
		
			windowFeaturesEnhancers = []
			
			print("Making features for enhancers")
			currentChr = None
			count = 0
			for pair in sortedPairs:
				
				pairFeatures = np.zeros(int(window*2/binSize+1))
			
				sv = pair[1:len(pair)]
				
				#only get a new chr subset if the chr changes
				if sv[0] != currentChr:
					#Get the subset of eQTls on this chromosome
					chr1Subset = self.enhancerData[np.where(self.enhancerData[:,0] == sv[0])]
					print(chr1Subset.shape)
					currentChr = sv[0]
				
				#The only eQTLs we consider are the ones that are starting or ending within the window, and are not within the SV.
				#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
				startMatches = (chr1Subset[:,1] <= sv[1]) * (chr1Subset[:,1] >= (sv[1] - window)) * (chr1Subset[:,2] <= sv[1])
						
				#The reverse for the end of the SV.
				endMatches = (chr1Subset[:,2] >= sv[5]) * (chr1Subset[:,2] <= (sv[5] + window)) * (chr1Subset[:,1] >= sv[5])
				
				matchingStart = chr1Subset[startMatches] # must be either matching on the left or right.
				matchingEnd = chr1Subset[endMatches] # must be either matching on the left or right.
				
				#For every position, determine the position relative to the SV and add it at the right place in the feature vector.
				for match in matchingStart:
					#the position in the feature vector is 2 MB - abs(the SV start - element pos)
					featurePos = window - abs(match[1] - int(sv[1]))
					
					correctBin = binMap[featurePos]
					pairFeatures[correctBin] = 1
				for match in matchingEnd:
					#Here, we do + 2MB
					featurePos = window + abs(match[1] - int(sv[5]))
					correctBin = binMap[featurePos]
					pairFeatures[correctBin] = 1
				
				windowFeaturesEnhancers.append(pairFeatures)
				
		#for now just return
		#return np.dstack((np.array(windowFeaturesEnhancers), np.array(windowFeaturesTads), np.array(windowFeaturesGenes)))
		return np.concatenate((np.array(windowFeaturesEnhancers), np.array(windowFeaturesTads), np.array(windowFeaturesGenes)), axis=1) #stitch columns together
		
		
	
	def setFeatureData(self):
		"""
			Get the feature data from all files. Re-use the original settings file from rule-based for now
		"""
		
		geneData = []
		if settings.general['genes'] == True:
			print("Getting genes")
			causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
			nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
			
			#Combine the genes into one set. 
			self.geneData = np.concatenate((causalGenes, nonCausalGenes), axis=0)
	
		tadData = []
		if settings.general['tads'] == True:
			print("Getting TADs")
			tadFile = settings.files['tadFile']
			self.tadData = InputParser().getTADsFromFile(tadFile)
		
		eQTLData = [] #Keep empty by default in case we do not use eQTLs
		if settings.general['eQTLs'] == True: #Always check if eQTLs are enabled in the settings
			
			eQTLFile = settings.files['eQTLFile']
			print("getting eQTLs")
			eQTLData = InputParser().getEQTLsFromFile(eQTLFile)

		#3. Get enhancers
		
		if settings.general['enhancers'] == True:
			print("getting enhancers")
			self.enhancerData = InputParser().getEnhancersFromFile(settings.files['enhancerFile'])
			
		#4. Get promoters
		
		if settings.general['promoters'] == True:
			print("getting promoters")
			promoterData = InputParser().getPromotersFromFile(settings.files['promoterFile'])
			
		#5. Get CpG islands
		if settings.general['cpgIslands'] == True:
			print("Getting cpg islands")
			cpgData = InputParser().getCpgIslandsFromFile(settings.files['cpgFile'])
		
		#6. Get Transcription factors
		if settings.general['transcriptionFactors'] == True:
			print("Getting transcription factors")

			tfData = InputParser().getTranscriptionFactorsFromFile(settings.files['tfFile'])
	
		
		#7. Get Hi-C data
		if settings.general['hiC'] == True:
			print("Getting Hi-C data")
			hicData = InputParser().getHiCInteractionsFromFile(settings.files['hicFile'])
			
		#8. Get histone marks
		if settings.general['histones'] == True:
			print("Getting histone marks")
			files = [settings.files['h3k9me3'], settings.files['h3k4me3'], settings.files['h3k27ac'], settings.files['h3k27me3'],
					 settings.files['h3k4me1'], settings.files['h3k36me3']]
			types = ['h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
			for histoneFileInd in range(0, len(files)):
				histoneData = InputParser().getHistoneMarksFromFile(files[histoneFileInd], types[histoneFileInd])
			
		#9. Get DNAse I hypersensitivty sites
		if settings.general['dnaseI'] == True:
			print("Getting DNAse I hypersensitivity sites")
			
			dnaseIData = InputParser().getDNAseIFromFile(settings.files['dnaseIFile'])
			