"""
	The main goal of this class is to take the genes (with neighborhood and dertivative neighborhoods mapped) and visualize the distributions (later to be used as channels) of the regulatory data within
	the nearest TAD boundaries. 

"""
import matplotlib.pyplot as plt
import numpy as np
import math

import settings

class ChannelVisualizer:
	
	
	
	def __init__(self, genes, mode):
		
		#Load the list of known cancer genes
		#Then determine if the known cancer genes have losses/gains
		#how do the counts compare to the other random genes?
		#for the genes that have losses/gains in the causal gene set, visualize the channels
		#Does that look different from a random gene?
		
		
		knownBcGenes = self.loadBCGenes("../../data/Genes/breastCancerCausalGenes.txt")
		causalGenes = self.loadCausalGenes('../../data/Genes/Census_allTue Apr 10 14_56_44 2018.tsv')
		
		#self.reportEQTLChanges(genes, knownBcGenes)
		self.reportEQTLChanges(genes, causalGenes)
		#self.reportSVOverlap(genes, causalGenes)
		#self.reportSNVOverlap(genes, causalGenes)
		#self.visualizeReference(genes)
		self.visualizeDerivative(genes)
	
	def loadCausalGenes(self, causalGeneFile):
		
		genes = []
		
		with open(causalGeneFile, 'r') as f:
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				genes.append(splitLine[0])
				
			
		
		return genes
		
	
	def loadBCGenes(self, causalGeneFile):

		genes = []
		
		with open(causalGeneFile, 'r') as f:
			
			for line in f:
				
				line = line.strip()
				
				genes.append(line)

		return genes
	
	def reportSVOverlap(self, genes, knownBcGenes):
		
		bcGenesWithOverlap = 0
		otherGenesWithOverlap = 0
		bcGenesWithSNVSVOverlap = 0
		otherGenesWithSNVSVOverlap = 0
		for gene in genes:
			
			if gene.SVs is not None:
				if len(gene.SVs) > 0:	
					if gene.name in knownBcGenes:
						
						bcGenesWithOverlap += 1
					else:
						otherGenesWithOverlap += 1	
			
			if (gene.SVs is not None and len(gene.SVs) > 0) or (gene.SNVs is not None and len(gene.SNVs) > 0): 
					if gene.name in knownBcGenes:
						bcGenesWithSNVSVOverlap += 1
					else:
						otherGenesWithSNVSVOverlap += 1
			
		print bcGenesWithOverlap
		print otherGenesWithOverlap
		
		print "bc with both overlap: ", bcGenesWithSNVSVOverlap
		print "other with both overlap: ", otherGenesWithSNVSVOverlap
				
		
		return

	def reportSNVOverlap(self, genes, knownBcGenes):
		
		bcGenesWithOverlap = 0
		otherGenesWithOverlap = 0
		for gene in genes:
			
			if gene.SNVs is not None:
				if len(gene.SNVs) > 0:	
					if gene.name in knownBcGenes:
						
						bcGenesWithOverlap += 1
					else:
						otherGenesWithOverlap += 1	
			
		print bcGenesWithOverlap
		print otherGenesWithOverlap
				
		exit()
		return

	
	
	def reportEQTLChanges(self, genes, knownBcGenes):
		
		#1. Determine the number of lost and gained eQTLs for the known BC genes compared to all other genes
		
		bcLossCount = 0
		bcGainCount = 0
		otherLossCount = 0
		otherGainCount = 0
		for gene in genes:
			# 
			# if gene.name == "APOBEC3B":
			# 	print gene.lostEQTLs
			# 	print gene.gainedEQTLs
			# 	exit()
			
			#The losses and gains are per sample. Do we average across the samples?
			if len(gene.lostEQTLs.keys()) < 1 and len(gene.gainedEQTLs.keys()) > 0:
				samples = gene.gainedEQTLs.keys() 
			elif len(gene.gainedEQTLs.keys()) < 1 and len(gene.lostEQTLs.keys()) > 0:
				samples = gene.lostEQTLs.keys()
			else: #combine the samples
				samples = gene.lostEQTLs.keys() + gene.gainedEQTLs.keys()
			
			#First aggregate across the samples? Count how often we see a gain or loss in total across all samples.
			
			totalLosses = []
			totalGains = []
			for sample in samples:
				
				if sample in gene.lostEQTLs:
					losses = len(gene.lostEQTLs[sample])
					totalLosses.append(losses)
				
				if sample in gene.gainedEQTLs:
					gains = len(gene.gainedEQTLs[sample])
					totalGains.append(gains)
			
			if len(totalLosses) > 0:
				if gene.name in knownBcGenes:
					print "gene ", gene.name, " loses eQTLs"
					bcLossCount += 1 
				else:
					otherLossCount += 1
				
			if len(totalGains) > 0:
				if gene.name in knownBcGenes:
					print "gene ", gene.name, " gains eQTLs"
					bcGainCount += 1
				else:
					otherGainCount += 1
				
		
		#Show the distribution
		
		print "bc loss: ", bcLossCount
		print "bc gain: ", bcGainCount
		print "other loss: ", otherLossCount
		print "other gain: ", otherGainCount
		
		
		print len(genes)
		print len(knownBcGenes)
		
		# 
		
		return
		
		
		
	def visualizeReference(self, genes):
		"""
			For every gene, if there are lost or gained eQTLs for the derivative, show where all the eQTLs normally are within the nearest TAD boundaries
		"""
		distances = []
		for gene in genes:
			
			if len(gene.gainedEQTLs) < 1 and len(gene.lostEQTLs) < 1:
				
				continue
			
			#Check if the total distance is something that we can easily plot or if we need to make bins
			
			if gene.leftTAD == None or gene.rightTAD == None: #Also skip if there is no TAD on the left/right for now
				continue
			
			distance = abs(gene.leftTAD.end - gene.rightTAD.start) #left comes before right
			#at each eQTL position, add a +1.
			tadRange = range(gene.leftTAD.end, gene.rightTAD.start)
			plotRange = range(0, distance)
			
			eQTLPositions = np.zeros([distance,1])
			for eQTL in gene.eQTLs:
				eQTLPos = eQTL.start
				if eQTLPos < gene.leftTAD.end or eQTLPos > gene.rightTAD.start: #skip the eQTL if it is outside of the TAD boundary
					continue
				
				#The position in the array depends on the positions within the actual tad
				tadInd = tadRange.index(eQTLPos)	
				eQTLPositions[tadInd] = 1
			
			plt.plot(plotRange, eQTLPositions)
			plt.ylim(-1,2)
			plt.show()
				
			return	
			
			#distances.append(distance)
		
		#plt.boxplot(distances)
		#plt.show()
	
		
		
		1+1
		
	def visualizeDerivative(self, genes):
		"""
			For each of the genes, determine if there are lost or gained eQTLs (these should already have been limited to within the nearest TAD boundaries)
		"""
		distances = []
		gainedEQTLs = []
		lostEQTLs = []
		for gene in genes:
			#print gene.name
			if gene.name != "CNTRL":
				continue
			else:
				print gene.name
			
			print gene.gainedEQTLs
			print gene.lostEQTLs
			
			if len(gene.gainedEQTLs) < 1 and len(gene.lostEQTLs) < 1:
				
				continue
			
			# if len(gene.lostEQTLs) > 0:
			# 	print gene.name
			# 	continue

			#Check if the total distance is something that we can easily plot or if we need to make bins
			
			if gene.leftTAD == None or gene.rightTAD == None: #Also skip if there is no TAD on the left/right for now
				continue
			
			#If the gene is within a TAD, the lost eQTLs are between the left tad end, and right tad end.
			#If the gene is not within TADs, the lost eQTLs are between the end of the left tad, and start of the right tad. 
			
			distance = abs(gene.leftTAD.end - gene.rightTAD.start) #left comes before right
			#at each eQTL position, add a +1.
			
			if gene.rightTAD.end != gene.leftTAD.end:
				tadRange = range(gene.leftTAD.end, gene.rightTAD.end)
				distance = abs(gene.leftTAD.end - gene.rightTAD.end)
			else:
				tadRange = range(gene.leftTAD.end, gene.rightTAD.start)
				distance = abs(gene.leftTAD.end - gene.rightTAD.start)
			
			plotRange = range(0, distance)
			
			
			
			geneGainedEQTLs = 0
			geneLostEQTLs = 0
			
			
			#Cover for the case when there are no gains, but losses and vv
			if len(gene.lostEQTLs.keys()) < 1 and len(gene.gainedEQTLs.keys()) > 0:
				samples = gene.gainedEQTLs.keys() 
			elif len(gene.gainedEQTLs.keys()) < 1 and len(gene.lostEQTLs.keys()) > 0:
				samples = gene.lostEQTLs.keys()
			else: #combine the samples
				samples = gene.lostEQTLs.keys() + gene.gainedEQTLs.keys()
			
			if len(samples) > 1 and len(samples) < 10:
				
				startPlot = int(str(int(math.ceil(len(samples) / 2.0))) + "21")
				print startPlot
			elif len(samples) > 9:
				startPlot = int(str(int(math.ceil(len(samples) / 2.0))) + "31")
				print startPlot
				
			else:
				startPlot = int(str(int(math.ceil(len(samples) / 2.0))) + "11")
			
			# if len(gene.lostEQTLs.keys()) > 7: #look for interesting genes with recurrent mutations
			# 	print gene.name
			# 	
			# else:
			# 	continue
			
			for sample in samples:
				
				
				
				if sample in gene.lostEQTLs:
					
					print "sample ", sample, "lost: ", len(gene.lostEQTLs[sample])
				if sample in gene.gainedEQTLs:
					
				
					print "sample ", sample, "gained: ", len(gene.gainedEQTLs[sample])
					
					#Get the TAD that the SV ends in. 
					otherTAD = gene.gainedEQTLs[sample][1]
					#The TAD range is the entire length of the TAD for now, but eventually we also need to take into account that inside this TAD some eQTLs are also deleted and are thus probably not gained. 
					
					print otherTAD
					gainedTadRange = range(otherTAD.start, otherTAD.end)
					gainedTadDistance = abs(otherTAD.start - otherTAD.end)
				
					gainedPlotRange = range(0, gainedTadDistance)
					gainedEQTLPositions = np.zeros(gainedTadDistance)
				
					for eQTL in gene.gainedEQTLs[sample][0]: #The TAD is the last element
				
						eQTLPos = eQTL.start
						
						geneGainedEQTLs += 1
						
						#The position in the array depends on the positions within the actual tad
						tadInd = gainedTadRange.index(eQTLPos)	
						gainedEQTLPositions[tadInd] = 1
				
					#Plot the gains per sample per gene
					print "Plotting gains for sample ", sample, " and gene ", gene.name
					plt.subplot(startPlot)
					plt.subplots_adjust(hspace=2)
					startPlot += 1
					plt.title(sample)
					plt.plot(gainedPlotRange, gainedEQTLPositions)
					plt.ylim(-2,2)
					label = str(otherTAD.start) + "-" + str(otherTAD.end)
					plt.xlabel(label)
					
					#plt.show()
				
				# if sample in gene.lostEQTLs:
				# 	eQTLPositions = np.zeros(distance)
				# 
				# 	geneStartPos = tadRange.index(gene.start)
				# 	geneEndPos = tadRange.index(gene.end)
				# 	eQTLPositions[geneStartPos:geneEndPos] = 1
				# 	print eQTLPositions
				# 	
				# 	for eQTL in gene.lostEQTLs[sample]:
				# 		
				# 		eQTLPos = eQTL.start
				# 		
				# 		if eQTLPos < gene.leftTAD.end or eQTLPos > gene.rightTAD.start: #skip the eQTL if it is outside of the TAD boundary
				# 			continue
				# 		
				# 		geneLostEQTLs += 1
				# 		
				# 		#The position in the array depends on the positions within the actual tad
				# 		tadInd = tadRange.index(eQTLPos)
				# 		
				# 		eQTLPositions[tadInd] = -1
				# 
				# 	 
				# 	
				# 
				# 	#Plot the gains and losses for each gene per sample
				# 	print "Plotting losses for sample ", sample, " and gene ", gene.name
				# 	print startPlot
				# 	plt.subplot(startPlot)
				# 	plt.subplots_adjust(hspace=2)
				# 	startPlot += 1
				# 	
				# 	if startPlot % 10 == 0:
				# 		startPlot += 1
				# 	
				# 	plt.title(sample)
				# 	plt.plot(plotRange, eQTLPositions)
				# 	plt.ylim(-2,2)
				# 	label = str(tadRange[0]) + "-" + str(tadRange[len(tadRange)-1])
				# 	plt.xlabel(label)
				# 	plt.show()
				# 	exit()
				
			plt.show()
			exit()
			# plt.plot(plotRange, eQTLPositions)
			# plt.ylim(-2,2)
			# plt.show()
			# 	
			#return
		# print "plotting"
		# plt.boxplot(gainedEQTLs)
		# plt.show()
		# plt.boxplot(lostEQTLs)
		# plt.show()
		#return
		
		
		
	def visualizeDelta(self, genes):
		
		1+1
	
	