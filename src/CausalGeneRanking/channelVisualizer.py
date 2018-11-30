"""
	The main goal of this class is to take the genes (with neighborhood and dertivative neighborhoods mapped) and visualize the distributions (later to be used as channels) of the regulatory data within
	the nearest TAD boundaries. 

"""
import matplotlib.pyplot as plt
import numpy as np

class ChannelVisualizer:
	
	
	
	def __init__(self, genes, mode):
		
		
		#self.visualizeReference(genes)
		self.visualizeDerivative(genes)
		
		
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
			
			if gene.name != "EPS15":
			 	continue
			
			
			if len(gene.gainedEQTLs) < 1 and len(gene.lostEQTLs) < 1:
				
				continue
			
			#Check if the total distance is something that we can easily plot or if we need to make bins
			
			if gene.leftTAD == None or gene.rightTAD == None: #Also skip if there is no TAD on the left/right for now
				continue
			
			#If the gene is within a TAD, the lost eQTLs are between the left tad end, and right tad end.
			#If the gene is not within TADs, the lost eQTLs are between the end of the left tad, and start of the right tad. 
			
			distance = abs(gene.leftTAD.end - gene.rightTAD.start) #left comes before right
			#at each eQTL position, add a +1.
			print "TAD: ", gene.leftTAD.start, gene.leftTAD.end, gene.rightTAD.start, gene.rightTAD.end
			if gene.rightTAD.end != gene.leftTAD.end:
				tadRange = range(gene.leftTAD.end, gene.rightTAD.end)
				distance = abs(gene.leftTAD.end - gene.rightTAD.end)
			else:
				tadRange = range(gene.leftTAD.end, gene.rightTAD.start)
				distance = abs(gene.leftTAD.end - gene.rightTAD.start)
			#print "tad range: ", min(tadRange), max(tadRange)
			
			plotRange = range(0, distance)
			
			eQTLPositions = np.zeros([distance,1])
			geneGainedEQTLs = 0
			geneLostEQTLs = 0
			
			
			samples = gene.lostEQTLs.keys() #if there is no loss, perhaps in gain
			 
			for sample in samples:
				
				
				if sample in gene.lostEQTLs:
					
					# if len(gene.lostEQTLs[sample]) > 0:
					# 	print gene.name
					# 	exit()
					print "sample ", sample, "lost: ", len(gene.lostEQTLs[sample])
				if sample in gene.gainedEQTLs:
				
					print "sample ", sample, "gained: ", len(gene.gainedEQTLs[sample])
		
					for eQTL in gene.gainedEQTLs[sample]:
						#eQTLPos = eQTL.start
						
						geneGainedEQTLs += 1
						
						#The position in the array depends on the positions within the actual tad
						#tadInd = tadRange.index(eQTLPos)	
						#eQTLPositions[tadInd] = 1
				
				if sample in gene.lostEQTLs:
					
					for eQTL in gene.lostEQTLs[sample]:
						
						eQTLPos = eQTL.start
						print "pos of lost eQTL: ", eQTLPos
						
						if eQTLPos < gene.leftTAD.end or eQTLPos > gene.rightTAD.start: #skip the eQTL if it is outside of the TAD boundary
							continue
						
						geneLostEQTLs += 1
						
						#The position in the array depends on the positions within the actual tad
						tadInd = tadRange.index(eQTLPos)	
						eQTLPositions[tadInd] = -1
			
			
				print "sample: ", sample, " gains ", geneGainedEQTLs, " gene: ", gene.name
				print "sample: ", sample, " loses ", geneLostEQTLs, " gene: ", gene.name
			
			#gainedEQTLs.append(geneGainedEQTLs)
			#lostEQTLs.append(geneLostEQTLs)
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
		return
		
		1+1
		
	def visualizeDelta(self, genes):
		
		1+1
	
	