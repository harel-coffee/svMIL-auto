#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"


### Imports ###
from databaseConnector import DatabaseConnector
from utilities import writeToCsv
from utilities import writeToCsvManual

import time
import sys

### Code ###

class Annotator:
	"""	
		Input: genomic regions sorted by their chromosomes and positions (chr1 < chr2, s1 < e1, s2 < e2). 
		Output: null
		Functionality: Class for handling annotation of genomic regions. First, a connection to the database will be made (abstracted). From this database, a number of features are requested. The returned features
		are attached to the input regions and written to an output file, which is a matrix of variants by features (including the positional information of the variants). See details at each function. 
		
		Attributes:
			databaseConnector (object): the handler that we will use to talk to the database without knowing which database it is.
	"""
	
	databaseConnector = None
	
	
	
	def __init__(self):
		"""
			Constructor. A database connection will always need to be made, so do that here. 
		"""
		self.obtainDatabaseConnection()
	
	def obtainDatabaseConnection(self):
		"""
			Instantiate a database connection. It will connect us to the database defined in the settings. 
		"""
		self.databaseConnector = DatabaseConnector()
	
	
	def annotate(self, regions):
		"""
			Annotate the provided regions with features.
			
			Parameters:
				regions (numpy matrix): a Numpy matrix of dtype 'object' with on the columns chr 1, s1, e1, chr2, s2, e2, and the regions on the rows.
		"""
		
		#1. Collect all annotations from the database
		#Gene-related features
		#nearestGeneFeatures = self.databaseConnector.database.computeNearestGeneFeatures(regions)
		
		#TAD-related features (could potentially be grouped with Hi-C, but comes from an independent file atm)
		#tadFeatures = self.databaseConnector.database.computeTADFeatures(regions)
		
		#Hi-C interaction-based features
		hiCFeatures = self.databaseConnector.database.computeHiCFeatures(regions)
		
		#2. Combine all features into one big dictionary of annotations
		allAnnotations = dict(nearestGeneFeatures.items() + tadFeatures.items())
		
		#3. Link the annotations to the regions
		annotatedRegions = self.annotateRegions(regions, allAnnotations)
		
		#4. Write the annotations to an output file
		self.writeAnnotationsToFile(annotatedRegions)
		
	def annotateRegions(self, regions, annotations):
		"""
			Links the provided regions and annotations together.
			
			Parameters:
				regions (numpy matrix): a Numpy matrix of dtype 'object' with on the columns chr 1, s1, e1, chr2, s2, e2, and the regions on the rows.
				annotations (dict): a dictionary with all feature names as keys, and a list of values where each value is the feature of a region. The value must be a list, even if it has only 1 value. 
				
			Returns:
				annotatedRegions (dict): a dictionary with keys chr 1, s1, e1, chr2, s2, e2 for the regions, merged with the data in the annotations parameter. 
		"""
		
		regionsDict = dict()
		regionsDict['chr1'] = regions[:,0]
		regionsDict['s1'] = regions[:,1]
		regionsDict['e1'] = regions[:,2]
		regionsDict['chr2'] = regions[:,3]
		regionsDict['s2'] = regions[:,4]
		regionsDict['e2'] = regions[:,5]
		
		#merge the regions and annotations
		annotatedRegions = dict(regionsDict.items() + annotations.items())
		
		return annotatedRegions
		
	#Write the annotations to a file, such that each SV still has a chromosome, start, end, and behind that columns with the annotations that we are interested in. 
	def writeAnnotationsToFile(self, annotatedRegions):
		"""
			Write the annotations to a file, with in the columns chr 1, s1, e1, chr2, s2, e2, followed by the features. Each region will be on a new line.
			
			Parameters:
				annotatedRegions (dict): output from the annotateRegions function. 
		"""

		#The writeToCsv does not seem to work somehow, what if we do this by hand? Can we then write to file?
		writeToCsvManual(sys.argv[2], annotatedRegions)
			
		#write the merged dictionary to csv, the order of the annotations and regions should column-wise be the same. 
		#writeToCsv('test.csv', annotatedRegions, False)	
		
		
		
		