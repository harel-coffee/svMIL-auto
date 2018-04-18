#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""
	This file contains a number of functions that are commonly used and are shared between different classes, for easy access. 
"""

### Imports ###
import csv

### Code ###

def writeToTsv(file, annotatedRegions):
	"""
		Give a dictionary, we use the keys as headers for the tsv, and immediately dump the values in a tsv file. All keys must have the same number of annotations for this piece of code to run smoothly. 
		
		Parameters:
			file (str): str containing the location of the file to write to
			annotatedRegions (dict): dictionary with values to write to a tsv (key: list)
	"""
	
	#As ZIP does not appear to be working, another version is to loop through the keys and write the data to a csv in a more manual way, but still without relying on providing the keys hardcoded
	
	keys = annotatedRegions.keys()
	
	with open(file, 'wb') as f:
		
		tsvHeader = '\t'.join(keys)
		
		f.write(tsvHeader)
		f.write("\n")

		print "number of annotated regions: ", len(annotatedRegions['chr1'])
		for annotationInd in range(0, len(annotatedRegions['chr1'])): #use chr1 as key as it should always correspond to the number of annotations
			line = ""
			for key in keys: #remove dependency on keys, this should work with any combination of keys
				line += str(annotatedRegions[key][annotationInd]) + '\t'
		
			f.write(line)
			f.write("\n")
			
			
def writeToCsvManual(outFile, annotatedRegions):
	"""
		@DEPRECATED
	
		Give a dictionary, we use the keys as headers for the csv, and immediately dump the values in a csv file.
		This function exists because the above one is somehow not functioning. It is very quick and dirty. 
		
		Parameters:
			outFile (str): str containing the location of the file to write to
			annotatedRegions (dict): dictionary with values to write to a csv (key: list)
			
	"""
	
	#First write the header
	header = 'chr1\ts1\te1\tchr2\ts2\te2\tidentifier\tnoOfGenesInWindow\tpLI\tRVIS\toverlappingTadBoundaries\thiCDegree\thiCBetweenness\t\n'
	
	with open(outFile, "wb") as f:
		
		f.write(header)
		
		for annotationInd in range(0, len(annotatedRegions[annotatedRegions.keys()[0]])):
			
			#quick and dirty, very hardcoded way of writing to a file
			line = annotatedRegions['chr1'][annotationInd] + '\t'
			line += annotatedRegions['s1'][annotationInd] + '\t'
			line += annotatedRegions['e1'][annotationInd] + '\t'
			line += annotatedRegions['chr2'][annotationInd] + '\t'
			line += annotatedRegions['s2'][annotationInd] + '\t'
			line += annotatedRegions['e2'][annotationInd] + '\t'
			line += annotatedRegions['identifier'][annotationInd] + '\t'
			
			line += str(annotatedRegions['noOfGenesInWindow'][annotationInd]) + '\t'
			line += str(annotatedRegions['pLI'][annotationInd]) + '\t'
			line += str(annotatedRegions['RVIS'][annotationInd]) + '\t'
			#line += str(annotatedRegions['HPO'][annotationInd]) + '\t'
			line += str(annotatedRegions['overlappingTadBoundaries'][annotationInd]) + '\t'
			line += str(annotatedRegions['hiCDegree'][annotationInd]) + '\t'
			line += str(annotatedRegions['hiCBetweenness'][annotationInd]) + '\t'
			
			f.write(line)
			f.write('\n')
			
			