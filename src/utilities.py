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

def writeToCsv(file, information, append):
	"""
		Give a dictionary, we use the keys as headers for the csv, and immediately dump the values in a csv file.
		
		Parameters:
			file (str): str containing the location of the file to write to
			information (dict): dictionary with values to write to a csv (key: list)
			append (bool): True or False, should we append to an existing file or overwrite, respectively? 
	"""
	if append is False:
	
		with open(file, "wb") as outfile:
			writer = csv.writer(outfile)
			writer.writerow(information.keys())
			writer.writerows(zip(*information.values()))
	else:
		with open(file, "a") as outfile: #here we use append
			writer = csv.writer(outfile)
			writer.writerows(zip(*information.values()))
			
			
def writeToCsvManual(outFile, annotatedRegions):
	"""
		Give a dictionary, we use the keys as headers for the csv, and immediately dump the values in a csv file.
		This function exists because the above one is somehow not functioning. It is very quick and dirty. 
		
		Parameters:
			outFile (str): str containing the location of the file to write to
			annotatedRegions (dict): dictionary with values to write to a csv (key: list)
			
		TODO:
			- Fix notation: it does not write to a csv, but to a tsv.
			- Fix the writeToCsv function, then this one will be deprecated
	"""
	
	#First write the header
	header = 'chr1\ts1\te1\tchr2\ts2\te2\tnearestGeneDistance\tpLI\tRVIS\toverlappingTadCount\n'
	
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
			
			line += str(annotatedRegions['nearestGeneDistance'][annotationInd]) + '\t'
			line += str(annotatedRegions['pLI'][annotationInd]) + '\t'
			line += str(annotatedRegions['RVIS'][annotationInd]) + '\t'
			#line += str(annotatedRegions['HPO'][annotationInd]) + '\t'
			line += str(annotatedRegions['overlappingTadCount'][annotationInd]) + '\t'
			
			f.write(line)
			f.write('\n')
			
			