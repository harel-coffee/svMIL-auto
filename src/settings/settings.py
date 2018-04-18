#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""
	Script in which settings are defined. Settings should be grouped in an understandable way. Each group is a dictionary and can contain multiple settings, in which the keys are the setting names, and the values the actual
	setting values. 
"""

database = dict(
	
	database = 'flatFile' #type of database to use, idea is that it should be easy to switch between databases
)

#Connection information for Neo4J
neo4j = dict(
	
	httpPath = 'bolt://localhost:7687',
	user = 'neo4j',
	password = 'test'
	
)

#determine which features to enable. The idea is that features can be turned off if necessary without having to modify the code.
#NOTE: does not work well with gene-based features yet, (including overlapping gene count, pLi and RVIS)
features = dict(
	numberOfOverlappingGenes = True, 
	numberOfDisruptedTadBoundaries = True,
	pLi = True,
	RVIS = True,
	hiCDegree = True,
	hiCBetweenness = True
)


#Locations of the input files
inputFiles = dict(
	geneList = '../data/Genes/MP_Genelist_HGNC_v2.txt',
	tads = '../data/tads/tads.csv',
	hiCDegree = '../data/HiC/HUVEC_intrachromosomal_degree.csv',
	hiCBetweenness = '../data/HiC/HUVEC_intrachromosomal_betweenness.csv'
)

parameters = dict(
	geneDistance = 200000 #we look at genes within 2Mb of the SV
)