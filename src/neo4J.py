#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""
	!!! No longer maintained and is deprecated as Neo4J will likely not be used in the project. !!!

"""

import sys
sys.path.append('settings/') #Add settings path
import settings
from neo4j.v1 import GraphDatabase, basic_auth


class Neo4J:
	
	#handle = None #handle for interacting with the database
	handle = None
	
	def __init__(self):
		self.initiateConnection()
	
	#Initiates connection to the Neo4J database at the path provided in the settings file. Username and password are also there
	####! This does not work from a function, if it is set, and later the handle is accessed, the connection has been closed. Not very neat. 
	def initiateConnection(self):
		driver = GraphDatabase.driver(settings.neo4j['httpPath'], auth=basic_auth(settings.neo4j['user'], settings.neo4j['password']))
		self.handle = driver.session()
	
	#Executes the given query and returns the results. 
	def executeQuery(self, query):
		######! I need to put the connection code here for some reason, otherwise the connection is closed according to Neo4J. Why?
		driver = GraphDatabase.driver(settings.neo4j['httpPath'], auth=basic_auth(settings.neo4j['user'], settings.neo4j['password']))
		self.handle = driver.session()
		######
		
		with self.handle.begin_transaction() as tx:
			results = list(tx.run(query))
			
		return results
		
	#May need an additional function for bulk querying
	
	
	#This class knows how to format our queries to cypher format.
	#There are probably better and more general ways to set up specific types of query (like with the other driver), but it is restricting in what we can do and also time-consuming.
	#I will be generating the queries 'manually' here for now.
	
	#Compute the nearest tad distance feature
	def computeNearestTadDistance(self, region):
		#Generate the queries
		tadDistanceQueries = self.generateTadDistanceQueries(region)
		#Execute the queries
		results = dict()
		for query in tadDistanceQueries:
			results[query] = self.executeQuery(query)
			
		#Process the results and obtain the nearest TAD distance
		[nearestTadDistance, nearestDistanceQuery] = self.computeGeneralDistanceToAnnotation(results)
		
		return nearestTadDistance
		
	#Generate the queries
	def generateTadDistanceQueries(self, region):
		
		matchTad = '"TAD" in labels(annotationType)' #specific part for TADs
		tadDistanceQueries = self.generateGeneralDistanceQueries(region, matchTad, None)
		return tadDistanceQueries
	
	#Still lots of overlap with the nearest TAD distance function, but ok for now
	def computeNearestEnhancerDistance(self, region):
		#Generate the queries
		enhancerDistanceQueries = self.generateEnhancerDistanceQueries(region)
		#Execute the queries
		results = dict()
		for query in enhancerDistanceQueries:
			results[query] = self.executeQuery(query)
			
		#Process the results and obtain the nearest TAD distance
		[nearestEnhancerDistance, nearestDistanceQuery] = self.computeGeneralDistanceToAnnotation(results)
		
		return nearestEnhancerDistance
		
		
	def generateEnhancerDistanceQueries(self, region):
		matchEnhancer = '"Enhancer" in labels(annotationType)' #specific part for TADs
		enhancerDistanceQueries = self.generateGeneralDistanceQueries(region, matchEnhancer, None)
		return enhancerDistanceQueries
		
	def generateTadOverlapQuery():
		
		1+1
		
	def generateEnhancerOverlapQuery():
		
		1+1
		
	#Rather than having these separate queries, we use the settings file to determine which features we wish to include
	#Then based on if/else, we determine what the query should look like.
	def queryGeneFeatures(self, region):
		
		geneFeaturesQueries = self.generateGeneFeaturesQueries(region)
		results = dict()
		for query in geneFeaturesQueries:
			results[query] = self.executeQuery(query)
			
		#Process the results and obtain the nearest TAD distance
		[nearestGeneDistance, nearestDistanceQueryIndex] = self.computeGeneralDistanceToAnnotation(results)
		
		##### Check if the query matches, it may not be the right index
		#Also obtain pLi and RVIS (for the gene that is the nearest, also obtain which query belongs to this)
		nearestDistanceQuery = results.keys()[nearestDistanceQueryIndex]
		
		pLi = results[nearestDistanceQuery][0]["pLi"]
		RVIS = results[nearestDistanceQuery][0]["RVIS"]
		
		print pLi
		print RVIS
		
		geneFeatures = dict()
		geneFeatures["pLi"] = pLi
		geneFeatures["RVIS"] = RVIS
		geneFeatures["nearestGeneDistance"] = nearestGeneDistance
		
		return geneFeatures
		
	def generateGeneFeaturesQueries(self, region):
		
		#based on the settings, check which components need to be in the query
		#we need to query for a gene, and at the same time also request the pLi and RVIS scores
		#This means that the match is different (match Gene instead of TAD or enhancer), and the return statement is more detailed
		
		returnAnnotations = 'return difference' #The downside of this approach is that we make 4 queries for obtaining the nearest gene, but for each we return the features as well, which is not necessary. Is there a better way? 
		
		if settings.features["pLi"] is True:
			returnAnnotations += ', annotationType.pLi as pLi'
		if settings.features["RVIS"] is True:
			returnAnnotations += ', annotationType.RVIS as RVIS'
		#The distance does not need to be added, this is returned in each query anyway
		
		if returnAnnotations == 'return difference':
			returnAnnotations = None #use the default in the distance queries
		else:
			returnAnnotations += ' limit 1'
		
		matchType = matchGene = '"Gene" in labels(annotationType)' #specific part for TADs
		geneFeaturesQueries = self.generateGeneralDistanceQueries(region, matchType, returnAnnotations)
		
		return geneFeaturesQueries

	#Query the database for the pLi given a gene identifier
	def obtainPliScore(self, nearestGeneId):
		
		pLiQuery = 'match (gene:Gene) where ID(gene) = ' + str(nearestGeneId) + ' return gene.pLi as pLi;'
		result = self.executeQuery(pLiQuery)
		
		pLi = result[0]["pLi"] #here we could in principle also immediately obtain the RVIS query
		
		return pLi
	
	#Mergign RVIS with pLi is an option
	def obtainRvisScore(self, nearestGeneId):
		
		rvisQuery = 'match (gene:Gene) where ID(gene) = ' + str(nearestGeneId) + ' return gene.RVIS as RVIS;'
		result = self.executeQuery(rvisQuery)
		
		RVIS = result[0]["RVIS"] #here we could in principle also immediately obtain the RVIS query
		
		return RVIS
		
	#Provide the query output from the nearest distance to enhancer or tad, and compute the actual distance regardless of annotation type. 
	def computeGeneralDistanceToAnnotation(self, results):
		
		queries = results.keys()
		
		nearestDistanceFromStart = results[queries[0]][0]["difference"]
		nearestDistanceFromEnd = results[queries[1]][0]["difference"]
		nearestDistanceFromStartToEnd = results[queries[2]][0]["difference"]
		nearestDistanceFromEndToStart = results[queries[3]][0]["difference"]
		
		nearestDistance = min(nearestDistanceFromStart, nearestDistanceFromEnd, nearestDistanceFromStartToEnd, nearestDistanceFromEndToStart)
		#Determine which query resulted in the nearest distance, we can use this to obtain other properties of this matched annotation
		nearestDistanceQuery = [nearestDistanceFromStart, nearestDistanceFromEnd, nearestDistanceFromStartToEnd, nearestDistanceFromEndToStart].index(nearestDistance)
		
		return [nearestDistance, nearestDistanceQuery]
	
	def getNearestGeneId(self, results):
		queries = results.keys()
		
		#The way to obtain the values is a bit weird, previously using "difference" as key worked, but not anymore
		nearestGeneId = results[queries[0]][0]["id"]
		
		return nearestGeneId
	
	#there will be large overlap in the above queries apart from the TAD/enhancer part, keep this in a simple function that can provide these queries with a few variables that may be different. 
	def generateGeneralDistanceQueries(self, region, matchType, returnAnnotations):
		
		if returnAnnotations is None: #use the default if the return is not pre-specified
			#returnAnnotations = 'return difference, ID(annotationType) as id limit 1' #Also return the ID so that we can use this query for genes as well and obtain the nearest gene
			returnAnnotations = 'return difference limit 1'
		
		matchRegion = 'match (region:Region)-[:has]->(annotation)-[:type]->(annotationType) using index region:Region(chromosome)'
		queries = []
		
		matchChromosome = 'region.chromosome = "' + region['chromosome'] + '"'
		
		
		
		#Make this type of query if the two chromosomes are the same
		if region['chromosome'] == region['chromosome2']:
			
			matchWithStart = 'with region, abs(region.start - ' + region['start'] + ') as difference order by difference asc'	
			matchWithEnd = 'with region, abs(region.end - ' + region['end'] + ') as difference order by difference asc'
			
			matchWithStartToEnd = 'with region, abs(region.start - ' + region['end'] + ') as difference order by difference asc'
			matchWithEndToStart = 'with region, abs(region.end - ' + region['start'] + ') as difference order by difference asc'	
			
			fullQueryStart = matchRegion + ' where ' + matchChromosome + ' and ' + matchType + ' ' + matchWithStart + ' ' + returnAnnotations + ';'
			fullQueryEnd = matchRegion + ' where ' + matchChromosome + ' and ' + matchType + ' ' + matchWithEnd + ' ' + returnAnnotations + ';'
		
			fullQueryStartToEnd = matchRegion + ' where ' + matchChromosome + ' and ' + matchType + ' ' + matchWithStartToEnd + ' ' + returnAnnotations + ';'
			fullQueryEndToStart = matchRegion + ' where ' + matchChromosome + ' and ' + matchType + ' ' + matchWithEndToStart + ' ' + returnAnnotations + ';'
		
			
		else: #Make 4 queries for translocations, include the correct chromosome match, the start and end positions are conceptually different than for the queries above
			
			matchChromosome2 = 'region.chromosome = "' + region['chromosome2'] + '"'
			
			matchWithStart = 'with region, abs(region.start - ' + region['start'] + ') as difference order by difference asc'	
			matchWithEnd = 'with region, abs(region.end - ' + region['end'] + ') as difference order by difference asc'
			
			matchWithStartToEnd = 'with region, abs(region.start - ' + region['end'] + ') as difference order by difference asc'
			matchWithEndToStart = 'with region, abs(region.end - ' + region['start'] + ') as difference order by difference asc'	
			
			fullQueryStart = matchRegion + ' where ' + matchChromosome + ' and ' + matchType + ' ' + matchWithStart + ' ' + returnAnnotations + ';'
			fullQueryEnd = matchRegion + ' where ' + matchChromosome + ' and ' + matchType + ' ' + matchWithEnd + ' ' + returnAnnotations + ';'
		
			fullQueryStartToEnd = matchRegion + ' where ' + matchChromosome + ' and ' + matchType + ' ' + matchWithStartToEnd + ' ' + returnAnnotations + ';'
			fullQueryEndToStart = matchRegion + ' where ' + matchChromosome + ' and ' + matchType + ' ' + matchWithEndToStart + ' ' + returnAnnotations + ';'
			
		
		queries.append(fullQueryStart)
		queries.append(fullQueryEnd)
		queries.append(fullQueryStartToEnd)
		queries.append(fullQueryEndToStart)
		
		return queries
		
		
	def generateGeneralOverlapQuery():
		
		1+1
	
	
	