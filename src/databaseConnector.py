#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

import sys
sys.path.append('settings/') #Add settings path
import settings
from neo4J import Neo4J
from flatFileDb import FlatFileDB

class DatabaseConnector:
	"""
		Interface for connecting to the database that is specified in the settings. Switching to another database is kept simple in this setup. 
		
		Attributes:
			database (object): this object can be called from other objects and will allow interaction with the instantiated database class. 
	"""
	
	database = None
	
	def __init__(self):
		"""
			Read the settings file and connect to the right database handler object. 
		"""
		
		if settings.database['database'] == 'neo4j':
			self.database = Neo4J() #initialize connection to Neo4J
			
		if settings.database['database'] == 'flatFile':
			self.database = FlatFileDB()