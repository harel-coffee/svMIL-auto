#!/usr/bin/env python

######### Documentation #########

#This file is a template for the FeatureHandlers, your code will work with the rest of the pipeline if you follow the format defined here. 

#Example author information. Your author information goes here, allowing everyone to contact the right person when necessary. If you edit an existing file, make sure to also add your information! 
__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

### Imports ###

#Put all your imports here

### Code ###

class TemplateFeatureHandler:
	"""
		This is a docstring. It is the standard way in Python to write documentation. Docstrings and good documentation are required to make collaboration possible.
		The first line in a docstring usually contains basic information about what the class/function does.
		
		Attributes:
			attribute1 (type): In here, you describe all attributes that the class has, what their types are, and what they are/do
	"""
	
	def annotate(self, regions, enabledFeatures = None):
		"""
			The annotate function is mandatory. This function should be called from a database handler, such as the FlatFileDb class. It should have one return value, which contains annotations.
			
			Parameters:
				regions (numpy matrix): a Numpy matrix of dtype 'object' with on the columns chr 1, s1, e1, chr2, s2, e2, and the regions on the rows.
				enabledFeatures (list): optional. Not yet implemented. The idea is that a list can be provided that contains the names of features that are enabled in the settings. If necessary, the listed features
				can be excluded from computation to save computational time. It can be easily implemented by having an if statement checking for that in this function.
				
			Returns:
				annotations (dict): a dictionary with annotations. {'featureName' : [values]}, where the values are in the same order as the provided regions. The value must be a list to guarantee that data is
				written to the output file in the same manner, even if it has only 1 value.
				
			TODO:
				- write your TODO here, then others can see what still needs to be fixed or added. 
		"""
		
		#annotations = your own function that returns the annotations in the right format, or it can be formatted within this function
		
		return annotations
	
	#Define the rest of your functions here