#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""
	This script reads the original TP.txt and TN.txt files for SVs and converts them to an input format that can be recognized by the annotation pipeline prototype.
	
	The accepted input format is:
	chr1	s1	e1	chr2	s2	e2	somatic
	
	The somatic column will have a yes/no value to distinguish SNPs from SNVs. 

"""

### Imports ###


### Code ###