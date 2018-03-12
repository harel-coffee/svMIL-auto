#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""	
	Input: tab-delimited file with variants (or simply, regions). Required format: columns chr 1, s1, e1, chr2, s2, e2, somatic (yes/no) per region on the rows. 
	Output: null
	Functionality: starting point of the annotation pipeline. Sorts the input file and starts the annotation process.
	
	TODO:
		- Update input data format
"""

### Imports ###
import sys
import numpy as np
from inputDataParser import InputDataParser
from annotator import Annotator

### Code ###

#1. Initialize
annotator = Annotator()

#2. Read input data
inFile = sys.argv[1]
inputDataParser = InputDataParser()
inputData = inputDataParser.parseInputData(inFile)

#3. Prepare input data, including sorting for speed
preparedInputData = inputDataParser.prepareInputData(inputData)

#4. Do the annotation
annotator.annotate(preparedInputData)

