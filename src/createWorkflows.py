"""
	Create workflows automatically

	Each cancer type can be run with different regulatory elements, so we need workflows
	that call the right settings and output folders.

"""


import os
import glob

mainOutFolder = 'workflows/'