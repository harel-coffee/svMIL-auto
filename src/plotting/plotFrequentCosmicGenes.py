"""
	Make figure 1 showing the most frequently affected cosmic genes in different cancer types

"""

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os.path
import pandas as pd
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
import gzip

path = sys.argv[1]
sys.path.insert(1, path)
sys.path.insert(1, 'linkSVsGenes/')

import settings
from inputParser import InputParser

class DriverPlotter:

	#the cancer types to show the results for
	cancerTypes = ['HMF_Breast_hmec', 'HMF_Ovary_ov', 'HMF_Lung_luad', 'HMF_Colorectal_coad',
				   'HMF_UrinaryTract_urinaryTract', 'HMF_Prostate_prostate', 'HMF_Esophagus_esophagus', 'HMF_Skin_skin',
				   'HMF_Pancreas_pancreas', 'HMF_Uterus_uterus', 'HMF_Kidney_kidney', 'HMF_NervousSystem_nervousSystem']
	cancerTypes = ['HMF_Breast_CTCF', 'HMF_Colorectal_CTCF', 'HMF_Lung_CTCF']


	#because the cancer types in the metadata have slashes and spaces, we cannot use them, so use those
	#converted names here to read in the data.
	cancerTypeMetadataNames = {'HMF_Breast_hmec': 'Breast', 'HMF_Ovary_ov': 'Ovary', 'HMF_Lung_luad': 'Lung',
					   'HMF_Colorectal_coad': 'Colon/Rectum',
					   'HMF_UrinaryTract_urinaryTract': 'Urinary tract',
					   'HMF_Prostate_prostate': 'Prostate',
					   'HMF_Esophagus_esophagus': 'Esophagus',
					   'HMF_Skin_skin': 'Skin', 'HMF_Pancreas_pancreas': 'Pancreas',
					   'HMF_Uterus_uterus': 'Uterus', 'HMF_Kidney_kidney': 'Kidney',
					   'HMF_NervousSystem_nervousSystem': 'Nervous system',
					   'HMF_Breast_CTCF': 'Breast', 'HMF_Colorectal_CTCF': 'Colon/Rectum',
					   'HMF_Lung_CTCF': 'Lung'}


	#cancerTypes = ['HMF_Colorectal']

	#use names to search for the tissue in the cosmic cancer type list
	#this is using regex, so partial matches are ok if necessary
	cancerTypeNames = {'HMF_Breast_hmec': ['breast'], 'HMF_Ovary_ov': ['ovarian'],
					   'HMF_Lung_luad': ['lung'],'HMF_Colorectal_coad': ['colorectal'],
					   'HMF_UrinaryTract_urinaryTract': ['bladder'], 'HMF_Prostate_prostate': ['prostate'],
					   'HMF_Esophagus_esophagus': ['esophagus', 'esophageal'],
					   'HMF_Skin_skin': ['skin', 'melanoma'], 'HMF_Pancreas_pancreas': ['pancreas', 'pancreatic'],
					   'HMF_Uterus_uterus': ['uter'], 'HMF_Kidney_kidney': ['kidney', 'renal'],
					   'HMF_NervousSystem_nervousSystem': ['brain', 'nervous'],
					   'HMF_Breast_CTCF': ['breast'], 'HMF_Colorectal_CTCF': ['colorectal'],
					   'HMF_Lung_CTCF': ['lung']}

	outDirPrefix = 'output/'
	svTypes = ['DEL', 'DUP', 'INV', 'ITX']




	

			
	




	
	
	
#2. Make the plot
#DriverPlotter().plotPathogenicSVFrequency()
#DriverPlotter().plotAUC()
DriverPlotter().plotCosmicFrequencyScatter()
#DriverPlotter().plotSVContributionVenn()
#DriverPlotter().plotTadCtcfOverlap()