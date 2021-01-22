# svMIL2: predicting the pathogenic effect of somatic non-coding structural variants disrupting the 3D genome through Multiple Instance Learning

svMIL2 is the improved version of svMIL, tested across multiple cancer types.

For details on the previous version of svMIL, refer to: https://doi.org/10.1093/bioinformatics/btaa802

# How to use: setting up the data

All used data (SVs, SNVs, CNVs and RNA-seq data) needs to be requested from the HMF: https://www.hartwigmedicalfoundation.nl/en/appyling-for-data/.

Used regulatory data is provided in the data folder. For more details on data sources and used tissue types, please refer to Table S1 (TBA).

# How to use: preprocessing

All steps for preprocessing are listed in preprocess.sh

Provided data in the data folder has been pre-processed. Only the TMM normalization step of the RNA-seq data from the HMF dataset is required.

# How to use: generating all paper figures

plotting.sh lists all steps and instructions to generate the figures.

# Requirements

All code was tested using:
- Python 3.7.3
- Numpy 1.16.4
- Scipy 1.3.0
- Scikit-learn 0.22.2.postl
- Matplotlib 3.2.1
- Statsmodels 0.10.0


