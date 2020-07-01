Folder containing SV information for the PCAWG datasets.

icgc folder: https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_sv/final_consensus_sv_bedpe_passonly.icgc.public.tgz, version 10-10-2019
metadata to map SVs/SNVs/CNVs to rna-seq:
icgc_metadata.tsv: https://dcc.icgc.org/api/v1/download?fn=/PCAWG/transcriptome/metadata/rnaseq.extended.metadata.aliquot_id.V4.tsv.gz, version 10-10-2019

SVs specific to cancer types:
liver_pcawg_parsed.txt: Processed from the icgc folder to a format that our tool can read using DataProcessing/parsePCAWGSVs.py
ov_pcawg_parsed.txt: Processed from the icgc folder to a format that our tool can read using DataProcessing/parsePCAWGSVs.py

Germline SVs:

Processing steps are also in preprocessing.sh

These are the 'old' germline SVs from DGV that we used:
GRCh37_hg19_variants_2020-02-25_filtered.txt.txt: processed from http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2020-02-25.txt using filterDgvVariants.py

These are the new germline SVs from gnomAD that are now in the latest paper version:
GRCh37_hg19_variants_2020-02-25_filtered.txt: processed from https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz using filterGnomadVariants.py