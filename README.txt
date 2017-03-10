North Temperate Lakes-Microbial Observatory
-------------------------------------------
The Microbial Observatory project seeks to observe microbial communities over long time scales.
Our time series includes bacterial 16S data collected in 2005, 2007, 2008, and 2009, from up to 8 bog lakes in Vilas County, Wisconsin.

Contact:
Alexandra Linz <amlinz16@gmail.com>
Katherine McMahon <trina.mcmahon@wisc.edu>

https://mcmahonlab.wisc.edu
https://github.com/McMahonLab

Copyright(c) 2017, Katherine D. McMahon
All rights reserved
Feel free to download and use the data and code in this repo for non-commerical purposes. If you publish your analysis, please cite:
Linz A., Crary B., Shade A., Owens S., Gilbert J.A., Knight R., McMahon K.D. "Trends in Bacterial Community Composition and Dynamics Span Multiple Years in Freshwater Bog Lakes." Manuscript in prep, March 10, 2017. 

Last updated March 10, 2017

HOW TO USE THIS REPO
--------------------
If you are looking to access the NTL-MO dataset:
	- Download the files in Data/16S_data
	- Download raw sequencing data at http://metagenomics.anl.gov/?page=MetagenomeProject&project=127
	- Install the R package "OTUtable" (included in OTUtable-package and on CRAN)

If you are looking to access the code used in our manuscript:
	- Visit the Scripts+Workflows directory - Sequence_processing includes R scripts written as part of the mothur workflow, and Analysis_scripts contains analyses performed on the output of the mothur workflow
	- Functions used in these scripts are in the R package OTUtable. Source codes is included in OTUtable-package/OTUtable/R

If you have any questions or concerns:
	- Contact us!

Repo Structure
--------------
North_Temperate_Lakes-Microbial_Observatory
|- README.txt							#This file
|- Additional_analyses
||- Networks						#LSA output of networks split by lake and layer. minReads = 100, minLSA = 0.75, no lag allowed.
|- Data
||- 16S_data
|||- bog_repseqs_07Jul15.fasta					#Representative sequences for each OTU
|||- bogs_OTUtable_07Jan15.csv					#Relative abundance table of OTUs
|||- bogs_reclassified_29Dec15.csv				#Classifications of each OTU
||- deblurring_output
|||- map.bogs.txt.gz						#Mapping file
|||- bogs.min25.tar.gz						#Fasta file and biom table after deblurring
||- indicator_analysis		
|||- dimictic indicators					#Source of data in Table 4
|||- indicators_by_mixing					#Output of mixing_regime_indicator_table.R
|||- meromictic indicators					#Soure of data in Table 5
|||- polymictic indicators					#Source of data in Table 3
||- metadata
|||- 2005_to_2009_field_book_for_package.csv			#Environmental metadata from each sample collection site and date
|||- NTL-MO_sample_metadata.xlsx				#MIMARKS-like file with information about the collection and sequencing of each sample
||- mothur_output
|||- qc.bogs.clean.min25.phylip.an.0.02.cons.taxonomy		#Taxonomy output file
|||- qc.bogs.clean.min25.phylip.an.0.02.subsample.shared	#Rarefied OTU table
||- unclustered_sequences					#Versions of the OTU table and classifications using the unclustered sequences from deblurring output
|- OTUtable-package
|||- data							#Datasets for loading in package OTUtable
|||- man
|||- R								#Source code for OTUtable
|||- DESCRIPTION						#Package description
|||- NAMESPACE							#Package namespace
|||- OTUtable.Rproj						#Rstudio package building project
||- OTUtable
||- OTUtable_1.0.tar.gz						#Downloadable package for Macs/Linux
||- OTUtable_1.0.zip						#Downloadable package for Windows
|- Scripts+Workflows
||- Sequence_processing
|||- 16S_deblurred_table_mothur_workflow.docx			#Workflow from sequence data to OTU table
|||- bog_repseq_generation.R					#Generates representative sequence fasta file
|||- names_from_shared.R					#Generates a .names from a .shared file
|||- remove_seqs_from_shared.R					#Removes sequences from .shared file that were removed in the .names file (from chimera and chloroplast removal steps)
|||- shared_to_count.R						#Converts .shared to .count file
||- Analysis_scripts						#R code used for figures and indicator analysis
