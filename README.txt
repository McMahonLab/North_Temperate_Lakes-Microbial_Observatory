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
Feel free to download and use the data and code in this repo for non-commerical purposes. If you publish analyses using our data or code, cite:
Linz A.M., Crary B.C., Shade A., Owens S., Gilbert J.A., Knight R., McMahon K.D. "Bacterial community composition and dynamics spanning five years in freshwater bog lakes." Manuscript accepted to mSphere. Preprint available on bioRxiv at http://biorxiv.org/content/early/2017/04/12/127035. 

Last updated June 15, 2017

HOW TO USE THIS REPO
--------------------
If you are looking to access the NTL-MO dataset:
	- Download the files in Data/16S_data
	- Download the sequencing data at http://qiita.microbio.me
	- Install the R package "OTUtable" (included in this repo and on CRAN)

If you are looking to access the code used in our manuscript:
	- Visit the Scripts+Workflows directory - Sequence_processing includes R scripts written as part of the mothur workflow, and Analysis_scripts contains analyses performed on the output of the mothur workflow
	- The final code used for analyses and figures in the manuscript is located at /Scripts+Workflows/Analysis_scripts/manuscript_plots_accepted_2017-06-07.R
	- Functions used in these scripts are in the R package OTUtable. Source codes is included in OTUtable-package/OTUtable/R
	- The classification workflow used is available at https://github.com/McMahonLab/TaxAss (manuscript in prep)

If you have any questions or concerns:
	- Contact us!

Repo Structure
--------------
North_Temperate_Lakes-Microbial_Observatory
|- README.txt							#This file
|- annual_trends_in_otus.R					#Code accompanying Figure S7
|- Additional_analyses
||- Networks							#LSA output of networks split by lake and layer. minReads = 100, minLSA = 0.75, no lag allowed.
|- Data
||- 16S_data
|||- bog_repseqs_07Jul15.fasta					#Representative sequences for each OTU
|||- bogs_OTUtable_07Jan15.csv					#Relative abundance table of OTUs
|||- bogs_reclassified_11Mar16.csv				#Classifications of each OTU
|||- unclustered_sequences					#Same data as above, but unique sequences are not clustered into OTUs
||- deblurring_output
|||- map.bogs.txt.gz						#Mapping file
|||- bogs.min25.tar.gz						#Fasta file and biom table after deblurring
||- indicator_analysis		
|||- indicators_by_mixing_2017-01-17.csv			#Raw output of indicator analysis
|||- indicators_by_mixing_2017-04-06.csv			#Labelled output of indicator analysis - included as supplemental dataset
|||- indicators_by_times.csv					#Analysis of indicators by month rather than site. Not included in manuscript
||- metadata
|||- 2005_to_2009_field_book_for_package.csv			#Environmental metadata from each sample collection site and date
|||- NTL-MO_sample_metadata.xlsx				#MIMARKS-like file with information about the collection and sequencing of each sample
||- mothur_output
|||- otus.98.70.70.taxonomy					#16STaxAss taxonomy output file
|||- qc.bogs.clean.min25.phylip.an.0.02.cons.taxonomy		#Taxonomy output file
|||- qc.bogs.clean.min25.phylip.an.0.02.subsample.shared	#Rarefied OTU table
||- unclustered_sequences					#Versions of the OTU table and classifications using the unclustered sequences from deblurring output
|- Environmental_data						#Data on environmental parameters collected by other entities
||- LTER_data							#Buoy data, chem data, and chlorophyll measurements for available study sites
||- Minocqua_airport_data					#Solar radiation data collected by local airport
||- NTLMO_phytodata						#Phytoplankton counts collected by previous grad students
|- OTUtable-package
|||- data							#Datasets for loading in package OTUtable
|||- man
|||- R								#Source code for OTUtable
|||- DESCRIPTION						#Package description
|||- NAMESPACE							#Package namespace
|||- OTUtable.Rproj						#Rstudio package building project
||- OTUtable
||- OTUtable_1.0.0.tar.gz					#Downloadable package for Macs/Linux
||- OTUtable_1.0.0.zip						#Downloadable package for Windows
||- OTUtable_1.1.0.tar.gz					#Updated release
||- OTUtable_1.1.0.zip						#Updated release
|- Plots							#Intermediate and final files for figures in manuscript + supplemental	
|- Scripts+Workflows
||- Sequence_processing
|||- 16S_deblurred_table_mothur_workflow.docx			#Workflow from sequence data to OTU table
|||- bog_repseq_generation.R					#Generates representative sequence fasta file
|||- names_from_shared.R					#Generates a .names from a .shared file
|||- remove_seqs_from_shared.R					#Removes sequences from .shared file that were removed in the .names file (from chimera and chloroplast removal steps)
|||- shared_to_count.R						#Converts .shared to .count file
||- Analysis_scripts						#R code used for figures and indicator analysis
|||- manuscript_plots_accepted_2017-06-07.R			#Script for analyses and figures in the final version of the manuscript
|||- manuscript_plots_submitted_2017-05-15.R			#Script for analyses and figures in the submitted version of the manuscript
|||- distance_metric_tests_4-9-16.R				#Script for comparing beta diversity metrics
|||- indicators_by_time.R					#Performs indicator analysis by month instead of site
|||- LSA_calcuation.R						#Script for performing LSA (network analysis)
|||- mixing_regime_indicator_table.R				#Produces the supplemental document of the results of indicator analysis using mixing regime as the pre-defined group
|- Manuscript_materials						
||- bogtags_accepted_bioRxiv_maintext.pdf			#PDF of accepted manuscript for submission to bioRxiv
||- bogtags_accepted_bioRxiv_supplemental.pdf			#PDF of supplemental figures accompanying accepted manuscript for submission to bioRxiv
||- Linz_supplemental_2017-04-12.pdf				#Submitted supplemental materials (pre-review)
||- Linz_manuscript_2017-04-12.docx				#Submitted manuscript (pre-review)
||- Linz_response_to_reviewers_2017-06-02.docx			#Response to reviewer comments
||- Linz_markedup_manuscript_2017-06-02.docx			#Tracked changes after review
||- response2reviewers.R					#Code for figures included/referenced in the response to reviewers
