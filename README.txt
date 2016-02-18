North Temperate Lakes-Microbial Observatory
-------------------------------------------
The Microbial Observatory project seeks to observe microbial communities over long time scales.
Our time series includes bacterial 16S data collected in 2005, 2007, 2008, and 2009, from up to 8 bog lakes in Vilas County, Wisconsin.

Contact:
Alexandra Linz <amlinz16@gmail.com>
Katherine McMahon <trina.mcmahon@wisc.edu>

https://mcmahonlab.wisc.edu
https://github.com/McMahonLab

Copyright(c) 2015, Katherine D. McMahon
All rights reserved
Last updated Feb 18, 2016

If you use this data or the package OTUtable in your research, please cite:
Linz A., Crary B., Shade A., Herren C.M., Owens S., Gilbert J.A., Knight R., McMahon K.D. "Seasonal mixing as a barrier to the development of stable bacterial communities in bog lakes." Manuscript in prep, Feb 18 2016. 

HOW TO USE THIS REPO
--------------------
If you are looking to access the NTL-MO dataset:
	- Download the files in Data/16S_data
	- Download raw sequencing data at http://metagenomics.anl.gov/?page=MetagenomeProject&project=127
	- Install the R package "OTUtable" (included in OTUtable-package and on CRAN)

If you are looking to access the code used in "Seasonal mixing as a barrier to the development of stable bacterial communities in bog lakes":
	- Visit the Scripts+Workflows directory - Sequence_processing includes R scripts written as part of the mothur workflow, and Analysis_scripts contains analyses performed on the output of the mothur workflow
	- Code for figures is available in Scripts+Workflows/Analysis_scripts/make_pdf_figures_28Jan16.R and MO-NTL_paper_supplemental_09Feb16.RMD
	- Code for tables is availabe in the same directory as above - mixing_regime_indicator_table.R and quantify_phenology.R
	- Functions used in these scripts are in the R package OTUtable. Source codes is included in OTUtable-package/OTUtable/R

If you have any questions or concerns:
	- Contact us!

Repo Structure
--------------
North_Temperate_Lakes-Microbial_Observatory
|- README.txt							#This file
|- Additional_analyses
||- Networks
|||- *_network.txt						#LSA output of networks split by lake and layer. minReads = 100, minLSA = 0.75, no lag allowed.
||- UniFrac_analysis
|||- unifrac_batch.cmd						#Batch script (Windows) used to calculate UniFrac distances
|||- parse_unifrac_output.R					#Peforms analysis from Figure 5 of the main text using UniFrac instead of Sorenson Index
||- additional_supplmentary.Rd					#Figures that were considered for the supplementary material, but not included in the final draft
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
||- allsamples_network_28Jan15.txt				#Output of LSA_calculations.R on full dataset, used to measure connectivity over time in Figure 3. minreads = 1000, minLSA = 0.7, no lag allowed.
|- OTUtable-package
|||- data							#Datasets for loading in package OTUtable
||||- metadata.rda						#Equivalent to 2005_to_2009_field.book.csv		
||||- otu_table.rda						#Equivalent to bogs_OTUtable_07Jan15.csv
||||- taxonomy.rda						#Equivalent to bogs_reclassified_expanded_19Oct15.csv
|||- man
||||- (many .Rd files)						#Help files for each function in OTUtable
|||- R								#Source code for OTUtable
||||- alpha_biodiversity_metrics.R
||||- data_processing.R
||||- format_mothur_output.R
||||- water_column_plots.R
|||- .Rbuildignore						#Ignore file for package building
|||- .Rdata							#Part of R history
|||- .Rhistory							#History of workspace
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
||- Analysis_scripts
|||- LSA_calculation.R						#Calculates LSA on NTL-MO dataset by calling fast_LSA. Can also be run on subsets of the dataset. Output is in Additional_analyses/Networks.
|||- make_pdf_figures.R						#Generates pdfs for figures in manuscript. Output is in Figures/.
|||- mixing_regime_indicator_table				#Generates Tables 3-5 in manuscript. Output is in Data/.
|||- MO-NTL_paper_supplemental_09Feb16.pdf			#Markdown output of supplemental document
|||- MO-NTL_paper_supplemental_09Feb16.Rmd			#Script for generating supplemental figures
