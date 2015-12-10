North Temperate Lakes-Microbial Observatory
-------------------------------------------
NOTE: THIS DATA IS NOT YET IN ITS FINAL FORM. FURTHER CLASSIFICATION NEEDED. (12-10-15)
The Microbial Observatory project seeks to observe microbial communities over long time scales.
Our time series includes bacterial 16S data collected in 2005, 2007, 2008, and 2009, from up to 8 bog lakes in Vilas County, Wisconsin.

Contact:
Alexandra Linz <amlinz16@gmail.com>
Katherine McMahon <trina.mcmahon@wisc.edu>

https://mcmahonlab.wisc.edu
https://github.com/McMahonLab

Copyright(c) 2015, Katherine D. McMahon
All rights reserved

If you use this data or the package OTUtable in your research, please cite:
"manuscript in prep"


Repo Structure
--------------
North_Temperate_Lakes-Microbial_Observatory
|- README.txt					#This file
|- Data
||- 16S_data
|||- bog_repseqs_07Jul15.fasta			#Representative sequences for each OTU
|||- bogs_OTUtable_07Jan15.csv			#Relative abundance table of OTUs
|||- bogs_reclassified_expanded_19Oct15.csv	#Classifications of each OTU		
||- metadata
|||- 2005_to_2009_field_book.csv		#Environmental metadata from each sample collection site and date
|||- NTL-MO_sample_metadata.xlsx		#Information about the collection and sequencing of each sample
||- mothur_output
|||- qc.bogs.clean.min25.phylip.an.0.02.cons.taxonomy	#Taxonomy output file
|||- qc.bogs.clean.min25.phylip.an.0.02.subsample.shared	#Rarefied OTU table
|- OTUtable-package
|||- data						#Datasets for loading in package OTUtable
||||- metadata.rda					#Equivalent to 2005_to_2009_field.book.csv		
||||- otu_table.rda					#Equivalent to bogs_OTUtable_07Jan15.csv
||||- taxonomy.rda					#Equivalent to bogs_reclassified_expanded_19Oct15.csv
|||- man
||||- (many .Rd files)					#Help files for each function in OTUtable
|||- R							#Source code for OTUtable
||||- alpha_biodiversity_metrics.R
||||- data_processing.R
||||- format_mothur_output.R
||||- water_column_plots.R
|||- .Rbuildignore					#Ignore file for package building
|||- .Rdata						#Part of R history
|||- .Rhistory						#History of workspace
|||- DESCRIPTION					#Package description
|||- NAMESPACE						#Package namespace
|||- OTUtable.Rproj					#Rstudio package building project
||- OTUtable
||- OTUtable_1.0.tar.gz					#Downloadable package for Macs/Linux
||- OTUtable_1.0.zip					#Downloadable package for Windows
|- Scripts+Workflows