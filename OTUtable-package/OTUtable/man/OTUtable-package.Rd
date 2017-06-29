\name{OTUtable-package}
\alias{OTUtable-package}
\alias{OTUtable}
\docType{package}
\title{
OTU table analysis functions
}
\description{
Contains functions for the analysis of an OTU table generated from 16S rRNA amplicon sequencing. It also includes the data from the North Temperate Lakes-Microbial Observatory used in the paper Linz, et al. (2017) "Bacterial community composition and dynamics spanning five years in freshwater bog lakes." Please cite this paper if you use OTUtable. The data and code used in this paper are available at <https://github.com/McMahonLab/North_Temperate_Lakes-Microbial_Observatory>.
Three data files are included: otu_table, taxonomy, and metadata. Access these by calling them with data(). There is also a fasta file associated with this dataset that is not included in this package - it can be found on our GitHub page in Data/16S_data. This package does not include functionality for fasta files; if you need this for analyses such as calculating UniFrac distance, please see the R package "phyloseq".
}
\details{
\tabular{ll}{
Package: \tab OTUtable\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2017-06-28\cr
License: \tab GPL-3 \cr
}
Functions include:
bog_subset \cr
chao1 \cr
clean_shared \cr
clean_mothur_taxonomy \cr
clean_TaxAss_taxonomy \cr
combine_otus \cr
extract_date \cr
filter_taxa \cr
grab_group \cr
make_do_matrix \cr
make_temp_matrix \cr
obs_richness \cr
pielou \cr
plot_column \cr
reduce_names \cr
remove_reps \cr
rotate \cr
shannon \cr
strat_metric \cr
year_subset \cr
zscore \cr
}
\author{
Alexandra Linz <amlinz16@gmail.com>
}
\references{
Alexandra M. Linz, Benjamin C. Crary, Ashley Shade, Sarah Owens, Jack A. Gilbert, Rob Knight, Katherine D. McMahon.
Bacterial Community Composition and Dynamics Spanning Five Years in Freshwater Bog Lakes. mSphere Jun 2017, 2 (3)
e00169-17; DOI: 10.1128/mSphere.00169-17
}
\keyword{ package }

