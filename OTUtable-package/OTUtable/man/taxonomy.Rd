\name{taxonomy}
\alias{taxonomy}

\title{
Taxonomic assignments of OTUs
}
\description{
A dataset containing the taxonomy of each OTU in the otu_table
Produced from mothur_output using clean_taxonomy() followed by expand_taxa()
}
\usage{
data(taxonomy)
}
\format{
A dataframe with 7 columns (taxonomic levels) and 6,208 rows	(OTUs)
}


\details{
Classified using our Freshwater database, followed by Greengenes - for the full workflow, visit the McMahonLab Github page
Some OTUs are missing; these were removed by subsampling of the OTU table
The presence of both blank (__) assignments and "unclassified" assignments are the result of the dual classification.
}

\author{
Alexandra Linz <amlinz16@gmail.com>
}

