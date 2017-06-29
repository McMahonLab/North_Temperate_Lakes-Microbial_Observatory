library(OTUtable)	# You will need these three packages
library(ggplot2)
library(reshape2)
data(otu_table)		# Load the OTU table

# Write function to plot multiple years at once. 
annual_trends <- function(lake, otu){
  bog <- bog_subset(lake, otu_table)
  year1 <- year_subset("05", bog)
  year2 <- year_subset("07", bog)
  year3 <- year_subset("08", bog)
  year4 <- year_subset("09", bog)
  
  # Since sites have different years sampled, these if statements identify which years are present
  if(dim(year1)[2] > 0){
    # Once years present are identified, normalize and combine into a single table
    year1 <- zscore(year1)
    year2 <- zscore(year2)
    year3 <- zscore(year3)
    year4 <- zscore(year4)
    
    ztable <- cbind(year1, year2, year3, year4)
  }else if(dim(year1)[2] == 0 & dim(year3)[2] > 0){
    year2 <- zscore(year2)
    year3 <- zscore(year3)
    year4 <- zscore(year4)
    
    ztable <- cbind(year2, year3, year4)
  }else if(dim(year1)[2] == 0 & dim(year3)[2] == 0 & dim(year4)[2] > 0){
    year2 <- zscore(year2)
    year4 <- zscore(year4)
    
    ztable <- cbind(year2, year4)
  }else{
    ztable <- zscore(year2)
  }
  # Format the final table
  ztable <- melt(ztable)
  ztable$Year <- substr(ztable$Var2, start = 9, stop = 10)
  ztable$Day <- format(extract_date(ztable$Var2), format = "%j")
  
  # Save the results for plotting
  plot <- ggplot(data = ztable[which(ztable$Var1 == otu), ], aes(x = Day, y = value, group = Year, color = Year)) + geom_point() + geom_line() + theme_bw() + labs(title = paste(lake, otu, sep = ": "))
  return(plot)
}

# Example Usage - 3 letter site code includes 1st 2 for site (see Table 1) and letter 3 for layer (E = epilimnion, H = hypolimnion. OTU designation is case sensitive, and number must contain 4 digits.
plot_this <- annual_trends("TBE", "Otu0012")
plot_this
# You may get warning messages about points being removed. That means the OTU was not present in those points
# If all points were removed and no plot is produced, it was not present in that site
