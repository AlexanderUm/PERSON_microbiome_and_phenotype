#-----------------------------------------------------------------------------
# Function to make short names. 
#-----------------------------------------------------------------------------

phy_shorten_tax_names <- function(phy_obj) {
  
  tax.t <- as(tax_table(phy_obj), "matrix")
  
  tax.t[tax.t == ""] <- NA
  
  tax.out <- c()
  
  excle.names <- "_NA|uncultured|gut metagenome|metagenome|human_gut|uncultured_archaeon|<empty>"
    
  for (r in 1:nrow(tax.t)) {
    
    ind.tax <- paste0(substr(colnames(tax.t), 1, 1), "_", tax.t[r, ])
    
    l.tax <-  tail(ind.tax[!grepl(excle.names, 
                                  ind.tax)], 1)
    
    if(length(l.tax) == 0) {
      
      l.tax <- NA
      
    }
    
    tax.out <- c(tax.out, gsub(" ", ".", l.tax))
    
  }
  
  tax.out <- gsub("-", "_", tax.out)
  
  tax.out <- gsub("\\[|\\]", "", tax.out)
  
  tax.out <- gsub("_group", "", tax.out)
  
  return(tax.out)
}
