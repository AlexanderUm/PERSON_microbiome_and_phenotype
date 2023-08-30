#-------------------------------------------------------------------------------
# Function: Filter phyloseqs based on the prevalence and 
#           root the tree at midpoint.  
#-------------------------------------------------------------------------------
phy_taxa_filter <- function(phy_obj, 
                          prev_fraction = 0, 
                          abund = 0, 
                          abund_transform_to_prop = TRUE,
                          group_col = NULL) {
  
  # Required libraries 
  require(phyloseq)
  require(dplyr)
  
  #-----------------------------------------------------------------------------
  # Messages 
  #-----------------------------------------------------------------------------
  
  if (!is.null(group_col)) {
    
    message("Groups column is selected. Minimal prevalece and abundance are
              applyed per group!")
    
    gr.vec <- sample_data(phy_obj) %>% 
              as.matrix() %>% 
              .[, group_col]
    
  } else {gr.vec <- rep("A", nsamples(phy_obj))}
  
  
  # IF abundance is not specified as count transform to proportion 
  if(abund_transform_to_prop){
    
    ps <- transform_sample_counts(physeq = phy_obj, function(x){x/sum(x)})
    
  } else {ps <- phy_obj}
  
  # OTU table 
  otu <- as(otu_table(ps), "matrix")
  
  # Transpose if necessary
  if(taxa_are_rows(ps)){otu <- t(otu)}
  
  # Make binary OTU table
  otu.b <- otu
  
  otu.b[otu.b > 0] <- 1
  
  
  tax.to.keep <- c()
  
  for(i in unique(gr.vec)) {
    
    # Make otu table subset
    otu.unique.b <- otu.b[gr.vec == i, ] 
    
    otu.unique <- otu[gr.vec == i, ]
    
    # Remove taxa with prevalence less then cutoff
    taxa.prev.cut <- (colSums(otu.unique.b)/nrow(otu.unique.b) >= prev_fraction)
    
    # Subset 
    otu.prev <- otu.unique[, taxa.prev.cut] %>% as.data.frame()
    
    tax.set <- apply(otu.prev, 2, function(x){TRUE %in% (x >= abund)}) %>% 
                    which(isTRUE(.)) %>% 
                    names()
    
    tax.to.keep <- c(tax.to.keep, tax.set)
  }
  
  
  prune_taxa(unique(tax.to.keep), phy_obj) %>% 
    return()
  
} # END OF THE FUNCTION


  

