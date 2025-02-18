#-------------------------------------------------------------------------------
# Calculate alpha diversity 
#-------------------------------------------------------------------------------
phy_alpha <- function(phylo, 
                      measures = c("Observed", "Shannon", 
                                   "InvSimpson", "PhyloDiversity")) {
  
  require("phyloseq")
  require("picante")
  require("tidyverse")
  
  
  # Calculate diversity indexes 
  alpha.df <- estimate_richness(physeq = phylo, 
                                measures = measures[measures != "PhyloDiversity"])
  
  if("PhyloDiversity" %in% measures) {
    # Calculate Phylogeny Diversity (package "picante")
    alpha.df$PhyloDiversity <- picante::pd(samp = t(otu_table(phylo)), 
                                          tree = phy_tree(phylo), 
                                          include.root = FALSE)  %>% 
                                          select("PD")  %>%  
                                          unlist()
  }
  
  rownames(alpha.df) <- sample_names(phylo)
  
  return(alpha.df)
}