#-------------------------------------------------------------------------------
# Normalize TMM 
#-------------------------------------------------------------------------------
# The function is adopted from MaasLin2 package 
phy_TMMnorm <- function(phylo) {
  
  NormPhylo <- phylo
  
  if(!taxa_are_rows(phylo)) {
    
    phylo <- phyloseq::t(phylo)
    
  }
  
  # Convert to Matrix from Data Frame
  OtuDf <- phylo %>%
              otu_table(phylo) %>% 
              as.matrix() %>% 
              as.data.frame()
    
 
  LibSize <- edgeR::calcNormFactors(OtuDf, method = "TMM")
  
  EffLibSize <- colSums(OtuDf) * LibSize
  
  RefLibSize <- mean(EffLibSize)
  
  # Use the mean of the effective library sizes as a reference library size
  OtuDfNorm = sweep(OtuDf, MARGIN = 2, EffLibSize, "/") * RefLibSize
  
  otu_table(NormPhylo) <- otu_table(OtuDf, taxa_are_rows = TRUE) 

  # Return as list
  return(NormPhylo)
}
