#-------------------------------------------------------------------------------
# Normalize ASV count with CSS
#-------------------------------------------------------------------------------
phy_css_norm <- function(phylo) {
  
  require("phyloseq")
  require("tidyverse")
  require("metagenomeSeq")
  
  otu_table(phylo) <- phylo %>% 
                        phyloseq_to_metagenomeSeq(.) %>% 
                        cumNorm(., p=cumNormStatFast(.)) %>% 
                        MRcounts(., norm=TRUE, log=TRUE) %>% 
                        as.data.frame() %>% 
                        otu_table(., taxa_are_rows = TRUE)
  
  return(phylo)
  
}