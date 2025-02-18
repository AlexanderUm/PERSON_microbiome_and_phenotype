#-------------------------------------------------------------------------------
# Calculate several dissimilarities distances and return a list with them
#-------------------------------------------------------------------------------\
phy_dist_ls <- function(phylo, 
                        dists =  c("unifrac", "wunifrac", 
                                   "jaccard", "bray")) {
  
  require(phyloseq)
  require(ape)
  
  ls.dist <- list()
  
  # Resolve not binary trees
  if(!is.binary(phy_tree(phylo))) {
    
    phy_tree(phylo) <- multi2di(phy_tree(phylo))
    
  }
  
  
  for (i.dist in dists) {
    
    if(i.dist == "jaccard") {
      
      dist.inst <- phyloseq::distance(phylo, i.dist, binary = TRUE)
      
    } else {
      
      dist.inst <- phyloseq::distance(phylo, i.dist)
      
    }
    
    ls.dist[[i.dist]] <- dist.inst
  }
  
  return(ls.dist)
  
}
