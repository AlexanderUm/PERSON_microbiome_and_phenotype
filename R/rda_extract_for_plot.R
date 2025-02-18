#-------------------------------------------------------------------------------
# Extract data from RDA (built automatically) for plotting.
#-------------------------------------------------------------------------------
rda_extract_for_plot <- function(dists_ls, metadata, form) {
  
  require(vegan)
  require(ggvegan)
  
  rda.form <- paste0("dist ~ ", form)
  
  rda.vars <- gsub("\\(.*", "", form) %>% 
              gsub("\\*|\\+", "|", .)
  
  df.ls <- list()
  
  for(i in names(dists_ls)) {
    
    dist <- dists_ls[[i]]
    
    rda.obj <- dbrda(as.formula(rda.form), data = metadata)
    
    obj.sum <- summary(rda.obj)
    
    # Extract scores
    df.ls[[i]][["main"]] <- obj.sum$sites[, 1:2] %>%
                                as.data.frame() %>%
                                bind_cols(metadata) 
    
    # Explained variations
    df.ls[[i]][["var_expl"]] <- round(obj.sum$cont$importance[2, 1:2]*100, 2) %>% 
                                paste0(names(.), " [", ., "%]")
    
    # Extract centroids and vectors
    rda.p.auto <- autoplot(rda.obj) %>%
                          ggplot_build()
    
    df.ls[[i]][["centr"]] <- rda.p.auto$data[[4]] %>%
                                mutate(Label = gsub(rda.vars, "", .$label)) %>% 
                                select(x, y, label, Label)
                    
    df.ls[[i]][["vect"]] <- rda.p.auto$data[[3]] %>%
                                mutate(Label = gsub(rda.vars, "", label)) %>% 
                                select(x, y, label, Label)
    
    
  }
  
  return(df.ls)
  
}




#-----------------------------------------------------------------------------
# Extract data for plot
#-----------------------------------------------------------------------------


