#-------------------------------------------------------------------------------
# Extract data from RDA (built automatically) for plotting.
#-------------------------------------------------------------------------------
cca_extract_for_plot <- function(cca_obj, metadata, const_var) {
  
  require(vegan)
  require(ggvegan)
  
  df.ls <- list()

  obj.sum <- summary(cca_obj)
    
  # Extract scores
  df.ls[["main"]] <- obj.sum$sites[, 1:2] %>%
                                  as.data.frame() %>%
                                  bind_cols(metadata) 
    
  # Explained variations
  df.ls[["var_expl"]] <- round(obj.sum$cont$importance[2, 1:2]*100, 2) %>% 
                                  paste0(names(.), " [", ., "%]")
    
  # Extract centroids and vectors
  rda.p.auto <- autoplot(cca_obj) %>%
                            ggplot_build()
  
  # Center - make new labels 
  center.df <- rda.p.auto$data[[4]] %>%
                      mutate(Label = label) %>% 
                      select(x, y, label, Label)
  
  for(i in const_var) {
    
    center.df <- center.df %>% 
                  mutate(Label = gsub(i, paste0(i, ": "), Label))
    
  }
  
  df.ls[["centr"]] <- center.df
    
  df.ls[["vect"]] <- rda.p.auto$data[[3]] %>%
                          mutate(Label = label) %>% 
                          select(x, y, label, Label)
  
  return(df.ls)
  
}

