
# Function for ploting PCoA scatter plot with centroids and poligons 
phy_pcoa_plot <- function(phyloseq_obj, 
                          group_column = NULL, 
                          distances = c("jaccard", "bray", 
                                        "wunifrac", "unifrac")) {
  
  axis.score.all <- NULL
  
  hull.df.all <- NULL
  
  for (i in distances) {
    
    # Calculate distance
    dist.inst <- distance(phyloseq_obj, method = i, type = "samples")
    
    # Make PCoA object 
    a.pcoa <- ape::pcoa(dist.inst)
    
    # Extract group column  
    gr.col <- phyloseq_obj %>% 
      sample_data() %>% 
      as.data.frame() %>% 
      pull(group_column)
    
    # Extract percentage 
    var.c <- round((a.pcoa$values$Relative_eig[1:2]*100), 1)
    
    # Extract scores 
    axis.score <- a.pcoa$vectors[,c(1,2)] %>% 
      as.data.frame() %>% 
      mutate(!!group_column := gr.col) %>% 
      group_by(across(group_column)) %>% 
      mutate(Center_x = mean(Axis.1), 
             Center_y = mean(Axis.2), 
             Distance = i, 
             Vrariations = paste0("(x=", var.c[1], 
                                  "%; y=", var.c[2], "%)"))
    
    axis.score.all <- rbind(axis.score.all, axis.score)
    
    # Getting convex hull for polygons 
    hull.df <- axis.score %>%
      group_by(across(group_column)) %>%
      slice(chull(Axis.1, Axis.2)) 
    
    hull.df.all <- rbind(hull.df.all, hull.df)
  }
  
  
  if(is.null(group_column)) {
    
    p <- ggplot(axis.score.all) + 
      geom_point(aes_string(x = "Axis.1", 
                            y = "Axis.2")) + 
      theme_bw() + 
      theme(axis.title = element_blank()) +
      facet_wrap(c("Distance", "Vrariations"), ncol = 2, scales = "free") 
    
  } else {
  
  # ggplot itself 
  p <- ggplot(axis.score.all) + 
    geom_point(aes_string(x = "Axis.1", 
                          y = "Axis.2", 
                          color = group_column)) + 
    
    geom_segment(aes_string(x = "Axis.1", 
                            xend = "Center_x", 
                            y = "Axis.2", 
                            yend = "Center_y", 
                            color = group_column), 
                 alpha = 0.2) +
    
    geom_point(aes_string(x = "Center_x", 
                                       y = "Center_y", 
                                       fill = group_column), 
                            size = 3,
                            color = "black", 
                            shape = 21) + 
    
    geom_polygon(data = hull.df.all, 
                 aes_string(x = "Axis.1", 
                            y = "Axis.2", 
                            fill = group_column), 
                 alpha = 0.20) + 
    
    theme_bw() + 
    theme(axis.title = element_blank()) +
    facet_wrap(c("Distance", "Vrariations"), ncol = 2, scales = "free") 
  
  }
  
  return(p)
  
}
