

# Function for ploting PCoA scatter plot with centroids and poligons 
phy_dbrda_plot <- function(phyloseq_obj, 
                           color_column, 
                           rda_formula, 
                           distances = c("jaccard", "bray", 
                                        "wunifrac", "unifrac")) {
  
  #-----------------------------------------------------------------------------
  library(vegan)
  library(tidyverse)
  #-----------------------------------------------------------------------------
  
  meta <- phyloseq_obj %>% 
            sample_data() %>% 
            as.matrix() %>% 
            as.data.frame()
  
  axis.score.all <- NULL
  
  hull.df.all <- NULL
  
  for (i in distances) {
    
    # Calculate distance
    dist.inst <- distance(phyloseq_obj, method = i, type = "samples")
    
    # Make RDA object object 
    obj.rda <- dbrda(as.formula(paste0("dist.inst ~ ", rda_formula)), 
                     data = meta)
    
    obj.sum <- summary(obj.rda)
    
    # Extract group column  
    gr.col <- meta %>% 
                pull(color_column)
    
    
    # Extract percentage 
    var.c <- round(obj.sum$cont$importance[2, 1:2]*100, 2)
    
    # Extract scores 
    axis.score <- obj.sum$sites[, 1:2] %>% 
      as.data.frame() %>% 
      mutate(cl.col = gr.col) %>% 
      setNames(c("Axis.1", "Axis.2", color_column)) %>% 
      group_by(across(color_column)) %>% 
      mutate(Center_x = mean(Axis.1), 
             Center_y = mean(Axis.2), 
             Distance = i, 
             Vrariations = paste0("(x=", var.c[1], 
                                  "%; y=", var.c[2], "%)"))
    
    axis.score.all <- rbind(axis.score.all, axis.score)
    
    # Getting convex hull for polygons 
    hull.df <- axis.score %>%
      group_by(across(color_column)) %>%
      slice(chull(Axis.1, Axis.2)) 
    
    hull.df.all <- rbind(hull.df.all, hull.df)
  }
  
  # ggplot itself 
  p <- ggplot(axis.score.all) + 
    geom_point(aes_string(x = "Axis.1", 
                          y = "Axis.2", 
                          color = color_column)) + 
    geom_segment(aes_string(x = "Axis.1", 
                            xend = "Center_x", 
                            y = "Axis.2", 
                            yend = "Center_y", 
                            color = color_column), 
                 alpha = 0.2) +
    geom_point(aes_string(x = "Center_x", 
                          y = "Center_y", 
                          fill = color_column), 
               size = 3,
               color = "black", 
               shape = 21) + 
    geom_polygon(data = hull.df.all, 
                 aes_string(x = "Axis.1", 
                            y = "Axis.2", 
                            fill = color_column), 
                 alpha = 0.15) + 
    theme_bw() + 
    theme(axis.title = element_blank()) +
    facet_wrap(c("Distance", "Vrariations"), ncol = 2, scales = "free")
  
  return(p)
  
}