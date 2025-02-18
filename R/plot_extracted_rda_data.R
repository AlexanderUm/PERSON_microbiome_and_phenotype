#-------------------------------------------------------------------------------
# Plot dbRDA scatterpots 
#-------------------------------------------------------------------------------

plot_extracted_rda_data <- function(extracted_data, 
                                    group_col,
                                    add_elepses = FALSE, 
                                    add_stat_text = FALSE,
                                    sig_text_size = 2.5, 
                                    color_vec = NULL, 
                                    dists_named_vec = NULL) {
  
  # Plot individual plots 
  ind.plots <- list()
  
  for(i in names(extracted_data)) {
    
    # Plot title 
    if(is.null(dists_named_vec)) {
      
      p.title <- paste0("Distance: ", i)
      
    } else {
      
      p.title <- paste0("Distance: ", 
                        names(dists_named_vec)[dists_named_vec == i])
      
    }
      
    # Extract data 
    p.df.inst <- extracted_data[[i]]$main 
    
    # Base plot
    p.out <- ggplot(p.df.inst,
                     aes(x = dbRDA1, 
                         y = dbRDA2)) +
              ggtitle(p.title) +
              coord_fixed() +
              theme_bw() + 
              theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) + 
              xlab(extracted_data[[i]]$var_expl[1]) +
              ylab(extracted_data[[i]]$var_expl[2]) +
              theme(plot.title = element_text(size=11, 
                                              face="italic"))
    
      # Add ellipse layer
      if(add_elepses) {
      
        p.out  <- p.out +
                    stat_ellipse(aes(fill = .data[[group_col]]), 
                                 geom = "polygon", 
                                 alpha = 0.15, 
                                 level = 0.90, 
                                 size =2) 
      }
    
      # Add custom color
      if(!is.null(color_vec)) {
        
        p.out <- p.out + 
                  scale_color_manual(values = color_vec) +
                  scale_fill_manual(values = color_vec)
      }
      
      # Add points layer 
      p.out <- p.out + 
                geom_point(aes(color = .data[[group_col]]),
                           alpha = 1, 
                           size = 2)
      
      # Add statistical text if present 
      if(add_stat_text) {
        
        sig.text.df <- extracted_data[[i]][["stat_lable"]]
        
        p.out <- p.out + 
                  geom_text(data = sig.text.df, 
                             aes(x = x_coord, 
                                 y = y_coord, 
                                 label = lab_text), 
                             hjust = 0, parse = T, 
                             size = sig_text_size)
        
      }
      
      ind.plots[[i]] <- p.out
      
  }
  
  #-------------------------------------------------------------------------------
  # Combine plots 
  #-------------------------------------------------------------------------------
  ind.plots.f <- lapply(ind.plots,
                        function(x){x + theme(legend.position="none")})
  
  legd.inst <- get_legend(ind.plots[[1]])
  
  p.grid.1 <- plot_grid(plotlist = ind.plots.f, 
                        ncol = 2)
  
  p.grid.2 <- plot_grid(p.grid.1, 
                        legd.inst, 
                        ncol = 2, 
                        rel_widths = c(0.85, .15))
  
  res.out <- list("Comb" = p.grid.2, 
                  "Ind" = ind.plots)
  
  return(res.out)
  
}
