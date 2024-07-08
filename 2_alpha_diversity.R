#-------------------------------------------------------------------------------
# Load data
#------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")

set.seed(prm.ls[["General"]][["Seed"]])

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("R/phy_taxa_filter.R")
source("R/phy_alpha.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
# Data sets to use
ps.set <- prm.ls$Alpha$data_set_ps

tax.lvls <- prm.ls$Alpha$Tax_lvl

norms <- prm.ls$Alpha$Norm

# Columns to use
gr.col <- prm.ls$General$Group_col

id.col <- prm.ls$General$Part_id_col

# Alpha parameters
alpha.ind <- prm.ls$Alpha$measures

kw.alpha <- prm.ls$Alpha$alpha_cut_kw

p.text.size <- 3.1

out.dir <- prm.ls[["Alpha"]][["out_dir"]]

# Results object
alpha.res.ls <- list()


################################################################################
# Run alpha diversity analysis
################################################################################

for(i.set in ps.set) {
  
  # Extract data
  ps.inst <- data.ls[[i.set]][["PS"]][[tax.lvls]][[norms]]
  
  meta.inst <- data.ls[[i.set]][["meta"]]
  
  # Calculate alpha diversity
  div.inst <- phy_alpha(ps.inst, measures = alpha.ind) %>% 
                bind_cols(., meta.inst)
  
  # Elongate the diversity table 
  alpha.long <- div.inst %>% 
                  pivot_longer(cols = all_of(alpha.ind), names_to = "Index")
  
  
  # Create output directories
  dir.create(paste0(out.dir, "/", i.set, "/tabs/"), 
             recursive = TRUE, showWarnings = FALSE)
  
  dir.create(paste0(out.dir, "/", i.set, "/plots/"), 
             recursive = TRUE, showWarnings = FALSE)
  
  
  #-----------------------------------------------------------------------------
  # Statistics 
  #-----------------------------------------------------------------------------
  stat.res.comb <- NULL
  
  for(i.ind in alpha.ind) { 
    
    form.inst <- paste0(i.ind, "~", gr.col)
    
    # Wilcox test
    if(length(levels(div.inst[[gr.col]])) == 2) {
      
      res.test <- wilcox.test(as.formula(form.inst), div.inst) %>% 
                    tidy() 
    }
    
    # Kruskal-Wallis -> Dunn
    if(length(levels(div.inst[[gr.col]])) > 2) {
      
      kw.res <- kruskal.test(as.formula(form.inst), div.inst) 
      
      res.test <- dunnTest(as.formula(form.inst), div.inst)$res  %>% 
                      mutate(Test = "Dunn", 
                             p.value = P.adj) %>% 
                      add_row(Comparison = "Overall", 
                              Test = "Kruskal-Wallis", 
                              p.value = kw.res$p.value,
                              .before = 1)
      
    }
    
    stat.res.comb <- res.test %>% 
                      mutate(Index = i.ind) %>% 
                      rbind(stat.res.comb, .)
    
    write.csv(stat.res.comb, 
              paste0(out.dir, "/", i.set, "/tabs/", "test_res.csv"), 
              row.names = FALSE)
    
  }
  
  
  #-----------------------------------------------------------------------------
  # Plot results
  #-----------------------------------------------------------------------------
  # Base plot
  p.base <- ggplot(alpha.long, 
                      aes(y = value, x = .data[[gr.col]])) + 
                    geom_jitter(aes(colour = .data[[gr.col]]), 
                                height = 0, 
                                width = 0.15, 
                                alpha = 0.5) +
                    geom_violin(fill = NA) + 
                    facet_grid(c("Index"), scales = "free") + 
                    theme_bw() + 
                    theme(axis.line = element_line(color='black'),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), 
                          axis.title.y = element_blank())
  
  # Plot annotation dataframe
  min.max.df <- alpha.long %>% 
                  reframe(y.min = min(value), 
                          y.max = max(value), 
                          .by = Index)
  
  res.inst <- stat.res.comb %>% 
                    mutate(p.short = round(p.value, 3)) %>% 
                    mutate(p.text = ifelse(p.short == 0, 
                                           "P<0.001", 
                                           paste0("P=", sprintf("%.3f", 
                                                                p.short)))) %>% 
                    left_join(., min.max.df, by = "Index")
  
  # Plot in correspondence with number of groups 
  # If 2 groups add wilcoxon results 
  if(length(levels(div.inst[[gr.col]])) == 2) {
    
    sig.df <- res.inst %>% 
                mutate(Start = levels(meta.inst[[gr.col]])[1], 
                       End = levels(meta.inst[[gr.col]])[2], 
                       y.sig = y.max + (y.max - y.min) * 0.2, 
                       y.inv.point = y.max + (y.max - y.min) * 0.4)
    
    p.final <- p.base +  
                  geom_signif(data = sig.df,
                              aes(xmin = Start,
                                  xmax = End,
                                  annotations = p.text,
                                  y_position = y.sig),
                              textsize = p.text.size, 
                              vjust = -0.2,
                              manual = TRUE, 
                              margin_top = 1) + 
                  geom_point(data = sig.df,
                             aes(x = End, 
                                 y = y.inv.point), 
                             x=NA) + 
                  scale_color_manual(values = aes.ls$col[[i.set]])
  }
  
  # More than two levels (Kruskal-Wallis)
  if(length(levels(div.inst[[gr.col]])) > 2) {
    
    kw.anot <- res.inst %>% 
                  filter(Test == "Kruskal-Wallis") %>% 
                  mutate(Text = paste0(Test, ": ", p.text), 
                         x.text = levels(meta.inst[[gr.col]])[1], 
                         y.text = y.max + (y.max - y.min) * 0.25, 
                         y.inv.point = y.max + (y.max - y.min) * 0.4)
    
    p.final <-  p.base + 
                  geom_text(data = kw.anot, 
                            aes(x = x.text, y = y.text, label = Text), 
                            hjust = 0.1, 
                            size = p.text.size) + 
                  geom_point(data = sig.df,
                             aes(x = x.text, 
                                 y = y.inv.point), 
                             x=NA) + 
                  scale_color_manual(values = aes.ls$col[[i.set]])
    
    if(any(kw.anot[["p.value"]] <= kw.alpha)) {
      
      # Data frame for significance levels 
      sig.df.dunn <- res.inst %>% 
                      filter(Test == "Dunn", 
                             p.value <= kw.alpha) %>% 
                      mutate(Start = str_split(.$Comparison, 
                                               " - ", simplify = TRUE)[, 1], 
                             End = str_split(.$Comparison, 
                                             " - ", simplify = TRUE)[, 2]) %>% 
                      group_by(Index) %>% 
                      mutate(y.sig = y.max + ((y.max - y.min)*(1:n()*0.35))) %>% 
                      mutate(y.inv.point = y.max + ((y.max - y.min)*((n()+2)*0.35)), 
                             y.text = y.max + ((y.max - y.min)*((n()+1)*0.45))) %>% 
                      ungroup()
      
      # Adjust y text location taking into account significance bars 
      kw.anot <- left_join(kw.anot, sig.df.dunn[c("Index", "y.text")], 
                           by = "Index", 
                           multiple = "first") %>% 
                    mutate(y.text = ifelse(is.na(y.text.y), 
                                           y.text.x, 
                                           y.text.y))
      
      # Final plot 
      p.final <- p.base +  
                  geom_text(data = kw.anot, 
                            aes(x = x.text, y = y.text, label = Text), 
                            hjust = 0.1, 
                            size = p.text.size) +
                  geom_signif(data = sig.df.dunn,
                              aes(xmin = Start,
                                  xmax = End,
                                  annotations = p.text,
                                  y_position = y.sig),
                              textsize = p.text.size, 
                              vjust = -0.2,
                              manual = TRUE, 
                              margin_top = 1) + 
                  geom_point(data = sig.df.dunn,
                             aes(x = End, 
                                 y = y.inv.point), 
                             x=NA) + 
                  scale_color_manual(values = aes.ls$col[[i.set]])
    } }
  
  alpha.res.ls[[i.set]][["plot"]] <- p.final
  
  alpha.res.ls[[i.set]][["res"]] <- stat.res.comb
  
  # Save plot 
  ggsave(filename = paste0(out.dir, "/", i.set, "/plots/alpha.png"), 
         plot = p.final, 
         width = length(levels(div.inst[[gr.col]]))*1.75 + 1.5, 
         height = length(alpha.ind)*2 + 0.25, 
         dpi = 600)
  
}


#-------------------------------------------------------------------------------
save(list = c("alpha.res.ls"), 
     file = "out/supp/2_alpha.Rdata")

# Clean environment
rm(list = ls())
gc()

