#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")

set.seed(prm.ls[["General"]][["Seed"]])


#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("R/phy_taxa_filter.R")
source("R/fix_taxa_names_for_plot.R")


#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
# Data sets to use
ps.set <- prm.ls$DA$data_set_ps

tax.lvls <- prm.ls$DA$Tax_lvl

norms <- prm.ls$DA$Norm

# Columns to use
gr.col <- prm.ls$General$Group_col

# Ancom parameters
min.prev <- prm.ls$DA$feature_min_prev

ancom.form <- prm.ls$DA$Formula

p.norm <- prm.ls$DA$Plot_norm

dir.out <- prm.ls[["DA"]][["out_dir"]]


# Results object
DA.res.ls <- list()

################################################################################
# Run DA analysis
################################################################################
set.gr <- expand.grid("Set" = ps.set, 
                       "Lvl" = tax.lvls, 
                       "Norm" = norms, 
                       stringsAsFactors = FALSE)

for(i.set in 1:nrow(set.gr)) {
  
  # Variables 
  set.inst <- set.gr[["Set"]][i.set]
  
  lvl.inst <- set.gr[["Lvl"]][i.set]
  
  norm.inst <- set.gr[["Norm"]][i.set]
  
  # Data sets 
  ps.inst <- data.ls[[set.inst]][["PS"]][[lvl.inst]][[norm.inst]]
  
  meta.inst <- data.ls[[set.inst]][["meta"]]
  
  # Output folders
  dir.create(paste0(dir.out, "/", set.inst, "/tabs"), 
             recursive = TRUE, showWarnings = FALSE)
  
  dir.create(paste0(dir.out, "/", set.inst, "/plots"), 
             recursive = TRUE, showWarnings = FALSE)
  
  # Filter based on prevalence
  ps.inst.f <- phy_taxa_filter(ps.inst, 
                               prev = min.prev, 
                               group_col = gr.col)
  
  anc.out <- ancombc(phyloseq = ps.inst.f, 
                      formula = ancom.form, 
                      p_adj_method = "fdr", 
                      # zero_cut = 0.90, 
                      # lib_cut = 5000, 
                      group = gr.col, 
                      struc_zero = FALSE, 
                      neg_lb = FALSE, 
                      tol = 1e-5, 
                      max_iter = 100, 
                      conserve = FALSE, 
                      alpha = 0.1,
                      global = TRUE)
  
  
  #-----------------------------------------------------------------------------
  # Plot results - ANCOM-BC 
  #-----------------------------------------------------------------------------
  for(i.p.set in p.norm) {
    
    # Singificant taxa 
    sig.tax <- anc.out$res$diff_abn %>% 
                    filter(.data[[grep(gr.col, colnames(.), value = TRUE)]]) %>% 
                    pull(taxon)
    
    DA.res.ls[["tabs"]][[set.inst]][[lvl.inst]][[norm.inst]][["sig_tax"]] <- sig.tax 
    
    
    # Extract and transform OTU table 
    otu.inst <- data.ls[[set.inst]][["PS"]][[lvl.inst]][[i.p.set]] %>% 
                   otu_table() %>% 
                   as.matrix() %>% 
                   as.data.frame() 
      
    otu.plot.sig <- otu.inst %>% 
                       mutate(feature = rownames(.)) %>% 
                       filter(feature %in% sig.tax) %>% 
                       pivot_longer(cols = - feature, 
                                    names_to = "SeqID", 
                                    values_to = "Abundance") %>% 
                       left_join(., meta.inst, by = "SeqID") %>% 
                       mutate(feature = fix_taxa_names_for_plot(feature))
  
    # Plot results
    dif.abund.tax_v2 <- ggplot(otu.plot.sig, 
                               aes(x = Abundance, 
                                   y = feature,
                                   color = .data[[gr.col]])) + 
                            geom_boxplot(outlier.alpha = 0.2) +
                            theme_bw() + 
                            xlab("Abundance (%)") + 
                            theme(axis.line = element_line(color='black'),
                                  plot.background = element_blank(),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), 
                                  axis.title.y = element_blank()) + 
                            xlim(c(0, 5)) + 
                            scale_color_manual(values = aes.ls$col[[i.set]])
    
    # Write results
    DA.res.ls[["plot"]][[set.inst]][[lvl.inst]][[norm.inst]] <- 
        list(p = dif.abund.tax_v2, 
             ntax = length(sig.tax))
    
    ggsave(filename = paste0(dir.out, "/", set.inst, "/plots/", 
                             "DA_", lvl.inst, "_", p.norm, ".png"), 
           plot = dif.abund.tax_v2, 
           width = 7, 
           height = length(sig.tax)*0.2 + 1, 
           dpi = 600)
    
    # Summary table
    ancom.res.df <- NULL
    
    for(i in names(anc.out$res)) {
      
      ancom.res.df <- anc.out$res[[i]] %>% 
                        select(-taxon) %>% 
                        setNames(paste0(i, "--", colnames(.))) %>% 
                        bind_cols(ancom.res.df, .)
    }
    
    ancom.res.df$feature <- fix_taxa_names_for_plot(anc.out$res[[1]][["taxon"]])
    
    dif.res.tab <- otu.plot.sig %>% 
                    group_by(across(all_of(c("feature", gr.col)))) %>% 
                    summarise(Mean = mean(Abundance), 
                              Median = median(Abundance), 
                              SD=sd(Abundance)) %>% 
                    pivot_wider(id_cols = feature, 
                                names_from = Phenotype, 
                                values_from = c(Mean, Median, SD)) %>% 
                    left_join(., ancom.res.df, by="feature") 
    
    DA.res.ls[["tabs"]][[set.inst]][[lvl.inst]][[norm.inst]][["summary"]] <- dif.res.tab
    
    write.csv(dif.res.tab, 
              paste0(dir.out, "/", set.inst, "/tabs/", 
                     "DA_", lvl.inst, "_", p.norm, ".png"))
                              
  }
  
}

#-------------------------------------------------------------------------------
# Combine plots 
#-------------------------------------------------------------------------------
w.all <- c(DA.res.ls$plot$phen$Genus$Raw$ntax, 
           DA.res.ls$plot$phen$ASV$Raw$ntax)

g.plot <- DA.res.ls$plot$phen$Genus$Raw$p + 
              theme(legend.position = "none")

g.grid <- plot_grid(g.plot, NULL, 
          ncol = 1, rel_heights = c(w.all[1]/sum(w.all)+ 0.2, 
                                    w.all[2]/sum(w.all)))

comb.plot <- plot_grid(NULL, g.grid, NULL,
                  DA.res.ls$plot$phen$ASV$Raw$p, 
                  rel_widths = c(0.02, 0.38, 0.02, 0.48), ncol = 4, 
                  labels = c("A", "", "B"))

ggsave(filename = paste0(dir.out, "/", set.inst, "/plots/", 
                         "comb_DA_", lvl.inst, "_", p.norm, ".png"), 
       plot = comb.plot, 
       width = 12, 
       height = max(w.all)*0.2 + 1, 
       dpi = 600)

save(file = "out/supp/DA_res.Rdata", list = "DA.res.ls")
