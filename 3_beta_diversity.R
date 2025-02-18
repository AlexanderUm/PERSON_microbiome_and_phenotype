#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")
load("out/supp/2_alpha.Rdata")

set.seed(prm.ls[["General"]][["Seed"]])

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("R/phy_taxa_filter.R")
source("R/phy_dist_ls.R")
source("R/rda_extract_for_plot.R")
source("R/plot_extracted_rda_data.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
MainVar <- prm.ls$Beta$main_var

CoVars <- prm.ls$Beta$fix_covar

TaxLvl <- prm.ls$Beta$Tax_lvl

PhyNorm <- prm.ls$Beta$Norm

Nperm <- prm.ls$Beta$n_perm

Dists <- prm.ls$Beta$distances

OutDir <- prm.ls$Beta$out_path

# Results object
beta.res.ls <- list()


################################################################################
# Run beta diversity analysis
################################################################################
for(i in prm.ls$Beta$data_set_ps) { 
  
  # Extract data
  ps.inst <- data.ls[[i]][["PS"]][[TaxLvl]][[PhyNorm]]
  
  meta.inst <- data.ls[[i]][["meta"]] %>% 
                  select(all_of(c(CoVars[[i]], 
                                  MainVar))) %>% 
                  mutate(across(where(is.numeric), 
                                function(x){replace(x, is.na(x), 
                                                    median(x, na.rm = TRUE))}))
  
  # Formula for dbRDA
  db.form.in <- paste0("dist.inst ~ ", 
                       paste(CoVars[[i]], collapse = " + "), 
                       " + ", MainVar)
  
  # Calculate distance 
  DistsLs <- phy_dist_ls(phylo = ps.inst, dists = Dists)
  
  names(DistsLs) <- names(Dists)
  
  # Create folders 
  dir.create(paste0(OutDir, "/", i, "/tabs"), 
             recursive = TRUE, showWarnings = FALSE)
  
  dir.create(paste0(OutDir, "/", i, "/plots"), 
             recursive = TRUE, showWarnings = FALSE)
  
  #-----------------------------------------------------------------------------
  # Run PERMANOVA
  #-----------------------------------------------------------------------------
  res.adon.df <- NULL
  
  for(j in names(DistsLs)) {
    
      dist.inst <- DistsLs[[j]]
    
      res.adon <- adonis2(formula = as.formula(db.form.in), 
                          data = meta.inst, 
                          by = "terms", 
                          permutations = Nperm, 
                          parallel = 4) %>% 
                      tidy() %>% 
                      mutate(Distance = j,
                             Data_set = i,
                             Formula = db.form.in) 
      
      res.adon.df <- bind_rows(res.adon.df, res.adon)
      
      beta.res.ls[[i]][["adonis2"]][[j]] <- res.adon
      
      write.csv(x = res.adon, 
                file = paste0(OutDir, "/", i, "/tabs", 
                             "/adonis_", i, "_", j, ".csv"))
  }
  
  
  #-----------------------------------------------------------------------------
  # Plot RDAs
  #-----------------------------------------------------------------------------
  # Extract data for plot
  rda.p.df <- rda_extract_for_plot(dists_ls = DistsLs, 
                                   metadata = meta.inst, 
                                   form = gsub(".*~", "", db.form.in))
  
  #-----------------------------------------------------------------------------
  # Add significance text data 
  #-----------------------------------------------------------------------------
  term.plot <- MainVar
  
  for(j in names(rda.p.df)) {
    
    # Extract coordinates 
    i.ax <- rda.p.df[[j]][["main"]]
    
    x.text <- min(i.ax[[1]]) 
    
    y.text <- (max(i.ax[[2]]) + (max(i.ax[[2]]) - min(i.ax[[2]]))*0.05)
    
    # Extract label text
    sig.text <- res.adon.df %>% 
                  filter(Distance == j, term == term.plot) %>% 
                  mutate(R2_text = ifelse(round(R2, 3) == 0, 
                                          "R^2<0.001", 
                                          paste0("R^2==", 
                                                 sprintf("%.3f", round(R2, 3)))), 
                         p_text = ifelse(round(p.value, 3) == 0, 
                                          "P<0.001", 
                                          paste0("P==", 
                                                 sprintf("%.3f", round(p.value, 3)))))
    
    # Combine into a data frame 
    text.df <- data.frame("Distance" = j, 
                           "x_coord" = as.numeric(x.text), 
                           "y_coord" = as.numeric(y.text), 
                           "lab_text" = paste0(sig.text[["p_text"]], 
                                               "~~(", sig.text[["R2_text"]], ")"))
    
    rda.p.df[[j]][["stat_lable"]] <- text.df
    
  }
  
  
  #-----------------------------------------------------------------------------
  # Plot RDAs
  #-----------------------------------------------------------------------------
  dbRDA.p.inst <- plot_extracted_rda_data(extracted_data = rda.p.df, 
                                           group_col = MainVar, 
                                           add_stat_text = TRUE, 
                                           sig_text_size = 3.1,
                                           add_elepses = TRUE, 
                                           color_vec = aes.ls$col[[i]])
  
  
  #-----------------------------------------------------------------------------
  # Combine alpha and beta plots
  #-----------------------------------------------------------------------------
  alpha.p <- alpha.res.ls[[i]][["plot"]] + 
                  theme(legend.position="none")
  
  alpha.gr <- plot_grid(NULL, 
                        alpha.p,
                        nrow = 2, 
                        rel_heights = c(0.03, 0.97))
  
  div.plot.comb <- plot_grid(NULL,
                             alpha.gr,
                             NULL,
                             dbRDA.p.inst$Comb, 
                             nrow = 1, 
                             rel_widths = c(0.01, 
                                            0.2 + (length(levels(meta.inst[[MainVar]]))-2)*0.05, 
                                            0.02, 0.725),
                             labels = c("A", "", "", "B"))
  
  save_plot(paste0(OutDir, "/", i, "/plots/", i, "_dbrda_beta.png"), 
            plot = plot(div.plot.comb), 
            base_width = 12, base_height = 7, units = "in", dpi = 600)
  
}

# Clean environment 
rm(list = ls())
gc()

