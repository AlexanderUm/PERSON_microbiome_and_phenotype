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
source("R/phy_alpha.R")
source("R/phy_dists_ls.R")
source("R/rda_extract_for_plot.R")
source("R/plot_extracted_rda_data.R")


#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
# Data sets to use
ps.set <- prm.ls$Beta$data_set_ps

tax.lvls <- prm.ls$Beta$Tax_lvl

norms <- prm.ls$Beta$Norm

# Columns to use
gr.col <- prm.ls$General$Group_col

id.col <- prm.ls$General$Part_id_col

# beta parameters
dists <- prm.ls$Beta$distances

n.perm <- prm.ls$Beta$n_perm

db.form <- prm.ls$Beta$Formula

out.path <- prm.ls$Beta$out_path

# Results object
beta.res.ls <- list()


################################################################################
# Run alpha diversity analysis
################################################################################
db.form.in <- paste0("dist.inst ~ ", db.form)

for(i.set in ps.set) { 
  
  # Extract data
  ps.inst <- data.ls[[i.set]][["PS"]][[tax.lvls]][[norms]]
  
  meta.inst <- data.ls[[i.set]][["meta"]] 
  
  # Calculate distance 
  dists.ls <- phy_dist_ls(phylo = ps.inst, dists = dists)
  
  names(dists.ls) <- names(dists)
  
  # Create folders 
  dir.create(paste0(out.path, "/", i.set, "/tabs"), 
             recursive = TRUE, showWarnings = FALSE)
  
  dir.create(paste0(out.path, "/", i.set, "/plots"), 
             recursive = TRUE, showWarnings = FALSE)
  
  #-----------------------------------------------------------------------------
  # Run PERMANOVA
  #-----------------------------------------------------------------------------
  res.adon.df <- NULL
  
  for(i.dist in names(dists.ls)) {
    
      dist.inst <- dists.ls[[i.dist]]
    
      res.adon <- adonis2(formula = as.formula(db.form.in), 
                          data = meta.inst, 
                          by = "terms", 
                          permutations = n.perm, 
                          parallel = 4) %>% 
                      tidy() %>% 
                      mutate(Distance = i.dist,
                             Data_set = i.set,
                             Formula = db.form.in) 
      
      res.adon.df <- bind_rows(res.adon.df, res.adon)
      
      beta.res.ls[[i.set]][["adonis2"]][[i.dist]] <- res.adon
      
      write.csv(x = res.adon, 
                file = paste0(out.path, "/", i.set, "/tabs", 
                             "/adonis_", i.set, "_", i.dist, ".csv"))
  }
  
  
  #-----------------------------------------------------------------------------
  # Plot RDAs
  #-----------------------------------------------------------------------------
  # Extract data for plot
  rda.p.df <- rda_extract_for_plot(dists_ls = dists.ls, 
                                   metadata = meta.inst, 
                                   form = db.form)
  
  #-----------------------------------------------------------------------------
  # Add significance text data 
  #-----------------------------------------------------------------------------
  term.plot <- gr.col
  
  for(i in names(rda.p.df)) {
    
    # Extract coordinates 
    i.ax <- rda.p.df[[i]][["main"]]
    
    x.text <- min(i.ax[[1]]) 
    
    y.text <- (max(i.ax[[2]]) + (max(i.ax[[2]]) - min(i.ax[[2]]))*0.05)
    
    # Extract label text
    sig.text <- res.adon.df %>% 
                  filter(Distance == i, term == term.plot) %>% 
                  mutate(R2_text = ifelse(round(R2, 3) == 0, 
                                          "R^2<0.001", 
                                          paste0("R^2==", 
                                                 sprintf("%.3f", round(R2, 3)))), 
                         p_text = ifelse(round(p.value, 3) == 0, 
                                          "P<0.001", 
                                          paste0("P==", 
                                                 sprintf("%.3f", round(p.value, 3)))))
    
    # Combine into a data frame 
    text.df <- data.frame("Distance" = i, 
                       "x_coord" = as.numeric(x.text), 
                       "y_coord" = as.numeric(y.text), 
                       "lab_text" = paste0(sig.text[["p_text"]], 
                                           "~~(", sig.text[["R2_text"]], ")"))
    
    rda.p.df[[i]][["stat_lable"]] <- text.df
    
  }
  
  
  #-----------------------------------------------------------------------------
  # Plot RDAs
  #-----------------------------------------------------------------------------
  dbRDA.p.inst <- plot_extracted_rda_data(extracted_data = rda.p.df, 
                                           group_col = gr.col, 
                                           add_stat_text = TRUE, 
                                           sig_text_size = 3.1,
                                           add_elepses = TRUE, 
                                           color_vec = aes.ls$col[[i.set]])
  
  
  #-----------------------------------------------------------------------------
  # Combine alpha and beta plots
  #-----------------------------------------------------------------------------
  alpha.p <- alpha.res.ls[[i.set]][["plot"]] + 
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
                                            0.2 + (length(levels(meta.inst[[gr.col]]))-2)*0.05, 
                                            0.02, 0.725),
                             labels = c("A", "", "", "B"))
  
  save_plot(paste0(out.path, "/", i.set, "/plots/", i.set, "_dbrda_beta.png"), 
            plot = plot(div.plot.comb), 
            base_width = 12, base_height = 7, units = "in", dpi = 600)
  
}

# Clean environment 
rm(list = ls())
gc()

