#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("out/supp/prm.Rdata")
load("out/supp/data.Rdata")

set.seed(prm.ls$General$Seed)

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("R/RF_supp_functions.R")
source("R/randomize_column.R")
source("R/fix_taxa_names_for_plot.R")
source("R/phy_taxa_filter.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
# Data sets to use
ps.set <- prm.ls$RF$data_set_ps
  
tax.lvls <- prm.ls$RF$tax_lvls
  
norms <- prm.ls$RF$ASV_norm
  
metabol.sets <- prm.ls$RF$metabol_tabs_indRF
  
# Columns to use
gr.col <- prm.ls$General$Group_col
  
id.col <- prm.ls$General$Part_id_col
  
# RF parameters
n.trees <- prm.ls$RF$n_trees
  
min.prev <- prm.ls$RF$feature_min_prev
  
n.perm <- prm.ls$RF$n_permut
  
n.rand <- prm.ls$RF$n_times_shufle_group_col
  
imp.max.p <- prm.ls$RF$max_pval_impotance
  
combRF.sets <- prm.ls$RF$sets_combRF
  
gr.size.prop <- prm.ls$RF$gr_size_prop

auto.tune.mtry <- prm.ls$RF$auto_tune_mtry
  
# Plotting parameters
roc.size <- c("h" = 3.5, "w" = 5)
  
sum.p.size <- c("h_prop" = 0.6, "w_coef" = 0.75, "w_add" = 1.75)
  
sig.p.size <- c("h_coef" = 0.2, "h_add" = 2, "w" = 12)
  
n.top.p <- 10
  
# Results object
out.dir <- prm.ls$RF$out_dir
  
rf.all.res.ls <- list()
  
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  
# Aestetics
accu.p.color <- c(aes.ls[["col"]]$phen, "Overall" = "black")
  
  
################################################################################
# Extract data for RF 
################################################################################
rf.data.ls <- list()
  
#-------------------------------------------------------------------------------
# OTU tables
#-------------------------------------------------------------------------------
ps.grid <- expand.grid("Samples_set" = ps.set, 
                         "Tax_lvl" = tax.lvls, 
                         "Count_norm" = norms, 
                         "Min_prevalence" = min.prev, 
                         stringsAsFactors = FALSE)
  
  
for(i.grid in 1:nrow(ps.grid)) {
    
    i.set <- ps.grid[i.grid, "Samples_set"]
    
    i.lvl <- ps.grid[i.grid, "Tax_lvl"]
    
    i.norm <- ps.grid[i.grid, "Count_norm"]
    
    i.prev <- ps.grid[i.grid, "Min_prevalence"]
    
    
    ps.inst <- data.ls[[i.set]][["PS"]][[i.lvl]][[i.norm]] 
    
    meta.inst <- data.ls[[i.set]][["meta"]]
    
    if(!is.null(ps.inst)) {
      
      i.name <- paste(ps.grid[i.grid, ], collapse = "--")
      
      # Extract data
      otu.tab.inst <- ps.inst %>% 
                          phy_taxa_filter(., prev_fraction = i.prev,  
                                          group_col = gr.col) %>% 
                          otu_table(.) %>% 
                          as.matrix() %>% 
                          t() %>% 
                          as.data.frame() %>% 
                          mutate(!!gr.col := meta.inst[[gr.col]]) %>% 
                          setNames(gsub("-|\\+", "_", names(.)))
      
      num.col <- grepl("^[0-9]", colnames(otu.tab.inst))
      
      names(otu.tab.inst)[num.col] <- paste0("p_", names(otu.tab.inst)[num.col])
      
      # Adjust row names to much rownames of metabolites tables
      rownames(otu.tab.inst) <- gsub("CIW1", "", rownames(otu.tab.inst))
      
      rf.data.ls[[i.name]] <- otu.tab.inst
      
  }
}
  
  
#-------------------------------------------------------------------------------
# Metabolites 
#-------------------------------------------------------------------------------
if (!is.null(metabol.sets)) {
    
    for(i.metabol in metabol.sets) {
      
      metabol.inst <- data.ls[["Metabol"]][[i.metabol]] %>% 
                        left_join(., 
                                  data.ls[[i.set]][["meta"]][, c(gr.col, id.col)], 
                                  by = id.col) %>% 
                        filter(!is.na(.[[gr.col]])) %>% 
                        column_to_rownames(var = id.col) %>% 
                        setNames(gsub("-|\\:|\\(|\\)|\\/| ", "_", names(.)))
      
      num.col <- grepl("[0-9]_", colnames(metabol.inst))
      
      names(metabol.inst)[num.col] <- paste0("met_", names(metabol.inst)[num.col])
      
      rf.data.ls[[paste0("metabol--",i.metabol)]] <- metabol.inst
  }
}
  
  
################################################################################
# Random forest model
################################################################################
#-------------------------------------------------------------------------------
# Run RF per dataset 
#-------------------------------------------------------------------------------
rf.formula <- as.formula(paste0(gr.col, "~."))
  
rf.res.ls <- list()
  
rf.res.df <- NULL
  
for(i.df in names(rf.data.ls)) {
  
  path.out <- paste0(out.dir, "/", i.df)
  
  dir.create(path.out, recursive = TRUE, showWarnings = FALSE)
    
  #-----------------------------------------------------------------------------
  # Full RF model 
  #-----------------------------------------------------------------------------
  data.inst <- rf.data.ls[[i.df]]
  
  rf.res.inst <- RF_run_with_par(data_in = data.inst, 
                                   RF_formula = rf.formula, 
                                   group_col = gr.col, 
                                   Mtry_auto_tune = auto.tune.mtry, 
                                   n_trees = n.trees, 
                                   gr_size_prop = gr.size.prop,
                                   n_perm = n.perm)
  
  rf.res.df <- rf.res.inst$Summary_df %>% 
                      mutate(Data_set = i.df, 
                             Feature_Set = "All features") %>% 
                      bind_rows(rf.res.df, .)
    
  write.csv(x = rf.res.inst$Summary_df, 
            file = paste0(path.out, "/sum_all__", i.df, ".csv"))
    
  ggsave(filename = paste0(path.out, "/roc_all__", i.df, ".png"), 
         plot = rf.res.inst$ROC_objects$plot, 
         width = roc.size[["w"]], 
         height = roc.size[["h"]], 
         dpi = 600)
  
  #-----------------------------------------------------------------------------
  # Features important for classification 
  #-----------------------------------------------------------------------------
  # Extract taxa
  imp.inst <- RF_extract_importance(rfPermut_results = rf.res.inst[["Full"]])
  
  write.csv(x = imp.inst, 
            file = paste0(path.out, "/features_all__", i.df, ".csv"))
  
  # Subset only significant taxa 
  imp.inst.sig <- imp.inst %>% 
                    select(-starts_with("MeanDecreaseGini")) %>% 
                    filter(if_any(ends_with("unscaled_pval"), ~ . <= imp.max.p)) 
  
  
  #-----------------------------------------------------------------------------                 
  # Plot all sig taxa 
  if(nrow(imp.inst.sig) > 0) {
      
     imp.inst.sig.p <- imp.inst.sig %>% 
                        mutate(feature = fix_taxa_names_for_plot(feature)) %>% 
                        mutate(feature = sub("_$", "", feature))
    
     imp.plot.sig <- list("p" = RF_plot_importance(imp.inst.sig.p, 
                                                   cols_to_plot = c("LIR_importance",
                                                                    "MIR_importance",
                                                                    "MeanDecreaseAccuracy_importance"), 
                                                   y_text_face_italic = FALSE), 
                          "h" = nrow(imp.inst.sig)*sig.p.size[["h_coef"]] + 
                                                   sig.p.size[["h_add"]], 
                          "w" = sig.p.size[["w"]])
     
     write.csv(x = imp.inst.sig.p, 
               file = paste0(path.out, "/features_sig__", i.df, ".csv"))
     
     ggsave(filename = paste0(path.out, "/features_sig__", i.df, ".png"), 
            plot = imp.plot.sig$p, 
            width = imp.plot.sig$w, 
            height = imp.plot.sig$h, 
            dpi = 600)
     
  } else { 
    
    imp.plot.sig <- NULL 
    
  }

  
  #-----------------------------------------------------------------------------
  # Plot top 20 taxa 
  if(nrow(imp.inst) < n.top.p) {
    
    top.cut <- nrow(imp.inst)
    
  } else { 
    
    top.cut <- n.top.p
    
  }
  
  imp.plot.topN.p <- imp.inst %>% 
                      arrange(desc(MeanDecreaseAccuracy_importance)) %>% 
                      slice(1:top.cut) %>% 
                      mutate(feature = fix_taxa_names_for_plot(feature)) %>% 
                      mutate(feature = sub("_$", "", feature)) %>% 
                      RF_plot_importance(., 
                                         cols_to_plot = c("LIR_importance",
                                                          "MIR_importance",
                                                          "MeanDecreaseAccuracy_importance"), 
                                         y_text_face_italic = FALSE)
  
  imp.plot.topN <- list("p" = imp.plot.topN.p, 
                        "h" = top.cut*sig.p.size[["h_coef"]] + sig.p.size[["h_add"]], 
                        "w" = sig.p.size[["w"]])
  
  ggsave(filename = paste0(path.out, "/top_", top.cut, "_features__", i.df, ".png"), 
         plot = imp.plot.topN$p, 
         width = imp.plot.topN$w, 
         height = imp.plot.topN$h, 
         dpi = 600)
  
  
  #-----------------------------------------------------------------------------
  # Top N significant 
  if(nrow(imp.inst.sig) > 0) {
    
    if(nrow(imp.inst.sig) < n.top.p) {
      
      top.cut.sig <- nrow(imp.inst.sig)
      
    } else {
      
      top.cut.sig <- n.top.p
    }
    
    imp.plot.sig.topN.p <- imp.inst.sig %>% 
                            arrange(desc(MeanDecreaseAccuracy_importance)) %>% 
                            slice(1:top.cut) %>% 
                            mutate(feature = fix_taxa_names_for_plot(feature)) %>% 
                            mutate(feature = sub("_$", "", feature)) %>% 
                            RF_plot_importance(., y_text_face_italic = FALSE, 
                                               cols_to_plot = c("LIR_importance",
                                                                "MIR_importance",
                                                                "MeanDecreaseAccuracy_importance"))
    
    imp.plot.sig.topN <- list("p" = imp.plot.sig.topN.p, 
                              "h" = top.cut.sig*sig.p.size[["h_coef"]] + 
                                                sig.p.size[["h_add"]], 
                              "w" = sig.p.size[["w"]])
    
    ggsave(filename = paste0(path.out, "/top_sig_", top.cut, "_features__", i.df, ".png"), 
           plot = imp.plot.sig.topN$p, 
           width = imp.plot.sig.topN$w, 
           height = imp.plot.sig.topN$h, 
           dpi = 600)
    
  } else {
    
    imp.plot.sig.topN <- NULL
    
  }
  
  
  #-----------------------------------------------------------------------------
  # RF with significant features
  #-----------------------------------------------------------------------------
  if(nrow(imp.inst.sig) > 0) {
    
  data.inst.sig <- data.inst %>% 
                    select(all_of(c(gr.col, imp.inst.sig[["feature"]])))
  
  rf.res.inst.sig <- RF_run_with_par(data_in = data.inst.sig, 
                                     RF_formula = rf.formula, 
                                     group_col = gr.col, 
                                     Mtry_auto_tune = auto.tune.mtry, 
                                     n_trees = n.trees, 
                                     gr_size_prop = gr.size.prop,
                                     n_perm = n.perm) 
  
  # Collect results   
  rf.res.df <- rf.res.inst.sig$Summary_df %>% 
                        mutate(Data_set = i.df, 
                               Feature_Set = "Significant features") %>% 
                        bind_rows(rf.res.df, .)
  
  write.csv(x = rf.res.inst.sig$Summary_df, 
            file = paste0(path.out, "/sum_sig__", i.df, ".csv"))
  
  ggsave(filename = paste0(path.out, "/roc_sig__", i.df, ".png"), 
         plot = rf.res.inst.sig$ROC_objects$plot, 
         width = roc.size[["w"]], 
         height = roc.size[["h"]], 
         dpi = 600)
  
  } else { 
    
    rf.res.inst.sig <- NULL
  
  }
  
  # Collect results   
  rf.res.ls[["RF"]][["All"]][[i.df]] <- rf.res.inst
  
  rf.res.ls[["RF"]][["Sig"]][[i.df]] <- rf.res.inst.sig
  
  rf.res.ls[["features"]][[i.df]][["Sig_tab"]] <- imp.inst.sig
  
  rf.res.ls[["features"]][[i.df]][["plot_sig_all"]] <- imp.plot.sig
  
  rf.res.ls[["features"]][[i.df]][["plot_topN"]] <- imp.plot.topN
  
  rf.res.ls[["features"]][[i.df]][["plot_sig_topN"]] <- imp.plot.sig.topN
  
}
  

#-------------------------------------------------------------------------------
# Combined data sets 
#-------------------------------------------------------------------------------
if(!is.null(combRF.sets)) {
  
  for(i.comb in names(combRF.sets)) { 
    
    #---------------------------------------------------------------------------
    # Prepare data
    names.inst <- combRF.sets[[i.comb]] %>% 
                    paste0(., collapse = "|") %>% 
                    grep(., names(rf.data.ls), value = TRUE)

    # Extract needed data sets
    data.inst <- rf.data.ls[names.inst] 
    
    # Find shared samples 
    samp.inst <- lapply(data.inst, function(x){rownames(x)}) %>% 
                  Reduce(intersect, .)
    
    # Extract significant features 
    sig.tax <- rf.res.ls$features[names.inst] %>% 
               lapply(., function(x){x$Sig_tab$feature}) %>% 
               unlist()
    
    # Subset significant taxa 
    data.inst.f <- lapply(data.inst, 
                           function(x){x[samp.inst, 
                                         colnames(x) %in% sig.tax]}) %>% 
                    bind_cols() %>% 
                    mutate(!!gr.col := data.inst[[1]][samp.inst, gr.col])
    
    
    #---------------------------------------------------------------------------
    # Random Forest 
    rf.res.inst.comb <- RF_run_with_par(data_in = data.inst.f, 
                                       RF_formula = rf.formula, 
                                       group_col = gr.col, 
                                       Mtry_auto_tune = auto.tune.mtry, 
                                       n_trees = n.trees, 
                                       gr_size_prop = gr.size.prop,
                                       n_perm = n.perm) 
    
    #---------------------------------------------------------------------------
    # Collect results   
    rf.res.df <- rf.res.inst.comb$Summary_df %>% 
                      mutate(Data_set = i.comb, 
                             Feature_Set = "Combined features") %>% 
                      bind_rows(rf.res.df, .)
    
    rf.res.ls[["RF"]][["Comb"]][[i.comb]] <- rf.res.inst.comb
    
    #---------------------------------------------------------------------------
    # Write out results 
    short.name.inst <- gsub(" ", "", i.comb)
    
    path.out <- paste0(out.dir, "/", short.name.inst)
    
    dir.create(path.out, recursive = TRUE, showWarnings = FALSE)
    
    
    write.csv(x = rf.res.inst.comb$Summary_df, 
              file = paste0(path.out, "/sum_sig__", 
                            short.name.inst, ".csv"))
    
    ggsave(filename = paste0(path.out, "/roc_sig__", 
                             short.name.inst, ".png"), 
           plot = rf.res.inst.comb$ROC_objects$plot, 
           width = roc.size[["w"]], 
           height = roc.size[["h"]], 
           dpi = 600)
    
    #-----------------------------------------------------------------------------                 
    # Features contribution 
    imp.inst <- RF_extract_importance(rfPermut_results = rf.res.inst.comb[["Full"]])
    
    write.csv(x = imp.inst, 
              file = paste0(path.out, "/features_all__", 
                            short.name.inst, ".csv"))
    
    # Subset only significant taxa 
    imp.inst.sig <- imp.inst %>% 
                      select(-starts_with("MeanDecreaseGini")) %>% 
                      filter(if_any(ends_with("unscaled_pval"), ~ . <= imp.max.p)) 
    
    if(nrow(imp.inst.sig) > 0) {
      
      imp.inst.sig.p <- imp.inst.sig %>% 
                          mutate(feature = fix_taxa_names_for_plot(feature)) %>% 
                          mutate(feature = sub("_$", "", feature))
      
      imp.plot.sig <- list("p" = RF_plot_importance(imp.inst.sig.p, 
                                                    cols_to_plot = c("LIR_importance",
                                                                     "MIR_importance",
                                                                     "MeanDecreaseAccuracy_importance"), 
                                                    y_text_face_italic = FALSE), 
                           "h" = nrow(imp.inst.sig)*sig.p.size[["h_coef"]] + 
                                                    sig.p.size[["h_add"]], 
                           "w" = sig.p.size[["w"]])
      
      write.csv(x = imp.inst.sig.p, 
                file = paste0(path.out, "/features_sig__", short.name.inst, ".csv"))
      
      ggsave(filename = paste0(path.out, "/features_sig__", short.name.inst, ".png"), 
             plot = imp.plot.sig$p, 
             width = imp.plot.sig$w, 
             height = imp.plot.sig$h, 
             dpi = 600) 
    }
  }
}
  

################################################################################
# Plots 
#-------------------------------------------------------------------------------
# Accuracy plot 
#-------------------------------------------------------------------------------
# Adjust Data frame for plotting 
rf.res.df <- rf.res.df %>%
                mutate(Data_set = sub(".*?--", "", Data_set)) %>%
                mutate(Data_set = sub("--.*", "", Data_set)) %>%
                mutate(Data_set = sub("metab_all", "", Data_set)) %>%
                mutate(Data_set = sub("_", " ", Data_set)) %>%
                mutate(Data_set = factor(Data_set, unique(Data_set)), 
                       Feature_Set = factor(Feature_Set, 
                                            levels = c("All features", 
                                                       "Significant features", 
                                                       "Combined features"))) 

# Plot with a custom ggplot2 function             
accu.p <- RF_plot_accuracy(RF_res_main = rf.res.df, 
                           x_col_name = "Data_set", 
                           color_col = "Group") + 
              scale_color_manual(values = accu.p.color) + 
              facet_grid(~ Feature_Set, scales = "free_x") +
              ylim(c(25, 90)) + 
              theme(panel.grid.major.x = element_blank(), 
                    panel.grid.minor.x = element_blank(), 
                    axis.text.x = element_text(angle = 25))

# Write results 
write.csv(x = rf.res.df, 
          file = paste0(out.dir, "/res_sum_comb_all.csv"))

ggsave(paste0(out.dir, "/accuracy_plot.png"), 
          accu.p, 
          width = roc.size[["w"]]*3, 
          height = roc.size[["h"]])


#-------------------------------------------------------------------------------
# AUC plots 
# Prepare data 
roc.df <- NULL

for(i1.res in names(rf.res.ls$RF)) {

  for(i2.res in names(rf.res.ls$RF[[i1.res]])) {
    
    roc.inst <- rf.res.ls$RF[[i1.res]][[i2.res]]$ROC_objects$roc
    
    rf.res.ls$RF[[i1.res]][[i2.res]]$ROC_objects$AUC
    
    roc.df <- data.frame(sensitivity = roc.inst$sensitivities, 
                         specificity = roc.inst$specificities, 
                         Feature_Set = i1.res, 
                         Data_set = i2.res, 
                         AUC = rf.res.ls$RF[[i1.res]][[i2.res]]$ROC_objects$AUC) %>% 
                 arrange(-row_number()) %>% 
                 bind_rows(roc.df, .)
    
  }
}

# Adjust for ploting 
roc.df <- roc.df %>% 
            mutate(Data_set = sub(".*?--", "", Data_set)) %>%
            mutate(Data_set = sub("--.*", "", Data_set)) %>%
            mutate(Data_set = sub("metab_all", "", Data_set)) %>%
            mutate(Data_set = sub("_", " ", Data_set)) %>%
            mutate(Data_set = factor(Data_set, unique(Data_set)),
                   Feature_Set = recode(Feature_Set, 
                                         All = "All features", 
                                         Sig = "Significant features", 
                                         Comb = "Combined features")) %>% 
            mutate(Feature_Set = factor(Feature_Set, 
                                        levels = c("All features", 
                                                   "Significant features", 
                                                   "Combined features")))


roc.comb.p <- ggplot(roc.df, aes(x = specificity, 
                                 y = sensitivity, 
                                 color = Data_set, alpha = AUC)) + 
                        geom_line(aes(group = Data_set), size = 1) + 
                        scale_x_reverse() + 
                        facet_grid(~Feature_Set) + 
                        geom_abline(intercept=1, 
                                    slope=1, 
                                    linetype="dashed") + 
                        theme_bw() + 
                        scale_color_brewer(name = "Features", 
                                           palette = "Dark2") + 
                        theme(panel.grid = element_blank(), 
                              axis.text.x = element_text(angle = 90), 
                              legend.key.height = unit(0.3, 'cm'))

ggsave(paste0(out.dir, "/roc_comb_plot.png"), 
          roc.comb.p, 
          width = roc.size[["w"]]*3, 
          height = roc.size[["h"]])

#-------------------------------------------------------------------------------
# Significantly contributing plots 
sig.plots.ls <- list()

sig.tax.df <- NULL

for(i1.res in names(rf.res.ls$features)) {
  
  name.inst <- sub(".*?--", "", i1.res) %>% 
               sub("--.*|_metab_all", "", .) 
  
  sig.plots.ls[[name.inst]] <- rf.res.ls$features[[i1.res]]$plot_sig_topN$p + 
                                theme(legend.position="none") + 
                                scale_y_discrete(position = "right") + 
                                theme(axis.text.x = element_text(angle = 90), 
                                      strip.text = element_text(size = 8))
}

#-------------------------------------------------------------------------------
# Combine plots into panel 
# Features
contr.plots <- plot_grid(sig.plots.ls$plasma,
                         sig.plots.ls$ASV, 
                         sig.plots.ls$fecal_FA, 
                         sig.plots.ls$Genus, 
                  ncol = 2, labels = c("C.1", "C.2", "C.3", "C.4"), 
                  align = "v", axis = "lr", 
                  vjust = 1, 
                  scale = 0.925)

contr.legend <- get_legend(rf.res.ls$features$`phen--Genus--CSS--0.25`$plot_sig_topN$p)

contr.plots.l <- plot_grid(contr.plots, contr.legend, rel_widths = c(0.95, 0.05))


# AUC and accuracy 
auc.grid <- plot_grid(accu.p, roc.comb.p, 
                      labels = c("A", "B"), scale = 0.925, 
                      align = "hv", axis = "l")

comb.plot <- plot_grid(auc.grid, contr.plots.l,
              ncol = 1,
              align = "v", axis = "l", rel_heights = c(0.35, 0.65))

save_plot(paste0(out.dir, "/comb_plot.png"), 
          comb.plot, base_width = 14, 
          base_height = 9)

#-------------------------------------------------------------------------------
save(list = c("rf.all.res.ls"),
       file = "out/supp/RF.Rdata")

