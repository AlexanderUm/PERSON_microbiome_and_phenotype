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
source("R/fix_taxa_names_for_plot.R")
source("R/phy_taxa_filter.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
# Data sets to use
PsSets <- prm.ls$RF$data_set_ps
  
TaxaLvs <- prm.ls$RF$tax_lvls
  
TaxaNorm <- prm.ls$RF$ASV_norm
  
MetabolSets <- prm.ls$RF$metabol_tabs_indRF

CombSets <- prm.ls$RF$sets_combRF

AllFeatSets <- prm.ls$RF$all_feat_sets

SigFeatSets <- prm.ls$RF$sig_feat_sets

TransFun <- prm.ls$RF$metabol_count_trans
  
# Columns to use
MainGrCol <- prm.ls$General$Group_col
  
ParticID <- prm.ls$General$Part_id_col
  
# RF parameters
RfNtrees <- prm.ls$RF$n_trees
  
TaxaMinPrev <- prm.ls$RF$feature_min_prev
  
RfNprem <- prm.ls$RF$n_permut
  
# ??? n.rand <- prm.ls$RF$n_times_shufle_group_col
  
RfImpMaxPval <- prm.ls$RF$max_pval_impotance
  
RfGrSizeBalance <- prm.ls$RF$gr_size_prop

RfAutoTuneMtry <- prm.ls$RF$auto_tune_mtry
  
# Plotting parameters
RfPlotRocSize <- prm.ls$RF$plot_roc_size
  
RfPlotSignifSize <- prm.ls$RF$plot_signif_size

RfTopNSign <- prm.ls$RF$plot_top_n_signif
  
# Results object
DirOut <- prm.ls$RF$out_dir

# Colors 
PlotColors <- c(aes.ls[["col"]]$phen, "Overall" = "black")

# Metadata 
Meta <- data.ls[[PsSets]][["meta"]]


################################################################################
# Extract data for RF 
################################################################################
RfDataFormLs <- list()

RfDataPrmDf <- NULL

#-------------------------------------------------------------------------------
# OTU tables
#-------------------------------------------------------------------------------
TaxaDataGrid <- expand.grid("Tax_lvl" = TaxaLvs, 
                            "Count_norm" = TaxaNorm, 
                            "Min_prevalence" = TaxaMinPrev, 
                            stringsAsFactors = FALSE)
  
if(nrow(TaxaDataGrid) > 0) {
  
  for(i in 1:nrow(TaxaDataGrid)) {
    
    # Variables 
    iLvls <- TaxaDataGrid[i, "Tax_lvl"]
    
    iTaxaNorm <- TaxaDataGrid[i, "Count_norm"]
    
    iTaxaMinPrev <- TaxaDataGrid[i, "Min_prevalence"]
    
    # Objects 
    iName <- paste(TaxaDataGrid[i, ], collapse = "--")
    
    iPs <- data.ls[[PsSets]][["PS"]][[iLvls]][[iTaxaNorm]] 
      
    # Extract data
    iTaxaTab <- iPs %>% 
                  phy_taxa_filter(prev_fraction = iTaxaMinPrev,  
                                    group_col = MainGrCol) %>% 
                    otu_table() %>% 
                    as.matrix() %>% 
                    t() %>% 
                    as.data.frame() %>% 
                    mutate(!!MainGrCol := Meta[[MainGrCol]]) 
      
    # Adjust row names to much rownames of metabolites tables
    rownames(iTaxaTab) <- gsub("CIW1", "", rownames(iTaxaTab))
      
    # Collect data
    RfDataFormLs[["Taxa"]][[iName]] <- iTaxaTab
    
    # Add to data frame 
    RfDataPrmDf <- data.frame(set_type = "Taxa", 
                              set = iName) %>% 
                      rbind(RfDataPrmDf, .)
    
  }
  
}


#-------------------------------------------------------------------------------
# Metabolites 
#-------------------------------------------------------------------------------
VecIdNameChange <- setNames(paste0("met_", ParticID), ParticID)

if (!is.null(MetabolSets)) {
    
    for(i in MetabolSets) { 
      
      iMetabolTab <- data.ls[["Metabol"]][[i]] %>% 
                        setNames(paste0("met_", colnames(.))) %>% 
                        rename(all_of(VecIdNameChange)) %>% 
                        left_join(., Meta[, c(MainGrCol, ParticID)], 
                                  by = ParticID) %>% 
                        filter(!is.na(.[[MainGrCol]])) %>% 
                        column_to_rownames(var = ParticID) 
      
      
      if(!is.null(TransFun)) {
        
        iMetabolTab <- iMetabolTab %>% 
                            mutate(across(dplyr::where(is.numeric), 
                                          TransFun))
        
      }
      
      RfDataFormLs[["Metabol"]][[i]] <- iMetabolTab
      
      # Add to data frame 
      RfDataPrmDf <- data.frame(set_type = "Metabol", 
                                set = i) %>% 
                          rbind(RfDataPrmDf, .)
      
    }
  
}


################################################################################
# Random forest model
################################################################################
RfResLs <- list()

RfDataSigLs <- list()

RfFormula <- as.formula(paste0(MainGrCol, "~ ."))

#-------------------------------------------------------------------------------
# Individual sets 
#-------------------------------------------------------------------------------
for(i in 1:nrow(RfDataPrmDf)) { 
  
  iPrm <- RfDataPrmDf[i, ]
  
  # All features
  iData <- RfDataFormLs[[iPrm[["set_type"]]]][[iPrm[["set"]]]]
  
  # Run model
  iRfResAll <- RF_run_with_par(data_in = iData,
                            RF_formula = RfFormula,
                            group_col = MainGrCol,
                            Mtry_auto_tune = RfAutoTuneMtry,
                            n_trees = RfNtrees,
                            gr_size_prop = RfGrSizeBalance,
                            n_perm = RfNprem)
  
  # Collect data 
  RfResLs[["AllFeat"]][[iPrm[["set_type"]]]][[iPrm[["set"]]]] <- iRfResAll
  
  # Only significant taxa 
  iFeatTax <- RF_extract_importance(iRfResAll$Full) %>% 
                select(-starts_with("MeanDecreaseGini")) %>% 
                filter(if_any(ends_with("unscaled_pval"), ~ . <= RfImpMaxPval)) %>% 
                pull(feature)
  
  # Data with only significant taxa
  iDataSig <- iData[, c(MainGrCol, iFeatTax)]
  
  RfDataSigLs[[iPrm[["set_type"]]]][[iPrm[["set"]]]] <- iDataSig
  
  # Run model
  iRfResSig <- RF_run_with_par(data_in = iDataSig,
                                RF_formula = RfFormula,
                                group_col = MainGrCol,
                                Mtry_auto_tune = RfAutoTuneMtry,
                                n_trees = RfNtrees,
                                gr_size_prop = RfGrSizeBalance,
                                n_perm = RfNprem)
  
  # Collect data 
  RfResLs[["SigFeat"]][[iPrm[["set_type"]]]][[iPrm[["set"]]]] <- iRfResSig
  
}


#-------------------------------------------------------------------------------
# Combined sets 
#-------------------------------------------------------------------------------
for(i in names(CombSets)) { 
  
  # Combine data 
  iTabsToComb <- c(RfDataSigLs$Taxa[CombSets[[i]]$Taxa], 
                   RfDataSigLs$Metabol[CombSets[[i]]$Metabol])
  
  iOverSamp <- iTabsToComb %>% 
                    lapply(function(x){rownames(x)}) %>% 
                    Reduce(intersect, .) %>% 
                    sort()
  
  iTabsToCombFilt <- iTabsToComb %>% 
                        lapply(function(x){x[iOverSamp, ]})
  
  iCombDf <- bind_cols(iTabsToCombFilt) %>% 
                mutate(!!MainGrCol := select(., starts_with(MainGrCol))[[1]]) %>% 
                select(-starts_with(paste0(MainGrCol, "..."))) %>% 
                suppressMessages()
  
  # Run model
  iRfResComb <- RF_run_with_par(data_in = iCombDf,
                               RF_formula = RfFormula,
                               group_col = MainGrCol,
                               Mtry_auto_tune = RfAutoTuneMtry,
                               n_trees = RfNtrees,
                               gr_size_prop = RfGrSizeBalance,
                               n_perm = RfNprem)
  
  # Collect data 
  RfResLs[["SigFeat"]][["Combs"]][[i]] <- iRfResComb
  
}


################################################################################
# Visualize and write results 
################################################################################
AllSetsPrmDf <- bind_rows(bind_cols(Type = "AllFeat", RfDataPrmDf),
                    bind_cols(Type = "SigFeat", RfDataPrmDf), 
                    data.frame(set = names(CombSets), 
                               Type = "SigFeat", 
                               set_type = "Combs"))

RocDataDf <- NULL

ResSummDf <- NULL

ResImpPlotsLs <- list()

for(i in 1:nrow(AllSetsPrmDf)) {
  
  # Variables 
  iType <- AllSetsPrmDf[i, "Type"]
  
  iSetType <- AllSetsPrmDf[i, "set_type"]
  
  iSet <- AllSetsPrmDf[i, "set"]
  
  # Create directories 
  iPathOut <- paste0(DirOut, "/", iSetType, "/", gsub(" ", "", iSet))
  
  dir.create(iPathOut, recursive = TRUE, showWarnings = FALSE)
  
  # Extract results 
  iRfRes <- RfResLs[[iType]][[iSetType]][[iSet]]
  
  # Data frame for combined ROC curves 
  RocDataDf <-  data.frame(sensitivity = iRfRes$ROC_objects$roc$sensitivities, 
                           specificity = iRfRes$ROC_objects$roc$specificities, 
                           AUC = iRfRes$ROC_objects$AUC) %>% 
                    bind_cols(AllSetsPrmDf[i, ]) %>% 
                    arrange(-row_number()) %>% 
                    bind_rows(RocDataDf, .)
                    
  # Summary data 
  ResSummDf <- iRfRes$Summary_df %>% 
                  bind_cols(AllSetsPrmDf[i, ]) %>% 
                  bind_rows(ResSummDf, .)
  
  
  # Features important for classification 
  #-----------------------------------------------------------------------------
  # Extract taxa
  iFeatImp <- RF_extract_importance(rfPermut_results = iRfRes[["Full"]])
  
  # Subset only significant taxa 
  iFeatImpSig <- iFeatImp %>% 
                  select(-starts_with("MeanDecreaseGini")) %>% 
                  filter(if_any(ends_with("unscaled_pval"), ~ . <= RfImpMaxPval)) 
  
  # Sets of features to plot 
  iFeatToPlotLs <- list("AllSig" = iFeatImpSig$feature, 
                        "TopNSig" = head(iFeatImpSig$feature, RfTopNSign), 
                        "TopNAll" = head(iFeatImp$feature, RfTopNSign))
  
  # Plot contributing features 
  for(j in names(iFeatToPlotLs)) { 
    
    jFeatToPlot <- iFeatImp %>% 
                      filter(feature %in% iFeatToPlotLs[[j]]) %>% 
                      mutate(feature = fix_taxa_names_for_plot(feature)) %>% 
                      mutate(feature = gsub("_$|met_", "", feature)) 
    
    if(nrow(jFeatToPlot) > 0) {
      
      # Object to collect plot and dimensions 
      jFeatPlot <- list()
      
      # Create plot using a custom ggplot function
      jFeatPlot[["p"]] <- RF_plot_importance(
                            jFeatToPlot, 
                            cols_to_plot = c("LIR_importance",
                                             "MIR_importance",
                                             "MeanDecreaseAccuracy_importance"), 
                            y_text_face_italic = FALSE)
      
      # Plot dim
      jFeatPlot[["h"]] <- (nrow(jFeatToPlot)*
                             RfPlotSignifSize[["h_coef"]] + 
                             RfPlotSignifSize[["h_add"]])                
      
      jFeatPlot[["w"]] <- RfPlotSignifSize[["w"]]
      
      # Collect plots 
      ResImpPlotsLs[[iSetType]][[gsub(" ", "", iSet)]] <- jFeatPlot 
      
      ggsave(filename = paste0(iPathOut, "/", j, "--", 
                               gsub(" ", "", iSet), ".png"),
             plot = jFeatPlot$p,
             width = jFeatPlot$w,
             height = jFeatPlot$h,
             dpi = 600)
      
    }
    
  }
  
  # Write out results 
  #-----------------------------------------------------------------------------
  # ROC plot 
  ggsave(filename = paste0(iPathOut, "/ROC--", gsub(" ", "", iSet), ".png"), 
         plot = iRfRes$ROC_objects$plot, 
         width = RfPlotRocSize[["w"]], 
         height = RfPlotRocSize[["h"]], 
         dpi = 600)
  
  # Results summary 
  write.csv(x = iRfRes$Summary_df, 
            file = paste0(iPathOut, "/ConfMat--", gsub(" ", "", iSet), ".csv"))
  
  # Features importance table 
  write.csv(x = iFeatImp, 
            file = paste0(iPathOut, "/ImpFeat--", gsub(" ", "", iSet), ".csv"))
  
}


################################################################################
# Combine plots - manual 
################################################################################
SummPlotsLs <- list()

# Function to adjust data 
fit_plot_data_loc <- function(x) {
  
  x %>% 
    filter(!(Type == "SigFeat" & 
               set == "fecal_FA")) %>% 
      mutate(set = gsub("--.*|_metab_all", "", set)) %>% 
      mutate(set = gsub("_", " ", set), 
             Type = case_match(Type, 
                               "AllFeat" ~ "All features", 
                               "SigFeat" ~ "Only significant")) %>% 
      mutate(set = factor(set, levels = unique(set)))
  
}

#-------------------------------------------------------------------------------
# Overall accuracy 
#-------------------------------------------------------------------------------
# Adjust data 
ResSummPlotDf <- ResSummDf %>% 
                    fit_plot_data_loc()
         
# Plot overall accuracy
SummPlotsLs[["A"]] <- RF_plot_accuracy(RF_res_main = ResSummPlotDf, 
                                       x_col_name = "set", 
                                       color_col = "Group") + 
                              scale_color_manual(values = PlotColors) + 
                              facet_grid(. ~ Type, 
                                         scales = "free_x", 
                                         space = "free_x") +
                              ylim(c(25, 90)) + 
                              theme(panel.grid.major.x = element_blank(), 
                                    panel.grid.minor.x = element_blank(), 
                                    axis.text.x = element_text(angle = 45))

#-------------------------------------------------------------------------------
# ROC curves 
#-------------------------------------------------------------------------------
# Adjust for plotting 
# Adjust data 
RocDataPlotDf <- RocDataDf %>% 
                    fit_plot_data_loc() %>% 
                    mutate(set_type = case_match(set_type, 
                                                  "Metabol" ~ "Metabolites", 
                                                  "Combs" ~ "Combined", 
                                                  "Taxa" ~ "Taxa")) %>% 
                    mutate(set_type = factor(set_type, unique(set_type)), 
                           Type = factor(Type, c("Only significant", 
                                                 "All features")))

# Plot ROC Curves 
SummPlotsLs[["B"]] <- ggplot(RocDataPlotDf, 
                             aes(x = specificity, 
                                 y = sensitivity, 
                                 color = set, 
                                 alpha = AUC)) + 
                        geom_line(aes(group = set), size = 1) + 
                        scale_x_reverse() + 
                        facet_wrap(Type~set_type, ncol = 3) + 
                        geom_abline(intercept=1, 
                                    slope=1, 
                                    linetype="dashed") + 
                        theme_bw() + 
                        scale_color_brewer(name = "Features", 
                                           palette = "Dark2") + 
                        theme(panel.grid = element_blank(), 
                              axis.text.x = element_text(angle = 90), 
                              legend.key.height = unit(0.3, 'cm'))

#-------------------------------------------------------------------------------
# Significantly contributing plots 
#-------------------------------------------------------------------------------
ImpFlatLs <- unlist(ResImpPlotsLs, recursive = FALSE)

for(i in names(ImpFlatLs)) {
  
  ImpFlatLs[[i]] <- ImpFlatLs[[i]]$p + 
                          theme(legend.position="none") + 
                          scale_y_discrete(position = "right") + 
                          theme(axis.text.x = element_text(angle = 90), 
                                strip.text = element_text(size = 8))
  
}


#-------------------------------------------------------------------------------
# Combine plots 
#-------------------------------------------------------------------------------
contr.plots <- plot_grid(ImpFlatLs$Metabol.plasma_metab_all,
                         ImpFlatLs$`Taxa.ASV--CSS--0.25`, 
                         ImpFlatLs$Metabol.fecal_FA, 
                         ImpFlatLs$`Taxa.Genus--CSS--0.25`, 
                         ncol = 2, labels = c("C.1", "C.2", "C.3", "C.4"), 
                         align = "v", axis = "lr", 
                         vjust = 1, 
                         scale = 0.925)

contr.legend <- get_legend(ResImpPlotsLs$Combs$`ASV&fecalFA`$p)

contr.plots.l <- plot_grid(contr.plots, contr.legend, rel_widths = c(0.95, 0.05))

# AUC and accuracy 
auc.grid <- plot_grid(SummPlotsLs$A, SummPlotsLs$B, 
                      labels = c("A", "B"), scale = 0.925, 
                      align = "hv", axis = "l")

comb.plot <- plot_grid(auc.grid, contr.plots.l,
                       ncol = 1,
                       align = "v", axis = "l", rel_heights = c(0.35, 0.65))

save_plot(paste0(DirOut, "/comb_plot.png"), 
          comb.plot, base_width = 15, 
          base_height = 9)

# #-------------------------------------------------------------------------------
# save(list = c("rf.all.res.ls"),
#        file = "out/supp/RF.Rdata")
# 
# # Clean environment 
# rm(list = ls())
# gc()


